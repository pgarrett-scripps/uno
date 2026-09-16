//! Phase 3: match single-run features to theoretical peptides by precursor mass.
//!
//! For each feature we hypothesise that koth may have mis-assigned the
//! monoisotopic peak by a small number of C13 neutrons. For each isotope offset
//! `o`, the implied peptide mass is `feature.mass - o * NEUTRON`; we collect all
//! peptides within `ppm_tolerance` of that center. Each (feature, peptide,
//! offset) is one Peptide-Feature Match (PFM), the unit FDR is later applied to.

use crate::database::PeptideDb;
use crate::input::FeatureMatrix;
use crate::mass::NEUTRON;
use rayon::prelude::*;

#[derive(Debug, Clone, serde::Deserialize)]
#[serde(default)]
pub struct MatchConfig {
    /// Mass tolerance in ppm for the feature<->peptide match.
    pub ppm_tolerance: f64,
    /// Isotope-error offsets to try. `o` means "koth's mono is o neutrons heavier
    /// than the true mono". `[0]` disables isotope-error search.
    pub isotope_offsets: Vec<i8>,
    /// Skip features below this charge. Singly-charged features are usually
    /// non-peptide contaminants (and CCS predictors don't model them).
    pub min_charge: u8,
    /// Skip features built from fewer than this many isotopic peaks (0 = off).
    /// The binomial background `p` is set by FEATURE density x mass window, not
    /// by database size, so on an instrument with poor mass accuracy the only way
    /// to keep `p` in range is to search fewer, better-characterised features.
    /// On timsTOF (mass sigma ~3 ppm vs ~0.55 on Orbitrap) this is the difference
    /// between 0 and ~300 proteins: 114,659 features -> p 0.554 -> 0 proteins;
    /// gating to `min_isotopes = 4` + `min_combined_score = 0.85` leaves 17,587
    /// -> p 0.206 -> 199-368 proteins at 1% honest entrapment FDP. The isotope
    /// count is the load-bearing half — `min_combined_score` alone (36,753
    /// features) still yields 0. Default 0 keeps Orbitrap byte-identical.
    pub min_isotopes: u32,
    /// Skip features whose `combined_score` (isotope x coelution) is below this
    /// (0 = off). Use with `min_isotopes`; see that field for the measurements.
    pub min_combined_score: f64,
}

impl Default for MatchConfig {
    fn default() -> Self {
        MatchConfig {
            ppm_tolerance: 10.0,
            isotope_offsets: vec![-1, 0, 1],
            min_charge: 1,
            min_isotopes: 0,
            min_combined_score: 0.0,
        }
    }
}

/// One peptide-feature match.
///
/// Field order + widths are chosen to keep the struct at 24 bytes (down from 32):
/// the `f64` leads (8-byte aligned), the two indices are `u32` (feature and
/// peptide counts are well under 2^32 — 211k features, ~1.2M peptides on run10),
/// and the two 1-byte tails pack into the trailing padding. At ~17M reduced PFMs
/// this saves ~136 MB versus `usize` indices. Values are exact, so output is
/// unchanged.
#[derive(Debug, Clone)]
pub struct Pfm {
    /// Signed mass error (peptide - offset-corrected feature center), in ppm.
    pub mass_error_ppm: f64,
    pub feature_idx: u32,
    pub peptide_idx: u32,
    pub isotope_offset: i8,
    pub decoy: bool,
}

/// Feature-quality gate shared by both match paths.
pub(crate) fn passes(feature: &crate::input::Feature, cfg: &MatchConfig) -> bool {
    feature.charge >= cfg.min_charge
        && feature.combined_score >= cfg.min_combined_score
        && (cfg.min_isotopes == 0
            || feature.n_isotopes.is_none_or(|n| n >= cfg.min_isotopes as f64))
}

/// Search the union of member mass windows once per aligned group. Overlapping
/// windows are merged, so a peptide/offset appears once regardless of run count.
fn joint_matches(matrix: &FeatureMatrix, fi: usize, db: &PeptideDb, cfg: &MatchConfig) -> Vec<Pfm> {
    let group = &matrix.joint_evidence.as_ref().unwrap()[fi];
    let feature = &matrix.features[fi];
    let mut out = Vec::new();
    for &offset in &cfg.isotope_offsets {
        let mut windows: Vec<(f64, f64)> = group.iter().filter(|o| passes(&o.feature, cfg))
            .filter_map(|o| {
                let c = o.feature.mass - offset as f64 * NEUTRON;
                (c > 0.0).then(|| (c * (1.0 - cfg.ppm_tolerance * 1e-6), c * (1.0 + cfg.ppm_tolerance * 1e-6)))
            }).collect();
        windows.sort_by(|a, b| a.0.total_cmp(&b.0));
        let mut merged: Vec<(f64, f64)> = Vec::new();
        for (lo, hi) in windows {
            if let Some(last) = merged.last_mut().filter(|p| p.1 >= lo) { last.1 = last.1.max(hi); }
            else { merged.push((lo, hi)); }
        }
        let center = feature.mass - offset as f64 * NEUTRON;
        for (lo, hi) in merged {
            let start = db.peptides.partition_point(|p| p.mass < lo);
            for (j, pep) in db.peptides[start..].iter().enumerate() {
                if pep.mass > hi { break; }
                out.push(Pfm { feature_idx: fi as u32, peptide_idx: (start + j) as u32,
                    isotope_offset: offset, mass_error_ppm: (pep.mass - center) / center * 1e6, decoy: pep.decoy });
            }
        }
    }
    out
}

/// Generate all PFMs for the feature matrix against the database.
///
/// `db.peptides` must be mass-sorted (it is, after [`PeptideDb::build`]); we use
/// the original sorted index as `peptide_idx`, so downstream code can resolve a
/// PFM back to its peptide via `db.peptides[pfm.peptide_idx]`.
pub fn match_features(matrix: &FeatureMatrix, db: &PeptideDb, cfg: &MatchConfig) -> Vec<Pfm> {
    matrix
        .features
        .par_iter()
        .enumerate()
        .flat_map_iter(|(fi, feature)| {
            if matrix.joint_evidence.is_some() { return joint_matches(matrix, fi, db, cfg).into_iter(); }
            let mut out = Vec::new();
            if !passes(feature, cfg) {
                return out.into_iter();
            }
            for &offset in &cfg.isotope_offsets {
                let center = feature.mass - (offset as f64) * NEUTRON;
                if center <= 0.0 {
                    continue;
                }
                let tol = center * cfg.ppm_tolerance * 1e-6;
                let lo = center - tol;
                let hi = center + tol;
                let start = db.peptides.partition_point(|p| p.mass < lo);
                for (j, pep) in db.peptides[start..].iter().enumerate() {
                    if pep.mass > hi {
                        break;
                    }
                    out.push(Pfm {
                        feature_idx: fi as u32,
                        peptide_idx: (start + j) as u32,
                        isotope_offset: offset,
                        mass_error_ppm: (pep.mass - center) / center * 1e6,
                        decoy: pep.decoy,
                    });
                }
            }
            out.into_iter()
        })
        .collect()
}

/// Generate PFMs already reduced to the top-`k` mass matches per feature, split
/// by target/decoy, without ever materialising the full (unreduced) PFM vector.
///
/// This is a memory-streaming fusion of [`match_features`] followed by the
/// downstream `reduce_top_k`: per feature we collect the candidate PFMs (in the
/// same offset-then-peptide scan order as [`match_features`]), then stable-sort
/// each of the target/decoy buckets by `|mass_error_ppm|` and truncate to `k`.
/// The emitted order — for each feature, its target bucket then its decoy bucket
/// — is byte-for-byte the order `reduce_top_k` produces, so re-running
/// `reduce_top_k` on this set is idempotent and downstream output is unchanged.
/// The full ~41M-PFM intermediate (~1.3GB on run10) is never allocated.
pub fn match_features_topk(
    matrix: &FeatureMatrix,
    db: &PeptideDb,
    cfg: &MatchConfig,
    k: usize,
) -> Vec<Pfm> {
    matrix
        .features
        .par_iter()
        .enumerate()
        .flat_map_iter(|(fi, feature)| {
            let mut targets: Vec<Pfm> = Vec::new();
            let mut decoys: Vec<Pfm> = Vec::new();
            if matrix.joint_evidence.is_some() {
                for pfm in joint_matches(matrix, fi, db, cfg) {
                    if pfm.decoy { decoys.push(pfm); } else { targets.push(pfm); }
                }
            } else if passes(feature, cfg) {
                for &offset in &cfg.isotope_offsets {
                    let center = feature.mass - (offset as f64) * NEUTRON;
                    if center <= 0.0 {
                        continue;
                    }
                    let tol = center * cfg.ppm_tolerance * 1e-6;
                    let lo = center - tol;
                    let hi = center + tol;
                    let start = db.peptides.partition_point(|p| p.mass < lo);
                    for (j, pep) in db.peptides[start..].iter().enumerate() {
                        if pep.mass > hi {
                            break;
                        }
                        let pfm = Pfm {
                            feature_idx: fi as u32,
                            peptide_idx: (start + j) as u32,
                            isotope_offset: offset,
                            mass_error_ppm: (pep.mass - center) / center * 1e6,
                            decoy: pep.decoy,
                        };
                        if pep.decoy {
                            decoys.push(pfm);
                        } else {
                            targets.push(pfm);
                        }
                    }
                }
            }
            for bucket in [&mut targets, &mut decoys] {
                bucket.sort_by(|a, b| {
                    a.mass_error_ppm.abs().partial_cmp(&b.mass_error_ppm.abs()).unwrap()
                });
                bucket.truncate(k);
            }
            targets.into_iter().chain(decoys)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::database::DigestConfig;

    #[test]
    fn matches_real_data() {
        let (Some(features), Some(fasta)) = (
            crate::test_input("tests/data/09744.features.tsv"),
            crate::test_input("tests/data/entrap_human_yeast.fasta"),
        ) else {
            return;
        };

        let matrix = FeatureMatrix::load(&features).unwrap();
        let db = PeptideDb::build(&fasta, &DigestConfig::default()).unwrap();
        let cfg = MatchConfig::default();
        let pfms = match_features(&matrix, &db, &cfg);

        let n_target = pfms.iter().filter(|p| !p.decoy).count();
        let n_decoy = pfms.iter().filter(|p| p.decoy).count();
        let n_feat_matched = {
            let mut v: Vec<u32> = pfms.iter().map(|p| p.feature_idx).collect();
            v.sort_unstable();
            v.dedup();
            v.len()
        };
        eprintln!(
            "PFMs: {} ({} target, {} decoy); {}/{} features matched at least one peptide",
            pfms.len(),
            n_target,
            n_decoy,
            n_feat_matched,
            matrix.len()
        );
        assert!(n_target > 0 && n_decoy > 0);
        // All errors within tolerance. Necessary but weak: it holds for an empty
        // match set and for a set that matched the wrong peptides, so the checks
        // below assert the properties the search actually depends on.
        assert!(pfms.iter().all(|p| p.mass_error_ppm.abs() <= cfg.ppm_tolerance + 1e-6));

        // Targets and decoys must arrive in comparable numbers. The decoy database
        // is the null the whole FDR rests on, so a lopsided ratio means the null is
        // not exchangeable with the targets, whatever the mass errors look like.
        let ratio = n_decoy as f64 / n_target as f64;
        assert!(
            (0.5..=2.0).contains(&ratio),
            "target/decoy PFM ratio is lopsided: {n_target} target vs {n_decoy} decoy"
        );

        // Every PFM must point at a real feature and a real peptide, and its charge
        // must be one the matcher was allowed to consider.
        assert!(pfms.iter().all(|p| (p.feature_idx as usize) < matrix.len()));
        assert!(pfms.iter().all(|p| (p.peptide_idx as usize) < db.peptides.len()));

        // The isotope offset is applied on the NEUTRAL mass in carbon-13 units, so
        // only the offsets the config lists may appear.
        assert!(pfms.iter().all(|p| cfg.isotope_offsets.contains(&p.isotope_offset)));

        // Charges below the floor are skipped as non-peptide contaminants.
        assert!(pfms
            .iter()
            .all(|p| passes(&matrix.features[p.feature_idx as usize], &cfg)));

        // A real search matches a substantial fraction of features, not a handful.
        assert!(
            n_feat_matched * 4 > matrix.len(),
            "only {n_feat_matched}/{} features matched a peptide",
            matrix.len()
        );
    }
}
