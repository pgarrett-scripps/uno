//! Phase 3: match aligned features to theoretical peptides by precursor mass.
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
}

impl Default for MatchConfig {
    fn default() -> Self {
        MatchConfig {
            ppm_tolerance: 10.0,
            isotope_offsets: vec![-1, 0, 1],
            min_charge: 1,
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
            let mut out = Vec::new();
            if feature.charge < cfg.min_charge {
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
            if feature.charge >= cfg.min_charge {
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
    use std::path::PathBuf;

    #[test]
    fn matches_real_data() {
        let dir =
            PathBuf::from("/home/patrick-garrett/Repos/koth_rust/benchmark/data/koth/align_output");
        let fasta = PathBuf::from(
            "/home/patrick-garrett/Repos/koth_rust/benchmark/data/fasta/human_ecoli.fasta",
        );
        if !dir.join("intensity_matrix.tsv").exists() || !fasta.exists() {
            eprintln!("skipping: benchmark inputs not present");
            return;
        }

        let matrix = FeatureMatrix::load(
            dir.join("intensity_matrix.tsv"),
            Some(dir.join("qvalue_matrix.tsv")),
        )
        .unwrap();
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
        // All errors within tolerance
        assert!(pfms.iter().all(|p| p.mass_error_ppm.abs() <= cfg.ppm_tolerance + 1e-6));
    }
}
