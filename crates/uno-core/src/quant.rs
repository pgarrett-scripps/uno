//! Multi-sample quantitation over a koth-aligned matrix: correct the
//! peptide->protein feature assignment using the across-sample intensity pattern,
//! then emit a protein-group × sample abundance matrix.
//!
//! uno is an ID + quant tool, not a stats framework. It stops at the abundance
//! matrix; normalization, differential-expression testing, and interpretation are
//! left to downstream tools (MSstats/limma/etc.), which take this matrix plus the
//! experimental design. So there is deliberately NO design file, contrast, or FDR
//! here — only the pieces that make the abundance values *correct*:
//!
//! - [`resolve_feature_competition`] — one feature quantifies for one protein
//!   (parsimony), so a minority peptide can't measure a mass-coincident feature
//!   owned by an abundant background.
//! - [`profile_reassign`] — re-route a feature to the protein whose across-sample
//!   profile it tracks, fixing mass+RT ties the ID stage couldn't (design-free:
//!   uses variation, not condition labels).
//! - [`group_abundance`] — sum each protein group's assigned peptides per run.

use crate::database::PeptideDb;
use crate::input::FeatureMatrix;
use crate::score::PeptideHit;
use std::collections::{HashMap, HashSet};

/// DirectMS1Quant two-condition statistics (reference port) and the uno
/// per-run evidence adapter. Exposed as `uno quant`.
pub mod directms1;

/// Feature-level competition for quant: each consensus feature may quantify for
/// only ONE protein — the best-scoring one whose peptide claims it. uno's ID stage
/// already does this (parsimony: one peak → one protein), but the per-peptide
/// quant feature-pick discards it and lets any peptide grab any mass+RT match. On
/// a mixed proteome that means a minority-species peptide (e.g. E. coli) quantifies
/// on a mass-coincident feature owned by an abundant background protein, reading
/// its flat profile instead of the real signal. Replaying the claim here — a
/// feature goes to the peptide whose protein scores highest — drops those stolen
/// features, so the minority peptide falls back to its real feature or drops out.
/// `prot_score` maps (decoy, accession) → the protein's binomial score (higher =
/// more confident).
pub fn resolve_feature_competition(
    peptides: &[PeptideHit],
    db: &PeptideDb,
    prot_score: &std::collections::HashMap<(bool, &str), f64>,
) -> Vec<PeptideHit> {
    // Best backing-protein score per hit (shared peptides take their strongest
    // protein). Non-scored proteins → -inf, so they lose every contested feature.
    let scores: Vec<f64> = peptides
        .iter()
        .map(|h| {
            db.peptides[h.peptide_idx as usize]
                .proteins
                .iter()
                .filter_map(|p| prot_score.get(&(h.decoy, p.as_ref())).copied())
                .fold(f64::NEG_INFINITY, f64::max)
        })
        .collect();
    let features: Vec<u32> = peptides.iter().map(|h| h.feature_idx).collect();
    feature_winners(&features, &scores)
        .into_iter()
        .map(|i| peptides[i])
        .collect()
}

/// The one-feature-one-protein claim behind [`resolve_feature_competition`],
/// on plain `(feature, score)` columns so per-run quantitative evidence can
/// replay exactly the same rule. Returns the winning row indices in ascending
/// (input) order. A feature goes to its highest-scoring row; ties keep the
/// lower index; rows with a non-finite score never win.
pub fn feature_winners(features: &[u32], scores: &[f64]) -> Vec<usize> {
    assert_eq!(features.len(), scores.len());
    let mut winner: std::collections::HashMap<u32, usize> = std::collections::HashMap::new();
    for (i, &feature) in features.iter().enumerate() {
        if !scores[i].is_finite() {
            continue;
        }
        match winner.get(&feature) {
            Some(&j) if scores[j] >= scores[i] => {}
            _ => {
                winner.insert(feature, i);
            }
        }
    }
    // Return in the input's own order, not the HashMap's. The winning SET is
    // deterministic (ties keep the lower index, since the incumbent is always
    // seen first), but `winner.values()` is process-random, and this Vec becomes
    // `--out-lfq` row order and the summation order inside `group_abundance` —
    // so without this the LFQ rows shuffle and the abundances differ in the low
    // bits between two runs on identical input.
    let mut idx: Vec<usize> = winner.values().copied().collect();
    idx.sort_unstable();
    idx
}

/// Single-pass profile-guided reassignment (combines ID and quant).
///
/// The feature competition settles mass ties by *ID confidence*; this settles the
/// ones it can't by the *quant pattern*. Starting from `base` (the competition
/// result), it builds each protein's across-sample profile, then for every
/// candidate `(feature, peptide)` pair whose feature actually VARIES, re-routes the
/// feature to the mass-legitimate protein whose profile it best CORRELATES with.
/// A rising E. coli feature that mass-matching first handed to a flat background
/// protein moves to the E. coli protein it tracks — recovering signal instead of
/// discarding it.
///
/// Guardrails against circularity: (1) only mass-candidate proteins are considered
/// (`candidates` comes from the search's PFMs — a feature can't join a protein with
/// no peptide at its mass); (2) only VARYING features/proteins are re-routed
/// (`min_var`) — a flat feature carries no profile information, so it keeps the
/// competition's assignment; (3) decoys go through the identical pass. `candidates`
/// = `(feature_idx, peptide_idx, decoy)` for confident-protein peptides.
pub fn profile_reassign(
    matrix: &FeatureMatrix,
    db: &PeptideDb,
    base: &[PeptideHit],
    candidates: &[(u32, u32, bool)],
    min_corr: f64,
    min_var: f64,
) -> Vec<PeptideHit> {
    let n_runs = matrix.run_names.len();
    let (run_med, grand) = run_normalizers(matrix, base, n_runs);
    // per-run low floor for imputing missing cells (5th pct of detected)
    let mut per_run: Vec<Vec<f64>> = vec![Vec::new(); n_runs];
    for h in base {
        for (r, &v) in matrix.features[h.feature_idx as usize].intensities.iter().enumerate().take(n_runs) {
            if v > 0.0 { per_run[r].push(v.log2()); }
        }
    }
    let run_floor: Vec<f64> = per_run.iter_mut().map(|v| {
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        if v.is_empty() { 0.0 } else { v[v.len() / 20] }
    }).collect();
    let profile = |fi: usize| -> Vec<f64> {
        let ints = &matrix.features[fi].intensities;
        (0..n_runs).map(|r| {
            let v = ints.get(r).copied().unwrap_or(0.0);
            let l = if v > 0.0 { v.log2() } else { run_floor[r] };
            l - run_med[r] + grand
        }).collect()
    };
    let variance = |p: &[f64]| -> f64 {
        let m = p.iter().sum::<f64>() / p.len() as f64;
        p.iter().map(|x| (x - m).powi(2)).sum::<f64>() / p.len() as f64
    };

    // Protein profiles (per-run median over the protein's base features).
    let mut prot_feats: HashMap<(bool, &str), Vec<usize>> = HashMap::new();
    for h in base {
        for prot in &db.peptides[h.peptide_idx as usize].proteins {
            prot_feats.entry((h.decoy, prot.as_ref())).or_default().push(h.feature_idx as usize);
        }
    }
    let prot_profile: HashMap<(bool, &str), Vec<f64>> = prot_feats
        .iter()
        .map(|(k, feats)| {
            let mut cols = vec![Vec::new(); n_runs];
            for &fi in feats {
                let p = profile(fi);
                for r in 0..n_runs { cols[r].push(p[r]); }
            }
            let prof: Vec<f64> = cols.iter_mut().map(|c| { c.sort_by(|a,b| a.partial_cmp(b).unwrap()); c[c.len()/2] }).collect();
            (*k, prof)
        })
        .collect();

    // Start from the competition assignment (one peptide/protein per feature).
    let mut by_feat: HashMap<u32, PeptideHit> = base.iter().map(|h| (h.feature_idx, *h)).collect();

    // For each VARYING candidate feature, find the best-correlating candidate protein.
    let mut best: HashMap<u32, (f64, u32, bool)> = HashMap::new(); // feat -> (corr, peptide, decoy)
    let mut fp_cache: HashMap<u32, (Vec<f64>, bool)> = HashMap::new(); // feat -> (profile, varies)
    for &(fi, pi, decoy) in candidates {
        let (fp, varies) = fp_cache.entry(fi).or_insert_with(|| {
            let p = profile(fi as usize);
            let v = variance(&p) >= min_var;
            (p, v)
        });
        if !*varies { continue; }
        let fp = fp.clone();
        for prot in &db.peptides[pi as usize].proteins {
            if let Some(pp) = prot_profile.get(&(decoy, prot.as_ref())) {
                if variance(pp) < min_var { continue; }
                let c = pearson(&fp, pp);
                let e = best.entry(fi).or_insert((f64::NEG_INFINITY, pi, decoy));
                if c > e.0 { *e = (c, pi, decoy); }
            }
        }
    }
    // Apply reassignments that clear the correlation bar.
    let mut moved = 0usize;
    for (fi, (c, pi, decoy)) in best {
        if c >= min_corr {
            let prev = by_feat.get(&fi).map(|h| h.peptide_idx);
            if prev != Some(pi) { moved += 1; }
            by_feat.insert(fi, PeptideHit { peptide_idx: pi, feature_idx: fi, decoy });
        }
    }
    log::info!("profile reassignment: {} features re-routed by quant profile", moved);
    // Same reason as in `resolve_feature_competition`: emit in feature order, not
    // HashMap order. This is the opt-in `--profile-reassign` path, which WAS used
    // for the reported quantitation, so it is on the live path to `--out-lfq`.
    let mut hits: Vec<PeptideHit> = by_feat.into_values().collect();
    hits.sort_unstable_by_key(|h| (h.feature_idx, h.peptide_idx, h.decoy));
    hits
}

/// Per-run log2 normalization offsets from the INVARIANT backbone: features
/// detected in every run (complete cases). These are dominated by the unchanging
/// background, so normalizing to them removes the systematic drift that a plain
/// all-feature median suffers when a large fraction of the sample changes
/// (e.g. a heavy spike-in): there, the all-feature median tracks the changing
/// species and makes the constant background appear to drift. Returns
/// `(offset_per_run, grand_mean)`; apply as `log2(v) - offset[r] + grand`.
/// Falls back to all detected cells if too few complete cases.
fn run_normalizers(matrix: &FeatureMatrix, peptides: &[PeptideHit], n_runs: usize) -> (Vec<f64>, f64) {
    let mut complete: Vec<Vec<f64>> = vec![Vec::new(); n_runs];
    let mut all: Vec<Vec<f64>> = vec![Vec::new(); n_runs];
    for h in peptides {
        let ints = &matrix.features[h.feature_idx as usize].intensities;
        let full = ints.len() >= n_runs && ints[..n_runs].iter().all(|&v| v > 0.0);
        for (r, &v) in ints.iter().enumerate().take(n_runs) {
            if v > 0.0 {
                all[r].push(v.log2());
                if full {
                    complete[r].push(v.log2());
                }
            }
        }
    }
    let med = |v: &mut Vec<f64>| -> f64 {
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        if v.is_empty() { 0.0 } else { v[v.len() / 2] }
    };
    let use_complete = complete.iter().map(|v| v.len()).min().unwrap_or(0) >= 200;
    let src = if use_complete { &mut complete } else { &mut all };
    let offs: Vec<f64> = src.iter_mut().map(|v| med(v)).collect();
    let grand = offs.iter().sum::<f64>() / n_runs.max(1) as f64;
    (offs, grand)
}


fn pearson(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    if n < 2.0 {
        return 0.0;
    }
    let mx = x.iter().sum::<f64>() / n;
    let my = y.iter().sum::<f64>() / n;
    let mut sxy = 0.0;
    let mut sxx = 0.0;
    let mut syy = 0.0;
    for i in 0..x.len() {
        let dx = x[i] - mx;
        let dy = y[i] - my;
        sxy += dx * dy;
        sxx += dx * dx;
        syy += dy * dy;
    }
    if sxx <= 0.0 || syy <= 0.0 {
        0.0
    } else {
        sxy / (sxx * syy).sqrt()
    }
}


/// Protein-group abundance matrix — uno's quant deliverable (design-free).
///
/// For each confident TARGET group, sums the intensity of its assigned
/// (group-unique + razor) peptide-features per run. This is a raw protein × sample
/// intensity table: uno stops here (ID + quant); normalization, differential
/// testing, and interpretation are downstream tools' job. Returns
/// `leading accession -> (per-run intensities, distinct peptide count)`.
pub fn group_abundance(
    matrix: &FeatureMatrix,
    peptides: &[PeptideHit],
    pep_group: &HashMap<u32, (bool, String)>,
    n_runs: usize,
) -> HashMap<String, (Vec<f64>, HashSet<u32>)> {
    let mut out: HashMap<String, (Vec<f64>, HashSet<u32>)> = HashMap::new();
    for h in peptides {
        if let Some((decoy, leading)) = pep_group.get(&h.peptide_idx) {
            if *decoy {
                continue; // abundance is reported for target groups only
            }
            let ints = &matrix.features[h.feature_idx as usize].intensities;
            let e = out
                .entry(leading.clone())
                .or_insert_with(|| (vec![0.0; n_runs], HashSet::new()));
            for (r, slot) in e.0.iter_mut().enumerate() {
                *slot += ints.get(r).copied().unwrap_or(0.0);
            }
            e.1.insert(h.peptide_idx);
        }
    }
    out
}

#[cfg(test)]
mod determinism_tests {
    use super::*;
    use crate::database::{PeptideDb, TheoPeptide};
    use std::sync::Arc;

    fn db_of(specs: &[(&str, bool, &[&str])]) -> PeptideDb {
        let peptides: Vec<TheoPeptide> = specs
            .iter()
            .map(|(seq, decoy, prots)| TheoPeptide {
                sequence: Arc::from(*seq),
                mass: 1000.0,
                decoy: *decoy,
                missed_cleavages: 0,
                proteins: prots.iter().map(|p| Arc::from(*p)).collect(),
            })
            .collect();
        let n_targets = peptides.iter().filter(|p| !p.decoy).count();
        let n_decoys = peptides.len() - n_targets;
        PeptideDb { peptides, n_targets, n_decoys }
    }

    /// `resolve_feature_competition` must return a stable ORDER, not just a stable set.
    ///
    /// It resolves winners through a `HashMap<feature_idx, usize>` and used to emit
    /// `winner.values()` directly. Two `HashMap`s in one process get different seeds,
    /// so this catches the bug without spawning a second process. The returned Vec is
    /// `--out-lfq` row order and the summation order inside `group_abundance`, so a
    /// shuffle there means the LFQ rows move and the abundances differ in the low bits.
    #[test]
    fn feature_competition_output_order_is_stable() {
        let db = db_of(&[
            ("PEPTIDEA", false, &["P1"][..]),
            ("PEPTIDEB", false, &["P2"][..]),
            ("PEPTIDEC", false, &["P1"][..]),
            ("PEPTIDED", false, &["P3"][..]),
            ("PEPTIDEE", false, &["P2"][..]),
            ("PEPTIDEF", false, &["P3"][..]),
        ]);
        let mut prot_score: std::collections::HashMap<(bool, &str), f64> =
            std::collections::HashMap::new();
        prot_score.insert((false, "P1"), 9.0);
        prot_score.insert((false, "P2"), 5.0);
        prot_score.insert((false, "P3"), 1.0);

        // Several peptides contend for the same features, so the map is exercised.
        let hits: Vec<PeptideHit> = (0..6u32)
            .map(|i| PeptideHit { peptide_idx: i, feature_idx: i % 3, decoy: false })
            .collect();

        let first = resolve_feature_competition(&hits, &db, &prot_score);
        for _ in 0..12 {
            let again = resolve_feature_competition(&hits, &db, &prot_score);
            assert_eq!(first, again, "feature competition output order is not stable");
        }
        // And the output really is in feature order.
        assert_eq!(first.len(), 3);
        assert!(first.windows(2).all(|w| w[0].feature_idx < w[1].feature_idx));
    }
}
