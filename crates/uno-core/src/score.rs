//! Phase 5: discriminant scoring, bootstrap RT calibration, and target-decoy FDR.
//!
//! Pipeline:
//!   1. Reduce the (huge) PFM set to the top-K best-mass-error target and decoy
//!      matches per feature. MS1 mass-matching yields ~hundreds of candidates
//!      per feature; all but a few are hopeless, and keeping every one makes the
//!      discriminant intractable. Top-K preserves the RT-correct peptide while
//!      cutting the set ~100x.
//!   2. Bootstrap RT calibration: score once without RT, accept a provisional
//!      target set at low FDR, fit abstract->observed RT calibration on it.
//!   3. Re-score with the `rt_diff` feature active, recalibrating each round.
//!   4. A 3-fold CV gbdt discriminant (split BY FEATURE so a feature's PFMs are
//!      never used to train the model that scores them) avoids overfitting to
//!      decoys. Target-decoy competition per feature -> q-values.

use crate::calibrate::Calibration;
use crate::database::PeptideDb;
use crate::input::{Feature, FeatureMatrix};
use crate::matching::Pfm;
use crate::protein::ProteinId;
use crate::rt::RtModel;
use rayon::prelude::*;
use std::cmp::Ordering;

/// Order-invariant "which feature best represents a peptide" comparator, BEST
/// first. When one peptide matches several features (or koth's exhaustive mode
/// emits near-duplicates that tie on mass error / deviation), the representative
/// must be chosen by feature *quality*, not by whichever row happened to be seen
/// first — otherwise the result depends on the input's feature row order. Higher
/// `combined_score`, then higher apex intensity, then more isotopes; ties broken
/// by ascending mass/rt/mz so the order is total and reproducible.
fn feature_cmp(a: &Feature, b: &Feature) -> Ordering {
    // Total intensity (sum over runs); for single-run input this is just the one
    // value, for an aligned matrix it's the cross-run total — a stable summary in
    // both cases (unlike `.first()`, which is run 0 and may be an empty cell).
    let ai: f64 = a.intensities.iter().sum();
    let bi: f64 = b.intensities.iter().sum();
    b.combined_score
        .partial_cmp(&a.combined_score)
        .unwrap()
        .then(bi.partial_cmp(&ai).unwrap())
        .then(b.n_isotopes.unwrap_or(0.0).partial_cmp(&a.n_isotopes.unwrap_or(0.0)).unwrap())
        .then(a.mass.partial_cmp(&b.mass).unwrap())
        .then(a.rt.partial_cmp(&b.rt).unwrap())
        .then(a.mz.partial_cmp(&b.mz).unwrap())
}

#[derive(Debug, Clone, serde::Deserialize)]
#[serde(default)]
pub struct ScoreConfig {
    /// Keep this many best-mass-error targets and decoys per feature. This caps
    /// the candidate set the discriminant sees, so setting it too low silently
    /// drops real peptides before RT can rescue them. Entrapment (human+yeast on
    /// single-run HeLa) showed 25 identified only ~half of what the search
    /// honestly supports: 25→105 IDs, 50→208, 100→275, all at ~1% entrapment
    /// FDP. 50 is the default (≈2× IDs over 25 at negligible cost); raise toward
    /// 100 for maximum depth at a runtime cost.
    pub top_k: usize,
    pub n_folds: usize,
    /// FDR used to select the provisional set for RT calibration.
    pub calibration_fdr: f64,
    /// Number of re-score rounds after the first calibration.
    pub bootstrap_rounds: usize,
    pub iterations: usize,
    pub max_depth: u32,
    pub shrinkage: f32,
    /// Cap on training rows per fold-model (0 = no cap). Bounds training time.
    pub max_train_rows: usize,
    /// Seed features must appear in at least this fraction of samples. Lower for
    /// many-sample designs (e.g. dose-response) where a peptide need not be in
    /// every run; the RANSAC seed tolerates the extra contamination.
    pub seed_min_runs_frac: f64,
    /// Fallback seed for scoring's random components. 0 = the shipped default
    /// streams (see [`crate::seed_mix`]). `uno --seed` sets this value; explicit
    /// `rt_seed`, `cv_seed`, or `classifier_seed` values take precedence.
    pub seed: u64,
    /// Override only RANSAC RT-calibration pair sampling. None uses `seed`.
    pub rt_seed: Option<u64>,
    /// Override only cross-validation fold assignment. None uses `seed`.
    pub cv_seed: Option<u64>,
    /// Override classifier row/column sampling and LightGBM hyperparameter
    /// search. Applies to both Sequoia and LightGBM. None uses `seed`.
    pub classifier_seed: Option<u64>,
    /// Mass tolerance (ppm) the mass-anchor bootstrap starts at, and the ceiling it
    /// tightens toward as it recalibrates the per-run mass shift + sigma. Default 10
    /// (Orbitrap). For lower-accuracy instruments (e.g. timsTOF) set this wider —
    /// e.g. 20 — plus `[matching] ppm_tolerance = 20`; the anchor then seeds wide and
    /// tightens to single-digit ppm once the mass error is recalibrated on-data.
    pub mass_anchor_ppm: f64,
    /// Discriminant: "deviation" (robust orthogonal RT+IM consistency, ms1searchpy
    /// style — default) or "gbdt" (semi-supervised boosted trees). On the
    /// MS1-only/saturated-mass regime the gbdt washes out the weak RT+IM signal,
    /// so the deviation score is the default.
    pub discriminant: String,
    /// Refit the RT predictor on *this run's* observed RTs (ridge on composition
    /// features) using the confident set, each iteration. The DirectMS1 trick:
    /// the frozen DeepLC model only rank-correlates ~0.3 on an unseen gradient,
    /// but a model fit to the data makes RT a strong discriminator. Default on.
    pub ondata_rt: bool,
    /// Ridge penalty for the on-data RT fit.
    pub ondata_rt_lambda: f64,
    /// Fine-tune mode for the on-data RT fit (requires `ondata_rt` and an MLP RT
    /// model). Instead of fitting a fresh ridge on raw composition features (a
    /// weak linear model that discards the pretrained representation), refit the
    /// readout over the MLP's frozen penultimate embedding — keeping its learned
    /// nonlinear RT structure while adapting offset/scale/shape to this run.
    pub ondata_rt_finetune: bool,
    /// Bootstrap the RT calibration from MASS-ONLY protein identifications
    /// (ms1searchpy stage-1) instead of the per-peptide RANSAC seed. On dense
    /// single-run data the per-peptide seed is ~99% contaminated; confident
    /// proteins (mass multiplicity) give far cleaner RT anchors.
    pub mass_anchor: bool,
    /// The RT library is already calibrated to this run's RT (e.g. DeepLC
    /// `calibrate_preds`). uno then uses rt_diff = obs - pred directly and does
    /// NOT apply its own calibration on top (which would double-calibrate and
    /// distort the predictions).
    pub rt_calibrated: bool,
    /// CV-gbdt PFM peptide-quality ranking, fed through a 10-step cumulative
    /// qpred-bin ensemble protein score (ms1searchpy `-ml 1` + `final_iteration`).
    /// VALIDATED (2026-07-07, PXD015669 human+yeast entrapment, native biosaur2
    /// features with isoerror/isoerror2/rt_width): 223 real/4 yeast (deviation
    /// alone) -> 251 real/1 yeast with this on — deeper AND more honest. The
    /// per-PFM classifier itself still barely separates decoys from dim-real
    /// peptides on uno's pure-Rust gbdt (dim_mean ~= DECOY_mean); the gain comes
    /// from the ENSEMBLE STRUCTURE (10 cumulative thresholds, averaged, like
    /// ms1searchpy's `final_iteration`), not classifier quality. When a feature
    /// file lacks the quality columns it correctly no-ops (guarded by a strict
    /// better-than-baseline check).
    ///
    /// ON by default (see `ScoreConfig::default`), and active in every benchmark
    /// reported in the manuscript. It costs a re-run of the protein engine per
    /// inclusion band, so it trades latency for depth. Set false for the plain
    /// deviation path.
    pub ml_rescore: bool,
}

impl Default for ScoreConfig {
    fn default() -> Self {
        ScoreConfig {
            top_k: 50,
            n_folds: 3,
            calibration_fdr: 0.01,
            bootstrap_rounds: 2,
            iterations: 100,
            max_depth: 6,
            shrinkage: 0.1,
            max_train_rows: 200_000,
            seed_min_runs_frac: 0.5,
            seed: 0,
            rt_seed: None,
            cv_seed: None,
            classifier_seed: None,
            mass_anchor_ppm: 10.0,
            discriminant: "deviation".to_string(),
            // On by default IN FINE-TUNE MODE: refit only the MLP readout over its
            // frozen embedding on this run's confident hits. This keeps the
            // pretrained representation (unlike the from-scratch composition ridge,
            // which degrades depth) and is FDR-honest under two-proteome entrapment
            // (0 entrapment false positives, same as the frozen path). Falls back
            // to frozen RT automatically when the model has no embedding (gbdt).
            ondata_rt: true,
            ondata_rt_lambda: 1.0,
            ondata_rt_finetune: true,
            mass_anchor: true,
            rt_calibrated: false,
            // ms1searchpy's `-ml 1` CV-boosted PFM rescorer + 10-step protein
            // ensemble. ON by default: run10 honest yeast-entrapment 412->610 @1%
            // (deterministic; +90% to 781 once the koth cosine_score/isotope_score
            // reach the ML — see run10 notes). Auto-no-regression guarded
            // (`ml_n > best_n` before it swaps in), so it can only deepen the
            // deviation result. Reproducible only with the LightGBM seed pinning +
            // `items` sort below.
            ml_rescore: true,
        }
    }
}

impl ScoreConfig {
    pub fn effective_rt_seed(&self) -> u64 {
        self.rt_seed.unwrap_or(self.seed)
    }

    pub fn effective_cv_seed(&self) -> u64 {
        self.cv_seed.unwrap_or(self.seed)
    }

    pub fn effective_classifier_seed(&self) -> u64 {
        self.classifier_seed.unwrap_or(self.seed)
    }
}

/// Result of scoring: the confident protein set from the mass-anchor mass+RT
/// binomial identifier (ms1searchpy-style). This IS uno's deliverable.
#[derive(Debug, Clone)]
pub struct ScoreResult {
    pub proteins: Vec<crate::protein::ProteinId>,
    /// RT+mass-gated best feature per peptide (target and decoy), for the
    /// multi-sample quant rollup. Empty when the ML rescore stage did not run.
    pub peptides: Vec<PeptideHit>,
}

/// Per-PFM feature vector for the discriminant. `rt_diff`/`im_diff` are `None`
/// before their calibration is available, in which case those columns carry no
/// signal (ion mobility is absent entirely on non-IM data).
fn featurize(
    pfm: &Pfm,
    matrix: &FeatureMatrix,
    db: &PeptideDb,
    rt_diff: Option<f64>,
    im_diff: Option<f64>,
) -> Vec<f32> {
    let f = &matrix.features[pfm.feature_idx as usize];
    let p = &db.peptides[pfm.peptide_idx as usize];
    let n_detected = f.intensities.iter().filter(|&&x| x > 0.0).count();
    let sum_int: f64 = f.intensities.iter().sum();
    vec![
        pfm.mass_error_ppm.abs() as f32,
        pfm.mass_error_ppm as f32,
        pfm.isotope_offset.unsigned_abs() as f32,
        rt_diff.unwrap_or(0.0).abs() as f32,
        f.combined_score as f32,
        ((1.0 + f.n_contributing_runs as f64).ln()) as f32,
        f.charge as f32,
        ((1.0 + sum_int).log10()) as f32,
        p.sequence.len() as f32,
        p.missed_cleavages as f32,
        n_detected as f32,
        im_diff.unwrap_or(0.0).abs() as f32,
        // Per-feature shape (biosaur2 single-run inputs): a well-formed feature
        // (more isotopic peaks, more scans) is more trustworthy than a junk
        // feature, which contributes random target+decoy matches equally. 0 when
        // the feature file lacks the column, so this is inert there.
        (1.0 + f.n_isotopes.unwrap_or(0.0)).ln() as f32,
        (1.0 + f.n_scans.unwrap_or(0.0)).ln() as f32,
        // Single-run PFM-quality columns (present on biosaur2 feature files).
        // These are the exact signals ms1searchpy's LightGBM PFM filter uses to
        // drop the random-match rate p to ~0.135 (its @5% depth). 0.0 when absent
        // (aligned matrices), so inert there — active only on single-run inputs.
        f.cosine_score.unwrap_or(0.0) as f32,
        f.isotope_score.unwrap_or(0.0) as f32,
        f.ppm_error.unwrap_or(0.0) as f32,
        f.neutron_offset.unwrap_or(0.0) as f32,
    ]
}

/// Stable fold assignment for a feature index.
fn fold_of(feature_idx: usize, n_folds: usize, seed: u64) -> usize {
    let mut x = (feature_idx as u64) ^ crate::seed_mix(seed);
    x ^= x >> 30;
    x = x.wrapping_mul(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    (x % n_folds as u64) as usize
}

fn median(v: &mut [f64]) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    if n == 0 {
        return f64::NAN;
    }
    if n % 2 == 1 {
        v[n / 2]
    } else {
        0.5 * (v[n / 2 - 1] + v[n / 2])
    }
}

/// Per-peptide posterior P(correct | dev) from the target vs decoy `dev`
/// distributions, for the continuous (sweep-free) protein scorer. Decoys are the
/// null; at a given `dev` the fraction of targets in excess of decoys is the correct
/// fraction: w = clamp(1 - scale * n_decoy / n_target, 0, 1), estimated in
/// equal-count `dev` bins, then made non-increasing in `dev` by pool-adjacent-
/// violators (larger deviation can only be *less* likely correct). `scale` corrects
/// for an unequal target:decoy DB size (~1.0 for uno's 1:1 decoys).
fn dev_posteriors(
    pep_dev: &std::collections::HashMap<usize, (bool, f64, usize)>,
    scale: f64,
) -> rustc_hash::FxHashMap<usize, f64> {
    let mut items: Vec<(f64, bool, usize)> =
        pep_dev.iter().map(|(&pi, &(dc, dev, _))| (dev, dc, pi)).collect();
    items.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let n = items.len();
    if n == 0 {
        return rustc_hash::FxHashMap::default();
    }
    let nbins = 256.min(n);
    let mut bin_w = vec![0.0f64; nbins];
    let mut bin_start = vec![0usize; nbins];
    for b in 0..nbins {
        let lo = b * n / nbins;
        let hi = (((b + 1) * n / nbins).max(lo + 1)).min(n);
        bin_start[b] = lo;
        let (mut tgt, mut dec) = (0.0f64, 0.0f64);
        for it in &items[lo..hi] {
            if it.1 { dec += 1.0; } else { tgt += 1.0; }
        }
        bin_w[b] = if tgt > 0.0 { (1.0 - scale * dec / tgt).clamp(0.0, 1.0) } else { 0.0 };
    }
    // Isotonic (non-increasing) via pool-adjacent-violators: merge a bin into its
    // predecessor whenever the predecessor is SMALLER (an increase = a violation).
    let mut vals: Vec<f64> = Vec::new();
    let mut cnt: Vec<usize> = Vec::new();
    for &w in &bin_w {
        vals.push(w);
        cnt.push(1);
        while vals.len() > 1 && vals[vals.len() - 2] < vals[vals.len() - 1] {
            let (v2, c2) = (vals.pop().unwrap(), cnt.pop().unwrap());
            let (v1, c1) = (vals.pop().unwrap(), cnt.pop().unwrap());
            vals.push((v1 * c1 as f64 + v2 * c2 as f64) / (c1 + c2) as f64);
            cnt.push(c1 + c2);
        }
    }
    let mut mono: Vec<f64> = Vec::with_capacity(nbins);
    for (v, c) in vals.iter().zip(&cnt) {
        for _ in 0..*c {
            mono.push(*v);
        }
    }
    let mut out: rustc_hash::FxHashMap<usize, f64> = rustc_hash::FxHashMap::default();
    for b in 0..nbins {
        let lo = bin_start[b];
        let hi = if b + 1 < nbins { bin_start[b + 1] } else { n };
        let w = mono.get(b).copied().unwrap_or(0.0);
        for it in &items[lo..hi] {
            out.insert(it.2, w);
        }
    }
    out
}

/// Final scoring output and the peptide evidence from the same calibration
/// candidate. Keep these owned together until candidate selection is complete.
struct FinalizedCandidate {
    proteins: Vec<ProteinId>,
    n_confident: usize,
    peptides: Vec<(u32, u32, bool)>,
}

fn soft_depth(proteins: &[ProteinId]) -> f64 {
    proteins.iter().filter(|p| !p.decoy).map(|p| (1.0 - p.q_value / 0.05).clamp(0.0, 1.0)).sum()
}

/// Preserve scored evidence even when no cut passes the reporting threshold.
/// Only all-zero 1% counts use the existing 25%-FDR ML training set as a tie
/// breaker. A nonzero count always wins; equal nonzero counts keep the first cut.
fn select_deviation_cut<'a>(cuts: impl IntoIterator<Item = &'a [ProteinId]>, retain_zero: bool) -> Option<usize> {
    let mut best: Option<(usize, (usize, usize))> = None;
    for (i, proteins) in cuts.into_iter().enumerate() {
        let n = proteins.iter().filter(|p| !p.decoy && p.q_value <= 0.01).count();
        if n == 0 && !retain_zero { continue; }
        let training = if n == 0 {
            proteins.iter().filter(|p| !p.decoy && p.q_value <= 0.25).count()
        } else { 0 };
        let rank = (n, training);
        if best.as_ref().map_or(true, |(_, previous)| rank > *previous) {
            best = Some((i, rank));
        }
    }
    best.map(|(i, _)| i)
}

fn select_finalized_candidate(candidates: Vec<FinalizedCandidate>, count_only: bool, protect_discoveries: bool) -> (usize, FinalizedCandidate) {
    let mut best: Option<(usize, FinalizedCandidate)> = None;
    let mut best_rank = (false, f64::NEG_INFINITY);
    for (i, candidate) in candidates.into_iter().enumerate() {
        let score = if count_only { candidate.n_confident as f64 } else { soft_depth(&candidate.proteins) };
        // Evidence retained for training must not displace actual discoveries.
        // Before zero-discovery cuts retained their lists, their soft depth was
        // zero, so a positive-discovery candidate already always outranked them.
        let rank = (protect_discoveries && candidate.n_confident > 0, score);
        // Strict comparison deliberately preserves the first candidate on ties.
        if rank > best_rank {
            best_rank = rank;
            best = Some((i, candidate));
        }
    }
    best.expect("at least one calibration candidate")
}

/// Recovery is a fallback after the complete original search, never an extra
/// competitor against a successful ranking. Keep the original evidence if the
/// retry also fails to produce any reported target discoveries.
fn recover_zero_discoveries(original: ScoreResult, retry: impl FnOnce() -> ScoreResult) -> ScoreResult {
    let has_discoveries = |r: &ScoreResult| r.proteins.iter().any(|p| !p.decoy && p.q_value <= 0.01);
    if has_discoveries(&original) { return original; }
    let recovered = retry();
    if has_discoveries(&recovered) {
        log::info!("mass-anchor: evidence-retention recovery produced reported 1% targets");
        recovered
    } else {
        log::info!("mass-anchor: recovery has no reported 1% targets; preserving original result");
        original
    }
}

/// Run complete V3 first, then at most one broad pass. Timing is diagnostic
/// only and never participates in candidate acceptance or ranking.
fn run_recovery_bootstraps(has_rt: bool, mut run_pass: impl FnMut(bool, bool) -> ScoreResult) -> ScoreResult {
    if !has_rt {
        return ScoreResult { proteins: Vec::new(), peptides: Vec::new() };
    }
    let v3_started = std::time::Instant::now();
    let initial_v3 = run_pass(false, false);
    let original_v3 = recover_zero_discoveries(initial_v3, || {
        log::info!("mass-anchor: original search has no reported 1% targets; attempting evidence-retention recovery");
        run_pass(false, true)
    });
    log::info!("bounded-rt-recovery timing stage=complete_v3 elapsed_seconds={:.6} reported_q1_targets={}",
        v3_started.elapsed().as_secs_f64(),
        original_v3.proteins.iter().filter(|p| !p.decoy && p.q_value <= 0.01).count());
    recover_zero_discoveries(original_v3, || {
        log::info!("bounded-rt-recovery: complete V3 has no reported 1% targets; starting one broad pass retain_zero=false applied_broad_floor=0.0");
        let broad_started = std::time::Instant::now();
        let recovered = run_pass(true, false);
        log::info!("bounded-rt-recovery timing stage=single_broad elapsed_seconds={:.6} reported_q1_targets={} retain_zero=false applied_broad_floor=0.0",
            broad_started.elapsed().as_secs_f64(),
            recovered.proteins.iter().filter(|p| !p.decoy && p.q_value <= 0.01).count());
        recovered
    })
}

/// Mass-only clean-anchor RT bootstrap (ms1searchpy's stage-1 trick).
///
/// On dense single-run data the per-peptide RANSAC seed is ~99% contaminated, so
/// RT calibration is fit on garbage. Instead, identify confident PROTEINS by MASS
/// ALONE (multiplicity: real proteins accumulate several tight-mass peptides while
/// decoy proteins stay singletons), then take *their* peptides as RT-calibration
/// anchors — a far cleaner set, because a present protein's tight-mass peptides are
/// enriched for genuinely-present peptides. Returns (calibration, anchor positions,
/// mass_shift_ppm, mass_sigma_ppm).
fn mass_anchor_calibration(
    reduced: &[usize],
    pfms: &[Pfm],
    matrix: &FeatureMatrix,
    db: &PeptideDb,
    abstract_rt: &[Option<f64>],
    cfg: &ScoreConfig,
    retain_zero_discovery_evidence: bool,
    broad_rt_bootstrap: bool,
) -> (Calibration, Vec<usize>, f64, f64) {
    use std::collections::{HashMap, HashSet};
    // Best mono (offset==0) match per peptide: (mass_error_ppm, reduced position).
    // A peptide is "matched" if this best error falls inside the current mass
    // window — binary set membership, exactly like ms1searchpy's `set(seqs)`.
    // Restrict to high-isotope features (ms1searchpy's `ci`=4 calibration quality)
    // when shape descriptors are present. Junk low-isotope features contribute most
    // of the random mass-coincidence background; excluding them drops the binomial's
    // p from ~0.44 toward ms1searchpy's ~0.10, which is the bulk of the mass-stage
    // depth gap. Inert (no filter) when the input lacks n_isotopes.
    const MIN_ISOTOPES: f64 = 3.0;
    let has_iso = matrix.features.iter().any(|f| f.n_isotopes.is_some());
    let mut best: HashMap<usize, (f64, usize)> = HashMap::new();
    for (pos, &i) in reduced.iter().enumerate() {
        if pfms[i].isotope_offset != 0 {
            continue;
        }
        if has_iso && matrix.features[pfms[i].feature_idx as usize].n_isotopes.unwrap_or(0.0) < MIN_ISOTOPES {
            continue;
        }
        let cand_err = pfms[i].mass_error_ppm;
        let e = best.entry(pfms[i].peptide_idx as usize).or_insert((f64::INFINITY, pos));
        let better = if cand_err.abs() < e.0.abs() {
            true
        } else if cand_err.abs() == e.0.abs() && e.0.is_finite() {
            // Exact mass-error tie (koth's exhaustive mode emits near-duplicate
            // features at the same mass). Break it on feature quality so the
            // representative feature — and thus this peptide's anchor RT — does not
            // depend on feature row order. Without this, `best` picks the first-seen
            // feature, making RANSAC calibration and the whole result swing with the
            // input's row order.
            let cf = &matrix.features[pfms[i].feature_idx as usize];
            let ef = &matrix.features[pfms[reduced[e.1]].feature_idx as usize];
            feature_cmp(cf, ef) == Ordering::Less
        } else {
            false
        };
        if better {
            *e = (cand_err, pos);
        }
    }
    if best.len() < 100 {
        return (Calibration::default(), Vec::new(), 0.0, 1.0);
    }
    // ALL mono PFMs per peptide (every feature it matches), no isotope filter —
    // used by the RT stage. A present peptide is "matched" if ANY of its features
    // is consistent in BOTH mass and RT (ms1searchpy's set membership over all
    // PFMs). Using only the best-MASS feature misses present peptides whose real
    // feature is out-massed by a random coincidence; checking every feature finds
    // the real one (tight mass AND tight RT) — this is the 7x RT multiplier.
    // ALL mono PFMs per peptide (every feature it matches), for the single-pass
    // any-feature RT stage: a peptide is "seen" if ANY of its features is consistent
    // in both mass and RT (recovers low-isotope real peptides -> the RT multiplier).
    let mut pep_pfms: HashMap<usize, Vec<(f64, usize)>> = HashMap::new();
    for (pos, &i) in reduced.iter().enumerate() {
        if pfms[i].isotope_offset != 0 {
            continue;
        }
        pep_pfms.entry(pfms[i].peptide_idx as usize).or_default().push((pfms[i].mass_error_ppm, pos));
    }
    // Canonicalize each peptide's PFM list best-feature-first (see `feature_cmp`),
    // so every downstream "best PFM" selection with a first-seen / strict-`<`
    // tie-break resolves identically regardless of the input's feature row order.
    // Without this, koth's exhaustive near-duplicate features (equal mass_error or
    // equal mass+RT deviation) let row order pick a different representative PFM,
    // perturbing the ML feature vectors and swinging the final protein count.
    for v in pep_pfms.values_mut() {
        v.sort_by(|&(_, pa), &(_, pb)| {
            let fa = &matrix.features[pfms[reduced[pa]].feature_idx as usize];
            let fb = &matrix.features[pfms[reduced[pb]].feature_idx as usize];
            feature_cmp(fa, fb)
        });
    }
    // Theoretical peptide count per (decoy, accession), and the decoy totals used
    // for the binomial background rate p = mean(decoy matched)/mean(decoy theo).
    let mut theo: rustc_hash::FxHashMap<(bool, &str), u32> = rustc_hash::FxHashMap::default();
    for pep in &db.peptides {
        for prot in &pep.proteins {
            *theo.entry((pep.decoy, prot.as_ref())).or_insert(0) += 1;
        }
    }
    let decoy_theo_sum: u64 = theo.iter().filter(|((d, _), _)| *d).map(|(_, &c)| c as u64).sum();

    // Iterative bootstrap: start at the full search tolerance, identify confident
    // proteins by the size-adjusted binomial, recalibrate mass on THEIR peptides,
    // tighten the window to +/-3 sigma, and repeat. The discrimination comes from
    // the binomial detecting each real protein's ~n_present excess of matched
    // peptides over the p*n random background — it works at high p too.
    let (mut shift, mut sigma) = (0.0f64, 3.5f64);
    let mut anchor_pos: Vec<usize> = Vec::new();
    let (mut last_nt, mut last_p) = (0usize, 0.0f64);
    let mut mass_confident: HashSet<(bool, &str)> = HashSet::new();
    // The tol=3*sigma feedback loop is not guaranteed to converge monotonically: an
    // early round admitting a few extra random matches can inflate that round's
    // recalibrated sigma, which widens next round's tolerance, admitting MORE noise,
    // in a runaway spiral. Measured repeatedly (run10, multiple configs): round 1
    // finds the MOST confident proteins with the tightest sigma, then rounds 2-4
    // monotonically WORSEN (e.g. 174->85->87->70 confident, sigma 0.64->0.74) even
    // though the loop keeps running and unconditionally uses the LAST round's
    // result. Track the best round by confident-protein count and use IT, instead
    // of assuming the final round is the best one.
    #[allow(clippy::type_complexity)]
    let mut best_round: Option<(usize, f64, f64, Vec<usize>, f64, HashSet<(bool, &str)>)> = None;
    for round in 0..5 {
        let tol = if round == 0 { cfg.mass_anchor_ppm } else { (3.0 * sigma).min(cfg.mass_anchor_ppm) };
        // Matched peptides -> per-protein positions.
        let mut prot_peps: HashMap<(bool, &str), Vec<usize>> = HashMap::new();
        for (&pi, &(err, pos)) in &best {
            if (err - shift).abs() > tol {
                continue;
            }
            let pep = &db.peptides[pi];
            for prot in &pep.proteins {
                prot_peps.entry((pep.decoy, prot.as_ref())).or_default().push(pos);
            }
        }
        let decoy_matched: u64 = prot_peps.iter().filter(|((d, _), _)| *d).map(|(_, v)| v.len() as u64).sum();
        let p = if decoy_theo_sum == 0 { 0.5 } else { (decoy_matched as f64 / decoy_theo_sum as f64).clamp(1e-4, 0.95) };
        // Size-adjusted binomial score per protein.
        let mut scored: Vec<(f64, bool, (bool, &str))> = prot_peps
            .iter()
            .map(|(&key, peps)| {
                let n = theo.get(&key).copied().unwrap_or(peps.len() as u32);
                let s = -crate::protein::binom_sf_ge_k(peps.len() as u32, n, p).log10();
                (s, key.0, key)
            })
            .collect();
        scored.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap().then_with(|| a.2.cmp(&b.2)));
        // Target/decoy FDR walk -> confident target proteins (q <= 0.05).
        let (mut t, mut d) = (0usize, 0usize);
        let mut confident: HashSet<(bool, &str)> = HashSet::new();
        for (_, decoy, key) in &scored {
            if *decoy { d += 1; } else { t += 1; }
            let q = if t == 0 { 1.0 } else { d as f64 / t as f64 };
            if q <= 0.05 && !*decoy {
                confident.insert(*key);
            }
            if t > 50 && q > 0.5 {
                break;
            }
        }
        // Anchors + recalibration set = matched peptides of confident proteins.
        anchor_pos.clear();
        let mut cerr: Vec<f64> = Vec::new();
        // Sorted key order: `best` is a HashMap, and `cerr` is reduced by
        // `clipped_sigma` below, which is a float reduction.
        let mut best_keys: Vec<usize> = best.keys().copied().collect();
        best_keys.sort_unstable();
        for pi in best_keys {
            let (err, pos) = best[&pi];
            if (err - shift).abs() > tol {
                continue;
            }
            let pep = &db.peptides[pi];
            // DECOY GUARD. `confident` is keyed (is_decoy, accession) and this
            // asks only for the TARGET entry, but a decoy peptide carries the
            // same parent accessions as the target it was shuffled from, so
            // without this check a decoy peptide of a confident target protein
            // is admitted as a mass anchor and into the recalibration sample.
            // The RT path (`est_rt_sig`) has always guarded this; the mass path
            // did not, which is the asymmetry that gave it away.
            if pep.decoy {
                continue;
            }
            if pep.proteins.iter().any(|pr| confident.contains(&(false, pr.as_ref()))) {
                anchor_pos.push(pos);
                cerr.push(err);
            }
        }
        // Recalibrate mass shift + Gaussian-core sigma on the confident peptides.
        if cerr.len() >= 50 {
            shift = median(&mut cerr.clone());
            sigma = clipped_sigma(cerr.iter().map(|&e| e - shift)).max(0.1);
        }
        last_nt = confident.len();
        last_p = p;
        mass_confident = confident.clone();
        log::info!(
            "mass-anchor round {}: tol={:.2} ppm p={:.3} -> {} confident target proteins, {} anchor peptides (shift={:.3} sigma={:.3})",
            round, tol, p, confident.len(), anchor_pos.len(), shift, sigma
        );
        let is_better = best_round.as_ref().map_or(true, |(bn, ..)| last_nt > *bn);
        if is_better {
            best_round = Some((last_nt, shift, sigma, anchor_pos.clone(), last_p, mass_confident.clone()));
        }
        if confident.is_empty() {
            break;
        }
    }
    if let Some((bn, bshift, bsigma, banchor, bp, bconf)) = best_round {
        last_nt = bn;
        shift = bshift;
        sigma = bsigma;
        anchor_pos = banchor;
        last_p = bp;
        mass_confident = bconf;
    }
    anchor_pos.sort_unstable();
    anchor_pos.dedup();
    log::info!(
        "mass-anchor bootstrap: {} confident target proteins, {} anchor peptides, p={:.3}, mass shift={:.3} ppm sigma={:.3} ppm (best-of-rounds)",
        last_nt, anchor_pos.len(), last_p, shift, sigma
    );
    // Experimental broad evidence ignores target/decoy and protein labels. The
    // mass bootstrap is unchanged; only its shift/gate defines this larger pool.
    let broad = if broad_rt_bootstrap {
        let broad_gate = cfg.mass_anchor_ppm.min((3.0 * sigma).max(1.0));
        let broad_evidence: Vec<_> = reduced.iter().enumerate().filter_map(|(pos, &i)| {
            let p = &pfms[i];
            let f = &matrix.features[p.feature_idx as usize];
            crate::broad_rt::observation(pos,p,f,&db.peptides[p.peptide_idx as usize],
                abstract_rt[pos],shift,broad_gate,has_iso)
        }).collect();
        let result = crate::broad_rt::fit(&broad_evidence);
        drop(broad_evidence);
        result
    } else { None };
    let broad_floor = 0.0; // Bounded candidate: diagnostic ambiguity floor is not applied.
    let freeze_rt = broad.is_some();
    let (cal, clean) = if let Some(b) = broad {
        (b.calibration, b.clean)
    } else {
        if anchor_pos.is_empty() {
            return (Calibration::default(), Vec::new(), shift, sigma);
        }
        // 6. RANSAC RT calibration on the (abstract, observed) pairs of the anchors.
        let pair_positions: Vec<usize> = if broad_rt_bootstrap {
            anchor_pos.iter().copied().filter(|&pos| abstract_rt[pos].is_some()).collect()
        } else { Vec::new() };
        let pair_source = if broad_rt_bootstrap { &pair_positions } else { &anchor_pos };
        let pairs: Vec<(f64, f64)> = pair_source
            .iter()
            .filter_map(|&pos| abstract_rt[pos].map(|a| (a, matrix.features[pfms[reduced[pos]].feature_idx as usize].rt)))
            .collect();
        if pairs.len() < 50 {
            return (Calibration::default(), anchor_pos, shift, sigma);
        }
        let (mut cal, mut inl) = Calibration::fit_ransac_seeded(&pairs, 2.5, 4000, cfg.effective_rt_seed());
        for thresh in [1.0_f64, 0.5, 0.3] {
            let sub: Vec<(f64, f64)> = inl.iter().map(|&k| pairs[k]).collect();
            if sub.len() < 50 { break; }
            let (c, si) = Calibration::fit_ransac_seeded(&sub, thresh, 4000, cfg.effective_rt_seed());
            if si.len() < 50 { break; }
            inl = si.iter().map(|&j| inl[j]).collect();
            cal = c;
        }
        // keep only the RANSAC-inlier anchors (cleanest)
        let clean: Vec<usize> = inl.iter().map(|&k| if broad_rt_bootstrap {
            pair_positions[k]
        } else {
            // Preserve V3 bookkeeping as well as its ranking on the first path.
            anchor_pos[k.min(anchor_pos.len()-1)]
        }).collect();
        let mut resid: Vec<f64> = inl.iter().map(|&k| (pairs[k].1 - cal.apply(pairs[k].0)).abs()).collect();
        let med = if resid.is_empty() { f64::NAN } else { median(&mut resid) };
        log::info!(
            "mass-anchor RT calibration: {} RANSAC inliers, median |residual| = {:.2} min",
            clean.len(), med
        );

        (cal, clean)
    };

    // === Stage 2 (SINGLE deterministic pass; ms1searchpy-faithful) ===
    // The RT calibration `cal` is now FIXED (fit on the clean mass-stage anchors).
    // The previous iterative recalibration loop oscillated (calibration <-> confident
    // set feedback) and made the result DB-composition-dependent. Here we estimate
    // RT_sigma ONCE and do ONE any-feature mass+RT binomial pass -> deterministic
    // and stable. A peptide is "seen" if ANY of its features is consistent in both
    // mass and RT (the low-isotope-recovery that gives the RT multiplier).
    let rt_diff_at = |cal: &Calibration, pos: usize| -> Option<f64> {
        abstract_rt[pos].map(|a| matrix.features[pfms[reduced[pos]].feature_idx as usize].rt - cal.apply(a))
    };
    let mass_gate = (3.0_f64 * sigma).min(cfg.mass_anchor_ppm);

    // RT accuracy from confident-protein peptides' best-RT (real-feature) residuals,
    // Gaussian core, clamped to weight RT ~its true accuracy in the deviation.
    // Floor applied to every RT-sigma estimate after a stage-1 collapse (see the
    // low-input note in the stage loop): an estimate from <25 proteins is not
    // trustworthy and the deviation search under-gates. 0.0 = no floor.
    let rt_sig_floor: std::cell::Cell<f64> = std::cell::Cell::new(broad_floor);
    let set_rt_floor = |floor: f64| rt_sig_floor.set(floor.max(broad_floor));
    // Multiplier on every RT-sigma estimate (UNO_SMOOTH candidate runs). 1.0 = off.
    let rt_sig_mult: std::cell::Cell<f64> = std::cell::Cell::new(1.0);
    // UNO_SMOOTH: choose between calibration candidates / stages / ML by a soft
    // depth (sum over targets of clamp(1 - q/0.05, 0, 1)) instead of the count at
    // q<=1%, which is a step function on short lists (one decoy needs 100 targets).
    // Calibration-candidate selection (default level 2, shipped 2026-09-04).
    // The stage loop and the ML stage each pick winners by the count of proteins
    // at q<=1%, which is a step function on short lists (one decoy needs 100
    // targets above it), so at low input two near-identical runs can land at
    // 187 and 866 proteins. Level 2 runs the unchanged pipeline at two RT sigmas
    // (x1.0 and x1.5) and keeps the candidate whose FINAL list has the larger
    // soft depth, sum over targets of clamp(1 - q/0.05, 0, 1). Validated on the
    // 14-file PXD066701 input series + 6 primary files: 1 ng rep4 187 -> 768,
    // 1 ng rep2 688 -> 925 at 1% entrapment FDR, 16 of 20 files byte-identical,
    // wall time +50-80% (the ML stage runs twice).
    // UNO_SMOOTH=0 restores the single-candidate path; 1 = soft objective also
    // inside the stage loop and ML take (misaligned with entrapment depth,
    // 10 ng rep1 1448 -> 1200, do not use); 3 = candidates chosen by the final
    // internal q<=1% count (worse: internal FDR is a poor proxy, 100 ng rep3 -298).
    let smooth_level: u8 = std::env::var("UNO_SMOOTH").ok().and_then(|v| v.parse().ok()).unwrap_or(2);
    let smooth = smooth_level == 1;
    let multi_cand = smooth_level >= 1;
    let est_rt_sig = |cal: &Calibration, confident: &HashSet<String>, m_shift: f64, m_gate: f64| -> f64 {
        let mut rr: Vec<f64> = Vec::new();       // best rt_diff per peptide (uno legacy)
        let mut all_rd: Vec<f64> = Vec::new();   // ALL confident-peptide PFM rt_diffs (ms1's noisygaus input)
        // Sorted key order, not `&pep_pfms` directly: HashMap iteration is
        // process-random, and while the resulting SET of residuals is the same
        // either way, the reducers below are float reductions whose result
        // depends on the order they see. Deterministic here, and `clipped_sigma`
        // sorts defensively as well.
        let mut pep_keys: Vec<usize> = pep_pfms.keys().copied().collect();
        pep_keys.sort_unstable();
        for pi in pep_keys {
            let pfm_list = &pep_pfms[&pi];
            let pep = &db.peptides[pi];
            if pep.decoy || !pep.proteins.iter().any(|pr| confident.contains(pr.as_ref())) {
                continue;
            }
            let mut best_rd: Option<f64> = None;
            for &(err, pos) in pfm_list {
                if (err - m_shift).abs() > m_gate { continue; }
                if let Some(rd) = rt_diff_at(cal, pos) {
                    all_rd.push(rd);
                    if best_rd.map_or(true, |b: f64| rd.abs() < b.abs()) { best_rd = Some(rd); }
                }
            }
            if let Some(rd) = best_rd { rr.push(rd); }
        }
        if let Some(v) = std::env::var("UNO_RT_SIG").ok().and_then(|s| s.parse::<f64>().ok()) {
            return v;
        }
        if std::env::var("UNO_RT_SIG_GCORE").is_ok() {
            return gaussian_core_sigma(rr.into_iter(), 15.0).clamp(0.3, 5.0);
        }
        // ms1searchpy's `calibrate_RT_gaus`: fit a Gaussian+baseline (noisygaus) to the
        // rt_diff histogram of the CONFIDENT peptides' best matches (ms1's df1[top_peps]),
        // take the peak sigma = the REAL per-peptide RT accuracy. A proper Gaussian fit
        // (vs `clipped_sigma`'s 3.5-clamped whole-spread or the HWHM's over-tight peak).
        // NOTE: faithfully porting ms1's noisygaus here gives rt_sig~2.15 and 791 @1%
        // (WORSE than clipped's 904) — ms1's tight gate is calibrated for DeepLC's ~1min
        // RT accuracy; uno's CNN RT is ~3min so the tight gate rejects real peptides.
        // The algorithm is faithful; uno's RT MODEL is the limiting input. Default =
        // clipped (works best with the CNN). UNO_RT_NOISYGAUS for the faithful ms1 fit.
        let _ = &all_rd;
        if std::env::var("UNO_RT_NOISYGAUS").is_ok() {
            let (_shift, s) = noisygaus_fit(&rr, 15.0, 300);
            return s.clamp(0.3, 5.0);
        }
        clipped_sigma(rr.into_iter()).clamp(0.3, 3.5).max(rt_sig_floor.get()) * rt_sig_mult.get()
    };

    // Per-stage MASS recalibration on the confident set (ms1searchpy's biggest p
    // drop is mass recal: p 0.68->0.30 at stage 1, then RT/ML tighten to 0.136;
    // uno's bootstrap tightened mass ONCE then froze it, leaving p~0.44). Each stage,
    // refit the mass shift + Gaussian-core sigma on the CURRENT confident target
    // peptides' best RT-consistent match, and tighten the gate to 3*sigma (NEVER
    // widen -- monotone, so coincidental matches fall outside the shrinking window
    // while real ones stay in). Returns (shift, sigma, gate). Gated by UNO_STAGE_MASS
    // (default on); when off, stages reuse the frozen bootstrap mass (old behavior).
    let recal_mass = |conf: &HashSet<String>, cal: &Calibration, cur_shift: f64, cur_sigma: f64, cur_gate: f64| -> (f64, f64, f64) {
        let mut errs: Vec<f64> = Vec::new();
        for (&pi, pfm_list) in &pep_pfms {
            let pep = &db.peptides[pi];
            if pep.decoy || !pep.proteins.iter().any(|pr| conf.contains(pr.as_ref())) { continue; }
            // mass error of this peptide's best RT-consistent match (min |rt_diff|
            // within the current mass gate) -- the same match the deviation pass keys on.
            let mut best: Option<(f64, f64)> = None; // (|rt_diff|, mass_err)
            for &(err, pos) in pfm_list {
                if (err - cur_shift).abs() > cur_gate { continue; }
                if let Some(rd) = rt_diff_at(cal, pos) {
                    if best.map_or(true, |(b, _): (f64, f64)| rd.abs() < b) { best = Some((rd.abs(), err)); }
                }
            }
            if let Some((_, err)) = best { errs.push(err); }
        }
        if errs.len() < 50 { return (cur_shift, cur_sigma, cur_gate); }
        let shift = median(&mut errs.clone());
        let sigma = clipped_sigma(errs.iter().map(|&e| e - shift)).max(0.05);
        // never widen the gate: monotone tightening only.
        let gate = (3.0 * sigma).min(cur_gate);
        (shift, sigma, gate)
    };
    // Default OFF: measured negative (per-stage recal RE-ESTIMATES sigma LOOSER than
    // the bootstrap's tightest-round subset -- 0.63->0.91 ppm -- so it widens, not
    // tightens; CNN 904->705, DeepLC 723->687). The bootstrap already tightens mass
    // optimally; mass is not uno's p bottleneck. Kept env-gated for future work.
    let stage_mass = std::env::var("UNO_STAGE_MASS").map(|v| v == "1").unwrap_or(false);

    // One deviation pass (given a calibration + rt_sig): per peptide, ANY-FEATURE
    // match by MINIMUM COMBINED deviation dev = mass_z^2 + (rt_diff/rt_sig)^2 (finds
    // the real feature, tight in BOTH dims; rejects coincidences tight in one). Scan
    // the cutoff, protein binomial, target/decoy FDR. Returns the best result and
    // the confident (q<=1%) target accessions (for two-stage recalibration).
    let run_deviation = |cal: &Calibration, rt_sig: f64, m_shift: f64, m_sigma: f64, m_gate: f64| -> (usize, f64, Vec<ProteinId>, HashSet<String>) {
        let mut pep_dev: HashMap<usize, (bool, f64, usize)> = HashMap::new();
        // Per peptide, ALL (dev, feature) pairs across its PFMs (charge states /
        // scans), so feature-competition can keep a peptide alive until EVERY one of
        // its peaks is claimed (ms1searchpy's banned_dict), not just its single best.
        // Feature-quality gate (koth combined_score / cosine etc.): exclude peptides
        // backed by low-quality (junk) features, which contribute disproportionately
        // many coincidental decoy/wrong matches. biosaur's combined_score is const
        // 1.0 so this is inert there; koth's is a real isotope-fit quality. Env knob.
        let cs_gate: f64 = std::env::var("UNO_CS_GATE").ok().and_then(|s| s.parse().ok()).unwrap_or(0.0);
        // `pep_all_feats` (the per-PFM (dev, feature) list) is ONLY read by the
        // UNO_MULTI_FEAT cut-sweep branch; when that's off (default) building it is
        // pure waste (a Vec alloc + 17M pushes + a map insert per stage). Read the
        // flag up front and skip the work entirely in the common case.
        let want_feats = std::env::var("UNO_MULTI_FEAT").is_ok();
        // Per-peptide min-combined-deviation match over all its PFMs. Independent per
        // peptide -> parallel map, deterministic collect (results reinserted into the
        // maps; key order is irrelevant to correctness).
        let per_pep: Vec<(usize, (bool, f64, usize), Vec<(f64, usize)>)> = pep_pfms
            .par_iter()
            .filter_map(|(&pi, pfm_list)| {
                let decoy = db.peptides[pi].decoy;
                let mut best: Option<(f64, usize)> = None;
                let mut feats: Vec<(f64, usize)> = Vec::new();
                // Every PFM of a peptide shares the same abstract (predicted) RT, so the
                // calibration interpolation is identical across the list -- compute it
                // ONCE per peptide instead of per PFM (17M `cal.apply` -> ~1.1M).
                let pred = pfm_list.first().and_then(|&(_, pos)| abstract_rt[pos]).map(|a| cal.apply(a));
                for &(err, pos) in pfm_list {
                    let mz = (err - m_shift) / m_sigma;
                    if mz.abs() > m_gate / m_sigma { continue; }
                    if cs_gate > 0.0 && matrix.features[pfms[reduced[pos]].feature_idx as usize].combined_score < cs_gate { continue; }
                    if let Some(pr) = pred {
                        let rd = matrix.features[pfms[reduced[pos]].feature_idx as usize].rt - pr;
                        let dev = mz * mz + (rd / rt_sig).powi(2);
                        if best.map_or(true, |(b, _)| dev < b) { best = Some((dev, pos)); }
                        if want_feats { feats.push((dev, pfms[reduced[pos]].feature_idx as usize)); }
                    }
                }
                best.map(|(dev, pos)| (pi, (decoy, dev, pos), feats))
            })
            .collect();
        let mut pep_all_feats: HashMap<usize, Vec<(f64, usize)>> = HashMap::with_capacity(per_pep.len());
        for (pi, dev, feats) in per_pep {
            pep_dev.insert(pi, dev);
            pep_all_feats.insert(pi, feats);
        }
        // Continuous-posterior scorer (UNO_CONTINUOUS): replace the ~30-point
        // deviation-cut SWEEP with ONE weighted pass. Each peptide's `dev` becomes a
        // posterior P(correct) from the target/decoy dev distributions; proteins are
        // scored by weighted evidence. No hyperparameter scan, no threshold-selection
        // bias, and it uses the full `dev` gradient instead of a hard gate.
        if std::env::var("UNO_CONTINUOUS").is_ok() {
            let scale = db.n_targets as f64 / db.n_decoys.max(1) as f64;
            let weights = dev_posteriors(&pep_dev, scale);
            let mut observed: rustc_hash::FxHashMap<(bool, &str), rustc_hash::FxHashSet<usize>> = rustc_hash::FxHashMap::default();
            let mut pep_feat: rustc_hash::FxHashMap<usize, Vec<usize>> = rustc_hash::FxHashMap::default();
            for (&pi, &(decoy, _dev, pos)) in &pep_dev {
                if weights.get(&pi).copied().unwrap_or(0.0) <= 1e-3 {
                    continue;
                }
                pep_feat.insert(pi, vec![pfms[reduced[pos]].feature_idx as usize]);
                for prot in &db.peptides[pi].proteins {
                    observed.entry((decoy, prot.as_ref())).or_default().insert(pi);
                }
            }
            let scored = crate::protein::score_proteins_weighted(
                observed,
                &pep_feat,
                &weights,
                &theo,
                decoy_theo_sum,
            );
            let nt = scored.iter().filter(|p| !p.decoy && p.q_value <= 0.01).count();
            let conf: HashSet<String> = scored
                .iter()
                .filter(|p| !p.decoy && p.q_value <= 0.01)
                .map(|p| p.accession.clone())
                .collect();
            log::info!(
                "continuous-posterior scorer: {} target proteins @1% FDR ({} peptides weighted, no cut grid)",
                nt,
                pep_feat.len()
            );
            return (nt, 0.0, scored, conf);
        }
        if let Ok(path) = std::env::var("UNO_DUMP_PEPSET") {
            use std::io::Write;
            if let Ok(f) = std::fs::File::create(&path) {
                let mut w = std::io::BufWriter::new(f);
                // One row per (peptide, MS1 feature) match with mass/RT residuals,
                // so an external scorer/ML (ms1searchpy) gets the same per-PFM
                // evidence its own PFMs carry. feature = peak id (ids).
                let _ = writeln!(w, "peptide\tdecoy\tdev\tfeature\tmass_ppm\trt_diff\tproteins");
                for (&pi, &(decoy, _, _)) in &pep_dev {
                    let pep = &db.peptides[pi];
                    let prots: Vec<&str> = pep.proteins.iter().map(|p| p.as_ref()).collect();
                    let pstr = prots.join(";");
                    // recompute per-PFM mass/RT residuals against the winning calibration
                    if let Some(list) = pep_pfms.get(&pi) {
                        for &(err, pos) in list {
                            let mz = (err - m_shift) / m_sigma;
                            if mz.abs() > m_gate / m_sigma { continue; }
                            if let Some(rd) = rt_diff_at(cal, pos) {
                                let dev = mz * mz + (rd / rt_sig).powi(2);
                                let feat = pfms[reduced[pos]].feature_idx as usize;
                                let _ = writeln!(w, "{}\t{}\t{:.3}\t{}\t{:.4}\t{:.4}\t{}",
                                    pep.sequence.as_ref(), decoy as u8, dev, feat, err, rd, pstr);
                            }
                        }
                    }
                }
            }
        }
        let mut best_proteins: Vec<ProteinId> = Vec::new();
        let mut best_n = 0usize;
        let mut best_cut = 0.0f64;
        // ms1searchpy ENSEMBLE: final_iteration runs the greedy at 10 inclusion
        // thresholds and scores each protein by the MEAN of its per-pass binomial
        // scores, instead of picking one best threshold. Accumulate (score_sum,
        // max_obs, theo) per protein across the cut sweep. UNO_ENSEMBLE.
        let ensemble = std::env::var("UNO_ENSEMBLE").is_ok();
        let multi_feat = std::env::var("UNO_MULTI_FEAT").is_ok();
        let no_feat_compete = std::env::var("UNO_NO_FEATURE_COMPETE").is_ok();
        let log_pepcount = std::env::var("UNO_LOG_PEPCOUNT").is_ok();
        let mut acc: HashMap<(bool, String), (f64, usize, usize)> = HashMap::new();
        // Deviation-cut sweep extent. The sweep is a `par_iter` over cuts; the best
        // cut sits near ~6, and the high-dev tail cuts (>8) are BOTH the most expensive
        // tasks (they admit ~all peptides, so their build + greedy dominate the sweep
        // makespan) AND actively harmful: the raw best-cut-by-target-count selection
        // over-fits to a loose, decoy-heavy high cut. Capping at 8 drops them.
        // Measured: run10 BYTE-IDENTICAL to the old 25.0 cap (its best cut was <=8);
        // deephela 1954 -> 2044 honest @1% (+4.6%, the loose cut regularized away) with
        // the per-stage sweep roughly halved (~189s -> ~97s across 4 stages). Both
        // deterministic. `cut_step` (fine step over 1..12) and `cut_max` stay env-
        // tunable; on this bigger single-input the tail was pure cost + FDR noise.
        let cut_step: f64 = std::env::var("UNO_CUT_STEP").ok().and_then(|s| s.parse().ok()).unwrap_or(0.5);
        let cut_max: f64 = std::env::var("UNO_CUT_MAX").ok().and_then(|s| s.parse().ok()).unwrap_or(8.0);
        let cut_grid: Vec<f64> = {
            let mut g: Vec<f64> = Vec::new();
            let mut c = 1.0;
            while c <= 12.0 { g.push(c); c += cut_step; }
            g.extend_from_slice(&[13.0, 14.0, 16.0, 18.0, 20.0, 25.0]);
            g.retain(|&x| x <= cut_max);
            g
        };
        // Each deviation cut is an INDEPENDENT protein-inference pass: build the
        // matched-peptide sets at that cut, then run the greedy competition. This is
        // the dominant cost (a ~30-way greedy sweep). Fan the cuts across cores; the
        // `collect` preserves cut_grid order, so the downstream best-cut pick and the
        // ensemble accumulation stay byte-identical to the old serial loop.
        let _sweep_t0 = std::time::Instant::now();
        let mut per_cut: Vec<(f64, usize, usize, usize, Vec<ProteinId>)> = cut_grid
            .par_iter()
            .map(|&cut| {
                // Matched peptides at this deviation cut -> protein -> peptide set.
                // GREEDY PARSIMONY (ms1searchpy final_iteration): assign each peptide
                // to one protein and pair-delete target/decoy.
                let mut observed: rustc_hash::FxHashMap<(bool, &str), rustc_hash::FxHashSet<usize>> = rustc_hash::FxHashMap::default();
                let mut decoy_seen: rustc_hash::FxHashSet<usize> = rustc_hash::FxHashSet::default();
                // Peptide -> ALL MS1 features it matched at this cut (feature-level
                // competition: a peptide survives until every one of its peaks is
                // claimed, not just its single best PFM).
                let mut pep_feat: rustc_hash::FxHashMap<usize, Vec<usize>> = rustc_hash::FxHashMap::default();
                for (&pi, &(decoy, dev, pos)) in &pep_dev {
                    if dev <= cut {
                        if decoy {
                            decoy_seen.insert(pi);
                        }
                        // Multi-feature retention (UNO_MULTI_FEAT) HURT honest depth;
                        // default = single best PFM (the single-feature reduction acts
                        // as a useful selectivity filter).
                        let feats: Vec<usize> = if multi_feat {
                            let mut f: Vec<usize> = pep_all_feats
                                .get(&pi)
                                .map(|v| v.iter().filter(|(d, _)| *d <= cut).map(|(_, f)| *f).collect())
                                .unwrap_or_default();
                            f.sort_unstable();
                            f.dedup();
                            if f.is_empty() { f.push(pfms[reduced[pos]].feature_idx as usize); }
                            f
                        } else {
                            vec![pfms[reduced[pos]].feature_idx as usize]
                        };
                        pep_feat.insert(pi, feats);
                        for prot in &db.peptides[pi].proteins {
                            observed.entry((decoy, prot.as_ref())).or_default().insert(pi);
                        }
                    }
                }
                // ms1searchpy-style feature-level competition + dynamic background
                // rate is the default. UNO_NO_FEATURE_COMPETE=1 restores the old scorer.
                let scored = if no_feat_compete {
                    crate::protein::score_proteins_greedy(
                        observed,
                        &theo,
                        decoy_theo_sum,
                        decoy_seen.len(),
                        0.01,
                    )
                } else {
                    crate::protein::score_proteins_feature_compete(
                        observed,
                        &pep_feat,
                        &theo,
                        decoy_theo_sum,
                        0.01,
                    )
                };
                let nt = scored.iter().filter(|p| !p.decoy && p.q_value <= 0.01).count();
                let tgt = if log_pepcount {
                    pep_dev.iter().filter(|(_, &(d, dev, _))| !d && dev <= cut).count()
                } else {
                    0
                };
                (cut, nt, tgt, decoy_seen.len(), scored)
            })
            .collect();
        log::info!("run_deviation cut-sweep: {} cuts (max {:.0}) in {:.2}s", cut_grid.len(), cut_max, _sweep_t0.elapsed().as_secs_f64());
        // Deterministic serial reduction over the ordered per-cut results: accumulate
        // ensemble sums in cut order, then pick the FIRST (smallest) cut achieving the
        // max target count -- identical to the old `if nt > best_n` scan.
        if ensemble {
            for (_, _, _, _, scored) in &per_cut {
                for p in scored {
                    let e = acc.entry((p.decoy, p.accession.clone())).or_insert((0.0, 0, 0));
                    e.0 += p.binomial_score;
                    e.1 = e.1.max(p.observed_peptides);
                    e.2 = p.theoretical_peptides;
                }
            }
        }
        let best_idx = select_deviation_cut(per_cut.iter().map(|cut| cut.4.as_slice()), retain_zero_discovery_evidence);
        if let Some(i) = best_idx {
            let (cut, nt, tgt, dec) = (per_cut[i].0, per_cut[i].1, per_cut[i].2, per_cut[i].3);
            best_n = nt;
            best_cut = cut;
            if log_pepcount {
                log::info!("run_deviation cut={:.1}: {} proteins; distinct peptides target={} decoy={} (decoy/target={:.2})",
                    cut, nt, tgt, dec, dec as f64 / tgt.max(1) as f64);
            }
            if nt == 0 {
                let training = per_cut[i].4.iter().filter(|p| !p.decoy && p.q_value <= 0.25).count();
                log::info!("run_deviation: no 1% discoveries; retaining cut {:.1} with {} training targets at q<=25%", cut, training);
            }
            best_proteins = std::mem::take(&mut per_cut[i].4);
        }
        if ensemble {
            // Mean per-pass binomial score; pair-compete target/decoy; (d+1)/t FDR.
            let nc = cut_grid.len() as f64;
            let mut ep: Vec<ProteinId> = acc
                .iter()
                .map(|((decoy, a), &(s, obs, theo))| ProteinId {
                    accession: a.clone(),
                    decoy: *decoy,
                    observed_peptides: obs,
                    theoretical_peptides: theo,
                    coverage: if theo > 0 { obs as f64 / theo as f64 } else { 0.0 },
                    binomial_score: s / nc,
                    mean_rt_diff: 0.0,
                    mean_mass_z: 0.0,
                    mean_pep_score: 0.0,
                    q_value: 1.0,
                })
                .collect();
            ep.sort_by(|a, b| b.binomial_score.partial_cmp(&a.binomial_score).unwrap()
                .then_with(|| a.accession.cmp(&b.accession)).then(a.decoy.cmp(&b.decoy)));
            let score_by: HashMap<(bool, String), f64> =
                ep.iter().map(|p| ((p.decoy, p.accession.clone()), p.binomial_score)).collect();
            ep.retain(|p| match score_by.get(&(!p.decoy, p.accession.clone())) {
                Some(&cs) => p.binomial_score > cs || (p.binomial_score == cs && !p.decoy),
                None => true,
            });
            let (mut t, mut d) = (0usize, 0usize);
            for p in ep.iter_mut() {
                if p.decoy { d += 1; } else { t += 1; }
                p.q_value = if t == 0 { 1.0 } else { (d as f64 + 1.0) / t as f64 };
            }
            let mut run = f64::INFINITY;
            for i in (0..ep.len()).rev() { run = run.min(ep[i].q_value); ep[i].q_value = run.min(1.0); }
            let ne = ep.iter().filter(|p| !p.decoy && p.q_value <= 0.01).count();
            log::info!("run_deviation ENSEMBLE over {} cuts: {} proteins (best-cut was {})", cut_grid.len(), ne, best_n);
            best_n = ne;
            best_proteins = ep;
        }
        let conf: HashSet<String> = best_proteins
            .iter()
            .filter(|pp| !pp.decoy && pp.q_value <= 0.01)
            .map(|pp| pp.accession.clone())
            .collect();
        (best_n, best_cut, best_proteins, conf)
    };

    // Stage 1: deviation on the mass-bootstrap calibration.
    let mass_conf_accs: HashSet<String> =
        mass_confident.iter().filter(|(d, _)| !d).map(|(_, a)| a.to_string()).collect();

    let build_cal_pairs = |conf: &HashSet<String>, cal: &Calibration, m_shift: f64, m_gate: f64| -> Vec<(f64, f64)> {
        let mut cp: Vec<(f64, f64)> = Vec::new();
        for (&pi, pfm_list) in &pep_pfms {
            let pep = &db.peptides[pi];
            if pep.decoy || !pep.proteins.iter().any(|pr| conf.contains(pr.as_ref())) { continue; }
            let mut best: Option<(f64, usize)> = None;
            for &(err, pos) in pfm_list {
                if (err - m_shift).abs() > m_gate { continue; }
                if let Some(rd) = rt_diff_at(cal, pos) {
                    if best.map_or(true, |(b, _): (f64, usize)| rd.abs() < b.abs()) { best = Some((rd, pos)); }
                }
            }
            if let Some((_, pos)) = best {
                if let Some(a) = abstract_rt[pos] {
                    cp.push((a, matrix.features[pfms[reduced[pos]].feature_idx as usize].rt));
                }
            }
        }
        cp.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap().then(a.1.partial_cmp(&b.1).unwrap()));
        cp
    };
    let n_stages: usize = std::env::var("UNO_RT_STAGES").ok().and_then(|s| s.parse().ok()).unwrap_or(4);

    // Stage 1 + the iterative multi-stage RT calibration (ms1searchpy drives the
    // random-match rate p 0.68->0.30->0.28->0.12 across stages; uno used to do ONE
    // recalibration and stop at p~0.23-0.44). Each stage recalibrates RT on the
    // CURRENT confident set, tightens rt_sig, and re-searches; the confident set
    // purifies each round so the gate tightens and p falls. Deterministic:
    // cal_pairs are total-ordered before the fit. Keeps the best-scoring stage.
    // `floor` is the RT-sigma floor (0 = none) and `seed` extra proteins added to
    // the stage-2 confident set; both are 0/empty on the normal path.
    let run_stages = |floor: f64, mult: f64, seed: &HashSet<String>| {
        set_rt_floor(floor);
        rt_sig_mult.set(mult);
        let rt_sig1 = est_rt_sig(&cal, &mass_conf_accs, shift, mass_gate);
        let (n1, cut1, prot1, conf1) = run_deviation(&cal, rt_sig1, shift, sigma, mass_gate);
        let mut best_soft = soft_depth(&prot1);
        let (mut best_n, mut best_cut, mut best_proteins, mut best_cal) =
            (n1, cut1, prot1.clone(), cal.clone());
        let (mut best_shift, mut best_sigma, mut best_gate) = (shift, sigma, mass_gate);
        let mut cur_conf: HashSet<String> = conf1.union(seed).cloned().collect();
        let mut cur_cal = cal.clone();
        let (mut cur_shift, mut cur_sigma, mut cur_gate) = (shift, sigma, mass_gate);
        for stage in 2..=n_stages.max(2) {
            let (s_shift, s_sigma, s_gate) = if stage_mass {
                recal_mass(&cur_conf, &cur_cal, cur_shift, cur_sigma, cur_gate)
            } else { (shift, sigma, mass_gate) };
            let cal_pairs = build_cal_pairs(&cur_conf, &cur_cal, s_shift, s_gate);
            if cal_pairs.len() < 100 { break; }
            let cal_s = crate::broad_rt::stage_map(freeze_rt.then_some(&cal), || Calibration::fit_robust(&cal_pairs, 60, 200));
            let rt_sig_s = est_rt_sig(&cal_s, &cur_conf, s_shift, s_gate);
            let (ns, cuts, prots, confs) = run_deviation(&cal_s, rt_sig_s, s_shift, s_sigma, s_gate);
            log::info!("mass-anchor stage {}: {} confident (rt_sig {:.2}, mass {:.2}+/-{:.2}ppm gate {:.2}, {} conf-peptides)",
                stage, ns, rt_sig_s, s_shift, s_sigma, s_gate, cal_pairs.len());
            let soft_s = soft_depth(&prots);
            let better = if smooth { soft_s > best_soft } else { ns > best_n };
            if better {
                best_n = ns; best_cut = cuts; best_proteins = prots; best_cal = cal_s.clone();
                best_shift = s_shift; best_sigma = s_sigma; best_gate = s_gate; best_soft = soft_s;
            }
            cur_shift = s_shift; cur_sigma = s_sigma; cur_gate = s_gate;
            if smooth { if !better { break; } } else if confs.len() <= cur_conf.len() { break; }
            cur_cal = cal_s;
            cur_conf = confs;
        }
        (best_n, best_cut, best_proteins, best_cal, best_shift, best_sigma, best_gate, conf1, rt_sig1, best_soft)
    };
    // Candidate calibration runs (see the UNO_SMOOTH note above). Level 0: one
    // run at the estimated sigma, plus the seeded low-input retries when stage 1
    // collapses; the candidate with the most proteins at q<=1% wins and the ML
    // stage runs once on it. Levels >=1 add the 1.5x-sigma candidate, run the ML
    // stage on EVERY candidate, and pick by the final list.
    type StageRun = (usize, f64, Vec<ProteinId>, Calibration, f64, f64, f64, HashSet<String>, f64, f64);
    let empty: HashSet<String> = HashSet::new();
    let mut cands: Vec<(f64, f64)> = vec![(0.0, 1.0)];
    if multi_cand { cands.push((0.0, 1.5)); }
    let mut runs: Vec<(StageRun, f64, f64)> = Vec::new();
    for &(fl, mu) in &cands {
        let r = run_stages(fl, mu, &empty);
        if multi_cand {
            log::info!("mass-anchor candidate (sigma x{:.1}): {} confident, soft depth {:.1}", mu, r.0, r.9);
        }
        runs.push((r, fl, mu));
    }
    // Low-input collapse (PXD066701 200 pg / 1 ng files, 2026-09-04): when stage 1
    // finds almost nothing at q<=1%, the stage loop has no confident set to
    // recalibrate on and the run exits EMPTY although the mass-anchor bootstrap
    // had a usable (small) confident set and a good RANSAC RT fit. Retry with the
    // bootstrap proteins seeding stage 2, once at the estimated sigma (which is
    // unreliable from <25 proteins) and once floored at 1.0 min. Only fires when
    // stage 1 collapses; every other run is unchanged.
    let base_n = runs.iter().map(|(r, _, _)| r.0).max().unwrap_or(0);
    if base_n < 20 && !mass_conf_accs.is_empty() {
        for (fl, mu) in [(0.0_f64, 1.0_f64), (1.0, 1.0)] {
            let r = run_stages(fl, mu, &mass_conf_accs);
            log::info!(
                "mass-anchor stage 1 collapsed ({} confident); retry seeded from {} bootstrap proteins, RT sigma floor {:.1}: {} confident (soft {:.1})",
                base_n, mass_conf_accs.len(), fl, r.0, r.9
            );
            runs.push((r, fl, mu));
        }
    }
    // The ML stage as a function of one calibration candidate (see below).
    let run_ml = |best_proteins: &Vec<ProteinId>, best_cal: &Calibration, best_shift: f64, best_sigma: f64, best_gate: f64, best_n: usize| -> Option<FinalizedCandidate> {
        // ms1searchpy's semi-supervised positive set is peptides of 25%-FDR proteins
        // (top_25_targets), NOT 1% -- a loose gate gives far more positive training
        // signal. (Validated: 25%-label ML peptide-filter -> 899->1073 honest @1%.)
        let conf: HashSet<String> = best_proteins
            .iter()
            .filter(|p| !p.decoy && p.q_value <= 0.25)
            .map(|p| p.accession.clone())
            .collect();
        if cfg.ml_rescore && conf.len() >= (if smooth { 10 } else { 20 }) {
            // One PFM per peptide: the MIN combined-deviation (mass_z^2 + rt_z^2) mono
            // feature under the winning calibration. Scoring every PFM and taking the
            // per-peptide MAX inflates random peptides (max-of-many); ms1searchpy
            // applies the ML to mass+RT-gated PFMs, so we gate to the best one here.
            let rt_sig_ml = est_rt_sig(&best_cal, &conf, best_shift, best_gate).max(0.5);
            // combined_score gate on the ML CANDIDATE set (koth feature-quality ML
            // score). Junk low-quality features carry most of the coincidental yeast
            // matches that flood the mid-tier (@5%); a real human peptide keeps its
            // high-quality feature, a yeast peptide loses its lucky junk one. Applied
            // here (unlike UNO_CS_GATE which only touches run_deviation) because THIS
            // pep_best set feeds the qpred ensemble that produces the final count.
            let cs_gate_ml: f64 = std::env::var("UNO_CS_GATE_ML").ok().and_then(|s| s.parse().ok()).unwrap_or(0.0);
            let mut pep_best: HashMap<usize, (bool, usize, f64)> = HashMap::new();
            for (&pi, pfm_list) in &pep_pfms {
                let decoy = db.peptides[pi].decoy;
                let mut best: Option<(f64, usize)> = None;
                for &(err, pos) in pfm_list {
                    let mz = (err - best_shift) / best_sigma;
                    if mz.abs() > best_gate / best_sigma {
                        continue;
                    }
                    if cs_gate_ml > 0.0 && matrix.features[pfms[reduced[pos]].feature_idx as usize].combined_score < cs_gate_ml {
                        continue;
                    }
                    if let Some(rd) = rt_diff_at(&best_cal, pos) {
                        let dev = mz * mz + (rd / rt_sig_ml).powi(2);
                        if best.map_or(true, |(b, _)| dev < b) {
                            best = Some((dev, pos));
                        }
                    }
                }
                if let Some((dev, pos)) = best {
                    pep_best.insert(pi, (decoy, pos, dev));
                }
            }
            // Carry this candidate's RT+mass-gated peptide→feature map with its
            // result. Publishing here would leak the last evaluated candidate's
            // evidence into the winner's grouping and quantification outputs.
            let mut peptides: Vec<(u32, u32, bool)> = pep_best
                .iter()
                .map(|(&pi, &(decoy, pos, _))| {
                    (pi as u32, pfms[reduced[pos]].feature_idx, decoy)
                })
                .collect();
            peptides.sort_unstable();
            // Candidate set = mass+RT-gated peptides. The gate must be LOOSE: uno's CNN
            // RT is ~2 min MAE (only ~half of real peptides within 1 min), so a tight
            // dev gate drops exactly the loose-RT real peptides the ML is meant to
            // recover (ms1searchpy promotes these via quality/consistency features). A
            // tight gate pre-decides the answer; the held-out-decoy FDR controls the
            // looseness. Keep a generous ceiling and let the classifier rank.
            let mut items: Vec<(usize, bool, usize, f64)> = pep_best
                .iter()
                .filter(|(_, &(_, _, dev))| dev <= 100.0)
                .map(|(&pi, &(d, pos, dev))| (pi, d, pos, dev))
                .collect();
            // `pep_best` is a std HashMap (randomized per-process iteration order),
            // so the collected row order is nondeterministic. LightGBM bagging picks
            // rows by index, so an unsorted order makes the trained model — and thus
            // the honest depth — vary run-to-run (observed 466..692 @1% for the SAME
            // config). Sort by peptide index to make the ML rescoring reproducible.
            items.sort_unstable_by_key(|&(pi, _, _, _)| pi);
            // Per-peptide group consistency stats over ALL the peptide's PFM matches
            // (ms1searchpy's id_count/charge_count/im_count + within-group |dev|
            // normalization): a real peptide recurs across charges / scans with a
            // CONSISTENT mass+RT; random coincidences do not. This multi-feature
            // evidence is exactly what the single-best-PFM deviation score discards.
            let group_stats: HashMap<usize, (f64, f64, f64, f64, f32, f32)> = {
                let med = |v: &mut Vec<f64>| -> f64 {
                    if v.is_empty() { return 0.0; }
                    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
                    v[v.len() / 2]
                };
                let mut m: HashMap<usize, (f64, f64, f64, f64, f32, f32)> = HashMap::new();
                for &(pi, _, _, _) in &items {
                    if m.contains_key(&pi) { continue; }
                    let (mut rtds, mut mzs) = (Vec::new(), Vec::new());
                    let mut chs: HashSet<u8> = HashSet::new();
                    let mut ims: HashSet<i64> = HashSet::new();
                    if let Some(list) = pep_pfms.get(&pi) {
                        for &(err, pos) in list {
                            mzs.push(((err - shift) / sigma).abs());
                            if let Some(rd) = rt_diff_at(&best_cal, pos) { rtds.push(rd.abs()); }
                            let f = &matrix.features[pfms[reduced[pos]].feature_idx as usize];
                            chs.insert(f.charge);
                            ims.insert((f.im.unwrap_or(0.0) * 100.0) as i64);
                        }
                    }
                    let (rtd_sum, mz_sum): (f64, f64) = (rtds.iter().sum(), mzs.iter().sum());
                    m.insert(pi, (med(&mut rtds), rtd_sum, med(&mut mzs), mz_sum, chs.len() as f32, ims.len() as f32));
                }
                m
            };
            // Feature set = uno's PFM features + amino-acid COMPOSITION (20) +
            // proline-context dipeptides (DP/KP/RP) + the per-peptide group
            // consistency features. This is what ms1searchpy's `-ml 1` LightGBM uses.
            const AAS: &[u8] = b"GASPVTCLINDQKEMHFRYW";
            let featurize_ml = |pos: usize, pi: usize| -> Vec<f32> {
                let mut v = featurize(&pfms[reduced[pos]], matrix, db, rt_diff_at(&best_cal, pos), None);
                let seq = db.peptides[pfms[reduced[pos]].peptide_idx as usize].sequence.as_bytes();
                let mut comp = [0f32; 20];
                let len = seq.len().max(1) as f32;
                for &b in seq {
                    if let Some(i) = AAS.iter().position(|&a| a == b) { comp[i] += 1.0; }
                }
                for c in comp.iter() { v.push(*c / len); }
                let (mut dp, mut kp, mut rp) = (0f32, 0f32, 0f32);
                for w in seq.windows(2) {
                    match w { b"DP" => dp += 1.0, b"KP" => kp += 1.0, b"RP" => rp += 1.0, _ => {} }
                }
                v.push(dp);
                v.push(kp);
                v.push(rp);
                let (rtd_med, rtd_sum, mz_med, mz_sum, n_ch, n_im) =
                    group_stats.get(&pi).copied().unwrap_or((0.0, 0.0, 0.0, 0.0, 1.0, 1.0));
                let id_count = pep_pfms.get(&pi).map_or(1, |l| l.len());
                let rtd_abs = rt_diff_at(&best_cal, pos).map(|x| x.abs()).unwrap_or(0.0);
                let mz_abs = ((pfms[reduced[pos]].mass_error_ppm - shift) / sigma).abs();
                v.push((1.0 + id_count as f64).ln() as f32);
                v.push((rtd_abs - rtd_med) as f32);
                v.push((rtd_abs / (rtd_sum + 0.01)) as f32);
                v.push((mz_abs - mz_med) as f32);
                v.push((mz_abs / (mz_sum + 0.01)) as f32);
                v.push((1.0 + n_ch as f64).ln() as f32);
                v.push((1.0 + n_im as f64).ln() as f32);
                // biosaur2 discriminators (isotope-pattern fit error + peak width):
                // the real-vs-random signal absent from n_isotopes/n_scans. 0 when the
                // feature file lacks them (inert there).
                let f = &matrix.features[pfms[reduced[pos]].feature_idx as usize];
                v.push(f.isoerror.unwrap_or(0.0) as f32);
                v.push(f.isoerror.map(|x| x.abs()).unwrap_or(0.0) as f32);
                v.push(f.isoerror2.unwrap_or(0.0) as f32);
                v.push(f.rt_width.unwrap_or(0.0) as f32);
                v
            };
            let feats: Vec<Vec<f32>> = items.iter().map(|&(pi, _, pos, _)| featurize_ml(pos, pi)).collect();
            // CV folds are PEPTIDE-disjoint (same peptide -> same fold): a peptide's
            // PFMs must never split across train/test or the model memorizes it.
            let folds: Vec<usize> = items.iter().map(|&(pi, _, _, _)| fold_of(pi, cfg.n_folds, cfg.effective_cv_seed())).collect();
            // SEMI-SUPERVISED labels (ms1searchpy's honesty mechanism): positives =
            // peptides of confident target proteins; negatives = dim target peptides;
            // DECOYS ARE EXCLUDED FROM TRAINING (label None) and held out purely to
            // estimate FDR. Training target-vs-decoy (the old scheme) leaks: real-but-
            // absent proteins look like targets in quality features, so the model
            // promotes them -> inflated FDR. Learning "what real IDs look like" from
            // confident-vs-dim targets keeps decoys a pristine null.
            let train_label: Vec<Option<f64>> = items
                .iter()
                .map(|&(pi, decoy, _, _)| {
                    if decoy {
                        None
                    } else if db.peptides[pi].proteins.iter().any(|pr| conf.contains(pr.as_ref())) {
                        Some(1.0)
                    } else {
                        Some(0.0)
                    }
                })
                .collect();
            // Nonlinear CV-gbdt (ms1searchpy uses LightGBM): the composition + group
            // consistency features interact nonlinearly, which a linear LDA cannot
            // capture (it collapsed to ~0 separation). Predicts ALL rows incl. the
            // held-out decoys for the FDR step.
            let ml = cv_train_predict_lgbm(&feats, &folds, &train_label, cfg);
            {
                let (mut ps, mut pn, mut ng, mut nn, mut ds, mut dn) =
                    (0.0f64, 0usize, 0.0f64, 0usize, 0.0f64, 0usize);
                for (k, &(_, _, _, _)) in items.iter().enumerate() {
                    match train_label[k] {
                        Some(l) if l >= 0.5 => { ps += ml[k]; pn += 1; }
                        Some(_) => { ng += ml[k]; nn += 1; }
                        None => { ds += ml[k]; dn += 1; }
                    }
                }
                // Target-peptide depth at 1% peptide-level decoy-FDR (the number that
                // matters: how many target peptides outrank the decoys).
                let mut rk: Vec<(f64, bool)> = items.iter().enumerate()
                    .map(|(k, &(_, decoy, _, _))| (ml[k], decoy)).collect();
                rk.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
                let (mut t, mut d, mut depth) = (0usize, 0usize, 0usize);
                for &(_, decoy) in &rk {
                    if decoy { d += 1; } else { t += 1; }
                    if t > 0 && (d as f64 / t as f64) <= 0.01 { depth = t; }
                }
                log::info!("ML diag: conf_mean={:.3} dim_mean={:.3} DECOY_mean={:.3} (conf={} dim={} decoy={}); target peptides @1% pep-FDR = {}",
                    ps / pn.max(1) as f64, ng / nn.max(1) as f64, ds / dn.max(1) as f64, pn, nn, dn, depth);
                if let Ok(path) = std::env::var("UNO_DUMP_ML") {
                    use std::io::Write;
                    if let Ok(f) = std::fs::File::create(&path) {
                        let mut w = std::io::BufWriter::new(f);
                        // label: 1=conf-target, 0=dim-target, -1=decoy ; ml = CV prediction; then feature vector
                        let _ = writeln!(w, "label\tdecoy\tml\tpeptide\tproteins\t{}", (0..feats[0].len()).map(|i| format!("f{i}")).collect::<Vec<_>>().join("\t"));
                        for (k, &(pi, decoy, _, _)) in items.iter().enumerate() {
                            let lab = match train_label[k] { Some(l) if l >= 0.5 => 1, Some(_) => 0, None => -1 };
                            let pep = &db.peptides[pi];
                            let prots: Vec<&str> = pep.proteins.iter().map(|p| p.as_ref()).collect();
                            let _ = writeln!(w, "{}\t{}\t{:.5}\t{}\t{}\t{}", lab, decoy as u8, ml[k], pep.sequence.as_ref(), prots.join(";"), feats[k].iter().map(|x| format!("{x:.4}")).collect::<Vec<_>>().join("\t"));
                        }
                        log::info!("dumped {} ML rows to {}", items.len(), path);
                    }
                }
            }
            // The ML feature matrix is only needed to fit/score/diagnose the
            // qpred ensemble above; the rest of this stage ranks by `ml` scores and
            // `pfms`/`reduced`, so free it now (one Vec<f32> per candidate peptide).
            drop(feats);
            // ms1searchpy `final_iteration` (main.py:140, 1965-2010): NOT a single
            // peptide-quality cutoff. Rank peptides by ML score (best first), then
            // sweep 10 CUMULATIVE inclusion thresholds (top 10%, 20%, ..., 100% of
            // ranked peptides == ms1searchpy's `mass_koef` 0..9 on `qpreds`), run the
            // full protein engine (feature competition + dynamic p) at each threshold,
            // and score every protein by the MEAN of its binomial score across the 10
            // passes (missing passes count as 0 -- consistency across thresholds is
            // rewarded, exactly like the mass/RT-deviation UNO_ENSEMBLE below).
            // DECISIVE (2026-06-29, uno-ml-protein-scoring memory): collapsing this to
            // a SINGLE cutoff (the old peptide-filter here) throws away most of the
            // lever -- an isolated Python test went 1073 (10-step ensemble) -> 632
            // (single cutoff) on the identical ranked peptide list. This is that missing
            // ensemble, ported to Rust.
            let pep_feat_ml: rustc_hash::FxHashMap<usize, Vec<usize>> = items
                .iter()
                .map(|&(pi, _, pos, _)| (pi, vec![pfms[reduced[pos]].feature_idx as usize]))
                .collect();
            // One row per peptide, ranked best-first. Default: ML score desc.
            // UNO_ENS_BY=dev ranks by mass+RT DEVIATION ascending instead -- this
            // mirrors ms1searchpy's `mass_diff <= mass_koef` sweep. The point: a
            // yeast/decoy peptide matching a human feature is a mass COINCIDENCE, so
            // its deviation is large; a real peptide matches tight. Sweeping the
            // deviation threshold therefore preferentially drops luck-matched false
            // peptides (the ML score does NOT -- DECOY_mean ~= dim_mean), which is
            // exactly the mid-tier (@5%) yeast contamination lever. UNO_ENS_BY=mix
            // blends: rank by ML but break ties / co-rank with deviation.
            let ens_by = std::env::var("UNO_ENS_BY").unwrap_or_default();
            let mut rows: Vec<(f64, bool, usize)> = match ens_by.as_str() {
                // small dev = best; negate so the shared desc sort keeps it first.
                "dev" => items.iter().map(|&(pi, decoy, _, dev)| (-dev, decoy, pi)).collect(),
                // blend: standardize ml (desc) and dev (asc) to ranks, sum. Computed
                // below after we have both orderings; here just carry ml, dev via pi.
                "mix" => {
                    let n = items.len().max(1) as f64;
                    // rank by ml desc
                    let mut by_ml: Vec<usize> = (0..items.len()).collect();
                    by_ml.sort_by(|&a, &b| ml[b].partial_cmp(&ml[a]).unwrap());
                    let mut ml_rank = vec![0f64; items.len()];
                    for (r, &k) in by_ml.iter().enumerate() { ml_rank[k] = r as f64 / n; }
                    // rank by dev asc
                    let mut by_dev: Vec<usize> = (0..items.len()).collect();
                    by_dev.sort_by(|&a, &b| items[a].3.partial_cmp(&items[b].3).unwrap());
                    let mut dev_rank = vec![0f64; items.len()];
                    for (r, &k) in by_dev.iter().enumerate() { dev_rank[k] = r as f64 / n; }
                    items.iter().enumerate()
                        .map(|(k, &(pi, decoy, _, _))| (-(ml_rank[k] + dev_rank[k]), decoy, pi))
                        .collect()
                }
                _ => items.iter().enumerate()
                    .map(|(k, &(pi, decoy, _, _))| (ml[k], decoy, pi))
                    .collect(),
            };
            rows.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap().then(a.2.cmp(&b.2)));
            let n = rows.len();
            // Ensemble cutoff sweep. The binomial background `p` is a global decoy
            // fraction, so LOOSE steps (keep 90-100%) admit many decoys -> p high ->
            // real peptide-rich proteins barely beat chance (scores collapse ~2.5 at
            // the FDR margin). Restricting to TIGHTER steps lowers the averaged p.
            // Env-tunable: UNO_ENS_STEPS resolution, UNO_ENS_LO..UNO_ENS_HI the
            // included step band.
            //
            // DEFAULT = steps=40, lo=1, hi=14: average only the 14 TIGHTEST of 40
            // inclusion bands (roughly the top 2.5-35% of ML-ranked peptides), dropping
            // the loose junk-heavy passes that inflate the random-match background p.
            // Validated 2026-07-20 on 2 large HeLa runs (run10 775->1218, deephela-240min
            // 1585->1752; neutral on a small run), dual-null calibrated (internal decoy
            // FDR and yeast-entrapment agree) and byte-deterministic. Takes uno from
            // behind ms1searchpy-DeepLC to a tie. Set `UNO_ENS_HI=10` with
            // `UNO_ENS_STEPS=10` for the old full uniform sweep (pre-2026-07-20).
            //
            // NOTE: this band was SELECTED by maximizing yeast-entrapment depth on
            // run10 and deephela, which are development datasets from the same family
            // as the published benchmark. The manuscript discloses it as a tuned
            // default; do not present it as an independent result.
            let n_steps: usize = std::env::var("UNO_ENS_STEPS").ok().and_then(|s| s.parse().ok()).unwrap_or(40);
            let lo: usize = std::env::var("UNO_ENS_LO").ok().and_then(|s| s.parse().ok()).unwrap_or(1).max(1);
            let hi: usize = std::env::var("UNO_ENS_HI").ok().and_then(|s| s.parse().ok()).unwrap_or(14).min(n_steps);
            let mut n_used = (hi + 1).saturating_sub(lo).max(1);
            // ms1searchpy's peptide filter (UNO_ML_DECOYCUT): instead of averaging a
            // 10-step cutoff ensemble, keep peptides in the top ML-quality quantile
            // bins up to where the per-bin DECOY fraction crosses `ml_correction`
            // (~0.5 for a balanced target/decoy db), then run ONE protein-binomial
            // pass. This is a loose, decoy-CALIBRATED cutoff (ms1 feeds ~19.4k peptides
            // vs uno's ensemble ~15.7k). `rows` is already ml-desc, one row/peptide.
            let decoy_cut = std::env::var("UNO_ML_DECOYCUT").is_ok();
            let ml_corr: f64 = std::env::var("UNO_ML_CORR").ok().and_then(|s| s.parse().ok()).unwrap_or(0.5);
            let steps: Vec<usize> = if decoy_cut {
                const NB: usize = 50;
                let mut qval_ok = 0usize;
                for b in 0..NB {
                    let (lo_i, hi_i) = (b * n / NB, ((b + 1) * n / NB).max(b * n / NB + 1).min(n));
                    if lo_i >= n { break; }
                    let bin = &rows[lo_i..hi_i];
                    let dfrac = bin.iter().filter(|(_, d, _)| *d).count() as f64 / bin.len().max(1) as f64;
                    if dfrac < ml_corr { qval_ok = b; } else { break; }
                }
                let keep = ((qval_ok + 1) * n / NB).max(1).min(n);
                n_used = 1;
                vec![keep]
            } else {
                (lo..=hi).map(|step| (n * step / n_steps).max(1).min(n)).collect()
            };
            // The inclusion steps are INDEPENDENT protein-inference passes -> run them
            // in parallel (each the dominant per-pass greedy cost). Results collect in
            // `steps` order, and the `acc` accumulation is then done SERIALLY in that
            // same order, so the float sums stay bit-identical to the old serial loop.
            let per_step: Vec<Vec<ProteinId>> = steps
                .par_iter()
                .map(|&keep| {
                    let mut observed: rustc_hash::FxHashMap<(bool, &str), rustc_hash::FxHashSet<usize>> = rustc_hash::FxHashMap::default();
                    for &(_, decoy, pi) in &rows[..keep] {
                        for prot in &db.peptides[pi].proteins {
                            observed.entry((decoy, prot.as_ref())).or_default().insert(pi);
                        }
                    }
                    crate::protein::score_proteins_feature_compete(
                        observed, &pep_feat_ml, &theo, decoy_theo_sum, 0.01)
                })
                .collect();
            let mut acc: HashMap<(bool, String), (f64, usize, usize)> = HashMap::new();
            for scored in &per_step {
                for p in scored {
                    let e = acc.entry((p.decoy, p.accession.clone())).or_insert((0.0, 0, 0));
                    e.0 += p.binomial_score;
                    e.1 = e.1.max(p.observed_peptides);
                    e.2 = p.theoretical_peptides;
                }
            }
            // Average per-pass binomial score; pair-compete target/decoy; (d+1)/t FDR
            // -- identical post-processing to the UNO_ENSEMBLE deviation path.
            let mut ensembled: Vec<ProteinId> = acc
                .into_iter()
                .map(|((decoy, accession), (s, obs, theo))| ProteinId {
                    accession,
                    decoy,
                    observed_peptides: obs,
                    theoretical_peptides: theo,
                    coverage: if theo > 0 { obs as f64 / theo as f64 } else { 0.0 },
                    binomial_score: s / n_used as f64,
                    mean_rt_diff: 0.0,
                    mean_mass_z: 0.0,
                    mean_pep_score: 0.0,
                    q_value: 1.0,
                })
                .collect();
            ensembled.sort_by(|a, b| b.binomial_score.partial_cmp(&a.binomial_score).unwrap()
                .then_with(|| a.accession.cmp(&b.accession)).then(a.decoy.cmp(&b.decoy)));
            let score_by: HashMap<(bool, String), f64> = ensembled
                .iter().map(|p| ((p.decoy, p.accession.clone()), p.binomial_score)).collect();
            ensembled.retain(|p| match score_by.get(&(!p.decoy, p.accession.clone())) {
                Some(&cs) => p.binomial_score > cs || (p.binomial_score == cs && !p.decoy),
                None => true,
            });
            let (mut t, mut d) = (0usize, 0usize);
            for p in ensembled.iter_mut() {
                if p.decoy { d += 1; } else { t += 1; }
                p.q_value = if t == 0 { 1.0 } else { (d as f64 + 1.0) / t as f64 };
            }
            let mut run = f64::INFINITY;
            for i in (0..ensembled.len()).rev() { run = run.min(ensembled[i].q_value); ensembled[i].q_value = run.min(1.0); }
            let ml_n = ensembled.iter().filter(|p| !p.decoy && p.q_value <= 0.01).count();
            log::info!("mass-anchor ML qpred-ensemble: {} peptides ranked, {}-step cumulative sweep; {} confident proteins (vs {} deviation)",
                n, n_used, ml_n, best_n);
            // VALIDATED (2026-07-07, human+yeast entrapment, PXD015669 QEHF1_09744_JB
            // native biosaur2 features): 223 real/4 yeast (deviation) -> 251 real/1
            // yeast (this ensemble) -- deeper AND more honest, not just deeper. The
            // `ml_n > best_n` check is itself a safe no-regression guard (only swaps
            // in a strictly deeper result), so no separate opt-in env var is needed.
            // If a feature file lacks the quality columns (isoerror/isoerror2 etc.)
            // the gbdt has no real signal over decoys (dim_mean ~= DECOY_mean) and
            // this correctly no-ops.
            return Some(FinalizedCandidate { proteins: ensembled, n_confident: ml_n, peptides });
        }
        None
    };

    // Finalize: default picks the candidate with the most proteins at q<=1% (first
    // wins ties, as before) and runs ML once; UNO_SMOOTH finalizes every candidate.
    let finalize = |r: &StageRun, fl: f64, mu: f64| -> FinalizedCandidate {
        set_rt_floor(fl);
        rt_sig_mult.set(mu);
        match run_ml(&r.2, &r.3, r.4, r.5, r.6, r.0) {
            Some(mut result) => {
                let take = if smooth { soft_depth(&result.proteins) > soft_depth(&r.2) } else { result.n_confident > r.0 };
                if !take {
                    // Preserve the deviation list and this same calibration's
                    // peptide evidence when its ML rescore does not improve it.
                    result.proteins = r.2.clone();
                    result.n_confident = r.0;
                }
                result
            }
            None => FinalizedCandidate { proteins: r.2.clone(), n_confident: r.0, peptides: Vec::new() },
        }
    };
    let chosen_idx = if multi_cand {
        let mut finals = Vec::new();
        for (i, (r, fl, mu)) in runs.iter().enumerate() {
            let f = finalize(r, *fl, *mu);
            log::info!("mass-anchor candidate {} (floor {:.1}, sigma x{:.1}): final {} proteins at q<=1%, soft depth {:.1}", i, fl, mu, f.n_confident, soft_depth(&f.proteins));
            finals.push(f);
        }
        let (best_i, result) = select_finalized_candidate(finals, smooth_level == 3, retain_zero_discovery_evidence);
        (best_i, Some(result))
    } else {
        let mut best_i = 0usize;
        for (i, (r, _, _)) in runs.iter().enumerate() { if r.0 > runs[best_i].0.0 { best_i = i; } }
        (best_i, None)
    };
    let (best, fl, mu) = runs.swap_remove(chosen_idx.0);
    set_rt_floor(fl);
    rt_sig_mult.set(mu);
    let (best_n, best_cut, best_proteins, best_cal, best_shift, best_sigma, best_gate, conf1, _rt_sig1, _best_soft) = best;
    let FinalizedCandidate { proteins: best_proteins, n_confident: best_n, peptides } = match chosen_idx.1 { Some(f) => f, None => finalize(&(best_n, best_cut, best_proteins.clone(), best_cal.clone(), best_shift, best_sigma, best_gate, conf1.clone(), _rt_sig1, _best_soft), fl, mu) };
    log::info!(
        "mass-anchor: {} confident target proteins (q<=1%, min-deviation any-feature, dev_cut={:.0})",
        best_n, best_cut
    );
    // Publish only the selected candidate's evidence. Evaluated losers must not
    // change protein grouping, peptide intensities, or LFQ feature assignments.
    MASS_ANCHOR_PEPTIDES.with(|c| *c.borrow_mut() = peptides);

    MASS_ANCHOR_PROTEINS.with(|c| *c.borrow_mut() = best_proteins);
    (cal, if anchor_pos.len() >= 50 { anchor_pos } else { clean }, shift, sigma)
}

thread_local! {
    /// Confident proteins from the mass-anchor mass+RT binomial, handed to `score`
    /// out-of-band (the seed function shares its signature with the legacy path).
    static MASS_ANCHOR_PROTEINS: std::cell::RefCell<Vec<crate::protein::ProteinId>> =
        std::cell::RefCell::new(Vec::new());
    /// RT+mass-gated best feature per peptide (the `pep_best` selection), handed
    /// out-of-band for the multi-sample quant rollup. `(peptide_idx, feature_idx,
    /// decoy)`. This is the feature uno actually trusts for the peptide (mass AND
    /// RT consistent), NOT the best-mass feature — the distinction matters for
    /// quant, where a mass-coincident contaminant would corrupt the intensity.
    static MASS_ANCHOR_PEPTIDES: std::cell::RefCell<Vec<(u32, u32, bool)>> =
        std::cell::RefCell::new(Vec::new());
}

/// One RT+mass-gated peptide→feature hit, exposed for the quant rollup.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PeptideHit {
    pub peptide_idx: u32,
    pub feature_idx: u32,
    pub decoy: bool,
}

/// Reduce all PFMs to the top-K best-mass-error targets and decoys per feature.
fn reduce_top_k(pfms: &[Pfm], n_features: usize, k: usize) -> Vec<usize> {
    // best |mass error| per feature, separately for targets and decoys
    let mut by_feature: Vec<(Vec<usize>, Vec<usize>)> = vec![(Vec::new(), Vec::new()); n_features];
    for (i, p) in pfms.iter().enumerate() {
        let bucket = if p.decoy {
            &mut by_feature[p.feature_idx as usize].1
        } else {
            &mut by_feature[p.feature_idx as usize].0
        };
        bucket.push(i);
    }
    let mut kept = Vec::new();
    for (targets, decoys) in by_feature.iter_mut() {
        for bucket in [targets, decoys] {
            bucket.sort_by(|&a, &b| {
                pfms[a]
                    .mass_error_ppm
                    .abs()
                    .partial_cmp(&pfms[b].mass_error_ppm.abs())
                    .unwrap()
            });
            bucket.truncate(k);
            kept.extend_from_slice(bucket);
        }
    }
    kept
}

/// Same peptide-disjoint CV / semi-supervised-label plumbing as
/// [`cv_train_predict`], but trains a real LightGBM classifier (via
/// `uno-lightgbm-sys`) instead of the pure-Rust `gbdt` crate. uno's own gbdt
/// crate has repeatedly failed to separate decoys from dim-real peptides on
/// this feature set (`ML diag: dim_mean ~= DECOY_mean` in the mass-anchor log,
/// see `ML_RESCORE_PORT.md`); LightGBM natively supports {0,1} binary labels
/// (no gbdt's {-1,+1} degenerate-gradient workaround) and is what
/// ms1searchpy itself uses, so this removes uno's own classifier as a
/// confound when judging whether the qpred-ensemble architecture is the
/// remaining lever.
/// Train on all folds != `test_fold` (pos=conf, neg=dim capped 3x, decoys
/// excluded) and write held-out predictions into `scores` for the test-fold rows.
/// Optionally rank-normalizes the test-fold predictions to `[0,1]` (ms1searchpy's
/// `rankdata(...)/n` — makes the per-fold scores comparable and robust to the raw
/// LightGBM output scale, which is what lets strong quality features help instead
/// of destabilize).
fn lgbm_fold(
    feats: &[Vec<f32>], folds: &[usize], train_label: &[Option<f64>],
    params: &str, iters: usize, test_fold: usize, rank_norm: bool, scores: &mut [f64],
) {
    let (mut pos_rows, mut neg_rows) = (Vec::new(), Vec::new());
    for p in 0..feats.len() {
        if folds[p] == test_fold { continue; }
        match train_label[p] {
            Some(l) if l >= 0.5 => pos_rows.push(p),
            Some(_) => neg_rows.push(p),
            None => {}
        }
    }
    if pos_rows.is_empty() || neg_rows.is_empty() { return; }
    let neg_cap = (pos_rows.len() * 3).max(20_000);
    if neg_rows.len() > neg_cap {
        neg_rows.sort_by_key(|&p| { let mut x = (p as u64).wrapping_mul(0x9e3779b97f4a7c15); x ^= x >> 29; x });
        neg_rows.truncate(neg_cap);
    }
    let train_feats: Vec<Vec<f32>> = pos_rows.iter().chain(neg_rows.iter()).map(|&p| feats[p].clone()).collect();
    let train_labels: Vec<f32> = pos_rows.iter().map(|_| 1.0f32).chain(neg_rows.iter().map(|_| 0.0f32)).collect();
    let model = match uno_lightgbm_sys::LgbmBooster::train(&train_feats, &train_labels, params, iters) {
        Ok(m) => m,
        Err(e) => { log::warn!("LightGBM fold {test_fold} training failed: {e}"); return; }
    };
    let test_idx: Vec<usize> = (0..feats.len()).filter(|&p| folds[p] == test_fold).collect();
    if test_idx.is_empty() { return; }
    let test_feats: Vec<Vec<f32>> = test_idx.iter().map(|&p| feats[p].clone()).collect();
    let preds = match model.predict(&test_feats) {
        Ok(p) => p,
        Err(e) => { log::warn!("LightGBM fold {test_fold} prediction failed: {e}"); return; }
    };
    if rank_norm {
        // ordinal rank / n within the fold (deterministic tie-break by index).
        let n = preds.len();
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| preds[a].partial_cmp(&preds[b]).unwrap().then(a.cmp(&b)));
        let mut rn = vec![0.0f64; n];
        for (rank, &k) in order.iter().enumerate() { rn[k] = (rank + 1) as f64 / n as f64; }
        for (k, &p) in test_idx.iter().enumerate() { scores[p] = rn[k]; }
    } else {
        for (k, &p) in test_idx.iter().enumerate() { scores[p] = preds[k] as f64; }
    }
}

/// Pure-Rust twin of [`lgbm_fold`] using `sequoia-boost` (an XGBoost reimplementation)
/// instead of the C LightGBM. IDENTICAL CV split, negative balancing (3x cap, same
/// hash order) and hyperparameters (binary:logistic, hist, depth/lr/iters, 0.8 row +
/// col subsample, seed 42, requested nthread=1). Only the tree library
/// differs. Gated by `UNO_ML_SEQUOIA`.
fn sequoia_fold(
    feats: &[Vec<f32>], folds: &[usize], train_label: &[Option<f64>],
    depth: u32, lr: f64, iters: usize, test_fold: usize, rank_norm: bool, seed: u64, scores: &mut [f64],
) {
    use sequoia_boost::prelude::*;
    let (mut pos_rows, mut neg_rows) = (Vec::new(), Vec::new());
    for p in 0..feats.len() {
        if folds[p] == test_fold { continue; }
        match train_label[p] {
            Some(l) if l >= 0.5 => pos_rows.push(p),
            Some(_) => neg_rows.push(p),
            None => {}
        }
    }
    if pos_rows.is_empty() || neg_rows.is_empty() { return; }
    let neg_cap = (pos_rows.len() * 3).max(20_000);
    if neg_rows.len() > neg_cap {
        neg_rows.sort_by_key(|&p| { let mut x = (p as u64).wrapping_mul(0x9e3779b97f4a7c15); x ^= x >> 29; x });
        neg_rows.truncate(neg_cap);
    }
    let ncol = feats.first().map(|f| f.len()).unwrap_or(0);
    if ncol == 0 { return; }
    let train_rows: Vec<usize> = pos_rows.iter().chain(neg_rows.iter()).copied().collect();
    let mut flat: Vec<f32> = Vec::with_capacity(train_rows.len() * ncol);
    for &p in &train_rows { flat.extend_from_slice(&feats[p]); }
    let labels: Vec<f32> = pos_rows.iter().map(|_| 1.0f32).chain(neg_rows.iter().map(|_| 0.0f32)).collect();
    let dtrain = match DMatrix::from_dense(&flat, train_rows.len(), ncol).and_then(|d| d.with_labels(&labels)) {
        Ok(d) => d, Err(e) => { log::warn!("sequoia fold {test_fold} DMatrix failed: {e}"); return; }
    };
    // Env-tunable hyperparameters (for the tuning sweep). depth/lr/iters come from the
    // caller (UNO_ML_DEPTH/LR/ITERS); these are sequoia/XGBoost-specific knobs.
    let ev = |k: &str, d: f64| std::env::var(k).ok().and_then(|s| s.parse().ok()).unwrap_or(d);
    let max_bin = ev("UNO_SEQ_BINS", 256.0) as usize;
    let lambda = ev("UNO_SEQ_LAMBDA", 1.0);
    let mcw = ev("UNO_SEQ_MCW", 80.0);
    let alpha = ev("UNO_SEQ_ALPHA", 0.0);
    let gamma = ev("UNO_SEQ_GAMMA", 0.0);
    let subsample = ev("UNO_SEQ_SUB", 0.8);
    let colsample = ev("UNO_SEQ_COL", 0.8);
    let params = match TrainingParams::builder()
        .objective("binary:logistic")
        .tree_method(TreeMethod::Hist)
        .max_depth(depth as usize)
        .eta(lr)
        .subsample(subsample)
        .colsample_bytree(colsample)
        .min_child_weight(mcw)
        .lambda(lambda)
        .alpha(alpha)
        .gamma(gamma)
        .max_bin(max_bin)
        .nthread(1) // Current Sequoia uses the ambient Rayon pool despite this setting.
        .seed(42 ^ crate::seed_mix(seed))
        .build() {
        Ok(p) => p, Err(e) => { log::warn!("sequoia fold {test_fold} params failed: {e}"); return; }
    };
    let model = match train(&params, &dtrain, iters) {
        Ok(m) => m, Err(e) => { log::warn!("sequoia fold {test_fold} training failed: {e}"); return; }
    };
    let test_idx: Vec<usize> = (0..feats.len()).filter(|&p| folds[p] == test_fold).collect();
    if test_idx.is_empty() { return; }
    let mut tflat: Vec<f32> = Vec::with_capacity(test_idx.len() * ncol);
    for &p in &test_idx { tflat.extend_from_slice(&feats[p]); }
    let dtest = match DMatrix::from_dense(&tflat, test_idx.len(), ncol) {
        Ok(d) => d, Err(e) => { log::warn!("sequoia fold {test_fold} test DMatrix failed: {e}"); return; }
    };
    let preds = match model.predict(&dtest) {
        Ok(p) => p, Err(e) => { log::warn!("sequoia fold {test_fold} prediction failed: {e}"); return; }
    };
    if rank_norm {
        let n = preds.len();
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| preds[a].partial_cmp(&preds[b]).unwrap().then(a.cmp(&b)));
        let mut rn = vec![0.0f64; n];
        for (rank, &k) in order.iter().enumerate() { rn[k] = (rank + 1) as f64 / n as f64; }
        for (k, &p) in test_idx.iter().enumerate() { scores[p] = rn[k]; }
    } else {
        for (k, &p) in test_idx.iter().enumerate() { scores[p] = preds[k] as f64; }
    }
}

/// AUC (Mann-Whitney) of conf(label 1) vs dim(label 0) over labeled rows; decoys
/// (label None) excluded. Used to score a hyperparameter set during the search.
fn lgbm_auc(scores: &[f64], train_label: &[Option<f64>]) -> f64 {
    let mut pairs: Vec<(f64, bool)> = Vec::new();
    for i in 0..scores.len() {
        match train_label[i] {
            Some(l) if l >= 0.5 => pairs.push((scores[i], true)),
            Some(_) => pairs.push((scores[i], false)),
            None => {}
        }
    }
    pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let n = pairs.len();
    let (mut rank_sum_pos, mut n_pos, mut n_neg) = (0.0f64, 0usize, 0usize);
    let mut i = 0;
    while i < n {
        let mut j = i;
        while j < n && pairs[j].0 == pairs[i].0 { j += 1; }
        let avg_rank = ((i + 1 + j) as f64) / 2.0;
        for k in i..j {
            if pairs[k].1 { rank_sum_pos += avg_rank; n_pos += 1; } else { n_neg += 1; }
        }
        i = j;
    }
    if n_pos == 0 || n_neg == 0 { return 0.5; }
    (rank_sum_pos - (n_pos as f64) * (n_pos as f64 + 1.0) / 2.0) / (n_pos as f64 * n_neg as f64)
}

/// Deterministic random LightGBM param string for search trial `t` (ms1searchpy's
/// grid: num_leaves, learning_rate, min_data, feature/bagging fractions, L1/L2).
fn lgbm_sample_params(t: u64, seed: u64) -> (String, usize) {
    let lgb_seed = 42u64.wrapping_add(seed);
    let mut s = t.wrapping_mul(0x9e3779b97f4a7c15).wrapping_add(0x2545f4914f6cdd1d) ^ crate::seed_mix(seed);
    let mut u = || { s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); (s >> 33) as f64 / (1u64 << 31) as f64 };
    let num_leaves = (15.0 + u() * 485.0) as u32;            // 15..500
    let lr = 0.01 + u() * 0.19;                              // 0.01..0.20
    let min_data = (20.0 + u() * 480.0) as usize;            // 20..500
    let ff = 0.4 + u() * 0.6;                                // 0.4..1.0
    let bf = 0.5 + u() * 0.5;                                // 0.5..1.0
    let l1 = u();                                            // 0..1
    let l2 = u();                                            // 0..1
    let iters = 100usize;
    let p = format!(
        "objective=binary num_leaves={num_leaves} learning_rate={lr:.4} \
         verbosity=-1 min_data_in_leaf={min_data} feature_fraction={ff:.3} bagging_fraction={bf:.3} bagging_freq=1 \
         lambda_l1={l1:.3} lambda_l2={l2:.3} \
         num_threads=1 deterministic=true force_row_wise=true seed={lgb_seed} \
         feature_fraction_seed={lgb_seed} bagging_seed={lgb_seed} data_random_seed={lgb_seed}"
    );
    (p, iters)
}

fn cv_train_predict_lgbm(feats: &[Vec<f32>], folds: &[usize], train_label: &[Option<f64>], cfg: &ScoreConfig) -> Vec<f64> {
    let n = feats.len();
    // ms1searchpy's per-dataset hyperparameter random search + peptide-disjoint CV
    // + rank-normalized predictions. UNO_ML_TUNE = number of search trials (default
    // 12; 0 = fixed-param legacy path). Each trial is scored by held-out conf-vs-dim
    // AUC on fold 0; the winner is refit over the full CV with rank-normalized preds.
    // Default OFF: the ms1-style search+CV+rank-norm (UNO_ML_TUNE=N) STABILIZES the
    // ML on strong quality features (aug/cosine-isotope 551->728, no longer collapses)
    // but does NOT raise the ceiling (plain koth 731->618) — uno's conf-vs-dim AUC caps
    // at ~0.649 regardless of hyperparameters, so the gap to ms1 (~1300) is FEATURE/
    // LABEL-limited, not tuning-limited. Kept env-gated for the feature-engineering work.
    let n_trials: usize = std::env::var("UNO_ML_TUNE").ok().and_then(|s| s.parse().ok()).unwrap_or(0);
    // depth 7 was the empirical fixed-path best on run10 (705@d6 -> 731@d7).
    let depth = std::env::var("UNO_ML_DEPTH").ok().and_then(|s| s.parse().ok()).unwrap_or(cfg.max_depth.max(7));
    if n_trials == 0 {
        let iters = std::env::var("UNO_ML_ITERS").ok().and_then(|s| s.parse().ok()).unwrap_or(300usize);
        let num_leaves = std::env::var("UNO_ML_LEAVES").ok().and_then(|s| s.parse().ok())
            .unwrap_or((1u32 << depth.min(8)).saturating_sub(1).max(2));
        let lr: f64 = std::env::var("UNO_ML_LR").ok().and_then(|s| s.parse().ok()).unwrap_or(cfg.shrinkage.min(0.05) as f64);
        let rank = std::env::var("UNO_ML_RANK").ok().map(|s| s != "0").unwrap_or(false);
        // DEFAULT ML backend = pure-Rust sequoia-boost (XGBoost clone), tuned mcw=80.
        // Validated 2026-07-21 (honest @1%): run10 sequoia 1218 = LightGBM 1218;
        // deephela-240min sequoia 1954 > LightGBM 1752 (+11.5%). Pure-Rust (drops the
        // vendored C lib_lightgbm.so) and byte-deterministic. Set UNO_ML_LIGHTGBM=1 to
        // fall back to the C LightGBM trainer.
        if !std::env::var("UNO_ML_LIGHTGBM").is_ok() {
            // Same byte-identical fold parallelization as the LightGBM path below: each
            // fold is seed-pinned and writes only
            // its own held-out rows, so training the folds concurrently and merging by
            // row-ownership reproduces the serial scatter exactly.
            let fold_scores: Vec<Vec<f64>> = (0..cfg.n_folds)
                .into_par_iter()
                .map(|f| {
                    let mut s = vec![0.0f64; n];
                    sequoia_fold(feats, folds, train_label, depth, lr, iters, f, rank, cfg.effective_classifier_seed(), &mut s);
                    s
                })
                .collect();
            let mut scores = vec![0.0f64; n];
            for i in 0..n { scores[i] = fold_scores[folds[i]][i]; }
            return scores;
        }
        let lgb_seed = 42u64.wrapping_add(cfg.effective_classifier_seed());
        let params = format!(
            "objective=binary num_leaves={num_leaves} max_depth={depth} learning_rate={lr} \
             verbosity=-1 min_data_in_leaf=20 feature_fraction=0.8 bagging_fraction=0.8 bagging_freq=1 \
             num_threads=1 deterministic=true force_row_wise=true seed={lgb_seed} \
             feature_fraction_seed={lgb_seed} bagging_seed={lgb_seed} data_random_seed={lgb_seed}");
        // The folds are INDEPENDENT (each trains on its complement and writes only its
        // own held-out test rows) and each LightGBM booster is single-threaded &
        // seed-pinned, so training them concurrently is byte-identical to the serial
        // loop -- only much faster (3 sequential single-thread trainings -> 3 parallel).
        // Each fold writes into its own buffer; the merge is by ownership (row i is a
        // test row of exactly fold `folds[i]`), so it reproduces the serial scatter
        // exactly (no float-equality heuristics).
        let fold_scores: Vec<Vec<f64>> = (0..cfg.n_folds)
            .into_par_iter()
            .map(|f| {
                let mut s = vec![0.0f64; n];
                lgbm_fold(feats, folds, train_label, &params, iters, f, rank, &mut s);
                s
            })
            .collect();
        let mut scores = vec![0.0f64; n];
        for i in 0..n {
            scores[i] = fold_scores[folds[i]][i];
        }
        return scores;
    }
    // Search: evaluate each trial on fold 0 (train folds 1..k, predict fold 0), pick
    // the highest held-out AUC.
    let mut best: Option<(f64, String, usize)> = None;
    for t in 0..n_trials as u64 {
        let (params, iters) = lgbm_sample_params(t, cfg.effective_classifier_seed());
        let mut sc = vec![0.0f64; n];
        lgbm_fold(feats, folds, train_label, &params, iters, 0, false, &mut sc);
        let auc = lgbm_auc(&sc, &train_label.iter().enumerate()
            .map(|(i, l)| if folds[i] == 0 { *l } else { None }).collect::<Vec<_>>());
        if best.as_ref().map_or(true, |(b, _, _)| auc > *b) { best = Some((auc, params, iters)); }
    }
    let (best_auc, params, iters) = best.unwrap_or_else(|| {
        let lgb_seed = 42u64.wrapping_add(cfg.effective_classifier_seed());
        (0.0, format!("objective=binary num_leaves=127 learning_rate=0.05 verbosity=-1 min_data_in_leaf=20 \
             feature_fraction=0.8 bagging_fraction=0.8 bagging_freq=1 num_threads=1 deterministic=true \
             force_row_wise=true seed={lgb_seed} feature_fraction_seed={lgb_seed} bagging_seed={lgb_seed} data_random_seed={lgb_seed}"), 100)
    });
    log::info!("ML hyperparam search: {} trials, best fold-0 AUC = {:.4}", n_trials, best_auc);
    // Final: full CV with the winning params + rank-normalized per-fold predictions.
    let mut scores = vec![0.0f64; n];
    for f in 0..cfg.n_folds { lgbm_fold(feats, folds, train_label, &params, iters, f, true, &mut scores); }
    scores
}

/// Robust spread (sigma estimate) of a centered-at-zero deviation sample, via
/// the median absolute value scaled to a Gaussian sigma. Resistant to the heavy
/// contamination in the positive seed.
fn robust_sigma(vals: impl Iterator<Item = f64>) -> f64 {
    let mut v: Vec<f64> = vals.filter(|x| x.is_finite()).map(|x| x.abs()).collect();
    if v.is_empty() {
        return f64::INFINITY;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mad = v[v.len() / 2];
    // |x| median -> sigma; guard against a degenerate zero spread.
    (mad / 0.6745).max(1e-9)
}

/// Sigma of the Gaussian CORE of a centered sample, via iterative sigma-clipping.
/// Real mass/RT errors form a tight Gaussian on a flat random-match background;
/// MAD is inflated by that background. Restricting to |x| <= 2.5σ and recomputing
/// the std converges to the true peak sigma (ms1searchpy fits a Gaussian for this).
fn clipped_sigma(vals: impl Iterator<Item = f64>) -> f64 {
    let mut v: Vec<f64> = vals.filter(|x| x.is_finite()).collect();
    // Sort before reducing. Every caller builds its input by iterating a HashMap,
    // so the values arrive in process-random order, and the two `sum::<f64>()`
    // below are not associative in floating point: the same sample gives a
    // slightly different sigma run to run. That sigma feeds the `dev <= cut`
    // acceptance, so the low bits reach the accepted protein set. Sorting also
    // makes the summation more accurate, since the magnitudes are then adjacent.
    v.sort_unstable_by(f64::total_cmp);
    if v.len() < 20 {
        return robust_sigma(v.into_iter());
    }
    let mut sig = robust_sigma(v.iter().copied());
    for _ in 0..8 {
        let core: Vec<f64> = v.iter().copied().filter(|&x| x.abs() <= 2.5 * sig).collect();
        if core.len() < 20 {
            break;
        }
        let mean = core.iter().sum::<f64>() / core.len() as f64;
        let var = core.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / core.len() as f64;
        let new_sig = var.sqrt().max(1e-9);
        if (new_sig - sig).abs() < 1e-4 {
            sig = new_sig;
            break;
        }
        sig = new_sig;
    }
    sig
}

/// ms1searchpy's `noisygaus` calibration, faithfully: fit `a*exp(-(x-x0)^2/2s^2)+b`
/// (Gaussian peak on a flat random-match baseline) to a histogram of the residuals,
/// return (x0=shift, |s|=sigma). ms1 uses scipy `curve_fit` (Levenberg-Marquardt);
/// here a grid over (x0, s) with the closed-form least-squares (a,b) per grid point
/// (the model is linear in a,b), minimizing SSE — deterministic, no LM needed.
/// `vals` = the raw residuals (e.g. rt_diff of ALL matched PFMs, ms1's input — NOT
/// the best-per-peptide order statistic, which collapses the peak).
fn noisygaus_fit(vals: &[f64], span: f64, nbins: usize) -> (f64, f64) {
    let v: Vec<f64> = vals.iter().copied().filter(|x| x.is_finite() && x.abs() < span).collect();
    if v.len() < 40 { return (0.0, robust_sigma(vals.iter().copied())); }
    let mut m = v.clone();
    m.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let med = m[m.len() / 2];
    let bw = 2.0 * span / nbins as f64;
    let mut hist = vec![0.0f64; nbins];
    for &x in &v {
        let b = (((x - med + span) / bw) as usize).min(nbins - 1);
        hist[b] += 1.0;
    }
    let centers: Vec<f64> = (0..nbins).map(|i| med - span + (i as f64 + 0.5) * bw).collect();
    // grid over x0 near the median, s in [0.1, span/2]; fit (a,b) by least squares.
    let n = nbins as f64;
    let mut best = (f64::INFINITY, med, 1.0f64);
    let mut x0 = med - 3.0;
    while x0 <= med + 3.0 {
        let mut s = 0.15;
        while s <= span * 0.4 {
            // g_i = exp(-(c_i - x0)^2 / 2 s^2); model h_i ~ a*g_i + b. Solve 2x2 normal eqs.
            let (mut sg, mut sgg, mut sh, mut sgh) = (0.0, 0.0, 0.0, 0.0);
            for i in 0..nbins {
                let g = (-((centers[i] - x0).powi(2)) / (2.0 * s * s)).exp();
                sg += g; sgg += g * g; sh += hist[i]; sgh += g * hist[i];
            }
            let det = sgg * n - sg * sg;
            if det.abs() > 1e-9 {
                let a = (sgh * n - sg * sh) / det;
                let b = (sgg * sh - sg * sgh) / det;
                if a > 0.0 {
                    let mut sse = 0.0;
                    for i in 0..nbins {
                        let g = (-((centers[i] - x0).powi(2)) / (2.0 * s * s)).exp();
                        let r = hist[i] - (a * g + b);
                        sse += r * r;
                    }
                    if sse < best.0 { best = (sse, x0, s); }
                }
            }
            s += 0.1;
        }
        x0 += bw;
    }
    (best.1, best.2.max(0.05))
}

/// Gaussian-CORE sigma via histogram half-width-at-half-maximum. Real RT/mass
/// residuals form a narrow Gaussian on a WIDE flat random-match background; MAD /
/// sigma-clipping stay inflated by that background (uno's `clipped_sigma` pins at
/// the 3.5 clamp), giving a too-loose gate that admits decoy matches. This finds
/// the narrow peak directly (ms1searchpy's `noisygaus` fit): histogram, locate the
/// mode, subtract the flat baseline, and measure the half-max half-width.
fn gaussian_core_sigma(vals: impl Iterator<Item = f64>, span: f64) -> f64 {
    let v: Vec<f64> = vals.filter(|x| x.is_finite()).collect();
    if v.len() < 40 {
        return robust_sigma(v.into_iter());
    }
    // center on the median, histogram over [-span, span].
    let mut m = v.clone();
    m.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let med = m[m.len() / 2];
    const NB: usize = 201; // odd -> a center bin
    let bw = 2.0 * span / NB as f64;
    let mut hist = vec![0u32; NB];
    for &x in &v {
        let d = x - med;
        if d.abs() >= span { continue; }
        let b = (((d + span) / bw) as usize).min(NB - 1);
        hist[b] += 1;
    }
    // light 3-bin smoothing.
    let sm: Vec<f64> = (0..NB).map(|i| {
        let lo = i.saturating_sub(1);
        let hi = (i + 1).min(NB - 1);
        (hist[lo] + hist[i] + hist[hi]) as f64 / (hi - lo + 1) as f64
    }).collect();
    let peak_idx = (0..NB).max_by(|&a, &b| sm[a].partial_cmp(&sm[b]).unwrap()).unwrap();
    let peak = sm[peak_idx];
    // flat background = median of the outer 40% of bins (tails).
    let mut tails: Vec<f64> = Vec::new();
    for i in 0..NB { if i < NB * 3 / 10 || i > NB * 7 / 10 { tails.push(sm[i]); } }
    tails.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let baseline = tails[tails.len() / 2];
    let half = (peak + baseline) / 2.0;
    if peak <= baseline + 1.0 {
        return clipped_sigma(v.into_iter());
    }
    // walk out from the peak to the half-max crossing on each side.
    let mut l = peak_idx;
    while l > 0 && sm[l] > half { l -= 1; }
    let mut r = peak_idx;
    while r < NB - 1 && sm[r] > half { r += 1; }
    let hwhm = ((r - l) as f64 / 2.0) * bw;
    (hwhm / 1.1774).max(0.05) // FWHM = 2.355 sigma -> HWHM = 1.1774 sigma
}

/// Run the full scoring + FDR pipeline. Returns one [`Identification`] per
/// feature that matched at least one peptide.
pub fn score(
    matrix: &FeatureMatrix,
    db: &PeptideDb,
    pfms: &[Pfm],
    rt_model: Option<&RtModel>,
    _im_model: Option<&crate::im::ImModel>,
    rt_library: Option<&std::collections::HashMap<String, f64>>,
    cfg: &ScoreConfig,
) -> ScoreResult {
    let n_features = matrix.features.len();
    let reduced = reduce_top_k(pfms, n_features, cfg.top_k);
    log::info!("reduced {} PFMs to {} (top-{} per feature)", pfms.len(), reduced.len(), cfg.top_k);

    // Penultimate embeddings per unique peptide_idx, captured in the SAME pass that
    // predicts abstract RT below (the on-data RT fine-tune re-embeds exactly these
    // peptides, so a separate embed pass would repeat ~1.1M CNN forwards). Filled
    // only when the fine-tune will consume them; reused at the rt_feat build.
    // Abstract RT per reduced PFM. Prefer a precomputed DeepLC library (sequence-
    // aware, generalizes — the discriminating signal); fall back to the distilled
    // model. Either is then calibrated on-data to this gradient.
    let abstract_rt: Vec<Option<f64>> = if let Some(lib) = rt_library {
        let mut hits = 0usize;
        let v: Vec<Option<f64>> = reduced
            .iter()
            .map(|&i| {
                let seq = db.peptides[pfms[i].peptide_idx as usize].sequence.as_ref();
                let r = lib.get(seq).copied();
                if r.is_some() {
                    hits += 1;
                }
                r
            })
            .collect();
        log::info!(
            "abstract RT from DeepLC library: {}/{} reduced PFMs covered",
            hits,
            reduced.len()
        );
        v
    } else {
        match rt_model {
            Some(m) => {
                let t = std::time::Instant::now();
                // Predict once per UNIQUE peptide_idx, then scatter back. The reduced
                // set repeats each peptide across many features (~9M PFMs, ~1M unique),
                // so deduplicating cuts RT inference ~9x.
                let mut uniq: Vec<usize> = reduced.iter().map(|&i| pfms[i].peptide_idx as usize).collect();
                uniq.sort_unstable();
                uniq.dedup();
                let seqs: Vec<&str> =
                    uniq.iter().map(|&pi| db.peptides[pi].sequence.as_ref()).collect();
                let upreds = crate::rt_cache::predict_many_cached(m, &seqs);
                let pred_of: std::collections::HashMap<usize, f64> =
                    uniq.iter().copied().zip(upreds).collect();
                log::info!(
                    "predicted abstract RT for {} unique peptides ({} PFMs) ({:.1}s)",
                    uniq.len(),
                    reduced.len(),
                    t.elapsed().as_secs_f64()
                );
                reduced.par_iter().map(|&i| pred_of.get(&(pfms[i].peptide_idx as usize)).copied()).collect()
            }
            None => vec![None; reduced.len()],
        }
    };

    // uno's deliverable is the mass-anchor protein set. It needs only the reduced
    // PFMs + the abstract RT computed above: RT calibration, mass recalibration,
    // the iterative protein-bootstrap, and the per-protein binomial all run inside
    // `mass_anchor_calibration`, which fills the MASS_ANCHOR_PROTEINS thread-local.
    // RT is required — without it there is no protein evidence and no result.
    let has_rt = abstract_rt.iter().any(|a| a.is_some());
    let run_anchor = |broad_rt_bootstrap, retain_zero| {
        // Early calibration exits publish nothing; each pass starts without
        // evidence from an earlier pass or search on this thread.
        MASS_ANCHOR_PROTEINS.with(|c| c.borrow_mut().clear());
        MASS_ANCHOR_PEPTIDES.with(|c| c.borrow_mut().clear());
        mass_anchor_calibration(&reduced, pfms, matrix, db, &abstract_rt, cfg, retain_zero, broad_rt_bootstrap);
        let proteins = MASS_ANCHOR_PROTEINS.with(|c| std::mem::take(&mut *c.borrow_mut()));
        let peptides = MASS_ANCHOR_PEPTIDES.with(|c| {
            std::mem::take(&mut *c.borrow_mut())
                .into_iter()
                .map(|(peptide_idx, feature_idx, decoy)| PeptideHit {
                    peptide_idx,
                    feature_idx,
                    decoy,
                })
                .collect()
        });
        ScoreResult { proteins, peptides }
    };
    if !has_rt {
        log::warn!("no RT model or RT library available; uno requires retention time -- no proteins produced");
    }
    run_recovery_bootstraps(has_rt, run_anchor)
}

#[cfg(test)]
mod determinism_tests {
    use super::*;

    fn candidate(n: usize, feature: u32, q: f64) -> FinalizedCandidate {
        FinalizedCandidate {
            proteins: (0..n).map(|i| ProteinId {
                accession: format!("candidate_{feature}_{i}"),
                decoy: false,
                observed_peptides: 1,
                theoretical_peptides: 1,
                coverage: 1.0,
                binomial_score: 1.0,
                mean_rt_diff: 0.0,
                mean_mass_z: 0.0,
                mean_pep_score: 0.0,
                q_value: q,
            }).collect(),
            n_confident: if q <= 0.01 { n } else { 0 },
            peptides: vec![(7, feature, false)],
        }
    }

    #[test]
    fn zero_discovery_cuts_retain_usable_training_evidence() {
        let weak = candidate(5, 10, 0.20).proteins;
        let useful = candidate(30, 20, 0.05).proteins;
        let mut decoys = candidate(100, 30, 0.005).proteins;
        for p in &mut decoys { p.decoy = true; }
        let cuts = [Vec::new(), weak, useful, decoys];
        let i = select_deviation_cut(cuts.iter().map(Vec::as_slice), true).unwrap();
        assert_eq!(i, 2);
        assert_eq!(cuts[i].iter().filter(|p| !p.decoy && p.q_value <= 0.01).count(), 0);
        assert_eq!(cuts[i].iter().filter(|p| !p.decoy && p.q_value <= 0.25).count(), 30);
        assert_eq!(select_deviation_cut(std::iter::empty(), true), None);
        // The first pass retains the exact legacy empty-cut behavior.
        assert_eq!(select_deviation_cut(cuts.iter().map(Vec::as_slice), false), None);
    }

    #[test]
    fn discovery_cut_selection_preserves_positive_counts_and_first_ties() {
        let training_only = candidate(1000, 10, 0.20).proteins;
        let discoveries = candidate(100, 20, 0.01).proteins;
        let mut tied = discoveries.clone();
        tied.extend(candidate(500, 30, 0.05).proteins);
        let cuts = [training_only, discoveries, tied];
        assert_eq!(select_deviation_cut(cuts.iter().map(Vec::as_slice), true), Some(1));
        assert_eq!(select_deviation_cut(cuts.iter().map(Vec::as_slice), false), Some(1));
        let empty = [Vec::new(), Vec::new()];
        assert_eq!(select_deviation_cut(empty.iter().map(Vec::as_slice), true), Some(0));
    }

    /// A losing candidate evaluated later must never replace the winner's
    /// peptide-to-feature mapping used for grouping and quantification.
    #[test]
    fn finalized_candidate_retains_its_own_peptide_evidence() {
        for count_only in [false, true] {
            let (i, winner) = select_finalized_candidate(
                vec![candidate(3, 10, 0.005), candidate(1, 20, 0.005)], count_only, true,
            );
            assert_eq!(i, 0);
            assert_eq!(winner.n_confident, 3);
            assert_eq!(winner.proteins[0].accession, "candidate_10_0");
            assert_eq!(winner.peptides, vec![(7, 10, false)]);
        }
    }

    #[test]
    fn training_only_evidence_cannot_displace_confident_identifications() {
        // A broad training list has greater soft depth but no reportable hits.
        // This reproduces the low-input seed-3 regression in candidate selection.
        let training = candidate(300, 10, 0.02);
        let confident = candidate(117, 20, 0.01);
        assert!(soft_depth(&training.proteins) > soft_depth(&confident.proteins));
        let (i, winner) = select_finalized_candidate(vec![training, confident], false, true);
        assert_eq!(i, 1);
        assert_eq!(winner.n_confident, 117);
        assert_eq!(winner.peptides, vec![(7, 20, false)]);
        // When every candidate has zero discoveries, retain useful ranked
        // evidence; no confidence threshold is weakened to manufacture hits.
        let (i, winner) = select_finalized_candidate(
            vec![candidate(100, 10, 0.02), candidate(300, 20, 0.02)], false, true,
        );
        assert_eq!(i, 1);
        assert_eq!(winner.n_confident, 0);
    }

    #[test]
    fn original_pass_preserves_legacy_soft_depth_selection() {
        // Experimental modes can produce training-only legacy rankings. The
        // first pass must preserve their original selection semantics too.
        let (i, winner) = select_finalized_candidate(
            vec![candidate(300, 10, 0.02), candidate(117, 20, 0.01)], false, false,
        );
        assert_eq!(i, 0);
        assert_eq!(winner.n_confident, 0);
    }

    fn recovery_result(n: usize, feature: u32, q: f64) -> ScoreResult {
        let c = candidate(n, feature, q);
        ScoreResult {
            proteins: c.proteins,
            peptides: c.peptides.into_iter().map(|(peptide_idx, feature_idx, decoy)|
                PeptideHit { peptide_idx, feature_idx, decoy }).collect(),
        }
    }

    #[test]
    fn successful_original_search_never_invokes_recovery() {
        let result = recover_zero_discoveries(recovery_result(117, 10, 0.01), || {
            panic!("a successful original ranking must not compete against rescue candidates")
        });
        assert_eq!(result.proteins.len(), 117);
        assert_eq!(result.peptides[0].feature_idx, 10);
    }

    #[test]
    fn failed_recovery_preserves_original_lower_confidence_evidence() {
        let mut calls = 0;
        let result = recover_zero_discoveries(recovery_result(30, 10, 0.03), || {
            calls += 1;
            recovery_result(300, 20, 0.02)
        });
        assert_eq!(calls, 1);
        assert_eq!(result.proteins.len(), 30);
        assert_eq!(result.peptides[0].feature_idx, 10);
    }

    #[test]
    fn recovery_requires_actual_target_discoveries_and_keeps_their_map() {
        let mut original = recovery_result(30, 10, 0.005);
        for p in &mut original.proteins { p.decoy = true; }
        let result = recover_zero_discoveries(original, || recovery_result(100, 20, 0.01));
        assert_eq!(result.proteins.len(), 100);
        assert_eq!(result.peptides[0].feature_idx, 20);
        let mut decoys = recovery_result(300, 30, 0.005);
        for p in &mut decoys.proteins { p.decoy = true; }
        let result = recover_zero_discoveries(recovery_result(30, 10, 0.03), || decoys);
        assert_eq!(result.peptides[0].feature_idx, 10);
    }

    #[test]
    fn rt_recovery_preserves_successful_initial_v3_lazily() {
        let mut calls = Vec::new();
        let result = run_recovery_bootstraps(true, |broad, retain| {
            calls.push((broad, retain));
            assert_eq!((broad, retain), (false, false), "successful V3 must not run another pass");
            recovery_result(117, 10, 0.01)
        });
        assert_eq!(calls, [(false, false)]);
        assert_eq!(result.proteins.len(), 117);
        assert_eq!(result.peptides[0].feature_idx, 10);
    }

    #[test]
    fn rt_recovery_waits_for_complete_v3_retention_policy() {
        let mut calls = Vec::new();
        let result = run_recovery_bootstraps(true, |broad, retain| {
            assert!(!broad, "V3 retention success must prevent broad RT from running");
            calls.push((broad, retain));
            if retain { recovery_result(50, 20, 0.01) } else { recovery_result(30, 10, 0.02) }
        });
        assert_eq!(calls, [(false, false), (false, true)]);
        assert_eq!(result.proteins.len(), 50);
        assert_eq!(result.peptides[0].feature_idx, 20);
    }

    #[test]
    fn rt_recovery_requires_actual_targets_and_never_runs_broad_retention() {
        let mut calls = Vec::new();
        let result = run_recovery_bootstraps(true, |broad, retain| {
            calls.push((broad, retain));
            match (broad, retain) {
                (false, false) => {
                    let mut result = recovery_result(300, 10, 0.001);
                    for protein in &mut result.proteins { protein.decoy = true; }
                    result
                },
                (false, true) => recovery_result(400, 20, 0.02),
                (true, false) => recovery_result(25, 40, 0.01),
                (true, true) => panic!("bounded recovery must never run broad retention"),
            }
        });
        assert_eq!(calls, [(false, false), (false, true), (true, false)]);
        assert_eq!(result.proteins.len(), 25);
        assert_eq!(result.peptides[0].feature_idx, 40);
    }

    #[test]
    fn failed_broad_recovery_preserves_original_v3_proteins_and_peptides() {
        let mut calls = Vec::new();
        let result = run_recovery_bootstraps(true, |broad, retain| {
            calls.push((broad, retain));
            let feature = 10 + 10 * calls.len() as u32;
            let mut result = recovery_result(if broad { 1000 } else { 30 }, feature, if broad { 0.001 } else { 0.03 });
            if broad { for protein in &mut result.proteins { protein.decoy = true; } }
            result
        });
        assert_eq!(calls, [(false, false), (false, true), (true, false)]);
        assert_eq!(result.proteins.len(), 30);
        assert!(result.proteins.iter().all(|p| !p.decoy && p.q_value == 0.03));
        assert_eq!(result.proteins[0].accession, "candidate_20_0");
        assert_eq!(result.peptides[0].feature_idx, 20);
    }

    #[test]
    fn rt_recovery_without_predictions_never_runs_any_bootstrap() {
        let result = run_recovery_bootstraps(false, |_, _| panic!("no RT must remain an early empty result"));
        assert!(result.proteins.is_empty());
        assert!(result.peptides.is_empty());
    }

    #[test]
    fn bounded_recovery_policy_repeats_identically_despite_fresh_timing() {
        let mut expected = None;
        for _ in 0..16 {
            let mut calls = Vec::new();
            let result = run_recovery_bootstraps(true, |broad, retain| {
                calls.push((broad, retain));
                assert!(!(broad && retain));
                recovery_result(if broad { 11 } else { 30 }, if broad { 40 } else { 10 }, if broad { 0.01 } else { 0.02 })
            });
            let signature = (calls,
                result.proteins.iter().map(|p| (p.accession.clone(), p.q_value.to_bits(), p.decoy)).collect::<Vec<_>>(),
                result.peptides.iter().map(|p| (p.peptide_idx, p.feature_idx, p.decoy)).collect::<Vec<_>>());
            if let Some(ref prior) = expected { assert_eq!(&signature, prior); }
            else { expected = Some(signature); }
        }
    }

    /// Ties retain the first candidate's entire result, while a later winner
    /// carries its own map; scoring mode changes must not change ownership.
    #[test]
    fn finalized_candidate_ties_and_later_winners_keep_matching_evidence() {
        for count_only in [false, true] {
            let (i, tie) = select_finalized_candidate(
                vec![candidate(2, 10, 0.005), candidate(2, 20, 0.005)], count_only, true,
            );
            assert_eq!(i, 0);
            assert_eq!(tie.peptides, vec![(7, 10, false)]);
            let (i, later) = select_finalized_candidate(
                vec![candidate(1, 10, 0.005), candidate(3, 20, 0.005)], count_only, true,
            );
            assert_eq!(i, 1);
            assert_eq!(later.proteins[0].accession, "candidate_20_0");
            assert_eq!(later.peptides, vec![(7, 20, false)]);
        }
    }

    /// `clipped_sigma` must not depend on the order it sees its sample.
    ///
    /// Every caller builds the input by iterating a `HashMap`, so the order is
    /// process-random, and the function reduces with `sum::<f64>()`, which is not
    /// associative in floating point. The resulting sigma feeds the `dev <= cut`
    /// acceptance, so a low-bit difference reaches the accepted protein set. This
    /// is the regression test for that fix: shuffled input, bit-identical output.
    #[test]
    fn clipped_sigma_is_order_invariant() {
        // A tight Gaussian-ish core plus a flat background, i.e. the shape the
        // function exists to handle, and large enough to pass its n >= 20 gate.
        let base: Vec<f64> = (0..400)
            .map(|i| {
                let x = (i % 41) as f64 - 20.0;
                if i % 17 == 0 { x * 3.7 } else { x * 0.11 }
            })
            .collect();

        let want = clipped_sigma(base.iter().copied());

        // Deterministic shuffles (no RNG: the test must not need one).
        for stride in [7usize, 13, 29, 97] {
            let mut shuffled = Vec::with_capacity(base.len());
            let mut i = 0usize;
            for _ in 0..base.len() {
                shuffled.push(base[i]);
                i = (i + stride) % base.len();
            }
            shuffled.reverse();
            let got = clipped_sigma(shuffled.into_iter());
            assert_eq!(
                want.to_bits(),
                got.to_bits(),
                "clipped_sigma changed with input order (stride {stride}): {want} vs {got}"
            );
        }
    }

    /// The same property for the degenerate branch (n < 20 -> `robust_sigma`).
    #[test]
    fn clipped_sigma_is_order_invariant_below_the_gate() {
        let base: Vec<f64> = vec![0.4, -1.2, 0.03, 2.9, -0.7, 0.11, -0.02, 1.4, 0.8, -0.35];
        let want = clipped_sigma(base.iter().copied());
        let mut rev = base.clone();
        rev.reverse();
        assert_eq!(want.to_bits(), clipped_sigma(rev.into_iter()).to_bits());
    }
}

#[cfg(test)]
mod stage_seed_tests {
    use super::*;

    #[test]
    fn old_configs_preserve_every_scoring_stream() {
        for seed in [0, 1, 42, u64::MAX] {
            let cfg: ScoreConfig = serde_json::from_value(serde_json::json!({ "seed": seed })).unwrap();
            assert_eq!(cfg.effective_rt_seed(), seed);
            assert_eq!(cfg.effective_cv_seed(), seed);
            assert_eq!(cfg.effective_classifier_seed(), seed);
            let old_folds: Vec<_> = (0..100).map(|i| fold_of(i, 3, seed)).collect();
            let new_folds: Vec<_> = (0..100).map(|i| fold_of(i, 3, cfg.effective_cv_seed())).collect();
            assert_eq!(old_folds, new_folds);
            assert_eq!(lgbm_sample_params(3, seed), lgbm_sample_params(3, cfg.effective_classifier_seed()));
        }
    }

    #[test]
    fn explicit_zero_pins_stages_across_parent_seeds() {
        let mut cfg: ScoreConfig = serde_json::from_value(serde_json::json!({
            "seed": 9, "rt_seed": 0, "cv_seed": 0, "classifier_seed": 0,
        })).unwrap();
        for seed in [0, 1, 9, u64::MAX] {
            cfg.seed = seed;
            assert_eq!(cfg.effective_rt_seed(), 0);
            assert_eq!(cfg.effective_cv_seed(), 0);
            assert_eq!(cfg.effective_classifier_seed(), 0);
        }
    }

    #[test]
    fn stage_overrides_are_independent() {
        let mut cfg = ScoreConfig { seed: 7, ..ScoreConfig::default() };
        cfg.rt_seed = Some(11);
        assert_eq!((cfg.effective_rt_seed(), cfg.effective_cv_seed(), cfg.effective_classifier_seed()), (11, 7, 7));
        cfg.cv_seed = Some(13);
        assert_eq!((cfg.effective_rt_seed(), cfg.effective_cv_seed(), cfg.effective_classifier_seed()), (11, 13, 7));
        cfg.classifier_seed = Some(17);
        assert_eq!((cfg.effective_rt_seed(), cfg.effective_cv_seed(), cfg.effective_classifier_seed()), (11, 13, 17));
        let original_folds: Vec<_> = (0..100).map(|i| fold_of(i, 3, 7)).collect();
        let changed_folds: Vec<_> = (0..100).map(|i| fold_of(i, 3, cfg.effective_cv_seed())).collect();
        assert_ne!(original_folds, changed_folds);
        assert_ne!(lgbm_sample_params(3, 7), lgbm_sample_params(3, cfg.effective_classifier_seed()));
    }
}
