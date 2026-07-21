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
use crate::input::FeatureMatrix;
use crate::matching::Pfm;
use crate::protein::ProteinId;
use crate::rt::RtModel;
use gbdt::config::Config;
use gbdt::decision_tree::{Data, DataVec};
use gbdt::gradient_boost::GBDT;
use rayon::prelude::*;

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
    /// ms1searchpy's `final_iteration`), not classifier quality. On koth-ALIGNED
    /// matrices (no isoerror/isoerror2 — lost in cross-run alignment) this
    /// correctly no-ops (guarded by a strict better-than-baseline check). OFF by
    /// default because it re-runs the full protein engine 10x (real runtime cost);
    /// turn on for single-run biosaur2-native searches where depth matters more
    /// than latency.
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

/// Best match for a single peptide across all features it matched (not subject
/// to per-feature competition). Used as protein-level evidence: a peptide is
/// "observed" if its best RT-consistent match is good enough.
#[derive(Debug, Clone)]
pub struct PeptideMatch {
    pub peptide_idx: usize,
    pub decoy: bool,
    pub abs_rt_diff: f64,
    pub abs_im_diff: f64,
    /// Recalibrated mass deviation in sigma units (|mass_err - shift| / mass_sigma).
    /// Used to require tight mass at the protein-evidence level, which (with RT)
    /// crushes the random-match background — ms1searchpy's mass-recalibration step.
    pub abs_mass_z: f64,
    pub score: f64,
    /// Shape quality of the feature backing this match (n_isotopes; 0 when the
    /// input lacks shape descriptors). Real peptides sit on well-formed features;
    /// random matches scatter onto junk features, so gating protein evidence on
    /// this drops the random-match background asymmetrically.
    pub feat_quality: f64,
}

/// Result of scoring: per-feature identifications plus per-peptide best matches.
#[derive(Debug, Clone)]
pub struct ScoreResult {
    pub ids: Vec<Identification>,
    pub peptide_matches: Vec<PeptideMatch>,
    /// Confident proteins from the mass-anchor mass+RT binomial identifier
    /// (ms1searchpy-style). `Some` only when `mass_anchor` is enabled; when set,
    /// these ARE the protein result and the caller should use them directly
    /// instead of `infer_proteins`.
    pub proteins: Option<Vec<crate::protein::ProteinId>>,
    /// Deviations of EVERY reduced PFM (not just the per-peptide representative):
    /// `(peptide_idx, decoy, abs_mass_z, abs_rt_diff)`. Lets the caller reproduce
    /// ms1searchpy's "matched peptides for top proteins" count, which credits a
    /// peptide if ANY of its features falls in the mass+RT window (not only its
    /// tightest-mass match). `abs_rt_diff` is in minutes; `INFINITY` when RT was
    /// unavailable for that PFM.
    pub pfm_dev: Vec<(usize, bool, f64, f64)>,
    /// Winning on-data RT calibration (abstract RT -> observed RT) from the mass
    /// anchor, for the hill second pass to reuse. `None` if unavailable.
    pub rt_cal: Option<Calibration>,
    /// RT residual sigma (min) of the confident feature-pass hits.
    pub rt_sigma: f64,
    /// Recalibrated mass shift and Gaussian-core sigma (ppm) from the mass anchor.
    pub mass_shift: f64,
    pub mass_sigma: f64,
}

/// One scored, FDR-controlled identification: the best peptide for a feature.
#[derive(Debug, Clone)]
pub struct Identification {
    pub feature_idx: usize,
    pub peptide_idx: usize,
    pub decoy: bool,
    pub score: f64,
    pub q_value: f64,
    pub mass_error_ppm: f64,
    pub isotope_offset: i8,
    pub rt_diff: f64,
    pub im_diff: f64,
}

const N_FEATURES: usize = 14;

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
    let f = &matrix.features[pfm.feature_idx];
    let p = &db.peptides[pfm.peptide_idx];
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
        // absent (koth aligned matrices), so this is inert there.
        (1.0 + f.n_isotopes.unwrap_or(0.0)).ln() as f32,
        (1.0 + f.n_scans.unwrap_or(0.0)).ln() as f32,
        // koth single-run PFM-quality columns (dropped by cross-run alignment).
        // These are the exact signals ms1searchpy's LightGBM PFM filter uses to
        // drop the random-match rate p to ~0.135 (its @5% depth). 0.0 when absent
        // (aligned matrices), so inert there — active only on single-run inputs.
        f.cosine_score.unwrap_or(0.0) as f32,
        f.isotope_score.unwrap_or(0.0) as f32,
        f.ppm_error.unwrap_or(0.0) as f32,
        f.neutron_offset.unwrap_or(0.0) as f32,
    ]
}

fn make_config(cfg: &ScoreConfig, n_features: usize) -> Config {
    let mut c = Config::new();
    c.set_feature_size(n_features);
    // ms1searchpy-style boosting: more rounds + stochastic sub-sampling (bagging +
    // feature fraction) so the PFM rescorer actually extracts the weak quality signal
    // (a 100-round full-sample fit barely separates; cf. lightgbm reaches ~1020 proteins
    // vs the old 367). Env overrides for tuning experiments.
    let iters = std::env::var("UNO_ML_ITERS").ok().and_then(|s| s.parse().ok()).unwrap_or(300);
    let depth = std::env::var("UNO_ML_DEPTH").ok().and_then(|s| s.parse().ok()).unwrap_or(cfg.max_depth);
    c.set_max_depth(depth);
    c.set_iterations(iters);
    c.set_shrinkage(cfg.shrinkage.min(0.05));
    c.set_data_sample_ratio(0.8);
    c.set_feature_sample_ratio(0.8);
    c.set_loss("LogLikelyhood");
    c.set_training_optimization_level(2);
    c
}

/// Stable fold assignment for a feature index.
fn fold_of(feature_idx: usize, n_folds: usize) -> usize {
    let mut x = feature_idx as u64;
    x ^= x >> 30;
    x = x.wrapping_mul(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    (x % n_folds as u64) as usize
}

/// Symmetric pseudo-random tiebreak in [0, 1) for a (feature, peptide) pair, so
/// target/decoy ties in the per-feature competition resolve without bias.
fn tiebreak(feature_idx: usize, peptide_idx: usize) -> f64 {
    let mut x = (feature_idx as u64)
        .wrapping_mul(0x9e3779b97f4a7c15)
        ^ (peptide_idx as u64).wrapping_mul(0xc2b2ae3d27d4eb4f);
    x ^= x >> 31;
    x = x.wrapping_mul(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    (x >> 11) as f64 / (1u64 << 53) as f64
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

/// Reusable mass-anchor-style protein rollup. Given per-peptide evidence
/// `(peptide_idx, decoy, dev, feature_id)` — where `dev = mass_z² + rt_z²` and
/// `feature_id` is the peak the peptide sits on — scan a deviation-cut grid, run
/// the feature-competition dynamic-background binomial at each cut, and return
/// the protein set at the cut that maximizes target IDs at `fdr`.
///
/// This is the exact scorer the feature pass uses internally (626 proteins),
/// exposed so a second pass (e.g. hills) can be scored the SAME way instead of
/// through the weaker `infer_proteins` (177). The caller supplies one match per
/// peptide (its best cluster); feature competition over `feature_id` resolves
/// peptides that claim the same peak.
pub fn score_by_deviation(
    db: &PeptideDb,
    pep_matches: &[(usize, bool, f64, usize)],
    fdr: f64,
) -> Vec<ProteinId> {
    let mut theo: rustc_hash::FxHashMap<(bool, &str), u32> = rustc_hash::FxHashMap::default();
    for pep in &db.peptides {
        for prot in &pep.proteins {
            *theo.entry((pep.decoy, prot.as_ref())).or_insert(0) += 1;
        }
    }
    let decoy_theo_sum: u64 = theo.iter().filter(|((d, _), _)| *d).map(|(_, &c)| c as u64).sum();
    let cut_grid: Vec<f64> = {
        let mut g: Vec<f64> = Vec::new();
        let mut c = 1.0;
        while c <= 12.0 {
            g.push(c);
            c += 0.5;
        }
        g.extend_from_slice(&[13.0, 14.0, 16.0, 18.0, 20.0, 25.0]);
        g
    };
    let mut best_proteins: Vec<ProteinId> = Vec::new();
    let mut best_n = 0usize;
    for &cut in &cut_grid {
        let mut observed: rustc_hash::FxHashMap<(bool, &str), rustc_hash::FxHashSet<usize>> = rustc_hash::FxHashMap::default();
        let mut pep_feat: rustc_hash::FxHashMap<usize, Vec<usize>> = rustc_hash::FxHashMap::default();
        for &(pi, decoy, dev, feat) in pep_matches {
            if dev <= cut {
                pep_feat.entry(pi).or_default().push(feat);
                for prot in &db.peptides[pi].proteins {
                    observed.entry((decoy, prot.as_ref())).or_default().insert(pi);
                }
            }
        }
        let scored = crate::protein::score_proteins_feature_compete(
            observed,
            &pep_feat,
            &theo,
            decoy_theo_sum,
            fdr,
        );
        let nt = scored.iter().filter(|p| !p.decoy && p.q_value <= fdr).count();
        if nt > best_n {
            best_n = nt;
            best_proteins = scored;
        }
    }
    best_proteins
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
        if has_iso && matrix.features[pfms[i].feature_idx].n_isotopes.unwrap_or(0.0) < MIN_ISOTOPES {
            continue;
        }
        let e = best.entry(pfms[i].peptide_idx).or_insert((f64::INFINITY, pos));
        if pfms[i].mass_error_ppm.abs() < e.0.abs() {
            *e = (pfms[i].mass_error_ppm, pos);
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
        pep_pfms.entry(pfms[i].peptide_idx).or_default().push((pfms[i].mass_error_ppm, pos));
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
        let tol = if round == 0 { 10.0 } else { (3.0 * sigma).min(10.0) };
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
        for (&pi, &(err, pos)) in &best {
            if (err - shift).abs() > tol {
                continue;
            }
            let pep = &db.peptides[pi];
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
    if anchor_pos.is_empty() {
        return (Calibration::default(), Vec::new(), shift, sigma);
    }
    // 6. RANSAC RT calibration on the (abstract, observed) pairs of the anchors.
    let pairs: Vec<(f64, f64)> = anchor_pos
        .iter()
        .filter_map(|&pos| abstract_rt[pos].map(|a| (a, matrix.features[pfms[reduced[pos]].feature_idx].rt)))
        .collect();
    if pairs.len() < 50 {
        return (Calibration::default(), anchor_pos, shift, sigma);
    }
    let (mut cal, mut inl) = Calibration::fit_ransac(&pairs, 2.5, 4000);
    for thresh in [1.0_f64, 0.5, 0.3] {
        let sub: Vec<(f64, f64)> = inl.iter().map(|&k| pairs[k]).collect();
        if sub.len() < 50 { break; }
        let (c, si) = Calibration::fit_ransac(&sub, thresh, 4000);
        if si.len() < 50 { break; }
        inl = si.iter().map(|&j| inl[j]).collect();
        cal = c;
    }
    // keep only the RANSAC-inlier anchors (cleanest)
    let clean: Vec<usize> = inl.iter().map(|&k| anchor_pos[k.min(anchor_pos.len()-1)]).collect();
    let mut resid: Vec<f64> = inl.iter().map(|&k| (pairs[k].1 - cal.apply(pairs[k].0)).abs()).collect();
    let med = if resid.is_empty() { f64::NAN } else { median(&mut resid) };
    log::info!(
        "mass-anchor RT calibration: {} RANSAC inliers, median |residual| = {:.2} min",
        clean.len(), med
    );

    // === Stage 2 (SINGLE deterministic pass; ms1searchpy-faithful) ===
    // The RT calibration `cal` is now FIXED (fit on the clean mass-stage anchors).
    // The previous iterative recalibration loop oscillated (calibration <-> confident
    // set feedback) and made the result DB-composition-dependent. Here we estimate
    // RT_sigma ONCE and do ONE any-feature mass+RT binomial pass -> deterministic
    // and stable. A peptide is "seen" if ANY of its features is consistent in both
    // mass and RT (the low-isotope-recovery that gives the RT multiplier).
    let rt_diff_at = |cal: &Calibration, pos: usize| -> Option<f64> {
        abstract_rt[pos].map(|a| matrix.features[pfms[reduced[pos]].feature_idx].rt - cal.apply(a))
    };
    let mass_gate = (3.0_f64 * sigma).min(10.0);

    // RT accuracy from confident-protein peptides' best-RT (real-feature) residuals,
    // Gaussian core, clamped to weight RT ~its true accuracy in the deviation.
    let est_rt_sig = |cal: &Calibration, confident: &HashSet<String>, m_shift: f64, m_gate: f64| -> f64 {
        let mut rr: Vec<f64> = Vec::new();       // best rt_diff per peptide (uno legacy)
        let mut all_rd: Vec<f64> = Vec::new();   // ALL confident-peptide PFM rt_diffs (ms1's noisygaus input)
        for (&pi, pfm_list) in &pep_pfms {
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
        clipped_sigma(rr.into_iter()).clamp(0.3, 3.5)
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
                    if cs_gate > 0.0 && matrix.features[pfms[reduced[pos]].feature_idx].combined_score < cs_gate { continue; }
                    if let Some(pr) = pred {
                        let rd = matrix.features[pfms[reduced[pos]].feature_idx].rt - pr;
                        let dev = mz * mz + (rd / rt_sig).powi(2);
                        if best.map_or(true, |(b, _)| dev < b) { best = Some((dev, pos)); }
                        if want_feats { feats.push((dev, pfms[reduced[pos]].feature_idx)); }
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
                pep_feat.insert(pi, vec![pfms[reduced[pos]].feature_idx]);
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
                                let feat = pfms[reduced[pos]].feature_idx;
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
                            if f.is_empty() { f.push(pfms[reduced[pos]].feature_idx); }
                            f
                        } else {
                            vec![pfms[reduced[pos]].feature_idx]
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
        let mut best_idx: Option<usize> = None;
        for i in 0..per_cut.len() {
            let (cut, nt, tgt, dec) = (per_cut[i].0, per_cut[i].1, per_cut[i].2, per_cut[i].3);
            if nt > best_n {
                best_n = nt;
                best_cut = cut;
                best_idx = Some(i);
                if log_pepcount {
                    log::info!("run_deviation cut={:.1}: {} proteins; distinct peptides target={} decoy={} (decoy/target={:.2})",
                        cut, nt, tgt, dec, dec as f64 / tgt.max(1) as f64);
                }
            }
        }
        if let Some(i) = best_idx {
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
    let rt_sig1 = est_rt_sig(&cal, &mass_conf_accs, shift, mass_gate);
    let (n1, cut1, prot1, conf1) = run_deviation(&cal, rt_sig1, shift, sigma, mass_gate);

    // ITERATIVE multi-stage RT calibration (ms1searchpy drives the random-match rate
    // p 0.68->0.30->0.28->0.12 across stages; uno used to do ONE recalibration and
    // stop at p~0.23-0.44). Each stage recalibrates RT on the CURRENT confident set,
    // tightens rt_sig, and re-searches; the confident set purifies each round so the
    // gate tightens and p falls. Deterministic: cal_pairs are total-ordered before
    // the fit (pep_pfms is a std HashMap; fit_robust's stable x-sort would otherwise
    // let ties wobble the fit and swing depth). Keep the best-scoring stage.
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
                    cp.push((a, matrix.features[pfms[reduced[pos]].feature_idx].rt));
                }
            }
        }
        cp.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap().then(a.1.partial_cmp(&b.1).unwrap()));
        cp
    };
    let n_stages: usize = std::env::var("UNO_RT_STAGES").ok().and_then(|s| s.parse().ok()).unwrap_or(4);
    let (mut best_n, mut best_cut, mut best_proteins, mut best_cal) =
        (n1, cut1, prot1.clone(), cal.clone());
    // Track the mass calibration (shift/sigma/gate) alongside the RT calibration so
    // the winning stage's TIGHTENED mass carries into the ML/hill/store steps.
    let (mut best_shift, mut best_sigma, mut best_gate) = (shift, sigma, mass_gate);
    let (mut cur_cal, mut cur_conf) = (cal.clone(), conf1.clone());
    let (mut cur_shift, mut cur_sigma, mut cur_gate) = (shift, sigma, mass_gate);
    for stage in 2..=n_stages.max(2) {
        // Retighten MASS on the current confident set FIRST (ms1's big p-drop lever),
        // then recalibrate RT under the tightened mass, then re-search.
        let (s_shift, s_sigma, s_gate) = if stage_mass {
            recal_mass(&cur_conf, &cur_cal, cur_shift, cur_sigma, cur_gate)
        } else { (shift, sigma, mass_gate) };
        let cal_pairs = build_cal_pairs(&cur_conf, &cur_cal, s_shift, s_gate);
        if cal_pairs.len() < 100 { break; }
        let cal_s = Calibration::fit_robust(&cal_pairs, 60, 200);
        let rt_sig_s = est_rt_sig(&cal_s, &cur_conf, s_shift, s_gate);
        let (ns, cuts, prots, confs) = run_deviation(&cal_s, rt_sig_s, s_shift, s_sigma, s_gate);
        log::info!("mass-anchor stage {}: {} confident (rt_sig {:.2}, mass {:.2}+/-{:.2}ppm gate {:.2}, {} conf-peptides)",
            stage, ns, rt_sig_s, s_shift, s_sigma, s_gate, cal_pairs.len());
        if ns > best_n {
            best_n = ns; best_cut = cuts; best_proteins = prots; best_cal = cal_s.clone();
            best_shift = s_shift; best_sigma = s_sigma; best_gate = s_gate;
        }
        // stop when the confident set stops growing (converged) or shrinks.
        cur_shift = s_shift; cur_sigma = s_sigma; cur_gate = s_gate;
        if confs.len() <= cur_conf.len() { break; }
        cur_cal = cal_s;
        cur_conf = confs;
    }
    let mut best_n = best_n;
    // An on-data (MLP-embedding + ridge) RT refit was tried here (re-running the
    // SAME protein-multiplicity binomial with a tighter per-run RT model, kept
    // only if it beat stage1/2's result): on run10 at 400k-feature scale it
    // reproduces the exact iter-0 divergence seen elsewhere in this file (rt_sig
    // pinned at the 3.5 clamp ceiling, correctly auto-rejected, 957 unchanged);
    // on a small/clean input it gave a mixed result (nominal +12, honest @1% -22%,
    // @2%/@5% +30-40%) -- not a validated win anywhere. Matches this file's other
    // documented on-data-RT-at-this-stage attempt (712->333, auto-rejected).
    log::info!(
        "mass-anchor: {} confident target proteins (q<=1%, min-deviation any-feature, dev_cut={:.0})",
        best_n, best_cut
    );

    // Hand the winning RT calibration (+ RT sigma, mass shift/sigma) to the hill
    // second pass so it predicts peptide RT in the SAME on-data-calibrated space.
    {
        let rt_sig_final = est_rt_sig(&best_cal, &conf1, best_shift, best_gate).max(0.02);
        MASS_ANCHOR_CAL.with(|c| {
            *c.borrow_mut() = Some((best_cal.clone(), rt_sig_final, best_shift, best_sigma));
        });
    }

    // === Stage 3: CV-gbdt PFM rescoring (ms1searchpy `-ml 1`) ===
    // The min-deviation match uses only mass_z + rt_z. ms1searchpy trains a CV
    // gradient-boosted classifier on the FULL PFM feature set (mass, RT, n_isotopes,
    // n_scans, charge, intensity, length, missed cleavages) — the nonlinear model
    // separates real PFMs from random where 2 features can't, recovering ~2x depth at
    // honest FDR. Positives = mono PFMs of deviation-confident proteins; negatives =
    // decoy PFMs; fold-by-feature CV prevents leakage. A peptide is "observed" if its
    // best PFM passes a scanned ML q-value; greedy parsimony then scores proteins.
    // Kept only if it beats the deviation result. (uno's prior LINEAR LDA failed; the
    // GBM nonlinearity + per-PFM q-thresholding is what ms1searchpy actually uses.)
    {
        // ms1searchpy's semi-supervised positive set is peptides of 25%-FDR proteins
        // (top_25_targets), NOT 1% -- a loose gate gives far more positive training
        // signal. (Validated: 25%-label ML peptide-filter -> 899->1073 honest @1%.)
        let conf: HashSet<String> = best_proteins
            .iter()
            .filter(|p| !p.decoy && p.q_value <= 0.25)
            .map(|p| p.accession.clone())
            .collect();
        if cfg.ml_rescore && conf.len() >= 20 {
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
                    if cs_gate_ml > 0.0 && matrix.features[pfms[reduced[pos]].feature_idx].combined_score < cs_gate_ml {
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
                            let f = &matrix.features[pfms[reduced[pos]].feature_idx];
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
                let seq = db.peptides[pfms[reduced[pos]].peptide_idx].sequence.as_bytes();
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
                // matrix lacks them (inert on koth aligned inputs).
                let f = &matrix.features[pfms[reduced[pos]].feature_idx];
                v.push(f.isoerror.unwrap_or(0.0) as f32);
                v.push(f.isoerror.map(|x| x.abs()).unwrap_or(0.0) as f32);
                v.push(f.isoerror2.unwrap_or(0.0) as f32);
                v.push(f.rt_width.unwrap_or(0.0) as f32);
                v
            };
            let feats: Vec<Vec<f32>> = items.iter().map(|&(pi, _, pos, _)| featurize_ml(pos, pi)).collect();
            // CV folds are PEPTIDE-disjoint (same peptide -> same fold): a peptide's
            // PFMs must never split across train/test or the model memorizes it.
            let folds: Vec<usize> = items.iter().map(|&(pi, _, _, _)| fold_of(pi, cfg.n_folds)).collect();
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
                .map(|&(pi, _, pos, _)| (pi, vec![pfms[reduced[pos]].feature_idx]))
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
            // DEFAULT = "top-30%" (lo=1, hi=3): average only the 3 TIGHTEST of the 10
            // inclusion bands (top 10/20/30% of ML-ranked peptides), dropping the loose
            // junk-heavy passes (top 40-100%) that inflate the random-match background p.
            // Validated 2026-07-20 on 2 large HeLa runs (run10 775->1218, deephela-240min
            // 1585->1752; neutral on a small run), dual-null calibrated (internal decoy
            // FDR and yeast-entrapment agree) and byte-deterministic. Takes uno from
            // behind ms1searchpy-DeepLC to a tie. Set `UNO_ENS_HI=10` for the old full
            // 1..=10 uniform sweep (the pre-2026-07-20 default).
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
            // On koth-ALIGNED matrices (isoerror/isoerror2 absent, single-run
            // biosaur2 shape columns lost in alignment) the gbdt has no real signal
            // over decoys (dim_mean ~= DECOY_mean) and this correctly no-ops.
            if ml_n > best_n {
                best_n = ml_n;
                best_proteins = ensembled;
            }
        }
    }
    let _ = best_n;

    MASS_ANCHOR_PROTEINS.with(|c| *c.borrow_mut() = best_proteins);
    (cal, if anchor_pos.len() >= 50 { anchor_pos } else { clean }, shift, sigma)
}

thread_local! {
    /// Confident proteins from the mass-anchor mass+RT binomial, handed to `score`
    /// out-of-band (the seed function shares its signature with the legacy path).
    static MASS_ANCHOR_PROTEINS: std::cell::RefCell<Vec<crate::protein::ProteinId>> =
        std::cell::RefCell::new(Vec::new());
    /// Winning RT calibration + RT sigma + mass shift/sigma from the mass-anchor,
    /// handed out-of-band so the hill second pass can reuse the SAME on-data
    /// calibration (no refitting). `(cal, rt_sigma_min, mass_shift_ppm, mass_sigma_ppm)`.
    static MASS_ANCHOR_CAL: std::cell::RefCell<Option<(Calibration, f64, f64, f64)>> =
        std::cell::RefCell::new(None);
}

/// Seed RT calibration WITHOUT the classifier, from features that are very
/// likely real: a single strong target match with tight mass error and a
/// high-quality, cross-run-reproducible koth feature. Reversed decoys are mass
/// twins of their targets, so mass alone cannot bootstrap; this heuristic
/// enriches for true peptides, then MAD-trimming locks onto the RT diagonal.
fn seed_calibration(
    reduced: &[usize],
    pfms: &[Pfm],
    matrix: &FeatureMatrix,
    abstract_rt: &[Option<f64>],
    min_runs_frac: f64,
) -> (Calibration, Vec<usize>) {
    let n_samples = matrix.n_samples().max(1);
    // Floor of 2 enforces cross-run reproducibility on aligned consensus matrices;
    // for a single run every feature has n_contributing_runs == 1, so drop to 1.
    let floor = if n_samples <= 1 { 1 } else { 2 };
    let min_runs = ((n_samples as f64 * min_runs_frac) as u32).max(floor);

    // Seed candidates: best (tightest mass) mono target per high-quality,
    // reproducible feature. This set is only ~30-50% correct (the tightest mass
    // match is often a wrong peptide), but correct matches form a sharp RT
    // diagonal that a robust mode fit extracts.
    let mut best: std::collections::HashMap<usize, (f64, f64, usize)> =
        std::collections::HashMap::new();
    for (pos, &ri) in reduced.iter().enumerate() {
        let p = &pfms[ri];
        if p.decoy || p.isotope_offset != 0 || p.mass_error_ppm.abs() > 5.0 {
            continue;
        }
        let f = &matrix.features[p.feature_idx];
        if f.combined_score < 0.5 || f.n_contributing_runs < min_runs {
            continue;
        }
        let a = match abstract_rt[pos] {
            Some(a) => a,
            None => continue,
        };
        let err = p.mass_error_ppm.abs();
        let entry = best.entry(p.feature_idx).or_insert((err, a, pos));
        if err < entry.0 {
            *entry = (err, a, pos);
        }
    }
    log::info!("RT calibration seed: {} candidate features (best-mass per feature)", best.len());

    // The seed is heavily contaminated (~5-10% correct): the tightest-mass match
    // is usually a wrong peptide at this database density. RANSAC isolates the
    // true-match diagonal; its inliers seed both calibration and the initial
    // positive training set.
    // Deterministic order is essential: HashMap iteration order is process-random,
    // and RANSAC samples point-pairs positionally, so an unordered seed makes the
    // calibration (and everything downstream) differ run-to-run. Sort by feature.
    let mut seed_items: Vec<(usize, f64, f64, usize)> = best
        .iter()
        .map(|(&fi, &(_, a, pos))| (fi, a, matrix.features[pfms[reduced[pos]].feature_idx].rt, pos))
        .collect();
    seed_items.sort_unstable_by_key(|&(fi, ..)| fi);
    let items: Vec<(f64, f64, usize)> =
        seed_items.into_iter().map(|(_, a, o, pos)| (a, o, pos)).collect();
    let pairs: Vec<(f64, f64)> = items.iter().map(|&(a, o, _)| (a, o)).collect();
    let (mut cal, mut inlier_idx) = Calibration::fit_ransac(&pairs, 2.5, 4000);
    // Progressive tightening: re-RANSAC on the current inliers at ever tighter
    // thresholds. With accurate (DeepLC) predictions the true matches form a
    // sharp diagonal that this locks onto (residual -> ~0.1-0.3 min); a weak
    // predictor just keeps the loose set. This is the key to calibrating a good
    // RT model on a contaminated seed without a separate clean search stage.
    for thresh in [1.0_f64, 0.5, 0.3] {
        let sub: Vec<(f64, f64)> = inlier_idx.iter().map(|&k| pairs[k]).collect();
        if sub.len() < 100 {
            break;
        }
        let (c, sub_in) = Calibration::fit_ransac(&sub, thresh, 4000);
        if sub_in.len() < 100 {
            break;
        }
        inlier_idx = sub_in.iter().map(|&j| inlier_idx[j]).collect();
        cal = c;
    }
    let seed_positions: Vec<usize> = inlier_idx.iter().map(|&k| items[k].2).collect();

    let mut resid: Vec<f64> = inlier_idx
        .iter()
        .map(|&k| (items[k].1 - cal.apply(items[k].0)).abs())
        .collect();
    let med_resid = if resid.is_empty() { f64::NAN } else { median(&mut resid) };
    log::info!(
        "RT calibration (RANSAC): {} inliers ({:.1}% of seed), median |residual| = {:.2} min",
        seed_positions.len(),
        100.0 * seed_positions.len() as f64 / items.len().max(1) as f64,
        med_resid
    );
    (cal, seed_positions)
}

/// Reduce all PFMs to the top-K best-mass-error targets and decoys per feature.
fn reduce_top_k(pfms: &[Pfm], n_features: usize, k: usize) -> Vec<usize> {
    // best |mass error| per feature, separately for targets and decoys
    let mut by_feature: Vec<(Vec<usize>, Vec<usize>)> = vec![(Vec::new(), Vec::new()); n_features];
    for (i, p) in pfms.iter().enumerate() {
        let bucket = if p.decoy {
            &mut by_feature[p.feature_idx].1
        } else {
            &mut by_feature[p.feature_idx].0
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

/// Semi-supervised CV training (Percolator/mokapot style): train each fold-out
/// model on the *confident positives* (`train_label == Some(1.0)`) versus
/// decoys (`Some(0.0)`); positions with `None` are scored but not trained on.
/// Splitting by feature keeps a feature's PFMs out of the model that scores it.
fn cv_train_predict(
    reduced: &[usize],
    feats: &[Vec<f32>],
    folds: &[usize],
    train_label: &[Option<f64>],
    cfg: &ScoreConfig,
) -> Vec<f64> {
    let mut scores = vec![0.0f64; reduced.len()];

    for test_fold in 0..cfg.n_folds {
        // Collect training rows from other folds, capped per class.
        let mut pos_rows = Vec::new();
        let mut neg_rows = Vec::new();
        for p in 0..reduced.len() {
            if folds[p] == test_fold {
                continue;
            }
            match train_label[p] {
                Some(l) if l >= 0.5 => pos_rows.push(p),
                Some(_) => neg_rows.push(p),
                None => {}
            }
        }
        if pos_rows.is_empty() || neg_rows.is_empty() {
            continue;
        }
        // Balance: cap negatives at 3x positives (hash-ordered subsample).
        let neg_cap = (pos_rows.len() * 3).max(20_000);
        if neg_rows.len() > neg_cap {
            neg_rows.sort_by_key(|&p| {
                let mut x = (reduced[p] as u64).wrapping_mul(0x9e3779b97f4a7c15);
                x ^= x >> 29;
                x
            });
            neg_rows.truncate(neg_cap);
        }
        let mut train: DataVec = pos_rows
            .iter()
            .map(|&p| Data::new_training_data(feats[p].clone(), 1.0, 1.0, None))
            .chain(
                // gbdt's LogLikelyhood loss expects labels in {-1, +1}; a 0.0 negative
                // label degenerates the gradient (the model collapses to a constant).
                neg_rows
                    .iter()
                    .map(|&p| Data::new_training_data(feats[p].clone(), 1.0, -1.0, None)),
            )
            .collect();
        let nf = feats.first().map(|f| f.len()).unwrap_or(N_FEATURES);
        let mut model = GBDT::new(&make_config(cfg, nf));
        model.fit(&mut train);

        let test_idx: Vec<usize> = (0..reduced.len()).filter(|&p| folds[p] == test_fold).collect();
        let test_data: DataVec = test_idx
            .iter()
            .map(|&p| Data::new_test_data(feats[p].clone(), None))
            .collect();
        if test_data.is_empty() {
            continue;
        }
        let preds = model.predict(&test_data);
        for (k, &p) in test_idx.iter().enumerate() {
            scores[p] = preds[k] as f64;
        }
    }
    scores
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
/// col subsample, seed 42, single-thread for determinism). Only the tree library
/// differs. Gated by `UNO_ML_SEQUOIA`.
fn sequoia_fold(
    feats: &[Vec<f32>], folds: &[usize], train_label: &[Option<f64>],
    depth: u32, lr: f64, iters: usize, test_fold: usize, rank_norm: bool, scores: &mut [f64],
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
        .nthread(1)
        .seed(42)
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
fn lgbm_sample_params(t: u64) -> (String, usize) {
    let mut s = t.wrapping_mul(0x9e3779b97f4a7c15).wrapping_add(0x2545f4914f6cdd1d);
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
         num_threads=1 deterministic=true force_row_wise=true seed=42 \
         feature_fraction_seed=42 bagging_seed=42 data_random_seed=42"
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
            // fold is single-threaded (sequoia nthread=1) + seed-pinned and writes only
            // its own held-out rows, so training the folds concurrently and merging by
            // row-ownership reproduces the serial scatter exactly.
            let fold_scores: Vec<Vec<f64>> = (0..cfg.n_folds)
                .into_par_iter()
                .map(|f| {
                    let mut s = vec![0.0f64; n];
                    sequoia_fold(feats, folds, train_label, depth, lr, iters, f, rank, &mut s);
                    s
                })
                .collect();
            let mut scores = vec![0.0f64; n];
            for i in 0..n { scores[i] = fold_scores[folds[i]][i]; }
            return scores;
        }
        let params = format!(
            "objective=binary num_leaves={num_leaves} max_depth={depth} learning_rate={lr} \
             verbosity=-1 min_data_in_leaf=20 feature_fraction=0.8 bagging_fraction=0.8 bagging_freq=1 \
             num_threads=1 deterministic=true force_row_wise=true seed=42 \
             feature_fraction_seed=42 bagging_seed=42 data_random_seed=42");
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
        let (params, iters) = lgbm_sample_params(t);
        let mut sc = vec![0.0f64; n];
        lgbm_fold(feats, folds, train_label, &params, iters, 0, false, &mut sc);
        let auc = lgbm_auc(&sc, &train_label.iter().enumerate()
            .map(|(i, l)| if folds[i] == 0 { *l } else { None }).collect::<Vec<_>>());
        if best.as_ref().map_or(true, |(b, _, _)| auc > *b) { best = Some((auc, params, iters)); }
    }
    let (best_auc, params, iters) = best.unwrap_or_else(|| {
        (0.0, format!("objective=binary num_leaves=127 learning_rate=0.05 verbosity=-1 min_data_in_leaf=20 \
             feature_fraction=0.8 bagging_fraction=0.8 bagging_freq=1 num_threads=1 deterministic=true \
             force_row_wise=true seed=42 feature_fraction_seed=42 bagging_seed=42 data_random_seed=42"), 100)
    });
    log::info!("ML hyperparam search: {} trials, best fold-0 AUC = {:.4}", n_trials, best_auc);
    // Final: full CV with the winning params + rank-normalized per-fold predictions.
    let mut scores = vec![0.0f64; n];
    for f in 0..cfg.n_folds { lgbm_fold(feats, folds, train_label, &params, iters, f, true, &mut scores); }
    scores
}

/// Public Percolator-style rescorer: semi-supervised cross-validated Fisher LDA
/// over arbitrary per-row features. `train_label[i]` = `Some(1.0)` for a
/// confident positive, `Some(0.0)` for a decoy, `None` for unlabeled; `pep_idx`
/// gives each row's peptide so folds split BY PEPTIDE (no leakage). Returns the
/// discriminant per row (higher = more target-like). Used by the hill second
/// pass to combine mass/RT/envelope/coelution/isotope evidence into one score.
pub fn cv_lda_rescore(
    feats: &[Vec<f32>],
    pep_idx: &[usize],
    train_label: &[Option<f64>],
    n_folds: usize,
) -> Vec<f64> {
    let reduced: Vec<usize> = (0..feats.len()).collect();
    let folds: Vec<usize> = pep_idx.iter().map(|&pi| fold_of(pi, n_folds)).collect();
    cv_lda_predict(&reduced, feats, &folds, train_label, n_folds)
}

/// Semi-supervised CV LDA discriminant (closed-form, deterministic): a fast,
/// FDR-honest alternative to the gbdt path. Same fold-by-feature, confident-
/// positive-vs-decoy training as [`cv_train_predict`], but fits a Fisher linear
/// discriminant per held-out fold — standardize features on the training rows,
/// compute class means, a shrunk pooled covariance, then w = Sigma^-1
/// (mu_pos - mu_neg) — and projects the held-out fold. No RNG -> deterministic.
fn cv_lda_predict(
    reduced: &[usize],
    feats: &[Vec<f32>],
    folds: &[usize],
    train_label: &[Option<f64>],
    n_folds: usize,
) -> Vec<f64> {
    let d = feats.first().map(|f| f.len()).unwrap_or(0);
    let mut scores = vec![0.0f64; reduced.len()];
    if d == 0 {
        return scores;
    }
    for test_fold in 0..n_folds {
        let mut pos_rows = Vec::new();
        let mut neg_rows = Vec::new();
        for p in 0..reduced.len() {
            if folds[p] == test_fold {
                continue;
            }
            match train_label[p] {
                Some(l) if l >= 0.5 => pos_rows.push(p),
                Some(_) => neg_rows.push(p),
                None => {}
            }
        }
        if pos_rows.len() < d + 2 || neg_rows.len() < d + 2 {
            continue;
        }
        // Balance negatives to 3x positives (deterministic hash subsample).
        let neg_cap = (pos_rows.len() * 3).max(2000);
        if neg_rows.len() > neg_cap {
            neg_rows.sort_by_key(|&p| {
                let mut x = (reduced[p] as u64).wrapping_mul(0x9e3779b97f4a7c15);
                x ^= x >> 29;
                x
            });
            neg_rows.truncate(neg_cap);
        }
        // Standardization stats over all training rows.
        let train_rows: Vec<usize> = pos_rows.iter().chain(neg_rows.iter()).copied().collect();
        let nt = train_rows.len() as f64;
        let mut mean = vec![0.0f64; d];
        for &p in &train_rows {
            for j in 0..d {
                mean[j] += feats[p][j] as f64;
            }
        }
        for m in mean.iter_mut() {
            *m /= nt;
        }
        let mut sd = vec![0.0f64; d];
        for &p in &train_rows {
            for j in 0..d {
                let dx = feats[p][j] as f64 - mean[j];
                sd[j] += dx * dx;
            }
        }
        for s in sd.iter_mut() {
            *s = (*s / nt).sqrt().max(1e-9);
        }
        let z = |p: usize, j: usize| -> f64 { (feats[p][j] as f64 - mean[j]) / sd[j] };
        // Standardized class means.
        let mut mp = vec![0.0f64; d];
        let mut mn = vec![0.0f64; d];
        for &p in &pos_rows {
            for j in 0..d {
                mp[j] += z(p, j);
            }
        }
        for &p in &neg_rows {
            for j in 0..d {
                mn[j] += z(p, j);
            }
        }
        for j in 0..d {
            mp[j] /= pos_rows.len() as f64;
            mn[j] /= neg_rows.len() as f64;
        }
        // Shrunk pooled within-class covariance.
        let mut cov = vec![vec![0.0f64; d]; d];
        for &p in &pos_rows {
            for a in 0..d {
                let da = z(p, a) - mp[a];
                for b in 0..d {
                    cov[a][b] += da * (z(p, b) - mp[b]);
                }
            }
        }
        for &p in &neg_rows {
            for a in 0..d {
                let da = z(p, a) - mn[a];
                for b in 0..d {
                    cov[a][b] += da * (z(p, b) - mn[b]);
                }
            }
        }
        let denom = (pos_rows.len() + neg_rows.len()) as f64;
        for a in 0..d {
            for b in 0..d {
                cov[a][b] /= denom;
            }
            cov[a][a] += 0.3; // diagonal loading (features are ~unit variance) -> robust
        }
        let rhs: Vec<f64> = (0..d).map(|j| mp[j] - mn[j]).collect();
        let w = match solve_linear(cov, rhs) {
            Some(w) => w,
            None => continue,
        };
        for p in 0..reduced.len() {
            if folds[p] == test_fold {
                let mut s = 0.0;
                for j in 0..d {
                    s += w[j] * z(p, j);
                }
                scores[p] = s;
            }
        }
    }
    scores
}

/// Dispatch to the configured discriminant ("gbdt", "lda", or default
/// "deviation"). The learned paths featurize each reduced PFM, label confident
/// positives vs decoys, and cross-validate by feature; deviation is the robust
/// closed-form mass+RT(+IM) consistency score.
#[allow(clippy::too_many_arguments)]
fn discriminant_scores(
    cfg: &ScoreConfig,
    reduced: &[usize],
    pfms: &[crate::matching::Pfm],
    matrix: &FeatureMatrix,
    db: &PeptideDb,
    rt_diff: &[Option<f64>],
    im_diff: &[Option<f64>],
    mass_z: &[f64],
    positives: &[usize],
    folds: &[usize],
) -> Vec<f64> {
    if cfg.discriminant == "gbdt" || cfg.discriminant == "lda" {
        let mut train_label = vec![None; reduced.len()];
        for &pos in positives {
            train_label[pos] = Some(1.0);
        }
        for pos in 0..reduced.len() {
            if pfms[reduced[pos]].decoy {
                train_label[pos] = Some(0.0);
            }
        }
        let feats: Vec<Vec<f32>> = reduced
            .iter()
            .enumerate()
            .map(|(pos, &i)| featurize(&pfms[i], matrix, db, rt_diff[pos], im_diff[pos]))
            .collect();
        if cfg.discriminant == "lda" {
            cv_lda_predict(reduced, &feats, folds, &train_label, cfg.n_folds)
        } else {
            cv_train_predict(reduced, &feats, folds, &train_label, cfg)
        }
    } else {
        deviation_scores(reduced, rt_diff, im_diff, mass_z, positives)
    }
}

/// Solve the symmetric positive-definite system `a x = b` by Gaussian
/// elimination with partial pivoting. `a` is consumed. Small (p≈69) so dense
/// elimination is fine.
fn solve_linear(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    for col in 0..n {
        // partial pivot
        let mut piv = col;
        for r in col + 1..n {
            if a[r][col].abs() > a[piv][col].abs() {
                piv = r;
            }
        }
        if a[piv][col].abs() < 1e-12 {
            return None;
        }
        a.swap(col, piv);
        b.swap(col, piv);
        let d = a[col][col];
        for r in col + 1..n {
            let f = a[r][col] / d;
            if f == 0.0 {
                continue;
            }
            for c in col..n {
                a[r][c] -= f * a[col][c];
            }
            b[r] -= f * b[col];
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut s = b[i];
        for j in i + 1..n {
            s -= a[i][j] * x[j];
        }
        x[i] = s / a[i][i];
    }
    Some(x)
}

/// Ridge regression weights for `X w ≈ y`: `(XᵀX + λI) w = Xᵀy`. Used to fit an
/// **on-data** RT predictor (composition features → this run's observed RT),
/// the DirectMS1/ms1searchpy trick that makes RT a real discriminator on a
/// gradient the frozen model never saw. Returns None if the system is singular.
fn fit_ridge(x: &[Vec<f64>], y: &[f64], lambda: f64) -> Option<Vec<f64>> {
    let p = x[0].len();
    let mut ata = vec![vec![0.0; p]; p];
    let mut aty = vec![0.0; p];
    for (row, &yi) in x.iter().zip(y) {
        for i in 0..p {
            aty[i] += row[i] * yi;
            let ri = row[i];
            for j in i..p {
                ata[i][j] += ri * row[j];
            }
        }
    }
    // symmetrize + ridge
    for i in 0..p {
        for j in i + 1..p {
            ata[j][i] = ata[i][j];
        }
        ata[i][i] += lambda;
    }
    solve_linear(ata, aty)
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
    let v: Vec<f64> = vals.filter(|x| x.is_finite()).collect();
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

/// Deviation discriminant (ms1searchpy/blitzms1 style): score each PFM by how
/// consistent its calibrated RT and IM are with the prediction, relative to the
/// spread of the confident-positive set. Lower combined z-deviation -> higher
/// score. On the saturated-mass MS1 regime this beats a gbdt because the only
/// information is the orthogonal RT/IM agreement, and the gbdt drowns it under
/// per-feature-constant features (mass, charge, intensity).
///
/// NOTE: a heteroscedastic (RT-binned) sigma and a 2D mass-recalibration surface
/// were both tested here (2026-06-24) and REJECTED — both LOWERED depth (mass
/// surface 218->168, RT-binning 218->213) without improving FDR. The koth
/// features are already in a well-calibrated reference-run space, so there is no
/// systematic (m/z, RT) ppm drift to remove (the poly only overfits noise), and
/// the post-bootstrap RT residual is genuinely homoscedastic (~0.025 min). A
/// single global mass shift + global RT sigma is optimal. See memory.
fn deviation_scores(
    reduced: &[usize],
    rt_diff: &[Option<f64>],
    im_diff: &[Option<f64>],
    mass_z: &[f64],
    positives: &[usize],
) -> Vec<f64> {
    let rt_sig = robust_sigma(positives.iter().map(|&p| rt_diff[p].unwrap_or(f64::NAN)));
    let im_sig = robust_sigma(positives.iter().map(|&p| im_diff[p].unwrap_or(f64::NAN)));
    log::info!(
        "deviation discriminant: rt_sigma={:.3} min, im_sigma={:.4} 1/K0 (mass z-scored)",
        rt_sig,
        im_sig
    );
    // Order-preserving parallel map over immutable inputs -> byte-identical to the
    // serial version, but fans the 17M-row deviation computation across cores.
    (0..reduced.len())
        .into_par_iter()
        .map(|pos| {
            // mass_z is already in sigma units (recalibrated), so it enters directly.
            let mut d2 = mass_z[pos].powi(2);
            if let Some(r) = rt_diff[pos] {
                d2 += (r / rt_sig).powi(2);
            }
            if let Some(i) = im_diff[pos] {
                d2 += (i / im_sig).powi(2);
            }
            -d2
        })
        .collect()
}

/// Best-scoring reduced PFM per feature, with target-decoy q-values.
fn compete_and_qvalue(
    reduced: &[usize],
    pfms: &[Pfm],
    scores: &[f64],
    n_features: usize,
) -> Vec<(usize, f64, f64)> {
    // Effective score adds a tiny symmetric tiebreak so target/decoy ties
    // (common before RT separates them) resolve without bias toward targets.
    let eff = |pos: usize| -> f64 {
        let p = &pfms[reduced[pos]];
        scores[pos] + 1e-9 * (tiebreak(p.feature_idx, p.peptide_idx) - 0.5)
    };

    // feature -> (reduced position, effective score)
    let mut best: Vec<Option<(usize, f64)>> = vec![None; n_features];
    for pos in 0..reduced.len() {
        let fi = pfms[reduced[pos]].feature_idx;
        let e = eff(pos);
        if best[fi].map_or(true, |(_, s)| e > s) {
            best[fi] = Some((pos, e));
        }
    }
    let mut winners: Vec<(usize, f64)> = best.into_iter().flatten().collect();
    // sort by effective score descending
    winners.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

    // walk down, count targets/decoys, FDR = decoys/targets
    let mut q = vec![0.0; winners.len()];
    let (mut t, mut d) = (0usize, 0usize);
    for (i, &(pos, _)) in winners.iter().enumerate() {
        if pfms[reduced[pos]].decoy {
            d += 1;
        } else {
            t += 1;
        }
        q[i] = if t == 0 { 1.0 } else { d as f64 / t as f64 };
    }
    // enforce monotone non-decreasing q from the bottom up
    let mut run = f64::INFINITY;
    for i in (0..q.len()).rev() {
        run = run.min(q[i]);
        q[i] = run.min(1.0);
    }
    winners
        .into_iter()
        .zip(q)
        .map(|((pos, score), qv)| (pos, score, qv))
        .collect()
}

/// Run the full scoring + FDR pipeline. Returns one [`Identification`] per
/// feature that matched at least one peptide.
pub fn score(
    matrix: &FeatureMatrix,
    db: &PeptideDb,
    pfms: &[Pfm],
    rt_model: Option<&RtModel>,
    im_model: Option<&crate::im::ImModel>,
    rt_library: Option<&std::collections::HashMap<String, f64>>,
    cfg: &ScoreConfig,
) -> ScoreResult {
    let n_features = matrix.features.len();
    let reduced = reduce_top_k(pfms, n_features, cfg.top_k);
    log::info!("reduced {} PFMs to {} (top-{} per feature)", pfms.len(), reduced.len(), cfg.top_k);

    // Abstract RT per reduced PFM. Prefer a precomputed DeepLC library (sequence-
    // aware, generalizes — the discriminating signal); fall back to the distilled
    // model. Either is then calibrated on-data to this gradient.
    let abstract_rt: Vec<Option<f64>> = if let Some(lib) = rt_library {
        let mut hits = 0usize;
        let v: Vec<Option<f64>> = reduced
            .iter()
            .map(|&i| {
                let seq = db.peptides[pfms[i].peptide_idx].sequence.as_ref();
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
                let mut uniq: Vec<usize> = reduced.iter().map(|&i| pfms[i].peptide_idx).collect();
                uniq.sort_unstable();
                uniq.dedup();
                let seqs: Vec<&str> =
                    uniq.iter().map(|&pi| db.peptides[pi].sequence.as_ref()).collect();
                let upreds = m.predict_many(&seqs);
                let pred_of: std::collections::HashMap<usize, f64> =
                    uniq.iter().copied().zip(upreds).collect();
                log::info!(
                    "predicted abstract RT for {} unique peptides ({} PFMs) ({:.1}s)",
                    uniq.len(),
                    reduced.len(),
                    t.elapsed().as_secs_f64()
                );
                reduced.par_iter().map(|&i| pred_of.get(&pfms[i].peptide_idx).copied()).collect()
            }
            None => vec![None; reduced.len()],
        }
    };

    // Predict 1/K0 for each reduced PFM (charge-aware), if an IM model is given
    // and the data has ion mobility.
    let has_im = matrix.features.iter().any(|f| f.im.is_some());
    let predicted_im: Vec<Option<f64>> = match (im_model, has_im) {
        (Some(m), true) => {
            let t = std::time::Instant::now();
            let items: Vec<(&str, u8)> = reduced
                .iter()
                .map(|&i| {
                    let p = &pfms[i];
                    (db.peptides[p.peptide_idx].sequence.as_ref(), matrix.features[p.feature_idx].charge)
                })
                .collect();
            let preds = m.predict_many(&items);
            log::info!("predicted 1/K0 for {} reduced PFMs ({:.1}s)", preds.len(), t.elapsed().as_secs_f64());
            preds.into_iter().map(Some).collect()
        }
        _ => vec![None; reduced.len()],
    };

    // Diagnostic dump: every reduced candidate's predicted vs observed RT & IM,
    // to test alignment recoverability offline.
    if std::env::var("UNO_DUMP_REDUCED").is_ok() {
        use std::io::Write;
        let mut w = std::io::BufWriter::new(std::fs::File::create("/tmp/uno_reduced.tsv").unwrap());
        writeln!(w, "feature_idx\tdecoy\tcharge\tmass_err_ppm\tabstract_rt\tobs_rt\tpred_im\tobs_im\tpeptide").ok();
        for (pos, &i) in reduced.iter().enumerate() {
            let p = &pfms[i];
            let f = &matrix.features[p.feature_idx];
            writeln!(
                w,
                "{}\t{}\t{}\t{:.3}\t{}\t{:.5}\t{}\t{}\t{}",
                p.feature_idx,
                p.decoy as u8,
                f.charge,
                p.mass_error_ppm,
                abstract_rt[pos].map(|v| format!("{:.4}", v)).unwrap_or_default(),
                f.rt,
                predicted_im[pos].map(|v| format!("{:.5}", v)).unwrap_or_default(),
                f.im.map(|v| format!("{:.5}", v)).unwrap_or_default(),
                db.peptides[p.peptide_idx].sequence.as_ref(),
            ).ok();
        }
        log::info!("dumped {} reduced candidates to /tmp/uno_reduced.tsv", reduced.len());
    }

    // Seed RT calibration AND the initial positive training set from a
    // classifier-independent high-confidence subset. Reversed decoys are mass
    // twins, so mass alone cannot bootstrap; koth quality + cross-run
    // reproducibility + tight mass + RT-inlier gives a clean positive seed.
    let has_rt = abstract_rt.iter().any(|a| a.is_some());
    let mut mass_anchor_shift: Option<(f64, f64)> = None;
    let mut mass_anchor_proteins: Option<Vec<ProteinId>> = None;
    // Winning RT calibration bundle for the hill second pass, filled by the mass
    // anchor: (cal, rt_sigma, mass_shift_ppm, mass_sigma_ppm).
    let mut hill_cal: Option<(Calibration, f64, f64, f64)> = None;
    let (mut calibration, seed_positions) = if !has_rt {
        (Calibration::default(), Vec::new())
    } else if cfg.mass_anchor {
        // NOTE: an on-data RT fine-tune (refit the MLP readout over its frozen
        // embedding on the confident peptides, re-predict, re-run) was tried here and
        // does NOT help this path — even with an RT-inlier seed it gave 712 -> 333
        // (auto-rejected). The robust two-stage calibration already captures the
        // run-level RT shift, and a 129-d readout refit on the confident subset
        // overfits its composition and generalizes worse than the pretrained posMLP.
        // The DeepLC-library path (rt_model = None) never ran it anyway.
        let (cal, anchors, shift, sigma) =
            mass_anchor_calibration(&reduced, pfms, matrix, db, &abstract_rt, cfg);
        mass_anchor_shift = Some((shift, sigma));
        mass_anchor_proteins = Some(MASS_ANCHOR_PROTEINS.with(|c| std::mem::take(&mut *c.borrow_mut())));
        hill_cal = MASS_ANCHOR_CAL.with(|c| c.borrow_mut().take());
        (cal, anchors)
    } else {
        seed_calibration(&reduced, pfms, matrix, &abstract_rt, cfg.seed_min_runs_frac)
    };

    // Fit IM calibration (im2deep reference -> instrument 1/K0) PER CHARGE on the
    // high-confidence RT seed — the offset between predicted and observed 1/K0 is
    // charge-dependent, so one global curve leaves higher charges biased. Charges
    // with too few seed points fall back to a pooled calibration.
    let im_calibration: std::collections::HashMap<u8, Calibration> =
        if !seed_positions.is_empty() && im_model.is_some() {
            let mut by_charge: std::collections::HashMap<u8, Vec<(f64, f64)>> =
                std::collections::HashMap::new();
            let mut pooled: Vec<(f64, f64)> = Vec::new();
            for &pos in &seed_positions {
                let f = &matrix.features[pfms[reduced[pos]].feature_idx];
                if let (Some(im_obs), Some(p)) = (f.im, predicted_im[pos]) {
                    by_charge.entry(f.charge).or_default().push((p, im_obs));
                    pooled.push((p, im_obs));
                }
            }
            let (pooled_cal, _) = Calibration::fit_ransac(&pooled, 0.03, 2000);
            let mut cals = std::collections::HashMap::new();
            for (z, pairs) in by_charge {
                let cal = if pairs.len() >= 30 {
                    let (c, inl) = Calibration::fit_ransac(&pairs, 0.02, 2000);
                    log::info!("IM calibration z={}: {} inliers / {} seed pairs", z, inl.len(), pairs.len());
                    c
                } else {
                    pooled_cal.clone()
                };
                cals.insert(z, cal);
            }
            cals
        } else {
            std::collections::HashMap::new()
        };
    let im_diff: Vec<Option<f64>> = reduced
        .iter()
        .enumerate()
        .map(|(pos, &i)| {
            let f = &matrix.features[pfms[i].feature_idx];
            let cal = im_calibration.get(&f.charge)?;
            let im_obs = f.im?;
            predicted_im[pos].map(|p| im_obs - cal.apply(p))
        })
        .collect();

    let folds: Vec<usize> = reduced
        .par_iter()
        .map(|&i| fold_of(pfms[i].feature_idx, cfg.n_folds))
        .collect();

    // rt_diff per reduced PFM under the current calibration.
    let compute_rt_diff = |calibration: &Calibration| -> Vec<Option<f64>> {
        // Order-preserving parallel indexed map (each row independent) -> byte-identical
        // to the serial version; parallelizes the 17M `calibration.apply` calls.
        reduced
            .par_iter()
            .enumerate()
            .map(|(pos, &i)| {
                let obs = matrix.features[pfms[i].feature_idx].rt;
                if cfg.rt_calibrated {
                    // Library is already in this run's RT space: use it directly.
                    return abstract_rt[pos].map(|a| obs - a);
                }
                if !calibration.is_fit() {
                    return None;
                }
                abstract_rt[pos].map(|a| obs - calibration.apply(a))
            })
            .collect()
    };

    let mut rt_diff = compute_rt_diff(&calibration);

    // On-data RT features (composition) for each reduced PFM, computed once.
    // Used to refit a per-run RT predictor (ridge) on the confident set.
    // On-data RT feature basis, DEDUPLICATED by peptide. The reduced set repeats each
    // peptide across ~26 features, so building the basis per-PFM re-embeds the same
    // ~1.2M unique rows ~26x (the dominant non-RT cost — 570s on run10) and would
    // materialize a ~16 GB matrix. Embed each UNIQUE peptide once; `rt_feat_idx[p]`
    // maps a reduced position to its unique row.
    let mut rt_feat_idx: Vec<u32> = Vec::new();
    let rt_feat_uniq: Vec<Vec<f64>> = if cfg.ondata_rt && rt_model.is_some() {
        let mut pi_to_uniq: rustc_hash::FxHashMap<usize, u32> = rustc_hash::FxHashMap::default();
        let mut uniq_pids: Vec<usize> = Vec::new();
        rt_feat_idx.reserve(reduced.len());
        for &i in &reduced {
            let pi = pfms[i].peptide_idx;
            let u = *pi_to_uniq.entry(pi).or_insert_with(|| {
                let n = uniq_pids.len() as u32;
                uniq_pids.push(pi);
                n
            });
            rt_feat_idx.push(u);
        }
        // Fine-tune mode: the MLP's frozen penultimate embedding as the feature basis
        // (keeps the pretrained representation). Falls back to raw composition features
        // if the model has no embedding (gbdt).
        if cfg.ondata_rt_finetune {
            let seqs: Vec<&str> = uniq_pids
                .iter()
                .map(|&pi| db.peptides[pi].sequence.as_ref())
                .collect();
            match rt_model.and_then(|m| m.embed_many(&seqs)) {
                Some(e) => {
                    log::info!("on-data RT fine-tune: MLP embedding basis ({}-d, {} unique peptides / {} PFMs)",
                        e.first().map_or(0, |v| v.len()), e.len(), reduced.len());
                    e
                }
                // Fine-tune requested but the model exposes no embedding (gbdt
                // student). Fall back to the FROZEN RT path (empty basis), NOT the
                // from-scratch composition ridge, which degrades depth.
                None => {
                    log::warn!("on-data RT fine-tune requested but RT model has no embedding; using frozen RT");
                    rt_feat_idx.clear();
                    Vec::new()
                }
            }
        } else {
            uniq_pids
                .par_iter()
                .map(|&pi| crate::rt::rt_features(db.peptides[pi].sequence.as_ref()))
                .collect()
        }
    } else {
        Vec::new()
    };
    // Refit the on-data RT model on `pos` features (observed RT as target) and
    // return rt_diff = obs_rt - on_data_pred for every reduced PFM. Falls back to
    // the existing rt_diff if the fit is degenerate or too few points.
    let obs_rt = |p: usize| -> f64 { matrix.features[pfms[reduced[p]].feature_idx].rt };
    let predict_rt = |w: &[f64], p: usize| -> f64 {
        rt_feat_uniq[rt_feat_idx[p] as usize].iter().zip(w).map(|(f, c)| f * c).sum()
    };
    // Iterative two-stage RT (ms1searchpy-style): fit a per-run composition->RT
    // model, reject RT outliers (contamination in the seed), refit on the clean
    // inliers, repeat. Each round tightens the residual sigma — the lever that
    // takes RT from a weak ~0.7 min to a sharp ~0.1-0.3 min discriminator.
    let refit_ondata_rt = |pos: &[usize], fallback: &[Option<f64>]| -> Vec<Option<f64>> {
        if rt_feat_uniq.is_empty() || pos.len() < 50 {
            return fallback.to_vec();
        }
        let mut keep: Vec<usize> = pos.to_vec();
        let mut weights: Option<Vec<f64>> = None;
        for _round in 0..5 {
            if keep.len() < 30 {
                break;
            }
            let x: Vec<Vec<f64>> = keep.iter().map(|&p| rt_feat_uniq[rt_feat_idx[p] as usize].clone()).collect();
            let y: Vec<f64> = keep.iter().map(|&p| obs_rt(p)).collect();
            let w = match fit_ridge(&x, &y, cfg.ondata_rt_lambda) {
                Some(w) => w,
                None => break,
            };
            let resids: Vec<f64> = keep.iter().map(|&p| obs_rt(p) - predict_rt(&w, p)).collect();
            let sigma = robust_sigma(resids.iter().copied()).max(0.01);
            let new_keep: Vec<usize> = keep
                .iter()
                .zip(&resids)
                .filter(|(_, r)| r.abs() <= 2.5 * sigma)
                .map(|(&p, _)| p)
                .collect();
            let converged = new_keep.len() == keep.len();
            weights = Some(w);
            if new_keep.len() >= 30 {
                keep = new_keep;
            }
            if converged {
                break;
            }
        }
        match weights {
            Some(w) => {
                let sig =
                    robust_sigma(keep.iter().map(|&p| obs_rt(p) - predict_rt(&w, p)));
                log::info!(
                    "on-data RT: {} seed -> {} inliers, residual sigma = {:.3} min",
                    pos.len(),
                    keep.len(),
                    sig
                );
                (0..reduced.len()).map(|i| Some(obs_rt(i) - predict_rt(&w, i))).collect()
            }
            None => fallback.to_vec(),
        }
    };

    // Initial positives = the high-confidence RT seed, further purified by IM
    // consistency when available. On fast gradients RT is weak but IM is sharp,
    // so the RT∩IM intersection is a much cleaner positive set to train on.
    let mut positives: Vec<usize> = if im_model.is_some() && im_diff.iter().any(|d| d.is_some()) {
        let before = seed_positions.len();
        let kept: Vec<usize> = seed_positions
            .into_iter()
            .filter(|&pos| im_diff[pos].map_or(false, |d| d.abs() <= 0.02))
            .collect();
        log::info!("IM-purified positives: {} / {} RT-seed features", kept.len(), before);
        kept
    } else {
        seed_positions
    };

    // Mass recalibration: estimate the systematic mass shift + sigma from the
    // confident seed, then express every PFM's mass error as a z-score. A tight,
    // recalibrated mass is ms1searchpy's single biggest background-suppression
    // lever (on the HeLa benchmark it cut p(random) from 0.023 to 0.004).
    let mut mass_z: Vec<f64> = if let Some((shift, sigma)) = mass_anchor_shift {
        // Mass shift/sigma already estimated cleanly in the mass-anchor bootstrap.
        log::info!("mass recalibration (mass-anchor): shift={:.3} ppm, sigma={:.3} ppm", shift, sigma);
        reduced.par_iter().map(|&i| (pfms[i].mass_error_ppm - shift) / sigma).collect()
    } else {
        let mut seed_mass: Vec<f64> = positives
            .iter()
            .map(|&p| pfms[reduced[p]].mass_error_ppm)
            .collect();
        if seed_mass.len() >= 20 {
            let shift = median(&mut seed_mass);
            // Gaussian-core sigma (sigma-clipped): the random-match tails inflate a
            // plain MAD, which keeps the mass window too loose to suppress background.
            let sigma =
                clipped_sigma(positives.iter().map(|&p| pfms[reduced[p]].mass_error_ppm - shift))
                    .max(0.05);
            log::info!(
                "mass recalibration: shift={:.3} ppm, sigma={:.3} ppm (from {} seed features)",
                shift,
                sigma,
                positives.len()
            );
            reduced
                .par_iter()
                .map(|&i| (pfms[i].mass_error_ppm - shift) / sigma)
                .collect()
        } else {
            log::info!("mass recalibration: too few seed ({}), mass not used", positives.len());
            vec![0.0; reduced.len()]
        }
    };

    // Seed the on-data RT model from the initial confident positives before the
    // first scoring round, so the deviation discriminant sees per-run RT.
    if cfg.ondata_rt && !positives.is_empty() {
        rt_diff = refit_ondata_rt(&positives, &rt_diff);
    }

    let n_iters = cfg.bootstrap_rounds.max(1) + 1;
    let mut scores: Vec<f64> = Vec::new();
    for iter in 0..n_iters {
        if positives.len() < 50 {
            log::warn!("iter {}: too few positives ({}) to train; stopping", iter, positives.len());
            break;
        }

        scores = discriminant_scores(
            cfg, &reduced, pfms, matrix, db, &rt_diff, &im_diff, &mass_z, &positives, &folds,
        );
        let winners = compete_and_qvalue(&reduced, pfms, &scores, n_features);

        // Next positives / RT-refit seed = the top mass-tight target winners by
        // score. This is far cleaner than the broad RANSAC seed (which is ~50%
        // random) and large enough to fit a sharp per-run RT model — the way
        // ms1searchpy bootstraps RT from confidently-scored peptides.
        let mut scored_targets: Vec<(usize, f64)> = winners
            .iter()
            .filter(|&&(pos, _, _)| !pfms[reduced[pos]].decoy && mass_z[pos].abs() <= 2.0)
            .map(|&(pos, s, _)| (pos, s))
            .collect();
        scored_targets.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        let n_keep = scored_targets.len().min(2000).max(0);
        let mut next: Vec<usize> = scored_targets[..n_keep].iter().map(|&(p, _)| p).collect();
        // Fall back to confident-FDR winners if the top-score set is too small.
        if next.len() < 50 {
            next = winners
                .iter()
                .filter(|&&(pos, _, q)| q <= cfg.calibration_fdr && !pfms[reduced[pos]].decoy)
                .map(|&(pos, _, _)| pos)
                .collect();
        }
        let cal_pairs: Vec<(f64, f64)> = next
            .iter()
            .filter_map(|&pos| {
                abstract_rt[pos].map(|a| (a, matrix.features[pfms[reduced[pos]].feature_idx].rt))
            })
            .collect();
        let n_id = winners
            .iter()
            .filter(|(pos, _, q)| *q <= 0.01 && !pfms[reduced[*pos]].decoy)
            .count();
        log::info!(
            "iter {}: trained on {} positives -> {} target IDs at 1% FDR",
            iter,
            positives.len(),
            n_id
        );

        if cfg.ondata_rt && next.len() >= 50 {
            // Refit the on-data RT model on this round's confident targets.
            rt_diff = refit_ondata_rt(&next, &rt_diff);
        } else if has_rt && cal_pairs.len() >= 50 {
            // Iteratively tighten the (fixed-prediction) calibration on the
            // cleaner confident set. With a DeepLC library the predictions don't
            // move, so this converges to the real diagonal without overfitting —
            // unlike the composition ridge.
            calibration = Calibration::fit(&cal_pairs);
            rt_diff = compute_rt_diff(&calibration);
        }
        positives = next;
    }

    // === Iterative protein-bootstrap RT calibration (ms1searchpy's key step) ===
    // Identify confident PROTEINS (>=2 RT/mass-consistent target peptides), then
    // recalibrate RT on the UNION of their peptides. Protein-level evidence is an
    // ORTHOGONAL voucher for peptide correctness (a peptide from a multi-peptide
    // protein is real even if marginal per-PFM), so this yields a large, clean,
    // RT-independent calibration set -> tight global RT -> keeps ~3x more real
    // peptides per protein. Without it uno gets ~2.2 vs ms1searchpy's ~7.7.
    if has_rt && !scores.is_empty() {
        for pround in 0..4 {
            // Best (highest-score) reduced position per peptide, with deviations.
            let mut best_pos: std::collections::HashMap<usize, (f64, usize)> =
                std::collections::HashMap::new();
            for pos in 0..reduced.len() {
                let pi = pfms[reduced[pos]].peptide_idx;
                let e = best_pos.entry(pi).or_insert((f64::NEG_INFINITY, pos));
                if scores[pos] > e.0 {
                    *e = (scores[pos], pos);
                }
            }
            // RT sigma from current confident positives to set the "consistent" gate.
            let rt_sig = robust_sigma(positives.iter().map(|&p| rt_diff[p].unwrap_or(f64::NAN)))
                .max(0.02);
            let rt_gate = (2.5 * rt_sig).min(1.0);
            // Count consistent target peptides per protein.
            let mut prot_peps: std::collections::HashMap<&str, Vec<usize>> =
                std::collections::HashMap::new();
            for (&pi, &(_, pos)) in &best_pos {
                if pfms[reduced[pos]].decoy {
                    continue;
                }
                let rt_ok = rt_diff[pos].map_or(false, |d| d.abs() <= rt_gate);
                let mass_ok = mass_z[pos].abs() <= 2.5;
                let im_ok = im_diff[pos].map_or(true, |d| d.abs() <= 0.025);
                if rt_ok && mass_ok && im_ok {
                    for prot in &db.peptides[pi].proteins {
                        prot_peps.entry(prot.as_ref()).or_default().push(pos);
                    }
                }
            }
            // Trusted calibration set = peptides of proteins with >=2 consistent peptides.
            let mut trusted: Vec<usize> = Vec::new();
            for (_, peps) in &prot_peps {
                if peps.len() >= 2 {
                    trusted.extend_from_slice(peps);
                }
            }
            trusted.sort_unstable();
            trusted.dedup();
            if trusted.len() < 50 {
                log::info!("protein-bootstrap round {}: only {} trusted peptides, stopping", pround, trusted.len());
                break;
            }
            // Recalibrate RT on the trusted set (RANSAC locks the diagonal cleanly).
            let pairs: Vec<(f64, f64)> = trusted
                .iter()
                .filter_map(|&pos| {
                    abstract_rt[pos].map(|a| (a, matrix.features[pfms[reduced[pos]].feature_idx].rt))
                })
                .collect();
            if pairs.len() < 50 {
                break;
            }
            let (cal, _) = Calibration::fit_ransac(&pairs, (3.0 * rt_sig).max(0.3), 4000);
            calibration = cal;
            rt_diff = compute_rt_diff(&calibration);
            // Re-estimate the mass recalibration from this clean trusted set: its
            // sigma is the true peak width (~0.3 ppm), not the contaminated-seed
            // 0.8 ppm — this tightens the mass window and suppresses background.
            {
                let mut tm: Vec<f64> = trusted.iter().map(|&p| pfms[reduced[p]].mass_error_ppm).collect();
                if tm.len() >= 20 {
                    let sh = median(&mut tm);
                    let sg = clipped_sigma(trusted.iter().map(|&p| pfms[reduced[p]].mass_error_ppm - sh)).max(0.05);
                    for i in 0..reduced.len() {
                        mass_z[i] = (pfms[reduced[i]].mass_error_ppm - sh) / sg;
                    }
                    if pround == 0 {
                        log::info!("mass re-estimate (trusted set): shift={:.3} ppm, sigma={:.3} ppm", sh, sg);
                    }
                }
            }
            positives = trusted.clone();
            scores = discriminant_scores(
                cfg, &reduced, pfms, matrix, db, &rt_diff, &im_diff, &mass_z, &positives, &folds,
            );
            log::info!(
                "protein-bootstrap round {}: {} trusted peptides from multi-peptide proteins, rt_gate={:.3} min",
                pround,
                trusted.len(),
                rt_gate
            );
        }
        // Dump the confident (peptide, observed RT) set for DeepLC's own
        // calibrate_preds() — ms1searchpy's two-stage RT. A monotone PWL on a
        // narrow seed doesn't generalize; DeepLC calibrated on this clean,
        // gradient-spanning set predicts every peptide tightly (sigma ~0.1).
        if let Ok(path) = std::env::var("UNO_DUMP_RT_CALIB") {
            use std::io::Write;
            if let Ok(f) = std::fs::File::create(&path) {
                let mut w = std::io::BufWriter::new(f);
                let mut seen = std::collections::HashSet::new();
                for &pos in &positives {
                    let p = &pfms[reduced[pos]];
                    let seq = db.peptides[p.peptide_idx].sequence.as_ref();
                    if seen.insert(seq.to_string()) {
                        let rt = matrix.features[p.feature_idx].rt;
                        let _ = writeln!(w, "{}\t{}", seq, rt);
                    }
                }
                log::info!("dumped {} RT-calibration peptides to {}", seen.len(), path);
            }
        }
    }

    // Final competition and identifications.
    if scores.is_empty() {
        log::warn!("no discriminant was trained; returning no identifications");
        let (rt_cal, rt_sigma, mass_shift, mass_sigma) = match hill_cal {
            Some((c, s, ms, msig)) => (Some(c), s, ms, msig),
            None => (None, 1.0, 0.0, 1.0),
        };
        return ScoreResult { ids: Vec::new(), peptide_matches: Vec::new(), proteins: mass_anchor_proteins, pfm_dev: Vec::new(), rt_cal, rt_sigma, mass_shift, mass_sigma };
    }
    // Whenever mass-anchor ran, its protein list is what main.rs ALWAYS reports
    // (unconditionally, regardless of this function's own discriminant result --
    // see main.rs's `match result.proteins.clone() { Some(p) => ..., None => ... }`).
    // But everything built below (ids, peptide_matches, pfm_dev, and in turn
    // main.rs's "matched peptides"/RT-gate diagnostic) used THIS function's own
    // discriminant-path `rt_diff`, which can be badly calibrated -- this file's
    // on-data iter-0 divergence, or later rounds that still don't generalize to the
    // full peptide population -- without that ever affecting the reported protein
    // counts (mass-anchor's own calibration is separate and validated). It just
    // silently corrupted these diagnostics (e.g. a nonsensical ~99 min "RT gate").
    // Recompute rt_diff under mass-anchor's OWN calibration before using it below.
    if let Some((ref mcal, _, _, _)) = hill_cal {
        rt_diff = compute_rt_diff(mcal);
    }
    let winners = compete_and_qvalue(&reduced, pfms, &scores, n_features);
    let n_id = winners
        .iter()
        .filter(|(pos, _, q)| *q <= 0.01 && !pfms[reduced[*pos]].decoy)
        .count();
    log::info!("final: {} target IDs at 1% FDR", n_id);

    let ids: Vec<Identification> = winners
        .iter()
        .map(|&(pos, score, qv)| {
            let ri = reduced[pos];
            let p = &pfms[ri];
            Identification {
                feature_idx: p.feature_idx,
                peptide_idx: p.peptide_idx,
                decoy: p.decoy,
                score,
                q_value: qv,
                mass_error_ppm: p.mass_error_ppm,
                isotope_offset: p.isotope_offset,
                rt_diff: rt_diff[pos].unwrap_or(f64::NAN),
                im_diff: im_diff[pos].unwrap_or(f64::NAN),
            }
        })
        .collect();

    // Best match per peptide (across all features, no competition): protein
    // evidence. Ranked by RT consistency; carries IM consistency alongside.
    // Rank a peptide's matches by the full deviation score (mass+rt+im), so the
    // representative is the most consistent overall — not just the closest in RT.
    // Select the representative match per peptide by TIGHTEST MASS (RT-independent),
    // NOT by the RT-inclusive discriminant score. A random peptide mass-coincides
    // with many features; picking the one with the smallest rt_diff is a
    // best-of-many cherry-pick that hands random matches an artificially tight RT,
    // collapsing RT's discriminating power to ~1x enrichment. Choosing by mass
    // leaves RT as an honest independent check (ms1searchpy does the same).
    // entry: (sel_key=-mass_z, score, abs_rt_diff, abs_im_diff, abs_mass_z, feat_q)
    let mut by_peptide: std::collections::HashMap<usize, (f64, f64, f64, f64, f64, f64)> =
        std::collections::HashMap::new();
    for pos in 0..reduced.len() {
        let p = &pfms[reduced[pos]];
        let ard = rt_diff[pos].map(|d| d.abs()).unwrap_or(f64::INFINITY);
        let aid = im_diff[pos].map(|d| d.abs()).unwrap_or(f64::INFINITY);
        let amz = mass_z[pos].abs();
        let fq = matrix.features[p.feature_idx].n_isotopes.unwrap_or(0.0);
        let sel = -amz; // tightest mass wins (RT-independent representative)
        let e = by_peptide.entry(p.peptide_idx).or_insert((
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
            f64::INFINITY,
            f64::INFINITY,
            f64::INFINITY,
            0.0,
        ));
        if sel > e.0 {
            *e = (sel, scores[pos], ard, aid, amz, fq);
        }
    }
    let peptide_matches: Vec<PeptideMatch> = by_peptide
        .into_iter()
        .map(|(peptide_idx, (_sel, score, abs_rt_diff, abs_im_diff, abs_mass_z, feat_quality))| PeptideMatch {
            peptide_idx,
            decoy: db.peptides[peptide_idx].decoy,
            abs_rt_diff,
            abs_im_diff,
            abs_mass_z,
            score,
            feat_quality,
        })
        .collect();

    if let Ok(path) = std::env::var("UNO_DUMP_IDS") {
        use std::io::Write;
        if let Ok(f) = std::fs::File::create(&path) {
            let mut w = std::io::BufWriter::new(f);
            let _ = writeln!(w, "decoy\tabs_rt_diff\tabs_mass_z\tq");
            for (pos, _s, qv) in winners.iter().map(|&(p, s, q)| (p, s, q)) {
                let amz = mass_z[pos].abs();
                let ard = rt_diff[pos].map(|d| d.abs()).unwrap_or(f64::INFINITY);
                let _ = writeln!(w, "{}\t{:.4}\t{:.4}\t{:.5}", pfms[reduced[pos]].decoy as u8, ard, amz, qv);
            }
            log::info!("dumped {} feature-winner ids to {}", winners.len(), path);
        }
    }
    if let Ok(path) = std::env::var("UNO_DUMP_PM") {
        use std::io::Write;
        if let Ok(f) = std::fs::File::create(&path) {
            let mut w = std::io::BufWriter::new(f);
            let _ = writeln!(w, "decoy\tabs_rt_diff\tabs_mass_z");
            for m in &peptide_matches {
                let _ = writeln!(w, "{}\t{:.4}\t{:.4}", m.decoy as u8, m.abs_rt_diff, m.abs_mass_z);
            }
            log::info!("dumped {} peptide_matches to {}", peptide_matches.len(), path);
        }
    }
    {
        // Stage-by-stage peptide-coverage diagnostics: where do real peptides go?
        let n_red_t = reduced.iter().filter(|&&i| !pfms[i].decoy).count();
        let n_red_d = reduced.len() - n_red_t;
        let mut dt = std::collections::HashSet::new();
        let mut dd = std::collections::HashSet::new();
        for &i in &reduced {
            if pfms[i].decoy { dd.insert(pfms[i].peptide_idx); } else { dt.insert(pfms[i].peptide_idx); }
        }
        let pm_t = peptide_matches.iter().filter(|m| !m.decoy).count();
        let pm_d = peptide_matches.iter().filter(|m| m.decoy).count();
        // how many distinct peptides have a TIGHT match (mass_z<=2 AND rt within 3 min)
        let tight_t = peptide_matches.iter().filter(|m| !m.decoy && m.abs_mass_z <= 2.0 && m.abs_rt_diff <= 3.0).count();
        let tight_d = peptide_matches.iter().filter(|m| m.decoy && m.abs_mass_z <= 2.0 && m.abs_rt_diff <= 3.0).count();
        log::info!(
            "PEPTIDE COVERAGE: reduced PFMs {} target / {} decoy; distinct peptides {} target / {} decoy; \
             peptide_matches {} target / {} decoy; tight(mass_z<=2,rt<=3min) {} target / {} decoy",
            n_red_t, n_red_d, dt.len(), dd.len(), pm_t, pm_d, tight_t, tight_d
        );
    }

    // Every reduced PFM's deviations, for the ms1searchpy-style any-feature
    // matched-peptide count in the caller.
    let pfm_dev: Vec<(usize, bool, f64, f64)> = (0..reduced.len())
        .map(|pos| {
            let p = &pfms[reduced[pos]];
            (
                p.peptide_idx,
                p.decoy,
                mass_z[pos].abs(),
                rt_diff[pos].map(|d| d.abs()).unwrap_or(f64::INFINITY),
            )
        })
        .collect();

    let (rt_cal, rt_sigma, mass_shift, mass_sigma) = match hill_cal {
        Some((c, s, ms, msig)) => (Some(c), s, ms, msig),
        None => (None, 1.0, 0.0, 1.0),
    };
    ScoreResult { ids, peptide_matches, proteins: mass_anchor_proteins, pfm_dev, rt_cal, rt_sigma, mass_shift, mass_sigma }
}
