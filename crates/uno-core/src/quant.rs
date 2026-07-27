//! Multi-sample protein quantitation over an aligned matrix (DirectMS1Quant-style).
//!
//! Peptide-level MS1 FDR is ~90%, so quant trusts nothing per-peptide. Instead it
//! rolls up to the PROTEIN: a real differentially-expressed protein has many
//! peptides that move the SAME direction between two conditions; a random-match
//! protein's peptides scatter. Per protein we count `n_sig` significant peptides
//! (|log2FC| past a decoy-calibrated threshold) and `k` = the largest co-directed
//! subset, then score co-direction with a binomial tail (P(≥k of n_sig share a
//! sign by chance), sign prob 0.5) and BH-correct across proteins. Decoy proteins
//! run through the identical pipeline as an independent FDR check.
//!
//! Intensities are median-normalized per run (the standard "most proteins are
//! unchanged" assumption). Abundance per condition is the linear sum of the
//! protein's peptides' mean condition intensity.

use anyhow::Context;

use crate::database::PeptideDb;
use crate::input::FeatureMatrix;
use crate::score::PeptideHit;
use std::collections::HashMap;

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
    let mut winner: std::collections::HashMap<u32, usize> = std::collections::HashMap::new();
    for (i, h) in peptides.iter().enumerate() {
        if !scores[i].is_finite() {
            continue;
        }
        match winner.get(&h.feature_idx) {
            Some(&j) if scores[j] >= scores[i] => {}
            _ => {
                winner.insert(h.feature_idx, i);
            }
        }
    }
    winner.values().map(|&i| peptides[i]).collect()
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
    by_feat.into_values().collect()
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

/// Per-protein differential-expression + abundance result.
#[derive(Debug, Clone)]
pub struct ProteinQuant {
    pub accession: String,
    pub decoy: bool,
    pub n_peptides: usize,
    pub n_significant: usize,
    /// Largest co-directed significant-peptide count.
    pub k: u32,
    pub log2_fc: f64,
    pub abundance_a: f64,
    pub abundance_b: f64,
    pub p_value: f64,
    pub q_value: f64,
}

/// Compute protein-level DE for the contrast `cond_a` vs `cond_b` (run-index
/// lists into `matrix.run_names`). `id_confident` is the ID-FDR-confident target
/// accession set; only those targets are tested (decoys are always kept as the
/// null). `min_peptides` drops proteins with too few quantified peptides.
pub fn differential(
    matrix: &FeatureMatrix,
    _db: &PeptideDb,
    peptides: &[PeptideHit],
    pep_group: &std::collections::HashMap<u32, (bool, String)>,
    cond_a: &[usize],
    cond_b: &[usize],
    min_peptides: usize,
) -> Vec<ProteinQuant> {
    let n_runs = matrix.run_names.len();
    // Normalize to the invariant backbone; a low per-run floor (5th pct of detected)
    // imputes missing cells so on/off peptides still yield a fold change.
    let (run_med, grand) = run_normalizers(matrix, peptides, n_runs);
    let mut per_run_logs: Vec<Vec<f64>> = vec![Vec::new(); n_runs];
    for h in peptides {
        let ints = &matrix.features[h.feature_idx as usize].intensities;
        for (r, &v) in ints.iter().enumerate().take(n_runs) {
            if v > 0.0 {
                per_run_logs[r].push(v.log2());
            }
        }
    }
    let pct5 = |v: &mut Vec<f64>| -> f64 {
        if v.is_empty() {
            return 0.0;
        }
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        v[v.len() / 20]
    };
    let run_floor: Vec<f64> = per_run_logs.iter_mut().map(|v| pct5(v)).collect();

    // Normalized log2 for run r of a feature's intensity (impute missing to floor).
    let norm = |ints: &[f64], r: usize| -> f64 {
        let v = ints.get(r).copied().unwrap_or(0.0);
        let l = if v > 0.0 { v.log2() } else { run_floor[r] };
        l - run_med[r] + grand
    };

    // Per-peptide log2FC + per-condition linear abundance (detected-only for abundance).
    struct Pep {
        fc: f64,
        ab_a: f64,
        ab_b: f64,
    }
    let cond_mean = |ints: &[f64], runs: &[usize]| -> f64 {
        let s: f64 = runs.iter().map(|&r| norm(ints, r)).sum();
        s / runs.len().max(1) as f64
    };
    let cond_lin = |ints: &[f64], runs: &[usize]| -> f64 {
        // linear abundance: mean of detected cells only (0 if none detected)
        let (mut s, mut n) = (0.0, 0u32);
        for &r in runs {
            let v = ints.get(r).copied().unwrap_or(0.0);
            if v > 0.0 {
                s += v;
                n += 1;
            }
        }
        if n == 0 {
            0.0
        } else {
            s / n as f64
        }
    };

    let mut pep_of: Vec<(usize, bool, Pep)> = Vec::with_capacity(peptides.len());
    for h in peptides {
        let ints = &matrix.features[h.feature_idx as usize].intensities;
        let a = cond_mean(ints, cond_a);
        let b = cond_mean(ints, cond_b);
        pep_of.push((
            h.peptide_idx as usize,
            h.decoy,
            Pep {
                fc: b - a,
                ab_a: cond_lin(ints, cond_a),
                ab_b: cond_lin(ints, cond_b),
            },
        ));
    }

    // Significance threshold: 2 SD of the decoy-peptide fold-change spread (the null).
    let dec_fc: Vec<f64> = pep_of.iter().filter(|(_, d, _)| *d).map(|(_, _, p)| p.fc).collect();
    let sd = std_dev(&dec_fc);
    let thr = 2.0 * sd;

    // Group peptides by PROTEIN GROUP, using only group-unique peptides (map
    // within a single group). Cross-group-shared peptides are absent from
    // `pep_group` and excluded — their intensity is a blend across groups.
    let mut by_prot: HashMap<(bool, &str), Vec<&Pep>> = HashMap::new();
    for (pi, _decoy, pep) in &pep_of {
        if let Some((d, leading)) = pep_group.get(&(*pi as u32)) {
            by_prot.entry((*d, leading.as_str())).or_default().push(pep);
        }
    }

    let mut out: Vec<ProteinQuant> = Vec::new();
    for ((decoy, acc), peps) in &by_prot {
        if peps.len() < min_peptides {
            continue;
        }
        let mut fcs: Vec<f64> = Vec::new();
        let mut sig: Vec<f64> = Vec::new();
        let (mut ab_a, mut ab_b) = (0.0, 0.0);
        for p in peps {
            ab_a += p.ab_a;
            ab_b += p.ab_b;
            fcs.push(p.fc);
            if p.fc.abs() >= thr {
                sig.push(p.fc);
            }
        }
        let kpos = sig.iter().filter(|&&x| x > 0.0).count() as u32;
        let kneg = sig.iter().filter(|&&x| x < 0.0).count() as u32;
        let k = kpos.max(kneg);
        let n_sig = sig.len() as u32;
        // t-test on ALL peptide fold changes (effect size + consistency), not the
        // co-direction count (which floors at 2·0.5^n and misses small strong
        // proteins). k/n_significant stay as informative columns.
        let p_value = ttest_p(&fcs, sd * sd * var_floor_mult());
        let log2_fc = median_of(&mut fcs.clone());
        out.push(ProteinQuant {
            accession: (*acc).to_string(),
            decoy: *decoy,
            n_peptides: peps.len(),
            n_significant: n_sig as usize,
            k,
            log2_fc,
            abundance_a: ab_a,
            abundance_b: ab_b,
            p_value,
            q_value: 1.0,
        });
    }

    benjamini_hochberg(&mut out);
    out.sort_by(|a, b| a.q_value.partial_cmp(&b.q_value).unwrap().then(a.p_value.partial_cmp(&b.p_value).unwrap()));
    out
}

/// BH-correct `p_value` -> `q_value` across TARGET proteins only (decoys keep q=1
/// but their p is comparable for an independent FDR check).
fn benjamini_hochberg(v: &mut [ProteinQuant]) {
    let mut idx: Vec<usize> = (0..v.len()).filter(|&i| !v[i].decoy).collect();
    idx.sort_by(|&a, &b| v[a].p_value.partial_cmp(&v[b].p_value).unwrap());
    let m = idx.len() as f64;
    let mut running = 1.0f64;
    for rank in (0..idx.len()).rev() {
        let i = idx[rank];
        let q = (v[i].p_value * m / (rank as f64 + 1.0)).min(1.0);
        running = running.min(q);
        v[i].q_value = running;
    }
}

/// Two-sided one-sample t-test p-value for H0: mean(v) = 0. Unlike a pure
/// co-direction count (which floors at 2·0.5^n and ignores effect size — so a
/// 7-peptide protein all shifted +2 log2 can't beat p=0.016), this uses the
/// magnitude AND consistency of the shift, so a small protein with a large,
/// tight fold change is correctly called. Returns 1.0 for n<2 or zero variance
/// with zero mean.
/// Multiplier on the population-noise variance floor (empirical-Bayes strength).
/// 1.0 = a protein can't look tighter than the baseline peptide scatter (safe,
/// conservative); lower = let genuinely-tight proteins keep more of their
/// confidence (more recall, more risk). Tunable via `UNO_VAR_FLOOR`.
fn var_floor_mult() -> f64 {
    // 0.25 chosen by ground-truth sweep (IonStar): the knee where recall is
    // recovered (10%->25% @5% FDR) at no precision cost, with margin against the
    // few-peptide degeneracy that a zero floor readmits.
    std::env::var("UNO_VAR_FLOOR").ok().and_then(|s| s.parse().ok()).unwrap_or(0.25)
}

fn ttest_p(v: &[f64], var_floor: f64) -> f64 {
    let n = v.len();
    if n < 2 {
        return 1.0;
    }
    let mean = v.iter().sum::<f64>() / n as f64;
    let var = v.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n as f64 - 1.0);
    // Floor the variance at the population per-peptide noise level. Without this a
    // protein with 2 near-identical peptides has se~0 -> spurious p~0. The floor is
    // an empirical-Bayes prior: a protein can't be measured more precisely than the
    // baseline peptide scatter.
    let var = var.max(var_floor);
    let se = (var / n as f64).sqrt();
    if se <= 0.0 {
        return 1.0;
    }
    let t = mean / se;
    let df = n as f64 - 1.0;
    // two-sided: p = I_{df/(df+t^2)}(df/2, 1/2)
    crate::protein::betai(df / 2.0, 0.5, df / (df + t * t)).clamp(1e-300, 1.0)
}

fn std_dev(v: &[f64]) -> f64 {
    if v.len() < 2 {
        return 1.0;
    }
    let m = v.iter().sum::<f64>() / v.len() as f64;
    (v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (v.len() - 1) as f64).sqrt()
}

fn median_of(v: &mut [f64]) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[v.len() / 2]
}

/// Dilution / monotonic-trend quant across ordered `levels` (each a condition's
/// run indices + its known level value, e.g. E. coli %). Stronger than a pairwise
/// contrast: instead of one fold change, each peptide's stat is the Spearman
/// correlation of its per-condition mean intensity with the level values, so a
/// protein whose peptides consistently track the dilution is caught even when no
/// single pair is extreme. Same co-direction binomial + BH rollup as
/// [`differential`]; `log2_fc` carries the median correlation, and abundance_a/b
/// are the lowest/highest level.
pub fn monotonic(
    matrix: &FeatureMatrix,
    _db: &PeptideDb,
    peptides: &[PeptideHit],
    pep_group: &std::collections::HashMap<u32, (bool, String)>,
    levels: &[(Vec<usize>, f64)],
    min_peptides: usize,
) -> Vec<ProteinQuant> {
    let n_runs = matrix.run_names.len();
    let (run_med, grand) = run_normalizers(matrix, peptides, n_runs);
    let level_vals: Vec<f64> = levels.iter().map(|(_, l)| *l).collect();

    // Per-peptide: Spearman(level, per-condition mean normalized-log2 of DETECTED
    // cells). Conditions with no detection are dropped; need >=3 points.
    struct Pep {
        rho: f64,
        first: f64,
        last: f64,
    }
    let cond_mean_det = |ints: &[f64], runs: &[usize]| -> Option<f64> {
        let (mut s, mut n) = (0.0, 0u32);
        for &r in runs {
            let v = ints.get(r).copied().unwrap_or(0.0);
            if v > 0.0 {
                s += v.log2() - run_med[r] + grand;
                n += 1;
            }
        }
        if n == 0 { None } else { Some(s / n as f64) }
    };
    let cond_lin = |ints: &[f64], runs: &[usize]| -> f64 {
        let (mut s, mut n) = (0.0, 0u32);
        for &r in runs {
            let v = ints.get(r).copied().unwrap_or(0.0);
            if v > 0.0 { s += v; n += 1; }
        }
        if n == 0 { 0.0 } else { s / n as f64 }
    };

    let mut pep_of: Vec<(usize, bool, Pep)> = Vec::with_capacity(peptides.len());
    for h in peptides {
        let ints = &matrix.features[h.feature_idx as usize].intensities;
        let mut lv = Vec::new();
        let mut mv = Vec::new();
        for (i, (runs, _)) in levels.iter().enumerate() {
            if let Some(m) = cond_mean_det(ints, runs) {
                lv.push(level_vals[i]);
                mv.push(m);
            }
        }
        if lv.len() < 3 {
            continue;
        }
        let rho = spearman(&lv, &mv);
        pep_of.push((
            h.peptide_idx as usize,
            h.decoy,
            Pep {
                rho,
                first: cond_lin(ints, &levels[0].0),
                last: cond_lin(ints, &levels[levels.len() - 1].0),
            },
        ));
    }

    let dec: Vec<f64> = pep_of.iter().filter(|(_, d, _)| *d).map(|(_, _, p)| p.rho).collect();
    // Significance = |rho| past the decoy-null spread, but CLAMPED to [0.4, 0.7]:
    // correlations are bounded [-1, 1], so a raw 2*SD threshold can exceed 1 and
    // flag nothing. A trending peptide needs a real monotonic correlation.
    let rho_sd = std_dev(&dec); let thr = (2.0 * rho_sd).clamp(0.4, 0.7);

    let mut by_prot: HashMap<(bool, &str), Vec<&Pep>> = HashMap::new();
    for (pi, _decoy, pep) in &pep_of {
        if let Some((d, leading)) = pep_group.get(&(*pi as u32)) {
            by_prot.entry((*d, leading.as_str())).or_default().push(pep);
        }
    }

    let mut out: Vec<ProteinQuant> = Vec::new();
    for ((decoy, acc), peps) in &by_prot {
        if peps.len() < min_peptides {
            continue;
        }
        let mut rhos: Vec<f64> = Vec::new();
        let mut sig: Vec<f64> = Vec::new();
        let (mut a, mut b) = (0.0, 0.0);
        for p in peps {
            a += p.first;
            b += p.last;
            rhos.push(p.rho);
            if p.rho.abs() >= thr {
                sig.push(p.rho);
            }
        }
        let kpos = sig.iter().filter(|&&x| x > 0.0).count() as u32;
        let kneg = sig.iter().filter(|&&x| x < 0.0).count() as u32;
        let k = kpos.max(kneg);
        let n_sig = sig.len() as u32;
        // t-test on all peptide correlations (H0: mean rho = 0) — magnitude-aware,
        // unlike the co-direction count.
        let p_value = ttest_p(&rhos, rho_sd * rho_sd * var_floor_mult());
        out.push(ProteinQuant {
            accession: (*acc).to_string(),
            decoy: *decoy,
            n_peptides: peps.len(),
            n_significant: n_sig as usize,
            k,
            log2_fc: median_of(&mut rhos.clone()),
            abundance_a: a,
            abundance_b: b,
            p_value,
            q_value: 1.0,
        });
    }
    benjamini_hochberg(&mut out);
    out.sort_by(|a, b| a.q_value.partial_cmp(&b.q_value).unwrap().then(a.p_value.partial_cmp(&b.p_value).unwrap()));
    out
}

/// Spearman rank correlation of two equal-length samples.
fn spearman(x: &[f64], y: &[f64]) -> f64 {
    let rank = |v: &[f64]| -> Vec<f64> {
        let mut idx: Vec<usize> = (0..v.len()).collect();
        idx.sort_by(|&a, &b| v[a].partial_cmp(&v[b]).unwrap());
        let mut r = vec![0.0; v.len()];
        let mut i = 0;
        while i < idx.len() {
            let mut j = i;
            while j + 1 < idx.len() && v[idx[j + 1]] == v[idx[i]] {
                j += 1;
            }
            let avg = (i + j) as f64 / 2.0 + 1.0; // average rank (1-based) for ties
            for &k in &idx[i..=j] {
                r[k] = avg;
            }
            i = j + 1;
        }
        r
    };
    pearson(&rank(x), &rank(y))
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

/// Parse a design file mapping runs to conditions and resolve a two-condition
/// contrast against `run_names`. Each line: `<run substring>\t<condition>`. A
/// matrix run column belongs to a condition if a design row's run field is a
/// substring of the column name. Returns (cond_a_runs, cond_b_runs) as indices.
pub fn resolve_contrast(
    run_names: &[String],
    design_tsv: &str,
    cond_a: &str,
    cond_b: &str,
) -> anyhow::Result<(Vec<usize>, Vec<usize>)> {
    let mut map: Vec<(String, String)> = Vec::new();
    for line in design_tsv.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut it = line.split('\t');
        if let (Some(run), Some(cond)) = (it.next(), it.next()) {
            map.push((run.trim().to_string(), cond.trim().to_string()));
        }
    }
    let cond_of = |rn: &str| -> Option<&str> {
        map.iter()
            .find(|(run, _)| rn == run || rn.contains(run.as_str()))
            .map(|(_, c)| c.as_str())
    };
    let (mut a, mut b) = (Vec::new(), Vec::new());
    for (i, rn) in run_names.iter().enumerate() {
        match cond_of(rn) {
            Some(c) if c == cond_a => a.push(i),
            Some(c) if c == cond_b => b.push(i),
            _ => {}
        }
    }
    if a.is_empty() || b.is_empty() {
        anyhow::bail!("contrast {cond_a} vs {cond_b}: matched {} / {} runs", a.len(), b.len());
    }
    Ok((a, b))
}

/// Run indices belonging to a single condition, per the design file (matched by
/// exact name or substring, same rule as [`resolve_contrast`]).
pub fn runs_for_condition(run_names: &[String], design_tsv: &str, cond: &str) -> Vec<usize> {
    let mut map: Vec<(String, String)> = Vec::new();
    for line in design_tsv.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut it = line.split('\t');
        if let (Some(run), Some(c)) = (it.next(), it.next()) {
            map.push((run.trim().to_string(), c.trim().to_string()));
        }
    }
    run_names
        .iter()
        .enumerate()
        .filter(|(_, rn)| {
            map.iter()
                .any(|(run, c)| c == cond && (rn.as_str() == run || rn.contains(run.as_str())))
        })
        .map(|(i, _)| i)
        .collect()
}

/// Resolve an ordered dilution spec `cond:level,cond:level,...` (e.g.
/// `A:3.3,B:10,C:33,D:67,E:90`) against the design, returning ordered
/// (run-indices, level) pairs. Conditions with no matched runs are an error.
pub fn resolve_levels(
    run_names: &[String],
    design_tsv: &str,
    spec: &str,
) -> anyhow::Result<Vec<(Vec<usize>, f64)>> {
    let mut out = Vec::new();
    for part in spec.split(',') {
        let (cond, level) = part
            .split_once(':')
            .with_context(|| format!("dilution spec part `{part}` must be cond:level"))?;
        let lv: f64 = level.trim().parse().with_context(|| format!("bad level `{level}`"))?;
        let runs = runs_for_condition(run_names, design_tsv, cond.trim());
        if runs.is_empty() {
            anyhow::bail!("dilution: no runs matched condition {}", cond.trim());
        }
        out.push((runs, lv));
    }
    Ok(out)
}
