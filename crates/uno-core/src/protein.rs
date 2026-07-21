//! Phase 6: protein-level aggregation + FDR.
//!
//! MS1 PSM-level scores are weak — random peptides matching in mass+RT are
//! hard to separate from real ones. But real peptides cluster into proteins
//! (several consistent peptides each), while random matches scatter as
//! singletons across random proteins. Aggregating to proteins and scoring each
//! by a binomial surprise (how unlikely its peptide count is by chance)
//! concentrates the signal. This is how DirectMS1/ms1searchpy reach depth.

use std::cmp::Ordering;
use rustc_hash::{FxHashMap, FxHashSet};
use std::collections::BinaryHeap;

/// One candidate in the greedy-parsimony max-heap: a protein `key`, its score,
/// and the peptide count `n_obs` the score was computed at. `n_obs` is a staleness
/// stamp — peptide counts only ever shrink, so an entry whose stored `n_obs` no
/// longer matches the protein's current count is outdated and skipped on pop.
/// Ordering reproduces the old `max_by` exactly: by score, ties broken by key.
struct GreedyEntry<'a> {
    score: f64,
    key: (bool, &'a str),
    n_obs: usize,
}
impl PartialEq for GreedyEntry<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}
impl Eq for GreedyEntry<'_> {}
impl Ord for GreedyEntry<'_> {
    fn cmp(&self, other: &Self) -> Ordering {
        // scores are finite (binom tail clamped to [1e-300,1]), so total_cmp
        // matches the old partial_cmp().unwrap(); key breaks ties identically.
        self.score.total_cmp(&other.score).then_with(|| self.key.cmp(&other.key))
    }
}
impl PartialOrd for GreedyEntry<'_> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Debug, Clone, serde::Deserialize)]
#[serde(default)]
pub struct ProteinConfig {
    /// A peptide counts as "observed" if its winning feature match has
    /// |rt_diff| within this window (minutes). Applied equally to targets and
    /// decoys, so the decoy-estimated hit rate stays calibrated.
    pub rt_window: f64,
    /// Additional ion-mobility consistency window (1/K0). When the data has IM,
    /// an observed peptide must also have |im_diff| within this window. IM is
    /// orthogonal to RT, so this cuts random coincidences sharply — the main
    /// reason ion mobility boosts MS1 depth. Ignored on non-IM data.
    pub im_window: f64,
    /// Recalibrated-mass window in sigma units. An observed peptide must match
    /// within this many sigma of the per-run mass error. Tight mass is the
    /// dominant background suppressor (ms1searchpy's mass recalibration).
    pub mass_z_window: f64,
    /// Target protein FDR the threshold optimization maximizes proteins at.
    pub fdr: f64,
    /// Optimize the RT/mass thresholds to maximize proteins at `fdr` (ms1searchpy's
    /// `worker` grid scan). When true, rt_window/mass_z_window are upper bounds for
    /// the scan; when false they are used as fixed cutoffs.
    pub optimize_thresholds: bool,
    /// Minimum feature shape quality (n_isotopes) for a peptide to count as
    /// observed. Real peptides sit on well-formed features; random matches scatter
    /// onto junk features, so this gate drops the random-match background
    /// asymmetrically. The scan tries [0, this] — 0 disables it. Inert when the
    /// input has no shape descriptors (feat_quality == 0 everywhere).
    pub max_feature_isotopes_gate: f64,
    /// Rescore proteins with a cross-validated Fisher LDA over protein-level
    /// features (binomial score, log#peptides, coverage, RT/mass consistency,
    /// mean peptide discriminant). Deterministic (closed-form, 3-fold CV by
    /// accession hash) so it does not reintroduce ms1searchpy-style ML variance.
    pub ml_rescore: bool,
}

impl Default for ProteinConfig {
    fn default() -> Self {
        ProteinConfig {
            rt_window: 3.0,
            im_window: 0.03,
            mass_z_window: 3.0,
            fdr: 0.01,
            optimize_thresholds: true,
            max_feature_isotopes_gate: 0.0,
            ml_rescore: false,
        }
    }
}

#[derive(Debug, Clone)]
pub struct ProteinId {
    pub accession: String,
    pub decoy: bool,
    pub observed_peptides: usize,
    pub theoretical_peptides: usize,
    pub coverage: f64,
    pub binomial_score: f64,
    /// Mean |rt_diff| of the protein's seen peptides (RT consistency). Real
    /// proteins cluster near 0; lucky decoy proteins sit near the window edge.
    pub mean_rt_diff: f64,
    /// Mean |mass_z| and mean per-peptide discriminant score of seen peptides
    /// (features for the optional LDA rescorer).
    pub mean_mass_z: f64,
    pub mean_pep_score: f64,
    pub q_value: f64,
}

/// Lanczos approximation to ln Γ(x).
fn ln_gamma(x: f64) -> f64 {
    const G: [f64; 8] = [
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7,
    ];
    if x < 0.5 {
        std::f64::consts::PI / ((std::f64::consts::PI * x).sin() * ln_gamma(1.0 - x).exp())
    } else {
        let x = x - 1.0;
        let mut a = 0.99999999999980993;
        let t = x + 7.5;
        for (i, &g) in G.iter().enumerate() {
            a += g / (x + (i as f64) + 1.0);
        }
        0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
    }
}

/// Continued fraction for the regularized incomplete beta (Numerical Recipes).
fn betacf(a: f64, b: f64, x: f64) -> f64 {
    let tiny = 1e-30;
    let qab = a + b;
    let qap = a + 1.0;
    let qam = a - 1.0;
    let mut c = 1.0;
    let mut d = 1.0 - qab * x / qap;
    if d.abs() < tiny {
        d = tiny;
    }
    d = 1.0 / d;
    let mut h = d;
    for m in 1..200 {
        let m = m as f64;
        let m2 = 2.0 * m;
        let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if d.abs() < tiny {
            d = tiny;
        }
        c = 1.0 + aa / c;
        if c.abs() < tiny {
            c = tiny;
        }
        d = 1.0 / d;
        h *= d * c;
        let aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if d.abs() < tiny {
            d = tiny;
        }
        c = 1.0 + aa / c;
        if c.abs() < tiny {
            c = tiny;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < 1e-12 {
            break;
        }
    }
    h
}

/// Regularized incomplete beta function I_x(a, b).
fn betai(a: f64, b: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let bt = (ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + a * x.ln() + b * (1.0 - x).ln()).exp();
    if x < (a + 1.0) / (a + b + 2.0) {
        bt * betacf(a, b, x) / a
    } else {
        1.0 - bt * betacf(b, a, 1.0 - x) / b
    }
}

/// Exact binomial upper tail P(X >= k) for X ~ Binomial(n, p), via the
/// regularized incomplete beta: P(X >= k) = I_p(k, n-k+1). This is ms1searchpy's
/// `binom.sf(k-1, n, p)`. The Poisson approximation it replaces is only valid for
/// small p; here p (the per-peptide random hit rate) is 5-20%, where Poisson is
/// ~1/(1-p) too wide and under-scores every protein.
pub(crate) fn binom_sf_ge_k(k: u32, n: u32, p: f64) -> f64 {
    if k == 0 {
        return 1.0;
    }
    if k as u64 > n as u64 {
        return 0.0;
    }
    betai(k as f64, (n - k + 1) as f64, p).clamp(1e-300, 1.0)
}

/// Continuous generalization of [`binom_sf_ge_k`]: the regularized incomplete beta
/// `I_p(k, n-k+1)` accepts a fractional `k`, so a WEIGHTED observed count (the sum of
/// per-peptide posterior probabilities) plugs straight into the same binomial tail
/// the discrete scorer uses. No Poisson approximation (which is too wide at the
/// p=5-20% random hit-rate). Reduces to `binom_sf_ge_k` exactly when `k` is integral.
pub(crate) fn binom_sf_ge_k_cont(k: f64, n: f64, p: f64) -> f64 {
    if k <= 0.0 {
        return 1.0;
    }
    if k >= n + 1.0 {
        return 0.0;
    }
    betai(k, (n - k + 1.0).max(1e-9), p).clamp(1e-300, 1.0)
}

/// Continuous-evidence protein scorer: the single-pass replacement for the
/// deviation-cut SWEEP. Instead of a hard `dev <= cut` gate (whose optimal cut is
/// found by re-running inference ~30x), each peptide carries a posterior weight
/// `w in [0,1]` = P(correct | dev) from the target/decoy `dev` distributions, and a
/// protein's evidence is the WEIGHTED sum of its (competed) peptides. Same machinery
/// as [`score_proteins_feature_compete`] otherwise: feature-level competition
/// (banning peaks), dynamic decoy background rate, greedy parsimony, target/decoy
/// pair competition + (d+1)/t FDR. `weight`: peptide index -> posterior.
pub(crate) fn score_proteins_weighted(
    observed: FxHashMap<(bool, &str), FxHashSet<usize>>,
    pep_feat: &FxHashMap<usize, Vec<usize>>,
    weight: &FxHashMap<usize, f64>,
    theoretical: &FxHashMap<(bool, &str), u32>,
    decoy_theo: u64,
) -> Vec<ProteinId> {
    let theo_of = |key: &(bool, &str)| theoretical.get(key).copied().unwrap_or(0) as usize;
    let w_of = |pi: usize| weight.get(&pi).copied().unwrap_or(0.0);
    let prot_peps = observed;
    let mut feat_peps: FxHashMap<usize, FxHashSet<usize>> = FxHashMap::default();
    let mut all_peps: FxHashSet<usize> = FxHashSet::default();
    for peps in prot_peps.values() {
        for &pi in peps {
            all_peps.insert(pi);
        }
    }
    let mut rem_feat: FxHashMap<usize, usize> = FxHashMap::default();
    for &pi in &all_peps {
        if let Some(fs) = pep_feat.get(&pi) {
            for &f in fs {
                feat_peps.entry(f).or_default().insert(pi);
            }
            rem_feat.insert(pi, fs.len().max(1));
        }
    }
    let mut pep_prots: FxHashMap<usize, Vec<(bool, &str)>> = FxHashMap::default();
    for (&key, peps) in &prot_peps {
        for &pi in peps {
            pep_prots.entry(pi).or_default().push(key);
        }
    }
    // Weighted observed evidence per protein (sum of member peptides' posteriors).
    let mut prot_peps = prot_peps;
    let mut obs_w: FxHashMap<(bool, &str), f64> =
        prot_peps.iter().map(|(&k, v)| (k, v.iter().map(|&pi| w_of(pi)).sum())).collect();
    let mut by_n: FxHashMap<usize, Vec<(bool, &str)>> = FxHashMap::default();
    for &k in prot_peps.keys() {
        let n = theo_of(&k);
        if n > 0 {
            by_n.entry(n).or_default().push(k);
        }
    }
    let decoy_theo_f = decoy_theo.max(1) as f64;
    let mut decoy_obs_sum: f64 = obs_w.iter().filter(|((d, _), _)| *d).map(|(_, &c)| c).sum();

    // Best (highest weighted evidence) protein of a given theoretical size; ties by
    // accession key for determinism.
    fn best_for_n<'a>(
        n: usize,
        by_n: &FxHashMap<usize, Vec<(bool, &'a str)>>,
        obs_w: &FxHashMap<(bool, &'a str), f64>,
    ) -> Option<(bool, &'a str)> {
        by_n.get(&n).and_then(|ks| {
            ks.iter()
                .filter(|k| obs_w.get(*k).copied().unwrap_or(0.0) > 1e-9)
                .max_by(|a, b| {
                    let (wa, wb) = (obs_w[*a], obs_w[*b]);
                    wa.partial_cmp(&wb).unwrap().then_with(|| b.cmp(a))
                })
                .copied()
        })
    }
    let mut best_per_n: FxHashMap<usize, (bool, &str)> = FxHashMap::default();
    for &n in by_n.keys() {
        if let Some(k) = best_for_n(n, &by_n, &obs_w) {
            best_per_n.insert(n, k);
        }
    }

    let mut proteins: Vec<ProteinId> = Vec::new();
    let (mut t, mut d) = (0usize, 0usize);
    let mut dead: FxHashSet<usize> = FxHashSet::default();
    let mut banned_feat: FxHashSet<usize> = FxHashSet::default();
    loop {
        let p = (decoy_obs_sum / decoy_theo_f).clamp(1e-9, 0.999);
        let mut best: Option<((bool, &str), f64, f64, usize)> = None;
        for (&n, &k) in &best_per_n {
            let kobs = obs_w.get(&k).copied().unwrap_or(0.0);
            if kobs <= 1e-9 {
                continue;
            }
            let s = -binom_sf_ge_k_cont(kobs, n as f64, p).log10();
            let better = match best {
                None => true,
                Some((bk, bs, _, _)) => s > bs || (s == bs && k < bk),
            };
            if better {
                best = Some((k, s, kobs, n));
            }
        }
        let (key, score, kobs, n) = match best {
            Some(x) => x,
            None => break,
        };
        let (decoy, acc) = key;
        if decoy {
            d += 1;
        } else {
            t += 1;
        }
        let q = if t == 0 { 1.0 } else { d as f64 / t as f64 };
        proteins.push(ProteinId {
            accession: acc.to_string(),
            decoy,
            observed_peptides: kobs.round() as usize,
            theoretical_peptides: n,
            coverage: if n > 0 { kobs / n as f64 } else { 0.0 },
            binomial_score: score,
            mean_rt_diff: 0.0,
            mean_mass_z: 0.0,
            mean_pep_score: 0.0,
            q_value: q,
        });
        if t > 20 && q > 0.5 {
            break;
        }
        let won_peps: Vec<usize> = prot_peps
            .get(&key)
            .map(|s| s.iter().copied().collect())
            .unwrap_or_default();
        let mut to_kill: Vec<usize> = Vec::new();
        for pi in &won_peps {
            match pep_feat.get(pi) {
                Some(fs) if !fs.is_empty() => {
                    for &f in fs {
                        if banned_feat.insert(f) {
                            if let Some(ps) = feat_peps.get(&f) {
                                for &pj in ps {
                                    if let Some(r) = rem_feat.get_mut(&pj) {
                                        if *r > 0 {
                                            *r -= 1;
                                            if *r == 0 {
                                                to_kill.push(pj);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                _ => to_kill.push(*pi),
            }
        }
        let mut dirty_n: FxHashSet<usize> = FxHashSet::default();
        for pj in to_kill {
            if !dead.insert(pj) {
                continue;
            }
            let wj = w_of(pj);
            if let Some(keys) = pep_prots.get(&pj) {
                for k2 in keys.clone() {
                    if let Some(s) = prot_peps.get_mut(&k2) {
                        if s.remove(&pj) {
                            if let Some(c) = obs_w.get_mut(&k2) {
                                *c -= wj;
                            }
                            if k2.0 {
                                decoy_obs_sum -= wj;
                            }
                            dirty_n.insert(theo_of(&k2));
                        }
                    }
                }
            }
        }
        obs_w.insert(key, 0.0);
        if let Some(s) = prot_peps.get_mut(&key) {
            s.clear();
        }
        dirty_n.insert(n);
        for dn in dirty_n {
            match best_for_n(dn, &by_n, &obs_w) {
                Some(k) => {
                    best_per_n.insert(dn, k);
                }
                None => {
                    best_per_n.remove(&dn);
                }
            }
        }
    }
    proteins.sort_by(|a, b| {
        b.binomial_score
            .partial_cmp(&a.binomial_score)
            .unwrap()
            .then_with(|| a.accession.cmp(&b.accession))
            .then(a.decoy.cmp(&b.decoy))
    });
    let score_by: FxHashMap<(bool, String), f64> = proteins
        .iter()
        .map(|p| ((p.decoy, p.accession.clone()), p.binomial_score))
        .collect();
    proteins.retain(|p| match score_by.get(&(!p.decoy, p.accession.clone())) {
        Some(&cs) => p.binomial_score > cs || (p.binomial_score == cs && !p.decoy),
        None => true,
    });
    let plus1 = if std::env::var("UNO_NO_FDR_PLUS1").is_ok() { 0.0 } else { 1.0 };
    let (mut t, mut d) = (0usize, 0usize);
    for p in proteins.iter_mut() {
        if p.decoy {
            d += 1;
        } else {
            t += 1;
        }
        p.q_value = if t == 0 { 1.0 } else { (d as f64 + plus1) / t as f64 };
    }
    let mut run = f64::INFINITY;
    for i in (0..proteins.len()).rev() {
        run = run.min(proteins[i].q_value);
        proteins[i].q_value = run.min(1.0);
    }
    proteins
}

/// ms1searchpy `final_iteration` core, reusable: greedy parsimony (assign each
/// peptide to its single best-scoring protein, removing it from competitors and
/// rescoring them) + target/decoy pair deletion + FDR walk. This is what keeps the
/// protein FDR honest — counting a shared/homologous peptide for every protein it
/// touches (as a plain per-protein binomial does) inflates spurious proteins (absent
/// organisms, decoys) and makes the target/decoy FDR optimistic.
///
/// `observed`: protein key (is_decoy, accession) -> set of matched peptide indices.
/// `theoretical`: protein key -> total tryptic peptide count (the binomial `n`).
/// `decoy_theo`: total decoy theoretical peptides; `decoy_seen`: distinct matched
/// decoy peptides (together set the random hit-rate `p`). `fdr`: target FDR (for the
/// early-stop). uno's decoys share their target's accession (differ only in the
/// decoy flag), so pair deletion compares (decoy, acc) against (!decoy, acc).
pub(crate) fn score_proteins_greedy(
    observed: FxHashMap<(bool, &str), FxHashSet<usize>>,
    theoretical: &FxHashMap<(bool, &str), u32>,
    decoy_theo: u64,
    decoy_seen: usize,
    fdr: f64,
) -> Vec<ProteinId> {
    // Decouple the two jobs (UNO_PROT_PERMUTE): decoys estimate ONLY the background
    // per-peptide hit rate p; the per-protein binomial gives a p-value
    // (= 10^-score = binom.sf(k-1, n, p)); Benjamini-Hochberg over those p-values
    // sets the FDR. No target/decoy competition in the FDR step — decoys feed in
    // solely through p. Valid as long as the binomial p-values are ~uniform under
    // the null (checked against the yeast entrapment).
    let permute = std::env::var("UNO_PROT_PERMUTE").is_ok();
    let hit_rate = if decoy_theo == 0 {
        1e-9
    } else {
        (decoy_seen as f64 / decoy_theo as f64).max(1e-9)
    };
    if std::env::var("UNO_LOG_HITRATE").is_ok() {
        log::info!("score_proteins_greedy: hit_rate(p)={:.4} decoy_seen={} decoy_theo={}", hit_rate, decoy_seen, decoy_theo);
    }
    let mut score_cache: FxHashMap<(u32, usize), f64> = FxHashMap::default();
    let mut score_of = |key: &(bool, &str), n_obs: usize| -> f64 {
        if n_obs == 0 {
            return 0.0;
        }
        let theo = theoretical.get(key).copied().unwrap_or(0);
        *score_cache
            .entry((theo, n_obs))
            .or_insert_with(|| -binom_sf_ge_k(n_obs as u32, theo, hit_rate).log10())
    };
    let mut prot_peps = observed;
    let mut pep_prots: FxHashMap<usize, Vec<(bool, &str)>> = FxHashMap::default();
    for (&key, peps) in &prot_peps {
        for &pi in peps {
            pep_prots.entry(pi).or_default().push(key);
        }
    }

    // ms1searchpy-style DYNAMIC background rate `p` (worker/final_iteration). The
    // static `hit_rate` above freezes the random-match probability at its initial
    // value; ms1searchpy instead RE-ESTIMATES p = (current decoy matches)/(decoy
    // theoretical) after every greedy pick. As real proteins consume their
    // peptides, decoy matches shrink, p drops, and the next real protein's binomial
    // becomes more significant -- a positive feedback that digs much deeper. This is
    // the engine difference behind ms1searchpy's ~2.8x depth -- BUT only when paired
    // with ms1searchpy's FEATURE-level competition (each MS1 peak supports one
    // protein). uno competes at the PEPTIDE level here (the peak->peptide identity is
    // already gone), so a real protein winning barely shrinks the decoy pool and the
    // virtuous cycle never fires: measured -12% honest depth (660 -> 580). Kept behind
    // UNO_DYNAMIC_P=1 for further work; OFF by default until feature-level competition
    // is wired through.
    if std::env::var("UNO_DYNAMIC_P").is_ok() {
        let theo_of = |key: &(bool, &str)| theoretical.get(key).copied().unwrap_or(0) as usize;
        // Current observed count per protein, and proteins grouped by theoretical
        // size n (constant). For a fixed (n, p) the binomial is monotincreasing in
        // the observed count k, so only the max-k protein per n can ever be the
        // global best -- we score one representative per n, like ms1searchpy.
        let mut obs: FxHashMap<(bool, &str), usize> =
            prot_peps.iter().map(|(&k, v)| (k, v.len())).collect();
        let mut by_n: FxHashMap<usize, Vec<(bool, &str)>> = FxHashMap::default();
        for &k in prot_peps.keys() {
            let n = theo_of(&k);
            if n > 0 {
                by_n.entry(n).or_default().push(k);
            }
        }
        let decoy_theo_f = decoy_theo.max(1) as f64;
        let mut decoy_obs_sum: f64 =
            obs.iter().filter(|((d, _), _)| *d).map(|(_, &c)| c as f64).sum();
        fn best_for_n<'a>(
            n: usize,
            by_n: &FxHashMap<usize, Vec<(bool, &'a str)>>,
            obs: &FxHashMap<(bool, &'a str), usize>,
        ) -> Option<(bool, &'a str)> {
            by_n.get(&n).and_then(|ks| {
                ks.iter()
                    .filter(|k| obs.get(*k).copied().unwrap_or(0) > 0)
                    .max_by_key(|k| (obs[*k], std::cmp::Reverse(**k)))
                    .copied()
            })
        }
        let mut best_per_n: FxHashMap<usize, (bool, &str)> = FxHashMap::default();
        for &n in by_n.keys() {
            if let Some(k) = best_for_n(n, &by_n, &obs) {
                best_per_n.insert(n, k);
            }
        }

        let mut proteins: Vec<ProteinId> = Vec::new();
        let (mut t, mut d) = (0usize, 0usize);
        loop {
            let p = (decoy_obs_sum / decoy_theo_f).clamp(1e-9, 0.999);
            // Pick the globally most-significant representative under the current p.
            let mut best: Option<((bool, &str), f64, usize, usize)> = None;
            for (&n, &k) in &best_per_n {
                let kobs = obs.get(&k).copied().unwrap_or(0);
                if kobs == 0 {
                    continue;
                }
                let s = -binom_sf_ge_k(kobs as u32, n as u32, p).log10();
                let better = match best {
                    None => true,
                    Some((bk, bs, _, _)) => s > bs || (s == bs && k < bk),
                };
                if better {
                    best = Some((k, s, kobs, n));
                }
            }
            let (key, score, kobs, n) = match best {
                Some(x) => x,
                None => break,
            };
            let (decoy, acc) = key;
            if decoy {
                d += 1;
            } else {
                t += 1;
            }
            let q = if t == 0 { 1.0 } else { d as f64 / t as f64 };
            proteins.push(ProteinId {
                accession: acc.to_string(),
                decoy,
                observed_peptides: kobs,
                theoretical_peptides: n,
                coverage: if n > 0 { kobs as f64 / n as f64 } else { 0.0 },
                binomial_score: score,
                mean_rt_diff: 0.0,
                mean_mass_z: 0.0,
                mean_pep_score: 0.0,
                q_value: q,
            });
            if t > 20 && q > (fdr * 25.0).max(0.5) {
                break;
            }
            // Consume the winner's peptides (parsimony): remove them from every
            // protein, decrement observed counts, and shrink the decoy pool so p
            // drops for the next pick. Track which n-groups changed.
            let won: Vec<usize> = prot_peps
                .get(&key)
                .map(|s| s.iter().copied().collect())
                .unwrap_or_default();
            let mut dirty_n: FxHashSet<usize> = FxHashSet::default();
            for pi in won {
                if let Some(keys) = pep_prots.get(&pi) {
                    for k2 in keys.clone() {
                        if let Some(s) = prot_peps.get_mut(&k2) {
                            if s.remove(&pi) {
                                if let Some(c) = obs.get_mut(&k2) {
                                    *c -= 1;
                                }
                                if k2.0 {
                                    decoy_obs_sum -= 1.0;
                                }
                                dirty_n.insert(theo_of(&k2));
                            }
                        }
                    }
                }
            }
            // The winner itself can never be picked again.
            obs.insert(key, 0);
            dirty_n.insert(n);
            for dn in dirty_n {
                match best_for_n(dn, &by_n, &obs) {
                    Some(k) => {
                        best_per_n.insert(dn, k);
                    }
                    None => {
                        best_per_n.remove(&dn);
                    }
                }
            }
        }
        proteins.sort_by(|a, b| {
            b.binomial_score
                .partial_cmp(&a.binomial_score)
                .unwrap()
                .then_with(|| a.accession.cmp(&b.accession))
                .then(a.decoy.cmp(&b.decoy))
        });
        // FDR: running decoy/target ratio over the score-sorted list, made monotone.
        let (mut t, mut d) = (0usize, 0usize);
        for p in proteins.iter_mut() {
            if p.decoy {
                d += 1;
            } else {
                t += 1;
            }
            p.q_value = if t == 0 { 1.0 } else { d as f64 / t as f64 };
        }
        let mut run = f64::INFINITY;
        for i in (0..proteins.len()).rev() {
            run = run.min(proteins[i].q_value);
            proteins[i].q_value = run.min(1.0);
        }
        return proteins;
    }

    let mut heap: BinaryHeap<GreedyEntry> = prot_peps
        .iter()
        .filter(|(_, peps)| !peps.is_empty())
        .map(|(&key, peps)| GreedyEntry { score: score_of(&key, peps.len()), key, n_obs: peps.len() })
        .collect();
    let mut proteins: Vec<ProteinId> = Vec::new();
    let (mut t, mut d) = (0usize, 0usize);
    loop {
        let (key, obs, score) = loop {
            let top = match heap.pop() {
                Some(e) => e,
                None => break (None, 0usize, 0.0f64),
            };
            match prot_peps.get(&top.key) {
                Some(peps) if !peps.is_empty() && peps.len() == top.n_obs => {
                    break (Some(top.key), top.n_obs, top.score);
                }
                _ => continue,
            }
        };
        let key = match key {
            Some(k) => k,
            None => break,
        };
        let (decoy, acc) = key;
        let theo = theoretical.get(&key).copied().unwrap_or(0) as usize;
        if decoy {
            d += 1;
        } else {
            t += 1;
        }
        let q = if t == 0 { 1.0 } else { d as f64 / t as f64 };
        proteins.push(ProteinId {
            accession: acc.to_string(),
            decoy,
            observed_peptides: obs,
            theoretical_peptides: theo,
            coverage: if theo > 0 { obs as f64 / theo as f64 } else { 0.0 },
            binomial_score: score,
            mean_rt_diff: 0.0,
            mean_mass_z: 0.0,
            mean_pep_score: 0.0,
            q_value: q,
        });
        if t > 20 && q > (fdr * 25.0).max(0.5) {
            break;
        }
        let won: Vec<usize> = prot_peps.get(&key).map(|s| s.iter().copied().collect()).unwrap_or_default();
        let mut affected: FxHashSet<(bool, &str)> = FxHashSet::default();
        for pi in won {
            if let Some(keys) = pep_prots.get(&pi) {
                for k in keys.clone() {
                    if let Some(s) = prot_peps.get_mut(&k) {
                        s.remove(&pi);
                        affected.insert(k);
                    }
                }
            }
        }
        for k in affected {
            if let Some(peps) = prot_peps.get(&k) {
                if !peps.is_empty() {
                    heap.push(GreedyEntry { score: score_of(&k, peps.len()), key: k, n_obs: peps.len() });
                }
            }
        }
    }
    proteins.sort_by(|a, b| {
        b.binomial_score
            .partial_cmp(&a.binomial_score)
            .unwrap()
            .then_with(|| a.accession.cmp(&b.accession))
            .then(a.decoy.cmp(&b.decoy))
    });
    if permute {
        // Decoy-free FDR: per-protein binomial p-value = 10^(-score); Benjamini-
        // Hochberg over the TARGET proteins (sorted by score desc = p-value asc).
        // Decoys are kept for diagnostics but get q=1 (they don't define the FDR).
        let m = proteins.iter().filter(|p| !p.decoy).count().max(1) as f64;
        let mut rank = 0usize;
        for p in proteins.iter_mut() {
            if p.decoy {
                p.q_value = 1.0;
                continue;
            }
            rank += 1;
            let pv = 10f64.powf(-p.binomial_score);
            p.q_value = (pv * m / rank as f64).min(1.0);
        }
        // Enforce monotone non-increasing q with rank (standard BH step-up).
        let mut run = 1.0f64;
        for p in proteins.iter_mut().rev() {
            if p.decoy {
                continue;
            }
            run = run.min(p.q_value);
            p.q_value = run;
        }
        return proteins;
    }
    // Target/decoy pair deletion (uno scheme: same accession, opposite decoy flag).
    let score_by: FxHashMap<(bool, String), f64> = proteins
        .iter()
        .map(|p| ((p.decoy, p.accession.clone()), p.binomial_score))
        .collect();
    proteins.retain(|p| match score_by.get(&(!p.decoy, p.accession.clone())) {
        Some(&cs) => p.binomial_score > cs || (p.binomial_score == cs && !p.decoy),
        None => true,
    });
    let (mut t, mut d) = (0usize, 0usize);
    for p in proteins.iter_mut() {
        if p.decoy {
            d += 1;
        } else {
            t += 1;
        }
        p.q_value = if t == 0 { 1.0 } else { d as f64 / t as f64 };
    }
    let mut run = f64::INFINITY;
    for i in (0..proteins.len()).rev() {
        run = run.min(proteins[i].q_value);
        proteins[i].q_value = run.min(1.0);
    }
    proteins
}

/// Faithful port of ms1searchpy's `worker`/`final_iteration` greedy protein
/// inference: FEATURE-level competition (each MS1 peak supports exactly one
/// protein) + a DYNAMICALLY re-estimated random-match probability `p`. `pep_feat`
/// maps each candidate peptide to the MS1 feature (peak id) backing its best match.
/// When a protein is accepted it claims those peaks; every other peptide sitting on
/// a claimed peak dies, which shrinks the decoy pool so `p` drops and the next real
/// protein scores higher — the feedback uno's peptide-level greedy can't produce.
pub(crate) fn score_proteins_feature_compete(
    observed: FxHashMap<(bool, &str), FxHashSet<usize>>,
    pep_feat: &FxHashMap<usize, Vec<usize>>,
    theoretical: &FxHashMap<(bool, &str), u32>,
    decoy_theo: u64,
    fdr: f64,
) -> Vec<ProteinId> {
    let theo_of = |key: &(bool, &str)| theoretical.get(key).copied().unwrap_or(0) as usize;
    let mut prot_peps = observed;
    // feature (peak id) -> peptides that matched it (a peptide can sit on several
    // peaks: charge states / scans). A peptide dies only when ALL its peaks are
    // claimed (ms1searchpy's banned_dict), so it survives competition the way ms1's
    // multi-charge evidence does -- retaining ~7 peptides/protein where the old
    // single-best-PFM mapping kept ~5.
    let mut feat_peps: FxHashMap<usize, FxHashSet<usize>> = FxHashMap::default();
    let mut all_peps: FxHashSet<usize> = FxHashSet::default();
    for peps in prot_peps.values() {
        for &pi in peps {
            all_peps.insert(pi);
        }
    }
    // Remaining (unbanned) peak count per peptide; peptide is dead at 0.
    let mut rem_feat: FxHashMap<usize, usize> = FxHashMap::default();
    for &pi in &all_peps {
        if let Some(fs) = pep_feat.get(&pi) {
            for &f in fs {
                feat_peps.entry(f).or_default().insert(pi);
            }
            rem_feat.insert(pi, fs.len().max(1));
        }
    }
    let mut pep_prots: FxHashMap<usize, Vec<(bool, &str)>> = FxHashMap::default();
    for (&key, peps) in &prot_peps {
        for &pi in peps {
            pep_prots.entry(pi).or_default().push(key);
        }
    }
    let mut obs: FxHashMap<(bool, &str), usize> =
        prot_peps.iter().map(|(&k, v)| (k, v.len())).collect();
    let decoy_theo_f = decoy_theo.max(1) as f64;
    let mut decoy_obs_sum: f64 =
        obs.iter().filter(|((d, _), _)| *d).map(|(_, &c)| c as f64).sum();

    // Per-theoretical-count (bucket) selection of the max-observed protein.
    //
    // `obs[k]` only ever DECREASES over the greedy (peptides get removed as their
    // peaks are claimed; the winner is zeroed). Exploit that monotonicity: keep,
    // per bucket `n`, a lazy max-heap of `(obs_snapshot, Reverse(key))`. The head
    // is the current best UNLESS its snapshot is stale (obs decremented since it
    // was pushed); on read we pop stale/zeroed heads, re-pushing decremented ones
    // at their current count, until the head's snapshot matches its live `obs`.
    // Because `obs` is non-increasing and every live key holds exactly one heap
    // entry with `snapshot >= live obs`, the first head with `snapshot == live`
    // is the exact `max_by_key(obs, Reverse(key))` the old full-bucket rescan
    // returned — byte-identical selection, but amortized O(log) per read instead
    // of O(bucket_size). (The old rescan cost O(iters * bucket_size) ~ billions of
    // scans; this bounds re-pushes by the total peptide-removal count.)
    type BucketHeap<'a> = std::collections::BinaryHeap<(usize, std::cmp::Reverse<(bool, &'a str)>)>;
    #[inline]
    fn peek_best<'a>(
        heap: &mut BucketHeap<'a>,
        obs: &FxHashMap<(bool, &'a str), usize>,
    ) -> Option<(usize, (bool, &'a str))> {
        while let Some(&(snap, std::cmp::Reverse(key))) = heap.peek() {
            let cur = obs.get(&key).copied().unwrap_or(0);
            if cur == 0 {
                heap.pop();
                continue;
            }
            if cur < snap {
                heap.pop();
                heap.push((cur, std::cmp::Reverse(key)));
                continue;
            }
            return Some((cur, key));
        }
        None
    }
    let mut heaps: FxHashMap<usize, BucketHeap> = FxHashMap::default();
    for (&k, &c) in &obs {
        let n = theo_of(&k);
        if n > 0 && c > 0 {
            heaps.entry(n).or_default().push((c, std::cmp::Reverse(k)));
        }
    }

    let mut proteins: Vec<ProteinId> = Vec::new();
    let (mut t, mut d) = (0usize, 0usize);
    let mut dead: FxHashSet<usize> = FxHashSet::default();
    let mut banned_feat: FxHashSet<usize> = FxHashSet::default();
    let mut empty_buckets: Vec<usize> = Vec::new();
    loop {
        let p = (decoy_obs_sum / decoy_theo_f).clamp(1e-9, 0.999);
        let mut best: Option<((bool, &str), f64, usize, usize)> = None;
        empty_buckets.clear();
        for (&n, heap) in heaps.iter_mut() {
            let (kobs, k) = match peek_best(heap, &obs) {
                Some(x) => x,
                None => {
                    empty_buckets.push(n);
                    continue;
                }
            };
            let s = -binom_sf_ge_k(kobs as u32, n as u32, p).log10();
            let better = match best {
                None => true,
                Some((bk, bs, _, _)) => s > bs || (s == bs && k < bk),
            };
            if better {
                best = Some((k, s, kobs, n));
            }
        }
        for n in empty_buckets.drain(..) {
            heaps.remove(&n);
        }
        let (key, score, kobs, n) = match best {
            Some(x) => x,
            None => break,
        };
        let (decoy, acc) = key;
        if decoy {
            d += 1;
        } else {
            t += 1;
        }
        let q = if t == 0 { 1.0 } else { d as f64 / t as f64 };
        proteins.push(ProteinId {
            accession: acc.to_string(),
            decoy,
            observed_peptides: kobs,
            theoretical_peptides: n,
            coverage: if n > 0 { kobs as f64 / n as f64 } else { 0.0 },
            binomial_score: score,
            mean_rt_diff: 0.0,
            mean_mass_z: 0.0,
            mean_pep_score: 0.0,
            q_value: q,
        });
        if t > 20 && q > (fdr * 25.0).max(0.5) {
            break;
        }
        // The winner claims the peaks of all its peptides; every peptide on a
        // claimed peak dies (removed from every protein), shrinking the pool.
        let won_peps: Vec<usize> = prot_peps
            .get(&key)
            .map(|s| s.iter().copied().collect())
            .unwrap_or_default();
        // The winner claims every peak under its peptides. Banning a peak decrements
        // each peptide on it; a peptide dies only when its LAST peak is claimed.
        let mut to_kill: Vec<usize> = Vec::new();
        for pi in &won_peps {
            match pep_feat.get(pi) {
                Some(fs) if !fs.is_empty() => {
                    for &f in fs {
                        if banned_feat.insert(f) {
                            if let Some(ps) = feat_peps.get(&f) {
                                for &pj in ps {
                                    if let Some(r) = rem_feat.get_mut(&pj) {
                                        if *r > 0 {
                                            *r -= 1;
                                            if *r == 0 {
                                                to_kill.push(pj);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                _ => to_kill.push(*pi),
            }
        }
        // Removing peptides only DECREMENTS `obs`; the per-bucket lazy heaps pick up
        // the new counts on their next read (no eager bucket recompute needed).
        for pj in to_kill {
            if !dead.insert(pj) {
                continue;
            }
            if let Some(keys) = pep_prots.get(&pj) {
                for k2 in keys.clone() {
                    if let Some(s) = prot_peps.get_mut(&k2) {
                        if s.remove(&pj) {
                            if let Some(c) = obs.get_mut(&k2) {
                                *c -= 1;
                            }
                            if k2.0 {
                                decoy_obs_sum -= 1.0;
                            }
                        }
                    }
                }
            }
        }
        // Winner can never be picked again (its peaks are all consumed).
        obs.insert(key, 0);
        if let Some(s) = prot_peps.get_mut(&key) {
            s.clear();
        }
    }
    proteins.sort_by(|a, b| {
        b.binomial_score
            .partial_cmp(&a.binomial_score)
            .unwrap()
            .then_with(|| a.accession.cmp(&b.accession))
            .then(a.decoy.cmp(&b.decoy))
    });
    // ms1searchpy's protein FDR (final_iteration:401-413): (1) target/decoy PAIR
    // competition -- for each accession keep only the better-scoring of X / DECOY_X;
    // (2) FDR with the conservative +1 decoy correction: q = (decoys+1)/targets
    // (pyteomics aux.filter correction=1). My first cut used a bare decoys/targets
    // with no pair deletion, which ran ~1.5x optimistic.
    let score_by: FxHashMap<(bool, String), f64> = proteins
        .iter()
        .map(|p| ((p.decoy, p.accession.clone()), p.binomial_score))
        .collect();
    proteins.retain(|p| match score_by.get(&(!p.decoy, p.accession.clone())) {
        Some(&cs) => p.binomial_score > cs || (p.binomial_score == cs && !p.decoy),
        None => true,
    });
    // Target/decoy ratio after pair competition, with ms1searchpy's conservative +1
    // phantom-decoy correction (pyteomics correction=1) ON by default. Dropping it
    // (UNO_NO_FDR_PLUS1=1) raised the reported count slightly but did NOT improve
    // honest entrapment depth -- it just reports nearer the (still-conservative)
    // decoy-estimated FDR. Keep the honest default.
    let plus1 = if std::env::var("UNO_NO_FDR_PLUS1").is_ok() { 0.0 } else { 1.0 };
    let (mut t, mut d) = (0usize, 0usize);
    for p in proteins.iter_mut() {
        if p.decoy {
            d += 1;
        } else {
            t += 1;
        }
        p.q_value = if t == 0 { 1.0 } else { (d as f64 + plus1) / t as f64 };
    }
    let mut run = f64::INFINITY;
    for i in (0..proteins.len()).rev() {
        run = run.min(proteins[i].q_value);
        proteins[i].q_value = run.min(1.0);
    }
    proteins
}
