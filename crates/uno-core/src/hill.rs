//! Peptide-centric second pass over koth **hills** (raw single-m/z XIC traces,
//! the layer below deisotoped features).
//!
//! The feature pass only sees peptides that survived deisotoping. Low-abundance
//! peptides whose isotopes fell below the feature finder's threshold never form
//! a feature even though their hills exist. This pass goes back to the hills and
//! asks, per peptide (targets AND uno's per-protein-pool decoys), "is there a
//! well-formed isotope cluster at my predicted m/z and predicted RT that looks
//! like ME?" — then target/decoy competition turns the fit scores into q-values.
//!
//! For MS1-only ID the discriminators are NOT "is there a hill" (decoys find
//! hills as readily as targets — mass alone separates nothing). They are:
//!   1. **predicted-RT constraint** — the cluster must sit at the peptide's
//!      on-data-calibrated predicted RT (from uno's own CNN);
//!   2. **isotope-ratio fit** — observed mono/+1/+2 intensities vs the averagine
//!      envelope for this mass (cosine similarity);
//!   3. **coelution** — the mono and +1 hills must rise and fall together
//!      (Pearson correlation of their per-scan intensity profiles);
//!   4. **charge-consistent spacing** — isotopes at exactly 1.00286/z.
//! A random mass coincidence rarely satisfies all four; a real peptide does.

use crate::database::PeptideDb;
use anyhow::{Context, Result};
use std::path::Path;

const PROTON: f64 = 1.007_276;
const NEUTRON: f64 = 1.002_866_5;
/// Mean number of 13C atoms for a peptide of neutral mass `m` (averagine:
/// 4.9384 C per 111.1254-Da residue, 13C natural abundance 0.0107).
fn averagine_lambda(m: f64) -> f64 {
    (m / 111.1254) * 4.9384 * 0.0107
}

/// One koth hill: a single-m/z intensity trace over retention time, with its
/// per-scan intensity profile for coelution scoring.
#[derive(Debug, Clone)]
pub struct Hill {
    pub mz: f64,
    pub rt: f64,
    pub rt_start: f64,
    pub rt_end: f64,
    pub scan_start: i64,
    pub scan_end: i64,
    pub intensity: f64,
    pub hill_score: f64,
    pub n_scans: f64,
    /// Per-scan intensities, indexed from `scan_start` to `scan_end`.
    pub profile: Vec<f64>,
}

impl Hill {
    /// Pearson correlation of this hill's profile with `other`'s over their
    /// overlapping absolute scan range. Returns 0 if they barely overlap.
    fn coelution(&self, other: &Hill) -> f64 {
        let lo = self.scan_start.max(other.scan_start);
        let hi = self.scan_end.min(other.scan_end);
        if hi - lo < 2 {
            return 0.0;
        }
        let mut a = Vec::new();
        let mut b = Vec::new();
        for s in lo..=hi {
            let ai = (s - self.scan_start) as usize;
            let bi = (s - other.scan_start) as usize;
            if ai < self.profile.len() && bi < other.profile.len() {
                a.push(self.profile[ai]);
                b.push(other.profile[bi]);
            }
        }
        pearson(&a, &b)
    }
}

fn pearson(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len().min(b.len());
    if n < 2 {
        return 0.0;
    }
    let ma = a[..n].iter().sum::<f64>() / n as f64;
    let mb = b[..n].iter().sum::<f64>() / n as f64;
    let mut num = 0.0;
    let mut da = 0.0;
    let mut db = 0.0;
    for i in 0..n {
        let x = a[i] - ma;
        let y = b[i] - mb;
        num += x * y;
        da += x * x;
        db += y * y;
    }
    if da <= 0.0 || db <= 0.0 {
        return 0.0;
    }
    num / (da.sqrt() * db.sqrt())
}

/// Number of isotope peaks in the theoretical envelope we predict + look for.
const MAX_ISO: usize = 5;

/// Theoretical averagine isotope envelope (relative intensities, L2-normalized)
/// for a peptide of neutral mass `mass`: Poisson(λ) on the averagine 13C count.
fn averagine_envelope(mass: f64) -> [f64; MAX_ISO] {
    let lambda = averagine_lambda(mass);
    let mut e = [0.0f64; MAX_ISO];
    let mut term = 1.0f64; // λ^k / k!
    for k in 0..MAX_ISO {
        if k > 0 {
            term *= lambda / k as f64;
        }
        e[k] = term;
    }
    let n: f64 = e.iter().map(|x| x * x).sum::<f64>().sqrt();
    if n > 0.0 {
        for x in e.iter_mut() {
            *x /= n;
        }
    }
    e
}

/// Cosine similarity between an observed and a theoretical envelope over the same
/// isotope positions. 1.0 = perfect match.
fn envelope_cosine(observed: &[f64], theo: &[f64]) -> f64 {
    let dot: f64 = observed.iter().zip(theo).map(|(o, e)| o * e).sum();
    let no: f64 = observed.iter().map(|o| o * o).sum::<f64>().sqrt();
    let ne: f64 = theo.iter().map(|e| e * e).sum::<f64>().sqrt();
    if no <= 0.0 || ne <= 0.0 {
        0.0
    } else {
        dot / (no * ne)
    }
}

/// Hills sorted by m/z, for binary-search lookup within a ppm window.
pub struct HillIndex {
    hills: Vec<Hill>,
    mz: Vec<f64>, // parallel, ascending
}

impl HillIndex {
    pub fn new(mut hills: Vec<Hill>) -> Self {
        hills.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap());
        let mz = hills.iter().map(|h| h.mz).collect();
        HillIndex { hills, mz }
    }
    pub fn len(&self) -> usize {
        self.hills.len()
    }
    pub fn is_empty(&self) -> bool {
        self.hills.is_empty()
    }
    fn window(&self, mz: f64, ppm: f64) -> std::ops::Range<usize> {
        let lo = mz * (1.0 - ppm * 1e-6);
        let hi = mz * (1.0 + ppm * 1e-6);
        let a = self.mz.partition_point(|&x| x < lo);
        let b = self.mz.partition_point(|&x| x <= hi);
        a..b
    }
    /// Best hill at `mz` (within ppm) that coelutes with `ref_hill`: highest
    /// coelution correlation. Returns (index, correlation).
    fn best_coeluting(&self, mz: f64, ppm: f64, ref_hill: &Hill) -> Option<(usize, f64)> {
        let mut best: Option<(usize, f64)> = None;
        for i in self.window(mz, ppm) {
            let r = self.hills[i].coelution(ref_hill);
            if best.map_or(true, |(_, br)| r > br) {
                best = Some((i, r));
            }
        }
        best
    }
}

/// Load koth `hills.tsv`. Requires `mz, rt, rt_start, rt_end, scan_start,
/// scan_end, intensity_max, hill_score, n_scans, intensity_profile`.
pub fn load_hills(path: impl AsRef<Path>) -> Result<Vec<Hill>> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path.as_ref())
        .with_context(|| format!("opening hills {}", path.as_ref().display()))?;
    let header = rdr.headers()?.clone();
    let col = |name: &str| -> Result<usize> {
        header
            .iter()
            .position(|h| h == name)
            .with_context(|| format!("hills file missing column `{name}`"))
    };
    let (c_mz, c_rt, c_rs, c_re, c_ss, c_se, c_int, c_sc, c_ns, c_prof) = (
        col("mz")?,
        col("rt")?,
        col("rt_start")?,
        col("rt_end")?,
        col("scan_start")?,
        col("scan_end")?,
        col("intensity_max")?,
        col("hill_score")?,
        col("n_scans")?,
        col("intensity_profile")?,
    );
    let mut out = Vec::new();
    for rec in rdr.records() {
        let rec = rec?;
        let gf = |i: usize| rec.get(i).unwrap_or("").parse::<f64>().unwrap_or(f64::NAN);
        let gi = |i: usize| rec.get(i).unwrap_or("").parse::<i64>().unwrap_or(0);
        let mz = gf(c_mz);
        if !mz.is_finite() {
            continue;
        }
        let profile: Vec<f64> = rec
            .get(c_prof)
            .unwrap_or("")
            .trim_matches(|c| c == '[' || c == ']')
            .split(',')
            .filter_map(|s| s.trim().parse::<f64>().ok())
            .collect();
        out.push(Hill {
            mz,
            rt: gf(c_rt),
            rt_start: gf(c_rs),
            rt_end: gf(c_re),
            scan_start: gi(c_ss),
            scan_end: gi(c_se),
            intensity: gf(c_int),
            hill_score: gf(c_sc),
            n_scans: gf(c_ns),
            profile,
        });
    }
    Ok(out)
}

#[derive(Debug, Clone)]
pub struct HillConfig {
    pub ppm: f64,
    pub charges: Vec<u8>,
    /// RT search half-window, in units of the confident-set RT sigma. Kept tight:
    /// a peptide is only sought at its predicted RT.
    pub rt_window_sigma: f64,
    /// Minimum coeluting isotopes (beyond the mono) required to count a peptide.
    pub min_isotopes: u32,
    /// Score weights for the envelope-fit cosine and the coelution correlation —
    /// the real MS1-only discriminators.
    pub w_envelope: f64,
    pub w_coelution: f64,
    /// Reject a cluster whose +1 coelution correlation is below this (a real
    /// isotope pair coelutes tightly; a random pair does not).
    pub min_coelution: f64,
    pub fdr: f64,
}

impl Default for HillConfig {
    fn default() -> Self {
        HillConfig {
            ppm: 10.0,
            charges: vec![2, 3, 4],
            rt_window_sigma: 3.0,
            min_isotopes: 1,
            w_envelope: 4.0,
            w_coelution: 4.0,
            min_coelution: 0.5,
            fdr: 0.01,
        }
    }
}

/// One peptide's best-fitting hill cluster: the per-peptide evidence handed to
/// uno's protein binomial (exactly like a feature-pass `PeptideMatch`). No
/// peptide-level FDR here — the target/decoy separation emerges at the PROTEIN
/// level, where a real protein accumulates many consistent peptides while decoy
/// proteins stay singletons.
#[derive(Debug, Clone)]
pub struct HillHit {
    pub peptide_idx: usize,
    pub decoy: bool,
    /// Fit score (envelope + coelution − mass_z² − rt_z²), for reference.
    pub score: f64,
    /// |predicted − observed| RT (min) of the best cluster.
    pub abs_rt_diff: f64,
    /// Mass deviation of the best cluster, in sigma units (ppm / (ppm_tol/3)).
    pub abs_mass_z: f64,
    /// Coeluting isotope count of the best cluster (shape quality gate).
    pub n_isotopes: f64,
    /// Averagine envelope-fit cosine of the best cluster, in [0, 1] (1 = perfect).
    pub envelope: f64,
    /// Mean coelution correlation of the best cluster's isotopes, in [0, 1].
    pub coelution: f64,
    /// Index of the mono hill this cluster is built on — the "peak" that
    /// peptides compete for (1 peak → 1 peptide).
    pub mono_hill: usize,
}

/// Score every peptide against the hills by ISOTOPE-CLUSTER FIT and return its
/// best cluster's evidence (RT residual, mass z, isotope count). The caller
/// feeds these to `infer_proteins` — the same binomial + razor parsimony the
/// feature pass uses.
///
/// `pred_rt[i]` is the on-data-calibrated predicted RT for `db.peptides[i]`.
pub fn hill_search(
    db: &PeptideDb,
    index: &HillIndex,
    pred_rt: &[Option<f64>],
    rt_sig: f64,
    cfg: &HillConfig,
) -> Vec<HillHit> {
    let rt_sig = rt_sig.max(0.02);
    let ppm_denom = (cfg.ppm / 3.0).max(1e-6);
    let rt_gate = cfg.rt_window_sigma * rt_sig;

    // ONE best cluster per peptide (across all charges). Multi-cluster (one per
    // charge) was tried and HURT (550→510, FDP 0.00%→0.20%): extra peaks keep a
    // peptide alive through competition longer, which keeps DECOY peptides alive
    // just as much — the same reason uno's feature-level UNO_MULTI_FEAT hurts.
    // Single-best acts as a selectivity filter. Peak competition (peptides
    // claiming the same mono hill) is resolved downstream by `score_by_deviation`.
    db.peptides
        .iter()
        .enumerate()
        .filter_map(|(pi, pep)| {
            let prt = pred_rt.get(pi).copied().flatten()?;
            let env_theo = averagine_envelope(pep.mass);
            // Best cluster: (score, rt_err, mass_z, niso, env, coel, mono_hill_idx).
            let mut best: Option<(f64, f64, f64, f64, f64, f64, usize)> = None;

            for &z in &cfg.charges {
                let z = z as f64;
                let mono_mz = (pep.mass + z * PROTON) / z;
                for i in index.window(mono_mz, cfg.ppm) {
                    let mono = &index.hills[i];
                    let rt_err = (mono.rt - prt).abs();
                    if rt_err > rt_gate {
                        continue;
                    }
                    // Match the predicted envelope over CONSECUTIVE isotopes from
                    // +1 (break at the first gap — requiring consecutive isotopes is
                    // the strong real-peptide gate; scanning all positions instead
                    // admitted decoys with a stray +2 and cost depth). "Consider
                    // multiple hills": every hill near each isotope m/z is a
                    // candidate via best_coeluting.
                    let mut obs = [0.0f64; MAX_ISO];
                    obs[0] = mono.intensity.max(0.0);
                    let mut coel_sum = 0.0;
                    let mut niso = 0u32;
                    for k in 1..MAX_ISO {
                        let iso_mz = mono_mz + k as f64 * NEUTRON / z;
                        match index.best_coeluting(iso_mz, cfg.ppm, mono) {
                            Some((j, r)) if r >= cfg.min_coelution => {
                                obs[k] = index.hills[j].intensity.max(0.0);
                                coel_sum += r;
                                niso += 1;
                            }
                            _ => break,
                        }
                    }
                    if niso < cfg.min_isotopes {
                        continue;
                    }
                    let mass_ppm = (mono.mz - mono_mz).abs() / mono_mz / 1e-6;
                    let mass_z = mass_ppm / ppm_denom;
                    let rt_z = rt_err / rt_sig;
                    // Envelope cosine over the found prefix only (mono..=last found),
                    // so a low-abundance peptide isn't penalized for isotopes below noise.
                    let env = envelope_cosine(&obs[..=niso as usize], &env_theo[..=niso as usize]);
                    let coel = coel_sum / niso as f64;
                    let score = cfg.w_envelope * env + cfg.w_coelution * coel
                        - mass_z * mass_z
                        - rt_z * rt_z;
                    if best.map_or(true, |(b, _, _, _, _, _, _)| score > b) {
                        best = Some((score, rt_err, mass_z, niso as f64, env, coel, i));
                    }
                }
            }
            best.map(|(score, abs_rt_diff, abs_mass_z, n_isotopes, envelope, coelution, mono_hill)| {
                HillHit {
                    peptide_idx: pi,
                    decoy: pep.decoy,
                    score,
                    abs_rt_diff,
                    abs_mass_z,
                    n_isotopes,
                    envelope,
                    coelution,
                    mono_hill,
                }
            })
        })
        .collect()
}
