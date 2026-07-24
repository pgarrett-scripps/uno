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

use crate::database::PeptideDb;
use crate::input::FeatureMatrix;
use crate::protein::binom_sf_ge_k;
use crate::score::PeptideHit;
use std::collections::HashMap;

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
    db: &PeptideDb,
    peptides: &[PeptideHit],
    id_confident: &std::collections::HashSet<String>,
    cond_a: &[usize],
    cond_b: &[usize],
    min_peptides: usize,
) -> Vec<ProteinQuant> {
    let n_runs = matrix.run_names.len();
    // Per-run median of log2 intensities (detected cells) + a low floor (5th pct)
    // for imputing missing cells so on/off peptides still yield a fold change.
    let mut per_run_logs: Vec<Vec<f64>> = vec![Vec::new(); n_runs];
    for h in peptides {
        let ints = &matrix.features[h.feature_idx as usize].intensities;
        for (r, &v) in ints.iter().enumerate().take(n_runs) {
            if v > 0.0 {
                per_run_logs[r].push(v.log2());
            }
        }
    }
    let median = |v: &mut [f64]| -> f64 {
        if v.is_empty() {
            return 0.0;
        }
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        v[v.len() / 2]
    };
    let pct5 = |v: &mut [f64]| -> f64 {
        if v.is_empty() {
            return 0.0;
        }
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        v[v.len() / 20]
    };
    let run_med: Vec<f64> = per_run_logs.iter_mut().map(|v| median(v)).collect();
    let grand: f64 = run_med.iter().sum::<f64>() / n_runs.max(1) as f64;
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

    // Group peptides by protein accession (targets restricted to ID-confident).
    let mut by_prot: HashMap<(bool, &str), Vec<&Pep>> = HashMap::new();
    for (pi, decoy, pep) in &pep_of {
        for prot in &db.peptides[*pi].proteins {
            let acc = prot.as_ref();
            if !*decoy && !id_confident.contains(acc) {
                continue;
            }
            by_prot.entry((*decoy, acc)).or_default().push(pep);
        }
    }

    let mut out: Vec<ProteinQuant> = Vec::new();
    for ((decoy, acc), peps) in &by_prot {
        if peps.len() < min_peptides {
            continue;
        }
        let mut sig: Vec<f64> = Vec::new();
        let (mut ab_a, mut ab_b) = (0.0, 0.0);
        for p in peps {
            ab_a += p.ab_a;
            ab_b += p.ab_b;
            if p.fc.abs() >= thr {
                sig.push(p.fc);
            }
        }
        let kpos = sig.iter().filter(|&&x| x > 0.0).count() as u32;
        let kneg = sig.iter().filter(|&&x| x < 0.0).count() as u32;
        let k = kpos.max(kneg);
        let n_sig = sig.len() as u32;
        // co-direction binomial: P(>= k of n_sig share a sign | p=0.5), two-sided.
        let p_value = if n_sig == 0 {
            1.0
        } else {
            (2.0 * binom_sf_ge_k(k, n_sig, 0.5)).min(1.0)
        };
        let log2_fc = median_of(&mut sig.clone());
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
