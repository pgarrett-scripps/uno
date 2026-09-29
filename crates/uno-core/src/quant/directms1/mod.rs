//! Rust reimplementation of pinned DirectMS1Quant contrast analysis, part of uno-core.
//! Adapted from Mark Ivanov's Apache-2.0 implementation; see THIRD_PARTY.md.
pub mod input;
pub mod order;
pub mod output;
pub mod rng;
pub mod stats;
pub mod uno_adapter;
use anyhow::{ensure, Context, Result};
use serde::{Deserialize, Serialize};
use stats::{mean, median};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Config {
    pub min_samples: usize,
    pub min_length: usize,
    pub normalization: u8,
    pub fold_change: f64,
    pub absolute_fold_change: bool,
    pub qvalue: f64,
    pub decoy_prefix: String,
    pub seed: u32,
}
impl Default for Config {
    fn default() -> Self {
        Self {
            min_samples: 0,
            min_length: 7,
            normalization: 2,
            fold_change: 2.,
            absolute_fold_change: false,
            qvalue: 0.05,
            decoy_prefix: "DECOY_".into(),
            seed: 42,
        }
    }
}
#[derive(Clone, Debug)]
pub struct InputPeptide {
    pub key: String,
    pub sequence: String,
    pub proteins: Vec<String>,
    pub values: Vec<Option<f64>>,
}
#[derive(Clone, Debug)]
pub struct Dataset {
    pub samples: Vec<String>,
    pub n_s1: usize,
    pub rows: Vec<InputPeptide>,
}
#[derive(Clone, Debug, Serialize)]
pub struct Peptide {
    pub peptide: String,
    pub sequence: String,
    pub protein: String,
    pub decoy: bool,
    pub intensities: Vec<f64>,
    pub missing: usize,
    pub observed_s1: usize,
    pub observed_s2: usize,
    pub intensity: f64,
    pub bin: usize,
    pub raw_fc: f64,
    pub fc: f64,
    pub corrected_fc: f64,
    pub pvalue: f64,
    pub running_fc: f64,
    pub better: usize,
    pub significant: bool,
    pub up: bool,
    pub down: bool,
}
#[derive(Clone, Debug, Serialize)]
pub struct Protein {
    pub name: String,
    pub group: usize,
    pub n_peptides: usize,
    pub n_changed: usize,
    pub score: f64,
    pub pvalue: f64,
    pub fc: f64,
    pub raw_fc: f64,
    pub all_fc: f64,
    pub all_raw_fc: f64,
    pub fc_pass: bool,
    pub fdr_pass: Option<bool>,
    pub bh_pass: bool,
}
#[derive(Clone, Debug, Serialize)]
pub struct Analysis {
    pub config: Config,
    pub samples: Vec<String>,
    pub n_s1: usize,
    pub normalization_factors: Vec<f64>,
    pub calibration_mean: f64,
    pub calibration_std: f64,
    pub fold_change_threshold: f64,
    pub decoy_probability: f64,
    pub n_decoy_sequences: usize,
    pub input_rows: usize,
    pub eligible_precursors: usize,
    pub peptides: Vec<Peptide>,
    pub proteins: Vec<Protein>,
}
fn observed(v: &[Option<f64>]) -> Vec<f64> {
    v.iter().flatten().copied().collect()
}
fn opt_mean(v: &[Option<f64>]) -> f64 {
    mean(&observed(v))
}
fn bounds(v: &[f64]) -> (f64, f64) {
    v.iter()
        .filter(|v| v.is_finite())
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &x| {
            (lo.min(x), hi.max(x))
        })
}
fn fill_absent(fc: &mut [f64], counts: &[(usize, usize)]) -> Result<()> {
    let (lo, hi) = bounds(fc);
    ensure!(
        lo.is_finite() && hi.is_finite(),
        "no peptides observed in both conditions"
    );
    for (f, &(a, b)) in fc.iter_mut().zip(counts) {
        if a == 0 {
            *f = hi;
        } else if b == 0 {
            *f = lo;
        }
    }
    Ok(())
}

pub fn analyze(data: &Dataset, cfg: &Config) -> Result<Analysis> {
    ensure!(
        data.n_s1 <= data.samples.len(),
        "condition size exceeds sample count"
    );
    ensure!(
        data.n_s1 >= 2 && data.samples.len() - data.n_s1 >= 2,
        "at least two replicates per condition are required"
    );
    ensure!(
        data.samples.iter().collect::<BTreeSet<_>>().len() == data.samples.len(),
        "duplicate sample names"
    );
    ensure!(
        cfg.normalization <= 2 && cfg.qvalue.is_finite() && cfg.qvalue > 0. && cfg.qvalue <= 1.,
        "invalid normalization or qvalue"
    );
    ensure!(
        cfg.fold_change.is_finite() && cfg.fold_change >= 0. && !cfg.decoy_prefix.is_empty(),
        "invalid fold change or decoy prefix"
    );
    let n = data.samples.len();
    let required = if cfg.min_samples == 0 {
        n / 2
    } else {
        cfg.min_samples
    };
    ensure!(
        required > 0 && required <= n,
        "min_samples exceeds sample count"
    );
    let mut seen = BTreeSet::new();
    for row in &data.rows {
        ensure!(seen.insert(&row.key), "duplicate precursor key {}", row.key);
        ensure!(
            row.values.len() == n && !row.proteins.is_empty(),
            "invalid row shape: {}",
            row.key
        );
        ensure!(
            row.values
                .iter()
                .flatten()
                .all(|v| v.is_finite() && *v > 0.),
            "intensities must be positive and finite; use missing for absent signal: {}",
            row.key
        );
    }
    let rows: Vec<_> = data
        .rows
        .iter()
        .filter(|r| r.values.iter().flatten().count() >= required)
        .collect();
    ensure!(
        rows.len() >= 6,
        "at least six eligible precursor rows needed for quintiles/calibration"
    );
    let counts: Vec<_> = rows
        .iter()
        .map(|r| {
            (
                r.values[..data.n_s1].iter().flatten().count(),
                r.values[data.n_s1..].iter().flatten().count(),
            )
        })
        .collect();
    let mut raw_fc: Vec<_> = rows
        .iter()
        .map(|r| (opt_mean(&r.values[data.n_s1..]) / opt_mean(&r.values[..data.n_s1])).log2())
        .collect();
    fill_absent(&mut raw_fc, &counts)?;
    let factors: Vec<f64> = (0..n)
        .map(|col| {
            let mut vals: Vec<_> = rows.iter().filter_map(|r| r.values[col]).collect();
            match cfg.normalization {
                0 => 1.,
                1 => median(&vals),
                _ => {
                    vals.sort_by(|a, b| b.total_cmp(a));
                    vals.iter().take(1000).sum()
                }
            }
        })
        .collect();
    ensure!(
        factors.iter().all(|x| x.is_finite() && *x > 0.),
        "a sample has no valid normalization factor"
    );
    let normalized: Vec<Vec<Option<f64>>> = rows
        .iter()
        .map(|r| {
            r.values
                .iter()
                .zip(&factors)
                .map(|(v, f)| v.map(|v| v / f))
                .collect()
        })
        .collect();
    let minima: Vec<_> = (0..n)
        .map(|i| {
            normalized
                .iter()
                .filter_map(|r| r[i])
                .fold(f64::INFINITY, f64::min)
        })
        .collect();
    ensure!(
        minima.iter().all(|x| x.is_finite() && *x > 0.),
        "a sample has no observed signal after filtering"
    );
    let means: Vec<_> = normalized
        .iter()
        .map(|r| (opt_mean(&r[..data.n_s1]), opt_mean(&r[data.n_s1..])))
        .collect();
    let intensity: Vec<_> = means
        .iter()
        .map(|&(a, b)| {
            if a.is_nan() {
                b
            } else if b.is_nan() {
                a
            } else {
                a.max(b)
            }
        })
        .collect();
    let edges: Vec<_> = (0..=5)
        .map(|i| stats::quantile(&intensity, i as f64 / 5.))
        .collect();
    ensure!(
        edges.windows(2).all(|w| w[0] < w[1]),
        "nonunique intensity quintile edges; upstream qcut is undefined for this panel"
    );
    let mut fc: Vec<_> = means.iter().map(|&(a, b)| (b / a).log2()).collect();
    let (lo, hi) = bounds(&fc);
    let calib: Vec<_> = fc
        .iter()
        .filter(|x| x.is_finite() && **x >= lo / 2. && **x <= hi / 2.)
        .copied()
        .collect();
    let (center, sigma) = stats::calibrate(&calib).context("fold-change calibration")?;
    fill_absent(&mut fc, &counts)?;
    let threshold = cfg.fold_change * if cfg.absolute_fold_change { 1. } else { sigma };
    let mut peptides = Vec::new();
    for (i, row) in rows.iter().enumerate() {
        let norm = &normalized[i];
        let completed: Vec<_> = norm
            .iter()
            .enumerate()
            .map(|(c, v)| v.unwrap_or(minima[c]))
            .collect();
        let log_a: Vec<_> = observed(&norm[..data.n_s1])
            .iter()
            .map(|v| v.log10())
            .collect();
        let log_b: Vec<_> = observed(&norm[data.n_s1..])
            .iter()
            .map(|v| v.log10())
            .collect();
        let mut p = stats::ttest(&log_a, &log_b);
        if !p.is_finite() {
            p = stats::ttest(
                &completed[..data.n_s1]
                    .iter()
                    .map(|x| x.log10())
                    .collect::<Vec<_>>(),
                &completed[data.n_s1..]
                    .iter()
                    .map(|x| x.log10())
                    .collect::<Vec<_>>(),
            );
        }
        if !p.is_finite() {
            p = 1.;
        }
        let bin = edges[1..].partition_point(|x| *x < intensity[i]).min(4);
        let decoy = row
            .proteins
            .iter()
            .all(|p| p.starts_with(&cfg.decoy_prefix));
        for protein in &row.proteins {
            peptides.push(Peptide {
                peptide: row.key.clone(),
                sequence: row.sequence.clone(),
                protein: protein.clone(),
                decoy,
                intensities: completed.clone(),
                missing: n - counts[i].0 - counts[i].1,
                observed_s1: counts[i].0,
                observed_s2: counts[i].1,
                intensity: intensity[i],
                bin,
                raw_fc: raw_fc[i],
                fc: fc[i],
                corrected_fc: fc[i] - center,
                pvalue: p,
                running_fc: 0.,
                better: 0,
                significant: false,
                up: false,
                down: false,
            });
        }
    }
    peptides.sort_by(|a, b| {
        a.missing
            .cmp(&b.missing)
            .then(b.intensity.total_cmp(&a.intensity))
    });
    let mut seen = BTreeSet::new();
    peptides.retain(|p| seen.insert((p.sequence.clone(), p.protein.clone())));
    let sorted = order::argsort(
        &peptides
            .iter()
            .map(|p| p.corrected_fc.abs())
            .collect::<Vec<_>>(),
        false,
    );
    peptides = sorted.into_iter().map(|i| peptides[i].clone()).collect();
    let mut groups: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (i, p) in peptides.iter().enumerate() {
        groups.entry(p.protein.clone()).or_default().push(i);
    }
    for indices in groups.values() {
        let mut sum = 0.;
        for (k, &i) in indices.iter().enumerate() {
            sum += peptides[i].corrected_fc;
            peptides[i].running_fc = sum / (k + 1) as f64;
        }
    }
    let mut pools = vec![Vec::new(); 5];
    for p in &peptides {
        if p.decoy {
            pools[p.bin].push((p.corrected_fc, p.pvalue));
        }
    }
    for (bin, pool) in pools.iter().enumerate() {
        ensure!(
            !pool.is_empty(),
            "no decoy observations in intensity quintile {bin}; cannot resample null"
        );
    }
    let mut rng = rng::PythonRandom::new(cfg.seed);
    for indices in groups.values() {
        for _ in 0..100 {
            let mut sum = 0.;
            for (k, &i) in indices.iter().enumerate() {
                let p = &mut peptides[i];
                let pair = pools[p.bin][rng.choice(pools[p.bin].len())];
                sum += pair.0;
                if p.running_fc.abs() >= (sum / (k + 1) as f64).abs() && p.pvalue <= pair.1 {
                    p.better += 1;
                }
            }
        }
        for &i in indices {
            let p = &mut peptides[i];
            p.significant = p.better >= 58;
            p.up = p.significant && p.corrected_fc > 0.;
            p.down = p.significant && p.corrected_fc < 0.;
        }
    }
    // The peptide export is before majority-direction replacement, as upstream.
    let exported = peptides.clone();
    for indices in groups.values() {
        let up = indices.iter().filter(|&&i| peptides[i].up).count();
        let down = indices.iter().filter(|&&i| peptides[i].down).count();
        if down > up {
            for &i in indices {
                peptides[i].up = peptides[i].down;
            }
        }
    }
    let mut by_support: Vec<_> = (0..peptides.len()).collect();
    by_support.sort_by(|&a, &b| {
        peptides[a]
            .missing
            .cmp(&peptides[b].missing)
            .then(peptides[b].intensity.total_cmp(&peptides[a].intensity))
    });
    let mut decoy_seen = BTreeSet::new();
    let mut decoy_changed = 0usize;
    for &i in &by_support {
        let p = &peptides[i];
        if p.decoy && decoy_seen.insert(&p.sequence) {
            decoy_changed += usize::from(p.up);
        }
    }
    ensure!(
        !decoy_seen.is_empty(),
        "no decoy sequences for protein calibration"
    );
    let p_up = decoy_changed as f64 / decoy_seen.len() as f64;
    let mut signatures: BTreeMap<String, usize> = BTreeMap::new();
    let mut proteins = Vec::new();
    for (name, indices) in &groups {
        let mut seqs: Vec<_> = indices
            .iter()
            .map(|&i| peptides[i].sequence.as_str())
            .collect();
        seqs.sort();
        let signature = seqs.concat();
        let next = signatures.len() + 1;
        let group = *signatures.entry(signature).or_insert(next);
        if name.starts_with(&cfg.decoy_prefix) {
            continue;
        }
        let selected: Vec<_> = indices
            .iter()
            .copied()
            .filter(|&i| peptides[i].up)
            .collect();
        let all_fc = median(
            &indices
                .iter()
                .map(|&i| peptides[i].corrected_fc)
                .collect::<Vec<_>>(),
        );
        let all_raw_fc = median(
            &indices
                .iter()
                .map(|&i| peptides[i].raw_fc)
                .collect::<Vec<_>>(),
        );
        let fc = if selected.is_empty() {
            all_fc
        } else {
            median(
                &selected
                    .iter()
                    .map(|&i| peptides[i].corrected_fc)
                    .collect::<Vec<_>>(),
            )
        };
        let raw_fc = if selected.is_empty() {
            all_raw_fc
        } else {
            median(
                &selected
                    .iter()
                    .map(|&i| peptides[i].raw_fc)
                    .collect::<Vec<_>>(),
            )
        };
        let score = stats::binomial_score(selected.len() as u64, indices.len() as u64, p_up);
        proteins.push(Protein {
            name: name.clone(),
            group,
            n_peptides: indices.len(),
            n_changed: selected.len(),
            score,
            pvalue: 10f64.powf(-score),
            fc,
            raw_fc,
            all_fc,
            all_raw_fc,
            fc_pass: fc.abs() >= threshold,
            fdr_pass: None,
            bh_pass: false,
        });
    }
    let sorted = order::argsort(&proteins.iter().map(|p| p.score).collect::<Vec<_>>(), true);
    proteins = sorted.into_iter().map(|i| proteins[i].clone()).collect();
    let n_groups = proteins
        .iter()
        .filter(|p| p.n_peptides >= 3 && p.fc_pass)
        .map(|p| p.group)
        .collect::<BTreeSet<_>>()
        .len();
    let mut ranked = BTreeSet::new();
    for p in &mut proteins {
        if p.n_peptides >= 3 && p.fc_pass {
            ranked.insert(p.group);
            let bh_threshold = -(ranked.len() as f64 * cfg.qvalue / n_groups as f64).log10();
            p.bh_pass = p.score >= bh_threshold;
            p.fdr_pass = Some(p.score >= -cfg.qvalue.log10());
        }
    }
    Ok(Analysis {
        config: cfg.clone(),
        samples: data.samples.clone(),
        n_s1: data.n_s1,
        normalization_factors: factors,
        calibration_mean: center,
        calibration_std: sigma,
        fold_change_threshold: threshold,
        decoy_probability: p_up,
        n_decoy_sequences: decoy_seen.len(),
        input_rows: data.rows.len(),
        eligible_precursors: rows.len(),
        peptides: exported,
        proteins,
    })
}
