//! Retention-time prediction: shared featurization + a gbdt model distilled
//! from DeepLC.
//!
//! Both `uno-rt-train` (offline distillation) and the search scorer use
//! [`rt_features`] so the model sees identical inputs at train and inference
//! time. [`rt_features`] emits the 99-d order-aware set described on
//! [`RT_N_FEATURES`] below. It REPLACED an earlier composition-only set in the
//! style of Klammer et al. (2007) (20 AA counts, 20 N-terminal-window counts, 20
//! C-terminal-window counts, length, ln(1+mass), intercept = 63 features), which
//! is no longer built here.

use crate::mass;
use anyhow::{Context, Result};
use gbdt::decision_tree::Data;
use gbdt::gradient_boost::GBDT;
use std::path::Path;

const AA: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";
const N_POS_BINS: usize = 3;

/// Number of features produced by [`rt_features`]. Order-aware "rich" set:
/// 20 AA composition + 6 terminal-residue identities (first 3, last 3) +
/// `N_POS_BINS`*20 position-binned composition + 13 scalar descriptors.
/// Validated to match DeepLC's MS1 protein recovery (vs the old 63-feature
/// composition-only set above, which trailed it), via order-aware features the
/// composition set lacked.
pub const RT_N_FEATURES: usize = 20 + 6 + N_POS_BINS * 20 + 13;

/// Per-residue physicochemical scales, indexed by position in `AA`
/// ("ACDEFGHIKLMNPQRSTVWY"). KD = Kyte-Doolittle hydropathy; RC =
/// Krokhin/SSRCalc-style retention coefficient; EI = Eisenberg consensus.
const KD: [f64; 20] = [
    1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -3.9, 3.8, 1.9, -3.5, -1.6, -3.5, -4.5, -0.8, -0.7,
    4.2, -0.9, -1.3,
];
const RC: [f64; 20] = [
    1.1, 0.5, 0.2, 1.1, 10.9, 0.2, -3.1, 8.0, -3.1, 9.6, 5.0, -1.2, 2.1, -2.0, -4.1, -2.1, 0.6,
    5.5, 11.1, 4.4,
];
const EI: [f64; 20] = [
    0.62, 0.29, -0.90, -0.74, 1.19, 0.48, -0.40, 1.38, -1.50, 1.06, 0.64, -0.78, 0.12, -0.85,
    -2.53, -0.18, -0.05, 1.08, 0.81, 0.26,
];

fn aa_index(c: u8) -> Option<usize> {
    AA.iter().position(|&a| a == c)
}

fn peptide_mass(seq: &[u8]) -> f64 {
    let mut m = mass::H2O;
    for &c in seq {
        m += mass::residue(c).unwrap_or(0.0);
    }
    m
}

/// `numpy.array_split(0..n, k)` boundary: cumulative size of the first `i` of
/// `k` contiguous chunks of `n` items (first `n % k` chunks get one extra).
fn split_bound(n: usize, k: usize, i: usize) -> usize {
    let (base, rem) = (n / k, n % k);
    i * base + i.min(rem)
}

/// Order-aware "rich" RT feature vector for a stripped peptide sequence.
/// Layout: AA composition (20), terminal residue identities (6), position-binned
/// composition (`N_POS_BINS`*20), then 13 scalar descriptors.
pub fn rt_features(seq: &str) -> Vec<f64> {
    let bytes = seq.as_bytes();
    let idx: Vec<usize> = bytes.iter().filter_map(|&c| aa_index(c)).collect();
    let n = idx.len();
    let mut v = vec![0.0f64; RT_N_FEATURES];

    // [0,20) AA composition.
    for &j in &idx {
        v[j] += 1.0;
    }
    // [20,26) terminal residue identities: first 3, last 3 (index value, -1 pad).
    for k in 0..3 {
        v[20 + k] = if k < n { idx[k] as f64 } else { -1.0 };
        v[23 + k] = if k < n { idx[n - 1 - k] as f64 } else { -1.0 };
    }
    // [26, 26+N_POS_BINS*20) position-binned composition (numpy.array_split order).
    let pb = 26;
    for b in 0..N_POS_BINS {
        let (lo, hi) = (
            split_bound(n, N_POS_BINS, b),
            split_bound(n, N_POS_BINS, b + 1),
        );
        for &j in &idx[lo..hi] {
            v[pb + b * 20 + j] += 1.0;
        }
    }
    // [pb+N_POS_BINS*20, +13) scalar descriptors.
    let s = pb + N_POS_BINS * 20;
    let nf = n.max(1) as f64;
    let (mut sum_kd, mut sum_rc, mut sum_ei) = (0.0, 0.0, 0.0);
    let (mut max_kd, mut min_kd) = (f64::NEG_INFINITY, f64::INFINITY);
    let (mut ssr, mut moment) = (0.0, 0.0);
    let (mut p_ct, mut arom, mut kr) = (0.0, 0.0, 0.0);
    for (k, &j) in idx.iter().enumerate() {
        sum_kd += KD[j];
        sum_rc += RC[j];
        sum_ei += EI[j];
        max_kd = max_kd.max(KD[j]);
        min_kd = min_kd.min(KD[j]);
        ssr += RC[j] * (1.0 + if k < 3 { 0.3 } else { 0.0 });
        moment += KD[j]
            * if n > 1 {
                k as f64 / (n - 1) as f64
            } else {
                0.0
            };
        match AA[j] {
            b'P' => p_ct += 1.0,
            b'F' | b'W' | b'Y' => arom += 1.0,
            b'K' | b'R' => kr += 1.0,
            _ => {}
        }
    }
    if n == 0 {
        max_kd = 0.0;
        min_kd = 0.0;
    }
    v[s] = n as f64;
    v[s + 1] = peptide_mass(bytes).ln_1p();
    v[s + 2] = sum_kd;
    v[s + 3] = sum_kd / nf;
    v[s + 4] = max_kd;
    v[s + 5] = min_kd;
    v[s + 6] = sum_rc;
    v[s + 7] = sum_ei;
    v[s + 8] = ssr;
    v[s + 9] = p_ct;
    v[s + 10] = arom;
    v[s + 11] = kr;
    v[s + 12] = moment;
    v
}

fn to_f32(v: &[f64]) -> Vec<f32> {
    v.iter().map(|&x| x as f32).collect()
}

/// Distilled RT backend. Predicts DeepLC abstract RT from sequence; the values
/// are in DeepLC's uncalibrated space (per-run calibration is applied downstream).
enum RtBackend {
    /// Distilled gbdt student (legacy; gbdt-rs serialized JSON).
    Gbdt(GBDT),
    /// Distilled MLP student (`UMLP` binary, see [`crate::rt_mlp`]).
    Mlp(crate::rt_mlp::RtMlp),
    /// Sequence-aware CNN trained on experimental iRT (`URTC` binary, see
    /// [`crate::rt_cnn`]).
    Cnn(crate::rt_cnn::RtCnn),
}

/// The supported pretrained CNN, embedded so installed binaries work anywhere.
pub const DEFAULT_RT_CNN: &[u8] = include_bytes!("../models/rt_cnn_chrono_global_v3_ch32.urtc");

/// A loaded RT model plus a content `fingerprint` (a stable hash of the model
/// file bytes). The fingerprint keys the on-disk prediction cache
/// ([`crate::rt_cache`]) so cached predictions are reused only for the exact
/// model that produced them.
pub struct RtModel {
    backend: RtBackend,
    fingerprint: u64,
}

impl RtModel {
    /// The bundled CNN, with no runtime file or Python dependency.
    pub fn bundled() -> Result<Self> {
        Ok(Self {
            backend: RtBackend::Cnn(crate::rt_cnn::RtCnn::from_bytes(DEFAULT_RT_CNN)?),
            fingerprint: fingerprint_bytes(DEFAULT_RT_CNN),
        })
    }

    /// CNN provenance when available; legacy backends have no CNN metadata.
    pub fn cnn_info(&self) -> Option<&serde_json::Value> {
        match &self.backend {
            RtBackend::Cnn(model) => Some(model.info()),
            _ => None,
        }
    }
    /// Load a distilled RT model, auto-detecting the backend by file magic:
    /// a leading `UMLP` => MLP, `URTC` => CNN, otherwise the gbdt-rs JSON format.
    pub fn load(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let bytes = std::fs::read(path)
            .map_err(|e| anyhow::anyhow!("opening rt model {}: {e}", path.display()))?;
        let fingerprint = fingerprint_bytes(&bytes);
        let magic = bytes.get(0..4).unwrap_or(&[]);
        let backend = if magic == b"UMLP" {
            RtBackend::Mlp(crate::rt_mlp::RtMlp::load(path)?)
        } else if magic == b"URTC" {
            RtBackend::Cnn(crate::rt_cnn::RtCnn::from_bytes(&bytes)?)
        } else {
            let model = GBDT::load_model(&path.to_string_lossy())
                .map_err(|e| anyhow::anyhow!("loading rt model {}: {e}", path.display()))?;
            RtBackend::Gbdt(model)
        };
        if !matches!(&backend, RtBackend::Cnn(_)) {
            log::warn!(
                "legacy RT backend: MLP/GBDT loading is deprecated; use a URTC CNN or --rt-library"
            );
        }
        Ok(RtModel {
            backend,
            fingerprint,
        })
    }

    /// Stable content hash of the model file, used as the prediction-cache key.
    pub fn fingerprint(&self) -> u64 {
        self.fingerprint
    }

    /// Predict abstract RT for one sequence.
    pub fn predict_seq(&self, seq: &str) -> f64 {
        match &self.backend {
            RtBackend::Gbdt(model) => {
                let data = vec![Data::new_test_data(to_f32(&rt_features(seq)), None)];
                model.predict(&data)[0] as f64
            }
            RtBackend::Mlp(m) => m.predict_seq(seq),
            RtBackend::Cnn(m) => m.predict_seq(seq),
        }
    }

    /// Penultimate-layer embeddings for fine-tuning the readout on-data. Only the
    /// MLP/CNN backends expose a usable embedding; the gbdt student returns `None`.
    pub fn embed_many<S: AsRef<str> + Sync>(&self, seqs: &[S]) -> Option<Vec<Vec<f64>>> {
        match &self.backend {
            RtBackend::Mlp(m) => Some(m.embed_many(seqs)),
            RtBackend::Cnn(m) => Some(m.embed_many(seqs)),
            RtBackend::Gbdt(_) => None,
        }
    }

    /// Predict abstract RT AND return the penultimate embedding per sequence in a
    /// single pass, when the on-data RT fine-tune will consume both. Returns
    /// `(preds, Some(embeds))` for backends with a usable embedding (MLP/CNN) and
    /// `(preds, None)` for the gbdt student. Byte-identical to calling
    /// `predict_many` and `embed_many` separately, but the CNN backend runs its
    /// forward pass once instead of twice over the same ~1.1M unique peptides.
    pub fn predict_embed_many<S: AsRef<str> + Sync>(
        &self,
        seqs: &[S],
    ) -> (Vec<f64>, Option<Vec<Vec<f64>>>) {
        match &self.backend {
            RtBackend::Cnn(m) => {
                let (p, e) = m.predict_embed_many(seqs);
                (p, Some(e))
            }
            // MLP: predict and embed remain separate passes (no fused kernel), but
            // still correct + byte-identical to the old two-call path.
            RtBackend::Mlp(m) => (m.predict_many(seqs), Some(m.embed_many(seqs))),
            RtBackend::Gbdt(_) => (self.predict_many(seqs), None),
        }
    }

    /// Predict abstract RT for many sequences in one batch.
    pub fn predict_many<S: AsRef<str> + Sync>(&self, seqs: &[S]) -> Vec<f64> {
        match &self.backend {
            RtBackend::Gbdt(model) => {
                let data: Vec<Data> = seqs
                    .iter()
                    .map(|s| Data::new_test_data(to_f32(&rt_features(s.as_ref())), None))
                    .collect();
                model.predict(&data).into_iter().map(|v| v as f64).collect()
            }
            RtBackend::Mlp(m) => m.predict_many(seqs),
            RtBackend::Cnn(m) => m.predict_many(seqs),
        }
    }
}

/// FNV-1a over the model bytes. Stable across processes/machines (unlike
/// `DefaultHasher`, which is randomized), so a cache file written by one run is
/// keyed identically by the next.
fn fingerprint_bytes(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bytes {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

/// Read a distillation label file: `SEQUENCE\tRT` per line (no header), as
/// produced by sherlock's DeepLC scripts. Returns (sequences, rts).
pub fn read_label_tsv(path: impl AsRef<Path>) -> Result<(Vec<String>, Vec<f64>)> {
    let path = path.as_ref();
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .flexible(true)
        .from_path(path)
        .with_context(|| format!("opening label tsv {}", path.display()))?;

    let mut seqs = Vec::new();
    let mut rts = Vec::new();
    for rec in rdr.records() {
        let rec = rec?;
        let seq = match rec.get(0) {
            Some(s) if !s.is_empty() => s,
            _ => continue,
        };
        let rt: f64 = match rec.get(1).and_then(|s| s.trim().parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        // skip a possible header row
        if seq.chars().all(|c| c.is_ascii_uppercase()) {
            seqs.push(seq.to_string());
            rts.push(rt);
        }
    }
    Ok((seqs, rts))
}

#[cfg(test)]
mod cnn_workflow_tests {
    use super::*;

    #[test]
    fn bundled_cnn_is_portable_and_matches_training_package() {
        assert_eq!(
            DEFAULT_RT_CNN,
            include_bytes!("../../uno-rt-train/python/uno_rt/data/default.urtc")
        );
        let model = RtModel::bundled().unwrap();
        let direct = crate::rt_cnn::RtCnn::from_bytes(DEFAULT_RT_CNN).unwrap();
        assert_eq!(
            model.predict_seq("PEPTIDEK"),
            direct.predict_seq("PEPTIDEK")
        );
        assert_eq!(model.fingerprint(), fingerprint_bytes(DEFAULT_RT_CNN));
    }

    #[test]
    fn cnn_rejects_invalid_dimensions_normalization_and_version() {
        for (offset, value) in [(4, 99u32), (12, 0), (16, 0), (24, 2), (28, u32::MAX)] {
            let mut bytes = DEFAULT_RT_CNN.to_vec();
            bytes[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
            assert!(crate::rt_cnn::RtCnn::from_bytes(&bytes).is_err());
        }
        assert!(crate::rt_cnn::RtCnn::from_bytes(&DEFAULT_RT_CNN[..80]).is_err());
    }
}
