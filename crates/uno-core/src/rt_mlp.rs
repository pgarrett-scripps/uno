//! Pure-Rust MLP retention-time model — a faster, smaller, more accurate student
//! than the distilled gbdt (see `rt.rs`). Trained offline (Keras) on the shared
//! [`crate::rt::rt_features`] and serialized to a compact binary (`models/rt_mlp.bin`).
//!
//! The student maps a peptide sequence -> DeepLC abstract RT via:
//!   rt_features (99-d) -> one-hot the 6 terminal-residue cols (-> 219-d)
//!   -> standardize -> ReLU MLP -> de-standardize.
//! Per-run calibration to observed RT is applied downstream (phase 5), exactly as
//! for the gbdt model.

use crate::rt::rt_features;
use anyhow::{bail, Context, Result};
use rayon::prelude::*;
use std::path::Path;

const MAGIC: &[u8; 4] = b"UMLP";

struct Layer {
    out_dim: usize,
    in_dim: usize,
    w: Vec<f32>, // row-major, out_dim x in_dim
    b: Vec<f32>, // out_dim
    relu: bool,
}

/// A distilled RT model backed by a small multilayer perceptron.
pub struct RtMlp {
    y_mean: f32,
    y_std: f32,
    x_mean: Vec<f32>, // in_dim
    x_std: Vec<f32>,  // in_dim
    in_dim: usize,
    /// Input featurization: 0 = legacy 99-d Klammer/order features (+ one-hot of
    /// the 6 terminal cols); 1 = full positional one-hot (MAXLEN*20). The
    /// positional encoding lets the MLP see sequence position, closing the
    /// student-teacher gap to DeepLC from ~12 min to ~4 min (the hand-crafted
    /// features collapse position into 3 bins and top out near 12 min RMSE).
    feature_kind: usize,
    layers: Vec<Layer>,
}

/// Max peptide length the positional one-hot covers (matches training).
const POS_MAXLEN: usize = 50;
const POS_NAA: usize = 20;
const POS_AA: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";

/// Little-endian byte cursor over the model file.
struct Cur<'a> {
    b: &'a [u8],
    p: usize,
}
impl<'a> Cur<'a> {
    fn u32(&mut self) -> Result<usize> {
        let e = self.p + 4;
        let v = self.b.get(self.p..e).context("truncated (u32)")?;
        self.p = e;
        Ok(u32::from_le_bytes(v.try_into().unwrap()) as usize)
    }
    fn f32(&mut self) -> Result<f32> {
        let e = self.p + 4;
        let v = self.b.get(self.p..e).context("truncated (f32)")?;
        self.p = e;
        Ok(f32::from_le_bytes(v.try_into().unwrap()))
    }
    fn f32s(&mut self, n: usize) -> Result<Vec<f32>> {
        let e = self.p + 4 * n;
        let v = self.b.get(self.p..e).context("truncated (f32 array)")?;
        self.p = e;
        Ok(v.chunks_exact(4)
            .map(|c| f32::from_le_bytes(c.try_into().unwrap()))
            .collect())
    }
}

impl RtMlp {
    pub fn load(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let bytes = std::fs::read(path)
            .with_context(|| format!("reading rt mlp model {}", path.display()))?;
        let mut c = Cur { b: &bytes, p: 0 };
        if c.b.get(0..4) != Some(&MAGIC[..]) {
            bail!("{}: not an UMLP model", path.display());
        }
        c.p = 4;
        let version = c.u32()?;
        // v2+ carries a feature_kind word right after the version; v1 had none.
        let feature_kind = if version >= 2 { c.u32()? } else { 0 };
        let y_mean = c.f32()?;
        let y_std = c.f32()?;
        let in_dim = c.u32()?;
        let x_mean = c.f32s(in_dim)?;
        let x_std = c.f32s(in_dim)?;
        let n_layers = c.u32()?;
        let mut layers = Vec::with_capacity(n_layers);
        let mut prev_in = in_dim;
        for li in 0..n_layers {
            let out_dim = c.u32()?;
            let in_d = c.u32()?;
            let relu = c.u32()? != 0;
            if in_d != prev_in {
                bail!("layer {li}: in_dim {in_d} != expected {prev_in}");
            }
            let w = c.f32s(out_dim * in_d)?;
            let b = c.f32s(out_dim)?;
            layers.push(Layer {
                out_dim,
                in_dim: in_d,
                w,
                b,
                relu,
            });
            prev_in = out_dim;
        }
        if layers.last().map(|l| l.out_dim) != Some(1) {
            bail!("final layer must have out_dim 1");
        }
        Ok(RtMlp {
            y_mean,
            y_std,
            x_mean,
            x_std,
            in_dim,
            feature_kind,
            layers,
        })
    }

    /// Featurize+standardize into `x` (len in_dim), dispatching on `feature_kind`.
    fn featurize_std(&self, seq: &str, x: &mut Vec<f32>) {
        x.clear();
        if self.feature_kind == 1 {
            // Positional one-hot: x[pos*20 + aa] = 1, padded/truncated to MAXLEN.
            x.resize(POS_MAXLEN * POS_NAA, 0.0);
            for (pos, &c) in seq.as_bytes().iter().take(POS_MAXLEN).enumerate() {
                if let Some(a) = POS_AA.iter().position(|&p| p == c) {
                    x[pos * POS_NAA + a] = 1.0;
                }
            }
            for i in 0..x.len() {
                x[i] = (x[i] - self.x_mean[i]) / self.x_std[i];
            }
            return;
        }
        let raw = rt_features(seq); // legacy 99-d
        for (j, &v) in raw.iter().enumerate() {
            if (20..26).contains(&j) {
                let vi = v as i32;
                let bin = if vi < 0 {
                    20usize
                } else {
                    (vi as usize).min(20)
                };
                for k in 0..21 {
                    x.push(if k == bin { 1.0 } else { 0.0 });
                }
            } else {
                x.push(v as f32);
            }
        }
        debug_assert_eq!(x.len(), self.in_dim);
        for i in 0..x.len() {
            x[i] = (x[i] - self.x_mean[i]) / self.x_std[i];
        }
    }

    /// Positional one-hot forward: the input is ~95% zeros (one active entry per
    /// residue), so the first layer is computed as a sparse sum of weight columns
    /// for the active indices `pos*20+aa` — collapsing a 1000x512 matmul to ~50x512
    /// adds — then the remaining layers run dense. Requires x_mean=0, x_std=1
    /// (the positional model is exported unstandardized); falls back to dense
    /// otherwise.
    fn forward_positional(&self, seq: &str, a: &mut Vec<f32>, b: &mut Vec<f32>) -> f32 {
        let l0 = &self.layers[0];
        a.clear();
        a.resize(l0.out_dim, 0.0);
        a.copy_from_slice(&l0.b[..l0.out_dim]);
        for (pos, &c) in seq.as_bytes().iter().take(POS_MAXLEN).enumerate() {
            if let Some(aa) = POS_AA.iter().position(|&p| p == c) {
                let idx = pos * POS_NAA + aa;
                for o in 0..l0.out_dim {
                    a[o] += l0.w[o * l0.in_dim + idx];
                }
            }
        }
        if l0.relu {
            for v in a.iter_mut() {
                *v = v.max(0.0);
            }
        }
        for layer in &self.layers[1..] {
            b.clear();
            b.resize(layer.out_dim, 0.0);
            for o in 0..layer.out_dim {
                let row = &layer.w[o * layer.in_dim..(o + 1) * layer.in_dim];
                let mut s = layer.b[o];
                for (wi, &ai) in row.iter().zip(a.iter()) {
                    s += wi * ai;
                }
                b[o] = if layer.relu { s.max(0.0) } else { s };
            }
            std::mem::swap(a, b);
        }
        a[0] * self.y_std + self.y_mean
    }

    /// True when the positional sparse path is valid (no input standardization).
    fn positional_sparse(&self) -> bool {
        self.feature_kind == 1
            && self.x_mean.iter().all(|&m| m == 0.0)
            && self.x_std.iter().all(|&s| s == 1.0)
    }

    /// Forward pass using two reusable activation buffers; returns RT in label units.
    fn forward(&self, x: &[f32], a: &mut Vec<f32>, b: &mut Vec<f32>) -> f32 {
        a.clear();
        a.extend_from_slice(x);
        for layer in &self.layers {
            b.clear();
            b.resize(layer.out_dim, 0.0);
            for o in 0..layer.out_dim {
                let row = &layer.w[o * layer.in_dim..(o + 1) * layer.in_dim];
                let mut s = layer.b[o];
                for (wi, &ai) in row.iter().zip(a.iter()) {
                    s += wi * ai;
                }
                b[o] = if layer.relu { s.max(0.0) } else { s };
            }
            std::mem::swap(a, b);
        }
        a[0] * self.y_std + self.y_mean
    }

    /// Penultimate-layer embedding for a sequence: the activations fed into the
    /// final linear readout. Fine-tuning refits a fresh linear map on THIS run's
    /// observed RT over these (frozen) features, keeping the MLP's learned
    /// nonlinear representation while adapting the readout to the gradient. A
    /// trailing 1.0 is appended as the ridge intercept term.
    pub fn embed(&self, seq: &str) -> Vec<f64> {
        let (mut x, mut a, mut b) = (Vec::with_capacity(self.in_dim), Vec::new(), Vec::new());
        self.featurize_std(seq, &mut x);
        a.extend_from_slice(&x);
        let last = self.layers.len() - 1; // stop before the final (readout) layer
        for layer in &self.layers[..last] {
            b.clear();
            b.resize(layer.out_dim, 0.0);
            for o in 0..layer.out_dim {
                let row = &layer.w[o * layer.in_dim..(o + 1) * layer.in_dim];
                let mut s = layer.b[o];
                for (wi, &ai) in row.iter().zip(a.iter()) {
                    s += wi * ai;
                }
                b[o] = if layer.relu { s.max(0.0) } else { s };
            }
            std::mem::swap(&mut a, &mut b);
        }
        let mut out: Vec<f64> = a.iter().map(|&v| v as f64).collect();
        out.push(1.0); // intercept
        out
    }

    /// Penultimate embeddings for many sequences (rayon, order-preserving).
    pub fn embed_many<S: AsRef<str> + Sync>(&self, seqs: &[S]) -> Vec<Vec<f64>> {
        seqs.par_iter().map(|s| self.embed(s.as_ref())).collect()
    }

    /// Predict abstract RT for one sequence.
    pub fn predict_seq(&self, seq: &str) -> f64 {
        let (mut a, mut b) = (Vec::new(), Vec::new());
        if self.positional_sparse() {
            return self.forward_positional(seq, &mut a, &mut b) as f64;
        }
        let mut x = Vec::with_capacity(self.in_dim);
        self.featurize_std(seq, &mut x);
        self.forward(&x, &mut a, &mut b) as f64
    }

    /// Predict abstract RT for many sequences (rayon, order-preserving).
    pub fn predict_many<S: AsRef<str> + Sync>(&self, seqs: &[S]) -> Vec<f64> {
        let sparse = self.positional_sparse();
        seqs.par_chunks(4096)
            .flat_map(|chunk| {
                let (mut x, mut a, mut b) =
                    (Vec::with_capacity(self.in_dim), Vec::new(), Vec::new());
                let mut out = Vec::with_capacity(chunk.len());
                for s in chunk {
                    if sparse {
                        out.push(self.forward_positional(s.as_ref(), &mut a, &mut b) as f64);
                    } else {
                        self.featurize_std(s.as_ref(), &mut x);
                        out.push(self.forward(&x, &mut a, &mut b) as f64);
                    }
                }
                out
            })
            .collect()
    }
}
