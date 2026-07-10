//! Pure-Rust inference for the `URTC` retention-time CNN — a sequence-aware
//! model (learned AA embedding -> residual dilated Conv1d stack -> masked
//! mean+max pool -> MLP head), trained offline on experimental iRT
//! (ProteomeTools/Prosit) by `crates/uno-rt-train/python/train_rt.py`.
//!
//! Replaces the positional one-hot MLP (`rt_mlp.rs`): a convolution shares
//! weights across positions and reads a learned embedding, so it generalizes
//! across length/composition far better than absolute-position one-hot.
//!
//! The trained model masks PAD positions after every conv, so a position-`t`
//! output never depends on positions >= len. Inference therefore runs over the
//! real peptide length L (no fixed-MAXLEN padding work), with input treated as
//! zero outside [0, L). Activations are laid out position-major (`x[t*ch + c]`)
//! and conv weights as `[out_c][k][in_c]` so the inner accumulation over input
//! channels is contiguous and autovectorizes.

use anyhow::{bail, Context, Result};
use rayon::prelude::*;
use std::path::Path;

const MAGIC: &[u8; 4] = b"URTC";
const MAXLEN: usize = 50;
const NAA: usize = 20;
const PAD: usize = NAA; // embedding pad index (21st row)
const AA: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";

/// Number of explicit global physicochemical features appended to the pooled
/// conv vector before the head (DeepLC-style: the CNN pooling infers these only
/// implicitly, so handing them to the net directly is what closes the RP-RT
/// accuracy gap). Must match `N_GLOBAL` / `global_feats` in
/// `crates/uno-rt-train/python/train_rt.py` EXACTLY (feature order + formulas).
/// URTC v1 models have none (n_global = 0); v2 models store the count in the header.
const N_GLOBAL: usize = 11;

/// Kyte-Doolittle hydropathy, the dominant RP-LC retention driver. Unknown
/// residues -> 0. Keep identical to `_KD` in train_rt.py.
fn kd(b: u8) -> f32 {
    match b {
        b'A' => 1.8,
        b'R' => -4.5,
        b'N' => -3.5,
        b'D' => -3.5,
        b'C' => 2.5,
        b'Q' => -3.5,
        b'E' => -3.5,
        b'G' => -0.4,
        b'H' => -3.2,
        b'I' => 4.5,
        b'L' => 3.8,
        b'K' => -3.9,
        b'M' => 1.9,
        b'F' => 2.8,
        b'P' => -1.6,
        b'S' => -0.8,
        b'T' => -0.7,
        b'W' => -0.9,
        b'Y' => -1.3,
        b'V' => 4.2,
        _ => 0.0,
    }
}

/// Global sequence features, computed identically in Rust (inference) and Python
/// (training). Order (v3, N_GLOBAL=11):
///   0 len/50            1 acidic(D+E)/len    2 basic(K+R)/len   3 His/len
///   4 hydrophobic(FWYLIVM)/len               5 KD hydropathy SUM /50
///   6 GRAVY (mean KD)   7 aromatic(F+W+Y)/len 8 Pro/len
///   9 N-term residue KD /4.5                 10 C-term residue KD /4.5
/// The SUM (5) and terminal descriptors (9,10) are what mean+max pooling destroys:
/// pooling yields the mean but not the total, and weights all positions equally so
/// terminal effects vanish. All are domain-invariant physics.
///
/// NOTE: with the OLD deviation-only scorer this v3 set scored lower than the
/// 5-feature set (248 vs 401 @1% on run10). That reversed once `ml_rescore` became
/// the default: deterministic run10 honest @1% is v3=610 > 5-feat=435 > DeepLC=561,
/// so v3 is the correct shipped default and RT accuracy is NOT the ceiling anymore.
/// See models/rt_cnn_README.md for the full A/B before changing this.
fn global_feats(bytes: &[u8]) -> [f32; N_GLOBAL] {
    let len = bytes.len().max(1) as f32;
    let (mut de, mut kr, mut h, mut phob) = (0.0f32, 0.0f32, 0.0f32, 0.0f32);
    let (mut kd_sum, mut arom, mut pro) = (0.0f32, 0.0f32, 0.0f32);
    for &b in bytes {
        match b {
            b'D' | b'E' => de += 1.0,
            b'K' | b'R' => kr += 1.0,
            b'H' => h += 1.0,
            b'P' => pro += 1.0,
            _ => {}
        }
        if matches!(b, b'F' | b'W' | b'Y' | b'L' | b'I' | b'V' | b'M') {
            phob += 1.0;
        }
        if matches!(b, b'F' | b'W' | b'Y') {
            arom += 1.0;
        }
        kd_sum += kd(b);
    }
    let nterm = bytes.first().map_or(0.0, |&b| kd(b));
    let cterm = bytes.last().map_or(0.0, |&b| kd(b));
    [
        bytes.len() as f32 / 50.0,
        de / len,
        kr / len,
        h / len,
        phob / len,
        kd_sum / 50.0,
        kd_sum / len,
        arom / len,
        pro / len,
        nterm / 4.5,
        cterm / 4.5,
    ]
}

/// 1-D convolution (cross-correlation), same-padding, dilation. Weights are
/// stored `[out_c][k][in_c]` (transposed from PyTorch's `[out_c][in_c][k]`) so
/// the inner loop over input channels is contiguous.
struct Conv {
    out_ch: usize,
    in_ch: usize,
    k: usize,
    dil: usize,
    w: Vec<f32>, // out_ch * k * in_ch
    b: Vec<f32>, // out_ch
}

impl Conv {
    /// `inp` is (L, in_ch) row-major; writes `out` (L, out_ch). Positions
    /// outside [0, L) are treated as zero (matches the trained masking).
    fn forward(&self, inp: &[f32], out: &mut [f32], len: usize) {
        let pad = (self.k - 1) / 2 * self.dil;
        let (kc, ic, oc) = (self.k, self.in_ch, self.out_ch);
        for t in 0..len {
            let obase = t * oc;
            for c in 0..oc {
                out[obase + c] = self.b[c];
            }
            for j in 0..kc {
                let p = t + j * self.dil;
                if p < pad {
                    continue;
                }
                let src = p - pad;
                if src >= len {
                    continue;
                }
                let in_row = &inp[src * ic..src * ic + ic];
                for c in 0..oc {
                    let wrow = &self.w[(c * kc + j) * ic..(c * kc + j) * ic + ic];
                    let mut s = 0.0f32;
                    for i in 0..ic {
                        s += wrow[i] * in_row[i];
                    }
                    out[obase + c] += s;
                }
            }
        }
    }
}

struct Linear {
    out_dim: usize,
    in_dim: usize,
    w: Vec<f32>, // out_dim * in_dim
    b: Vec<f32>,
}

struct Block {
    c1: Conv,
    c2: Conv,
}

/// RT CNN backed by the `URTC` binary.
pub struct RtCnn {
    y_mean: f32,
    y_std: f32,
    emb: usize,
    ch: usize,
    emb_w: Vec<f32>, // (NAA+1) * emb, row-major [token][feat]
    stem: Conv,
    blocks: Vec<Block>,
    head0: Linear,
    head2: Linear,
    n_global: usize, // 0 for URTC v1; N_GLOBAL for v2
}

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
        Ok(v.chunks_exact(4).map(|c| f32::from_le_bytes(c.try_into().unwrap())).collect())
    }
}

/// Load a Conv from PyTorch `[out_c][in_c][k]` layout, transposing to
/// `[out_c][k][in_c]`.
fn load_conv(c: &mut Cur, out_ch: usize, in_ch: usize, k: usize, dil: usize) -> Result<Conv> {
    let raw = c.f32s(out_ch * in_ch * k)?; // [o][i][j]
    let mut w = vec![0.0f32; out_ch * k * in_ch];
    for o in 0..out_ch {
        for i in 0..in_ch {
            for j in 0..k {
                w[(o * k + j) * in_ch + i] = raw[(o * in_ch + i) * k + j];
            }
        }
    }
    let b = c.f32s(out_ch)?;
    Ok(Conv { out_ch, in_ch, k, dil, w, b })
}

impl RtCnn {
    pub fn load(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let bytes = std::fs::read(path)
            .with_context(|| format!("reading rt cnn model {}", path.display()))?;
        let mut c = Cur { b: &bytes, p: 0 };
        if c.b.get(0..4) != Some(&MAGIC[..]) {
            bail!("{}: not a URTC model", path.display());
        }
        c.p = 4;
        let version = c.u32()?;
        if version != 1 && version != 2 {
            bail!("URTC version {version} unsupported");
        }
        let y_mean = c.f32()?;
        let y_std = c.f32()?;
        let emb = c.u32()?;
        let ch = c.u32()?;
        let k = c.u32()?;
        let n_blocks = c.u32()?;
        let dils: Vec<usize> = (0..n_blocks).map(|_| c.u32()).collect::<Result<_>>()?;
        let head = c.u32()?;
        // v2 stores the count of appended global features right after `head`.
        let n_global = if version >= 2 { c.u32()? } else { 0 };
        // `global_feats` computes the full N_GLOBAL-length vector; the 5-feature
        // (v2) set is a strict PREFIX of the 11-feature (v3) set (same order/
        // formulas), so any model with n_global <= N_GLOBAL loads by slicing the
        // computed features to the header's declared count (see `pool`).
        if n_global > N_GLOBAL {
            bail!("URTC declares {n_global} global features, this build computes at most {N_GLOBAL}");
        }
        let head_in = 2 * ch + n_global;

        let emb_w = c.f32s((NAA + 1) * emb)?;
        let stem = load_conv(&mut c, ch, emb, k, 1)?;
        let mut blocks = Vec::with_capacity(n_blocks);
        for &dil in &dils {
            let c1 = load_conv(&mut c, ch, ch, k, dil)?;
            let c2 = load_conv(&mut c, ch, ch, k, dil)?;
            blocks.push(Block { c1, c2 });
        }
        let head0 = Linear { out_dim: head, in_dim: head_in, w: c.f32s(head * head_in)?, b: c.f32s(head)? };
        let head2 = Linear { out_dim: 1, in_dim: head, w: c.f32s(head)?, b: c.f32s(1)? };
        if c.p != bytes.len() {
            bail!("URTC: {} trailing bytes", bytes.len() - c.p);
        }
        Ok(RtCnn { y_mean, y_std, emb, ch, emb_w, stem, blocks, head0, head2, n_global })
    }

    /// Pooled 2*ch feature vector for a sequence (masked mean ++ masked max),
    /// written into `buf.pooled`.
    fn pool(&self, seq: &str, buf: &mut Buffers) {
        let ch = self.ch;
        let bytes = seq.as_bytes();
        let len = bytes.len().min(MAXLEN).max(1);
        // tokens -> embedding (L, emb); PAD rows are zero by construction.
        let x = &mut buf.x;
        x.clear();
        x.resize(len * self.emb, 0.0);
        for (t, &byte) in bytes.iter().take(len).enumerate() {
            let tok = AA.iter().position(|&p| p == byte).unwrap_or(PAD);
            let row = &self.emb_w[tok * self.emb..tok * self.emb + self.emb];
            x[t * self.emb..t * self.emb + self.emb].copy_from_slice(row);
        }
        // stem -> relu (mask is implicit: we only compute t in 0..len)
        let a = &mut buf.a;
        a.resize(len * ch, 0.0);
        self.stem.forward(x, a, len);
        for v in a[..len * ch].iter_mut() {
            *v = v.max(0.0);
        }
        let h = &mut buf.h;
        let h2 = &mut buf.h2;
        h.resize(len * ch, 0.0);
        h2.resize(len * ch, 0.0);
        for blk in &self.blocks {
            blk.c1.forward(a, h, len);
            for v in h[..len * ch].iter_mut() {
                *v = v.max(0.0);
            }
            blk.c2.forward(h, h2, len);
            // relu(a + c2(h)) -> a
            for idx in 0..len * ch {
                a[idx] = (a[idx] + h2[idx]).max(0.0);
            }
        }
        // masked mean + max pool over valid positions -> 2*ch
        let out = &mut buf.pooled;
        out.clear();
        out.resize(2 * ch, 0.0);
        for c in 0..ch {
            let mut sum = 0.0f32;
            let mut mx = f32::NEG_INFINITY;
            for t in 0..len {
                let v = a[t * ch + c];
                sum += v;
                if v > mx {
                    mx = v;
                }
            }
            out[c] = sum / len as f32;
            out[ch + c] = if mx.is_finite() { mx } else { 0.0 };
        }
        // Append explicit global features (v2 models) so the head sees them
        // directly. Computed from the FULL sequence (not the MAXLEN-truncated
        // token view) to match training.
        if self.n_global > 0 {
            let g = global_feats(seq.as_bytes());
            // n_global is <= N_GLOBAL (checked at load); the 5-feature set is the
            // first 5 of the computed 11, so slice to the model's declared count.
            out.extend_from_slice(&g[..self.n_global]);
        }
    }

    fn head(&self, pooled: &[f32], hid: &mut Vec<f32>) -> f64 {
        hid.clear();
        hid.resize(self.head0.out_dim, 0.0);
        for o in 0..self.head0.out_dim {
            let row = &self.head0.w[o * self.head0.in_dim..(o + 1) * self.head0.in_dim];
            let mut s = self.head0.b[o];
            for (wi, &ai) in row.iter().zip(pooled.iter()) {
                s += wi * ai;
            }
            hid[o] = s.max(0.0);
        }
        let mut y = self.head2.b[0];
        for (wi, &ai) in self.head2.w.iter().zip(hid.iter()) {
            y += wi * ai;
        }
        (y * self.y_std + self.y_mean) as f64
    }

    /// Predict abstract RT (iRT units) for one sequence.
    pub fn predict_seq(&self, seq: &str) -> f64 {
        let mut buf = Buffers::default();
        self.pool(seq, &mut buf);
        self.head(&buf.pooled, &mut buf.hid)
    }

    /// Predict abstract RT for many sequences (rayon, order-preserving).
    pub fn predict_many<S: AsRef<str> + Sync>(&self, seqs: &[S]) -> Vec<f64> {
        seqs.par_chunks(2048)
            .flat_map(|chunk| {
                let mut buf = Buffers::default();
                chunk
                    .iter()
                    .map(|s| {
                        self.pool(s.as_ref(), &mut buf);
                        self.head(&buf.pooled, &mut buf.hid)
                    })
                    .collect::<Vec<_>>()
            })
            .collect()
    }

    /// Pooled embedding (+intercept) for on-data readout fine-tuning, mirroring
    /// `RtMlp::embed`: refit a linear map over these frozen pooled features.
    pub fn embed(&self, seq: &str) -> Vec<f64> {
        let mut buf = Buffers::default();
        self.pool(seq, &mut buf);
        let mut out: Vec<f64> = buf.pooled.iter().map(|&v| v as f64).collect();
        out.push(1.0);
        out
    }

    pub fn embed_many<S: AsRef<str> + Sync>(&self, seqs: &[S]) -> Vec<Vec<f64>> {
        seqs.par_iter().map(|s| self.embed(s.as_ref())).collect()
    }
}

/// Reusable per-thread scratch buffers (avoids per-sequence allocation).
#[derive(Default)]
struct Buffers {
    x: Vec<f32>,
    a: Vec<f32>,
    h: Vec<f32>,
    h2: Vec<f32>,
    pooled: Vec<f32>,
    hid: Vec<f32>,
}
