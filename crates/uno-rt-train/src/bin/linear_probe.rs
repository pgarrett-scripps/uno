//! linear_probe — experiment: can a ridge linear model match the distilled GBDT
//! on DeepLC RT, using the SAME `rt_features`?
//!
//! Reports held-out R^2/MAE/RMSE on the SAME deterministic val split as
//! `uno-rt-train` (so it's apples-to-apples), plus train time and inference
//! throughput. Tests two feature encodings:
//!   raw     — the 99-d rt_features as-is (terminal-residue cols are ordinal)
//!   onehot  — one-hot the 6 terminal-residue identity cols [20,26) -> width-21
//!             blocks (bin 20 = pad/-1), the linear-appropriate encoding.
//!
//! Usage: linear_probe [LABELS_TSV]   (default: sherlock ecoli_dda rt_deeplc.tsv)

use anyhow::Result;
use rayon::prelude::*;
use std::time::Instant;
use uno_core::rt::{read_label_tsv, rt_features};

const VAL_FRAC: f64 = 0.1;
const LAMBDAS: &[f64] = &[0.01, 0.1, 1.0, 10.0, 100.0, 1000.0];

/// Same val split as uno-rt-train (splitmix hash of the row index).
fn is_val(idx: usize) -> bool {
    let mut x = idx as u64;
    x ^= x >> 30;
    x = x.wrapping_mul(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    let u = (x >> 11) as f64 / (1u64 << 53) as f64;
    u < VAL_FRAC
}

/// The 13 scalar physicochemical descriptors live at raw indices [86,99).
const DESC: std::ops::Range<usize> = 86..99;

/// Expand a raw rt_features row into f64.
/// `onehot`: one-hot the 6 terminal-residue identity cols [20,26) -> width-21 blocks.
/// `quad`:  squares + pairwise products of the 13 scalar descriptors (2nd order).
/// `cross`: descriptor x AA-composition products (13*20) — the interaction class
///          quad misses (where the GBDT earns most of its edge).
/// `cube`:  all degree-3 monomials of the 13 descriptors (C(15,3)=455).
fn expand_into(row: &[f32], onehot: bool, quad: bool, cross: bool, cube: bool, out: &mut Vec<f64>) {
    out.clear();
    if !onehot {
        out.extend(row.iter().map(|&x| x as f64));
    } else {
        for (j, &x) in row.iter().enumerate() {
            if (20..26).contains(&j) {
                let v = x as i32;
                let bin = if v < 0 { 20usize } else { (v as usize).min(20) };
                let mut oh = [0.0f64; 21];
                oh[bin] = 1.0;
                out.extend_from_slice(&oh);
            } else {
                out.push(x as f64);
            }
        }
    }
    if quad {
        for a in DESC {
            let va = row[a] as f64;
            for b in a..DESC.end {
                out.push(va * row[b] as f64);
            }
        }
    }
    if cross {
        for a in DESC {
            let va = row[a] as f64;
            for j in 0..20 {
                out.push(va * row[j] as f64);
            }
        }
    }
    if cube {
        for a in DESC {
            let va = row[a] as f64;
            for b in a..DESC.end {
                let vab = va * row[b] as f64;
                for c in b..DESC.end {
                    out.push(vab * row[c] as f64);
                }
            }
        }
    }
}

fn expanded_dim(onehot: bool, quad: bool, cross: bool, cube: bool, raw_dim: usize) -> usize {
    let d = DESC.len();
    let mut p = if onehot { raw_dim - 6 + 6 * 21 } else { raw_dim };
    if quad {
        p += d * (d + 1) / 2;
    }
    if cross {
        p += d * 20;
    }
    if cube {
        p += d * (d + 1) * (d + 2) / 6; // multisets of size 3
    }
    p
}

/// Solve (A) x = b with A symmetric positive-definite (full p*p row-major, lower
/// triangle used). In-place Cholesky.
fn cholesky_solve(a: &mut [f64], b: &[f64], p: usize) -> Vec<f64> {
    for j in 0..p {
        let mut d = a[j * p + j];
        for k in 0..j {
            d -= a[j * p + k] * a[j * p + k];
        }
        let ljj = d.max(1e-10).sqrt();
        a[j * p + j] = ljj;
        for i in (j + 1)..p {
            let mut s = a[i * p + j];
            for k in 0..j {
                s -= a[i * p + k] * a[j * p + k];
            }
            a[i * p + j] = s / ljj;
        }
    }
    let mut y = vec![0.0; p];
    for i in 0..p {
        let mut s = b[i];
        for k in 0..i {
            s -= a[i * p + k] * y[k];
        }
        y[i] = s / a[i * p + i];
    }
    let mut x = vec![0.0; p];
    for i in (0..p).rev() {
        let mut s = y[i];
        for k in (i + 1)..p {
            s -= a[k * p + i] * x[k];
        }
        x[i] = s / a[i * p + i];
    }
    x
}

fn metrics(y: &[f64], pred: &[f64]) -> (f64, f64, f64) {
    let n = y.len() as f64;
    let mean = y.iter().sum::<f64>() / n;
    let (mut ss_res, mut ss_tot, mut abs) = (0.0, 0.0, 0.0);
    for (&yi, &pi) in y.iter().zip(pred) {
        let e = yi - pi;
        ss_res += e * e;
        ss_tot += (yi - mean) * (yi - mean);
        abs += e.abs();
    }
    (1.0 - ss_res / ss_tot, abs / n, (ss_res / n).sqrt())
}

fn run_variant(
    feats: &[Vec<f32>],
    ys: &[f64],
    train_idx_full: &[usize],
    val_idx: &[usize],
    onehot: bool,
    quad: bool,
    cross: bool,
    cube: bool,
    train_frac: f64,
    tag: &str,
) {
    // Strided subsample of the training rows (representative regardless of file
    // ordering) to trace a learning curve.
    let stride = (1.0 / train_frac).round().max(1.0) as usize;
    let train_idx: Vec<usize> = train_idx_full.iter().copied().step_by(stride).collect();
    let train_idx = train_idx.as_slice();
    let p = expanded_dim(onehot, quad, cross, cube, feats[0].len());
    let nt = train_idx.len() as f64;

    // Standardization params (mean/std per col over train, expanded).
    let t = Instant::now();
    let (sum, sumsq) = train_idx
        .par_iter()
        .fold(
            || (vec![0.0f64; p], vec![0.0f64; p], Vec::<f64>::new()),
            |(mut s, mut sq, mut buf), &i| {
                expand_into(&feats[i], onehot, quad, cross, cube, &mut buf);
                for j in 0..p {
                    s[j] += buf[j];
                    sq[j] += buf[j] * buf[j];
                }
                (s, sq, buf)
            },
        )
        .map(|(s, sq, _)| (s, sq))
        .reduce(
            || (vec![0.0f64; p], vec![0.0f64; p]),
            |(mut s1, mut q1), (s2, q2)| {
                for j in 0..p {
                    s1[j] += s2[j];
                    q1[j] += q2[j];
                }
                (s1, q1)
            },
        );
    let mean: Vec<f64> = (0..p).map(|j| sum[j] / nt).collect();
    let std: Vec<f64> = (0..p)
        .map(|j| {
            let var = sumsq[j] / nt - mean[j] * mean[j];
            if var > 1e-12 {
                var.sqrt()
            } else {
                1.0 // constant col -> standardized value 0, drops out
            }
        })
        .collect();
    let ybar = train_idx.iter().map(|&i| ys[i]).sum::<f64>() / nt;

    // Normal-equations Gram matrix G = Z^T Z (upper triangle) and b = Z^T (y-ybar).
    let (mut gram, bvec) = train_idx
        .par_chunks(16384)
        .map(|chunk| {
            let mut g = vec![0.0f64; p * p];
            let mut bb = vec![0.0f64; p];
            let mut z = vec![0.0f64; p];
            let mut buf = Vec::<f64>::new();
            for &i in chunk {
                expand_into(&feats[i], onehot, quad, cross, cube, &mut buf);
                for j in 0..p {
                    z[j] = (buf[j] - mean[j]) / std[j];
                }
                let yc = ys[i] - ybar;
                for a in 0..p {
                    let za = z[a];
                    if za == 0.0 {
                        continue;
                    }
                    bb[a] += za * yc;
                    let row = a * p;
                    for c in a..p {
                        g[row + c] += za * z[c];
                    }
                }
            }
            (g, bb)
        })
        .reduce(
            || (vec![0.0f64; p * p], vec![0.0f64; p]),
            |(mut g1, mut b1), (g2, b2)| {
                for k in 0..p * p {
                    g1[k] += g2[k];
                }
                for k in 0..p {
                    b1[k] += b2[k];
                }
                (g1, b1)
            },
        );
    // Mirror upper -> lower for the Cholesky.
    for i in 0..p {
        for j in 0..i {
            gram[i * p + j] = gram[j * p + i];
        }
    }
    let fit_time = t.elapsed().as_secs_f64();

    // Sweep lambda; pick best held-out R^2.
    let mut best = (f64::NEG_INFINITY, 0.0, 0.0, 0.0, Vec::new());
    for &lam in LAMBDAS {
        let mut a = gram.clone();
        for j in 0..p {
            a[j * p + j] += lam;
        }
        let w = cholesky_solve(&mut a, &bvec, p);
        let pred: Vec<f64> = val_idx
            .par_iter()
            .map(|&i| {
                let mut buf = Vec::<f64>::new();
                expand_into(&feats[i], onehot, quad, cross, cube, &mut buf);
                let mut acc = ybar;
                for j in 0..p {
                    acc += (buf[j] - mean[j]) / std[j] * w[j];
                }
                acc
            })
            .collect();
        let val_y: Vec<f64> = val_idx.iter().map(|&i| ys[i]).collect();
        let (r2, mae, rmse) = metrics(&val_y, &pred);
        if r2 > best.0 {
            best = (r2, mae, rmse, lam, w);
        }
    }

    println!(
        "[{tag}] n_train={nt:<7.0} p={p:<4}  best λ={lam:<6}  R²={r2:.4}  MAE={mae:.2}  RMSE={rmse:.2}  | fit {fit:.1}s",
        tag = tag,
        nt = nt,
        p = p,
        lam = best.3,
        r2 = best.0,
        mae = best.1,
        rmse = best.2,
        fit = fit_time,
    );
}

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let dump = args.get(1).map(|s| s == "dump").unwrap_or(false);
    let labels = args
        .get(if dump { 2 } else { 1 })
        .cloned()
        .unwrap_or_else(|| {
            "/home/patrick-garrett/Repos/sherlock/paper/benchmark/ecoli_dda/rt_deeplc.tsv".to_string()
        });

    let t = Instant::now();
    let (seqs, ys) = read_label_tsv(&labels)?;
    eprintln!("loaded {} (seq, rt) labels in {:.1}s", seqs.len(), t.elapsed().as_secs_f64());

    let t = Instant::now();
    let feats: Vec<Vec<f32>> = seqs
        .par_iter()
        .map(|s| rt_features(s).into_iter().map(|x| x as f32).collect())
        .collect();
    eprintln!("featurized in {:.1}s ({}-d)", t.elapsed().as_secs_f64(), feats[0].len());

    // Dump mode: write raw 99-d features + labels + val flag for the Python MLP
    // probe, so it trains on EXACTLY the same features and val split.
    if dump {
        use std::io::Write;
        let p = feats[0].len();
        let mut xw = std::io::BufWriter::new(std::fs::File::create("/tmp/rt_X.f32")?);
        let mut row = Vec::<u8>::with_capacity(p * 4);
        for f in &feats {
            row.clear();
            for &v in f {
                row.extend_from_slice(&v.to_le_bytes());
            }
            xw.write_all(&row)?;
        }
        xw.flush()?;
        let mut yw = std::io::BufWriter::new(std::fs::File::create("/tmp/rt_y.f32")?);
        for &v in &ys {
            yw.write_all(&(v as f32).to_le_bytes())?;
        }
        yw.flush()?;
        let mut vw = std::io::BufWriter::new(std::fs::File::create("/tmp/rt_val.u8")?);
        for i in 0..feats.len() {
            vw.write_all(&[is_val(i) as u8])?;
        }
        vw.flush()?;
        eprintln!("dumped {} rows x {} feats -> /tmp/rt_X.f32, rt_y.f32, rt_val.u8", feats.len(), p);
        return Ok(());
    }

    let train_idx: Vec<usize> = (0..feats.len()).filter(|&i| !is_val(i)).collect();
    let val_idx: Vec<usize> = (0..feats.len()).filter(|&i| is_val(i)).collect();
    eprintln!("train rows: {}, val rows: {}", train_idx.len(), val_idx.len());
    println!("--- ridge (student vs DeepLC teacher), same split as uno-rt-train ---");

    println!("--- learning curve: onehot+quad ridge, val R² vs training-set size (same val split) ---");
    for &frac in &[0.02f64, 0.05, 0.1, 0.25, 0.5, 1.0] {
        let tag = format!("quad frac={frac:<4}");
        run_variant(&feats, &ys, &train_idx, &val_idx, true, true, false, false, frac, &tag);
    }
    println!("--- capacity vs data: a bigger model (+cube) on LESS data beats quad on ALL data ---");
    run_variant(&feats, &ys, &train_idx, &val_idx, true, true, true, true, 0.1, "cube frac=0.1 ");
    run_variant(&feats, &ys, &train_idx, &val_idx, true, true, true, true, 1.0, "cube frac=1.0 ");
    Ok(())
}
