//! mlp_bench — parity-check the pure-Rust RtMlp against the Keras export, and
//! benchmark its inference throughput. Reports held-out R²/MAE/RMSE on the same
//! val split as uno-rt-train so it lines up with the GBDT/ridge numbers.
//!
//! Usage: mlp_bench [MODEL.bin] [LABELS.tsv]

use anyhow::Result;
use std::time::Instant;
use uno_core::rt::read_label_tsv;
use uno_core::RtMlp;

fn is_val(idx: usize) -> bool {
    let mut x = idx as u64;
    x ^= x >> 30;
    x = x.wrapping_mul(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    let u = (x >> 11) as f64 / (1u64 << 53) as f64;
    u < 0.1
}

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let model = args.get(1).map(|s| s.as_str()).unwrap_or("models/rt_mlp.bin");
    let labels = args
        .get(2)
        .map(|s| s.as_str())
        .unwrap_or("/home/patrick-garrett/Repos/sherlock/paper/benchmark/ecoli_dda/rt_deeplc.tsv");

    let m = RtMlp::load(model)?;
    let (seqs, ys) = read_label_tsv(labels)?;
    eprintln!("loaded model {} and {} labels", model, seqs.len());

    // Predict all (rayon) and time it.
    let t = Instant::now();
    let preds = m.predict_many(&seqs);
    let par_s = t.elapsed().as_secs_f64();

    // Single-thread timing on a 100k slice.
    let n_st = seqs.len().min(100_000);
    let t = Instant::now();
    let mut sink = 0.0;
    for s in &seqs[..n_st] {
        sink += m.predict_seq(s);
    }
    let st_s = t.elapsed().as_secs_f64();
    std::hint::black_box(sink);

    // Held-out metrics (same val split).
    let (mut ss_res, mut ss_tot, mut abs, mut n) = (0.0, 0.0, 0.0, 0.0);
    let mut mean = 0.0;
    for (i, &y) in ys.iter().enumerate() {
        if is_val(i) {
            mean += y;
            n += 1.0;
        }
    }
    mean /= n;
    for (i, (&y, &p)) in ys.iter().zip(&preds).enumerate() {
        if is_val(i) {
            let e = y - p;
            ss_res += e * e;
            ss_tot += (y - mean) * (y - mean);
            abs += e.abs();
        }
    }
    let r2 = 1.0 - ss_res / ss_tot;
    let mae = abs / n;
    let rmse = (ss_res / n).sqrt();

    let np = seqs.len() as f64;
    println!("Rust RtMlp parity (val, same split):  R²={r2:.4}  MAE={mae:.2}  RMSE={rmse:.2}  (n={n:.0})");
    println!(
        "inference: 1-thread {st:.0} pep/s  |  rayon {par:.0} pep/s ({alln} peptides in {pars:.2}s)",
        st = n_st as f64 / st_s,
        par = np / par_s,
        alln = seqs.len(),
        pars = par_s,
    );
    Ok(())
}
