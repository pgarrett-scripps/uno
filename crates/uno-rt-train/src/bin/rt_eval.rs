//! rt_eval — load any uno RT model (UMLP / URTC / gbdt, auto-detected by magic)
//! and report held-out accuracy + inference throughput on a `SEQUENCE\tRT` file.
//!
//! Reports raw metrics (model trained directly on the label scale) and
//! affine-aligned metrics (least-squares pred->label, for models on a different
//! abstract scale, e.g. a DeepLC distillate). dt95 = 95th pct |error|.
//!
//! Usage: rt_eval MODEL LABELS.tsv [--csv preds.tsv]

use anyhow::Result;
use std::time::Instant;
use uno_core::rt::read_label_tsv;
use uno_core::RtModel;

fn metrics(y: &[f64], p: &[f64]) -> (f64, f64, f64, f64) {
    let n = y.len() as f64;
    let mean = y.iter().sum::<f64>() / n;
    let (mut ss_res, mut ss_tot, mut abs) = (0.0, 0.0, 0.0);
    let mut ae: Vec<f64> = Vec::with_capacity(y.len());
    for (&yi, &pi) in y.iter().zip(p) {
        let e = yi - pi;
        ss_res += e * e;
        ss_tot += (yi - mean) * (yi - mean);
        abs += e.abs();
        ae.push(e.abs());
    }
    ae.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let dt95 = ae[((ae.len() as f64 * 0.95) as usize).min(ae.len() - 1)];
    (1.0 - ss_res / ss_tot, abs / n, (ss_res / n).sqrt(), dt95)
}

/// Least-squares affine fit p -> y; returns aligned predictions a*p+b.
fn affine_align(y: &[f64], p: &[f64]) -> Vec<f64> {
    let n = p.len() as f64;
    let (mut sp, mut sy, mut spp, mut spy) = (0.0, 0.0, 0.0, 0.0);
    for (&yi, &pi) in y.iter().zip(p) {
        sp += pi;
        sy += yi;
        spp += pi * pi;
        spy += pi * yi;
    }
    let denom = n * spp - sp * sp;
    let a = if denom.abs() < 1e-12 {
        1.0
    } else {
        (n * spy - sp * sy) / denom
    };
    let b = (sy - a * sp) / n;
    p.iter().map(|&pi| a * pi + b).collect()
}

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let model_path = args.get(1).expect("usage: rt_eval MODEL LABELS.tsv");
    let labels = args.get(2).expect("usage: rt_eval MODEL LABELS.tsv");

    let model = RtModel::load(model_path)?;
    let (seqs, ys) = read_label_tsv(labels)?;
    eprintln!("model {} | {} labels", model_path, seqs.len());

    // throughput: full parallel batch
    let t = Instant::now();
    let preds = model.predict_many(&seqs);
    let par = t.elapsed().as_secs_f64();

    // single-thread on up to 100k
    let n_st = seqs.len().min(100_000);
    let t = Instant::now();
    let mut sink = 0.0;
    for s in &seqs[..n_st] {
        sink += model.predict_seq(s);
    }
    let st = t.elapsed().as_secs_f64();
    std::hint::black_box(sink);

    if let Some(dump) = args.get(3) {
        use std::io::Write;
        let mut f = std::fs::File::create(dump)?;
        for (s, p) in seqs.iter().zip(&preds) {
            writeln!(f, "{}\t{}", s, p)?;
        }
        eprintln!("dumped {} preds -> {}", preds.len(), dump);
    }

    let (r2, mae, rmse, dt95) = metrics(&ys, &preds);
    let aligned = affine_align(&ys, &preds);
    let (ar2, amae, armse, adt95) = metrics(&ys, &aligned);

    println!(
        "raw     : R2={:.4} MAE={:.3} RMSE={:.3} dt95={:.3}",
        r2, mae, rmse, dt95
    );
    println!(
        "aligned : R2={:.4} MAE={:.3} RMSE={:.3} dt95={:.3}",
        ar2, amae, armse, adt95
    );
    println!(
        "speed   : {:.0} pred/s parallel ({:.2}s/{} seqs) | {:.0} pred/s single-thread",
        seqs.len() as f64 / par,
        par,
        seqs.len(),
        n_st as f64 / st
    );
    Ok(())
}
