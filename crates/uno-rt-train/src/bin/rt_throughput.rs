//! rt_throughput — steady-state inference throughput of a uno RT model as a
//! function of batch size, for the manuscript's predictor benchmark.
//!
//! Loads the model once, warms up the thread pool on the largest batch, then
//! for each requested batch size times `predict_many` on the first N peptides
//! of the library REPEATS times. Prints one TSV row per repeat:
//! `n\trepeat\tseconds`. The Python runner computes medians and peptides/s so
//! both engines are summarised by the same code.
//!
//! Usage: rt_throughput MODEL LIBRARY.tsv REPEATS N1,N2,...

use anyhow::{Context, Result};
use std::time::Instant;
use uno_core::RtModel;

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let usage = "usage: rt_throughput MODEL LIBRARY.tsv REPEATS N1,N2,...";
    let model_path = args.get(1).context(usage)?;
    let library = args.get(2).context(usage)?;
    let repeats: usize = args.get(3).context(usage)?.parse()?;
    let sizes: Vec<usize> = args
        .get(4)
        .context(usage)?
        .split(',')
        .map(|s| s.trim().parse::<usize>())
        .collect::<std::result::Result<_, _>>()?;

    let model = RtModel::load(model_path)?;
    let seqs: Vec<String> = std::fs::read_to_string(library)?
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty() && *l != "sequence")
        .map(|l| l.split('\t').next().unwrap_or("").to_string())
        .collect();
    let max = *sizes.iter().max().context("no sizes")?;
    anyhow::ensure!(seqs.len() >= max, "library has {} peptides, need {}", seqs.len(), max);

    // warm-up: thread pool creation and first-touch allocations
    std::hint::black_box(model.predict_many(&seqs[..max]));

    println!("n\trepeat\tseconds");
    for &n in &sizes {
        for r in 0..repeats {
            let t = Instant::now();
            let p = model.predict_many(&seqs[..n]);
            let s = t.elapsed().as_secs_f64();
            std::hint::black_box(p);
            println!("{}\t{}\t{:.6}", n, r, s);
        }
    }
    Ok(())
}
