//! uno-rt-train — Phase 4: distill a DeepLC RT model into a fast pure-Rust gbdt.
//!
//! Input is a `SEQUENCE\tRT` label file (DeepLC predictions, e.g. sherlock's
//! `rt_deeplc.tsv`). We featurize with the shared Klammer features, train a
//! gbdt regressor, report the held-out distillation residual (student vs
//! teacher), and serialize the model.

use anyhow::{bail, Result};
use clap::Parser;
use gbdt::config::Config;
use gbdt::decision_tree::{Data, DataVec};
use gbdt::gradient_boost::GBDT;
use std::path::PathBuf;
use std::time::Instant;
use uno_core::im::{im_features, read_im2deep_csv, IM_N_FEATURES};
use uno_core::rt::{read_label_tsv, rt_features, RT_N_FEATURES};

#[derive(Parser, Debug)]
#[command(name = "uno-rt-train", version)]
struct Cli {
    /// DeepLC label file: `SEQUENCE\tRT` per line.
    #[arg(long)]
    labels: PathBuf,

    /// IM (CCS) mode: train a charge-aware 1/K0 model from an im2deep `-i` output
    /// CSV (`modified_seq,charge,predicted IM`) instead of an RT label file.
    #[arg(long)]
    im_labels: Option<PathBuf>,
    /// Output model path.
    #[arg(long, default_value = "rt_model.gbdt")]
    out: PathBuf,
    /// Cap training rows (subsample for speed); 0 = use all.
    #[arg(long, default_value_t = 0)]
    max_train: usize,
    /// Fraction held out for validation.
    #[arg(long, default_value_t = 0.1)]
    val_frac: f64,
    #[arg(long, default_value_t = 300)]
    iterations: usize,
    #[arg(long, default_value_t = 5)]
    max_depth: u32,
    #[arg(long, default_value_t = 0.1)]
    shrinkage: f32,

    /// Predict mode: load `--out` as a model, read sequences (one per line) from
    /// this file, print `seq\tpred` to stdout. Skips training.
    #[arg(long)]
    predict_seqs: Option<PathBuf>,
}

/// Deterministic train/val split by a hash of the row index (no rng dep).
fn is_val(idx: usize, val_frac: f64) -> bool {
    // simple splitmix-ish hash -> [0,1)
    let mut x = idx as u64;
    x ^= x >> 30;
    x = x.wrapping_mul(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    let u = (x >> 11) as f64 / (1u64 << 53) as f64;
    u < val_frac
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    if let Some(seq_file) = &cli.predict_seqs {
        let model = uno_core::RtModel::load(&cli.out)?;
        let text = std::fs::read_to_string(seq_file)?;
        let seqs: Vec<&str> = text.lines().map(|l| l.trim()).filter(|s| !s.is_empty()).collect();
        let preds = model.predict_many(&seqs);
        let mut out = String::with_capacity(seqs.len() * 24);
        for (s, p) in seqs.iter().zip(&preds) {
            out.push_str(s);
            out.push('\t');
            out.push_str(&p.to_string());
            out.push('\n');
        }
        use std::io::Write;
        std::io::stdout().write_all(out.as_bytes())?;
        return Ok(());
    }

    // Load labels + featurize. IM mode is charge-aware (im2deep CSV); RT mode
    // is a plain SEQUENCE\tRT file.
    let t = Instant::now();
    let (feats, ys, feature_size): (Vec<Vec<f32>>, Vec<f64>, usize) = if let Some(p) = &cli.im_labels
    {
        let (seqs, charges, ims) = read_im2deep_csv(p)?;
        if seqs.is_empty() {
            bail!("no IM labels parsed from {}", p.display());
        }
        eprintln!("loaded {} (seq, charge, 1/K0) IM labels in {:.1}s", seqs.len(), t.elapsed().as_secs_f64());
        let feats = seqs
            .iter()
            .zip(&charges)
            .map(|(s, &z)| im_features(s, z).into_iter().map(|x| x as f32).collect())
            .collect();
        (feats, ims, IM_N_FEATURES)
    } else {
        let (seqs, rts) = read_label_tsv(&cli.labels)?;
        if seqs.is_empty() {
            bail!("no labels parsed from {}", cli.labels.display());
        }
        eprintln!("loaded {} (sequence, rt) labels in {:.1}s", seqs.len(), t.elapsed().as_secs_f64());
        let feats = seqs
            .iter()
            .map(|s| rt_features(s).into_iter().map(|x| x as f32).collect())
            .collect();
        (feats, rts, RT_N_FEATURES)
    };

    // Build train / val splits.
    let mut train: DataVec = Vec::new();
    let mut val_feat: Vec<Vec<f32>> = Vec::new();
    let mut val_y: Vec<f64> = Vec::new();
    for (i, (f, &y)) in feats.iter().zip(ys.iter()).enumerate() {
        if is_val(i, cli.val_frac) {
            val_feat.push(f.clone());
            val_y.push(y);
        } else {
            if cli.max_train != 0 && train.len() >= cli.max_train {
                continue;
            }
            train.push(Data::new_training_data(f.clone(), 1.0, y as f32, None));
        }
    }
    eprintln!("train rows: {}, val rows: {}", train.len(), val_y.len());

    let mut cfg = Config::new();
    cfg.set_feature_size(feature_size);
    cfg.set_max_depth(cli.max_depth);
    cfg.set_iterations(cli.iterations);
    cfg.set_shrinkage(cli.shrinkage);
    cfg.set_data_sample_ratio(1.0);
    cfg.set_feature_sample_ratio(1.0);
    cfg.set_loss("SquaredError");
    cfg.set_training_optimization_level(2);

    let t = Instant::now();
    let mut model = GBDT::new(&cfg);
    model.fit(&mut train);
    eprintln!("trained in {:.1}s", t.elapsed().as_secs_f64());

    // Validate: distillation residual vs the DeepLC teacher.
    let val_data: DataVec = val_feat
        .iter()
        .map(|f| Data::new_test_data(f.clone(), None))
        .collect();
    let preds = model.predict(&val_data);
    let (r2, mae, rmse) = regression_metrics(&val_y, &preds);
    eprintln!(
        "validation (student vs DeepLC):  R^2 = {:.4}   MAE = {:.2}   RMSE = {:.2}",
        r2, mae, rmse
    );

    model
        .save_model(&cli.out.to_string_lossy())
        .map_err(|e| anyhow::anyhow!("saving model: {e}"))?;
    eprintln!("saved model to {}", cli.out.display());
    Ok(())
}

fn regression_metrics(y: &[f64], pred: &[f32]) -> (f64, f64, f64) {
    let n = y.len() as f64;
    let mean = y.iter().sum::<f64>() / n;
    let mut ss_res = 0.0;
    let mut ss_tot = 0.0;
    let mut abs = 0.0;
    for (&yi, &pi) in y.iter().zip(pred.iter()) {
        let e = yi - pi as f64;
        ss_res += e * e;
        ss_tot += (yi - mean) * (yi - mean);
        abs += e.abs();
    }
    let r2 = 1.0 - ss_res / ss_tot;
    let mae = abs / n;
    let rmse = (ss_res / n).sqrt();
    (r2, mae, rmse)
}
