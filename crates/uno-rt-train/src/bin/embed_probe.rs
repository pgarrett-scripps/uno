//! embed_probe — experiment: can a NONLINEAR readout (gbdt) on the frozen CNN
//! embedding, fit on a run's confident (peptide, observed-RT) pairs, recover the
//! run's RT well? Tests whether the pretrained embedding CONTAINS run-specific RT
//! information that a deeper-than-linear on-data adapter could exploit.
//!
//! Usage: embed_probe <model.urtc> <train.tsv> <test_seqs.txt>
//!   train.tsv: `SEQUENCE\tRT` (this run's confident anchors)
//!   prints `seq\tpred` for each test seq to stdout (eval externally with PWL-MAE).

use anyhow::{bail, Result};
use gbdt::config::Config;
use gbdt::decision_tree::{Data, DataVec};
use gbdt::gradient_boost::GBDT;
use uno_core::rt::read_label_tsv;
use uno_core::RtModel;

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 4 {
        bail!("usage: embed_probe <model.urtc> <train.tsv> <test_seqs.txt>");
    }
    let model = RtModel::load(&args[1])?;
    let (tr_seqs, tr_rt) = read_label_tsv(std::path::Path::new(&args[2]))?;
    let tr_emb = model
        .embed_many(&tr_seqs)
        .ok_or_else(|| anyhow::anyhow!("model has no embedding"))?;
    let dim = tr_emb.first().map_or(0, |v| v.len());
    eprintln!("embedding dim = {}, train rows = {}", dim, tr_seqs.len());

    let mut train: DataVec = Vec::new();
    for (e, &y) in tr_emb.iter().zip(&tr_rt) {
        let f: Vec<f32> = e.iter().map(|&x| x as f32).collect();
        train.push(Data::new_training_data(f, 1.0, y as f32, None));
    }

    let mut cfg = Config::new();
    cfg.set_feature_size(dim);
    cfg.set_max_depth(5);
    cfg.set_iterations(300);
    cfg.set_shrinkage(0.1);
    cfg.set_data_sample_ratio(1.0);
    cfg.set_feature_sample_ratio(1.0);
    cfg.set_loss("SquaredError");
    cfg.set_training_optimization_level(2);
    let mut gb = GBDT::new(&cfg);
    gb.fit(&mut train);

    let test_txt = std::fs::read_to_string(&args[3])?;
    let test_seqs: Vec<&str> = test_txt.lines().map(|l| l.trim()).filter(|s| !s.is_empty()).collect();
    let te_emb = model.embed_many(&test_seqs).unwrap();
    let te_data: DataVec = te_emb
        .iter()
        .map(|e| Data::new_test_data(e.iter().map(|&x| x as f32).collect(), None))
        .collect();
    let preds = gb.predict(&te_data);
    let mut out = String::new();
    for (s, p) in test_seqs.iter().zip(&preds) {
        out.push_str(s);
        out.push('\t');
        out.push_str(&p.to_string());
        out.push('\n');
    }
    print!("{}", out);
    Ok(())
}
