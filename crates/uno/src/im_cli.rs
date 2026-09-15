//! IM uses the RT CNN/training engine, with an explicit charge feature and units.
use anyhow::{bail, ensure, Context, Result};
use clap::{Args, Subcommand, ValueEnum};
use serde_json::json;
use sha2::{Digest, Sha256};
use std::collections::{BTreeMap, HashSet};
use std::path::{Path, PathBuf};
use uno_core::im::ImModel;

#[derive(Debug, Subcommand)]
pub enum ImCommand {
    /// Train/adapt a charge-aware sequence CNN with the optional uno-rt package.
    Train(TrainArgs),
    /// Inspect the CNN and its training provenance without Python.
    Inspect { model: PathBuf },
    /// Evaluate independent observations with the Rust search predictor.
    Evaluate(EvaluateArgs),
    /// Predict a headered sequence/charge TSV with the Rust search predictor.
    Predict {
        #[arg(long)]
        model: PathBuf,
        #[arg(long)]
        data: PathBuf,
        #[arg(long)]
        out: PathBuf,
    },
}
#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum Units {
    InverseMobility,
    ProportionalInverseMobility,
}
impl Units {
    fn name(self) -> &'static str {
        match self {
            Self::InverseMobility => "inverse_mobility",
            Self::ProportionalInverseMobility => "proportional_inverse_mobility",
        }
    }
}
#[derive(Debug, Args)]
pub struct TrainArgs {
    #[arg(long)]
    data: PathBuf,
    #[arg(long)]
    out: PathBuf,
    #[arg(long, value_enum)]
    units: Units,
    #[arg(long)]
    method: String,
    #[arg(
        long,
        conflicts_with = "from_scratch",
        required_unless_present = "from_scratch"
    )]
    init: Option<PathBuf>,
    #[arg(long)]
    from_scratch: bool,
    #[arg(long, conflicts_with = "holdout_run")]
    validation: Option<PathBuf>,
    #[arg(long)]
    holdout_run: Option<String>,
    #[arg(long, default_value_t = 0)]
    seed: u32,
    #[arg(long, default_value_t = 40)]
    epochs: u32,
    #[arg(long)]
    config: Option<PathBuf>,
    #[arg(long, default_value = "python3")]
    python: PathBuf,
}
#[derive(Debug, Args)]
pub struct EvaluateArgs {
    #[arg(long)]
    model: PathBuf,
    #[arg(long)]
    data: PathBuf,
    #[arg(long, value_enum)]
    units: Units,
    /// Separate peptide-disjoint observations for per-charge multiplicative calibration.
    #[arg(long)]
    calibration_data: Option<PathBuf>,
    #[arg(long)]
    out_predictions: Option<PathBuf>,
}
struct Row {
    sequence: String,
    charge: u8,
    im: f64,
}
fn read(path: &Path, labels: bool) -> Result<Vec<Row>> {
    let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').from_path(path)?;
    let headers = reader.headers()?.clone();
    ensure!(
        headers.iter().collect::<HashSet<_>>().len() == headers.len(),
        "duplicate column headers"
    );
    let col = |name: &str| {
        headers
            .iter()
            .position(|h| h == name)
            .with_context(|| format!("missing {name} header"))
    };
    let (s, z) = (col("sequence")?, col("charge")?);
    let im = if labels { Some(col("im")?) } else { None };
    let mut grouped: BTreeMap<(String, u8), Vec<f64>> = BTreeMap::new();
    for record in reader.records() {
        let r = record?;
        let sequence = r[s].trim().to_string();
        ensure!(
            (1..=50).contains(&sequence.len())
                && sequence
                    .bytes()
                    .all(|b| b"ACDEFGHIKLMNPQRSTVWY".contains(&b)),
            "require 1–50 canonical residues with fixed CAM C"
        );
        let charge: u8 = r[z].parse()?;
        ensure!((1..=8).contains(&charge), "charge must be 1–8");
        let value = if let Some(i) = im {
            r[i].parse::<f64>()?
        } else {
            1.0
        };
        ensure!(
            value.is_finite() && value > 0.0,
            "IM must be finite and positive"
        );
        grouped.entry((sequence, charge)).or_default().push(value);
    }
    ensure!(!grouped.is_empty(), "no observations");
    Ok(grouped
        .into_iter()
        .map(|((sequence, charge), mut v)| Row {
            sequence,
            charge,
            im: median(&mut v),
        })
        .collect())
}
fn median(v: &mut [f64]) -> f64 {
    v.sort_by(f64::total_cmp);
    (v[(v.len() - 1) / 2] + v[v.len() / 2]) / 2.0
}
fn predictions(model: &ImModel, rows: &[Row]) -> Result<Vec<f64>> {
    let p = model.predict_many(
        &rows
            .iter()
            .map(|r| (r.sequence.as_str(), r.charge))
            .collect::<Vec<_>>(),
    );
    ensure!(
        p.iter().all(|v| v.is_finite() && *v > 0.0),
        "model produced invalid IM predictions"
    );
    Ok(p)
}
fn write(path: &Path, rows: &[Row], p: &[f64], labels: bool) -> Result<()> {
    let f = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)?;
    let mut w = csv::WriterBuilder::new().delimiter(b'\t').from_writer(f);
    if labels {
        w.write_record(["sequence", "charge", "observed_im", "predicted_im"])?;
    } else {
        w.write_record(["sequence", "charge", "predicted_im"])?;
    }
    for (r, p) in rows.iter().zip(p) {
        let mut values = vec![r.sequence.clone(), r.charge.to_string()];
        if labels {
            values.push(r.im.to_string());
        }
        values.push(p.to_string());
        w.write_record(values)?;
    }
    w.flush()?;
    Ok(())
}
fn cnn(path: &Path) -> Result<ImModel> {
    let m = ImModel::load(path)?;
    ensure!(m.cnn_info().is_some(), "IM commands require a UIMC CNN");
    Ok(m)
}

pub fn run(command: ImCommand) -> Result<()> {
    match command {
        ImCommand::Train(a) => {
            let mut c = std::process::Command::new(&a.python);
            c.args(["-m", "uno_rt.im", "--data"])
                .arg(a.data)
                .arg("--out")
                .arg(&a.out)
                .arg("--units")
                .arg(a.units.name())
                .arg("--method")
                .arg(a.method)
                .arg("--seed")
                .arg(a.seed.to_string())
                .arg("--epochs")
                .arg(a.epochs.to_string());
            if let Some(p) = a.init {
                c.arg("--init").arg(p);
            }
            if a.from_scratch {
                c.arg("--from-scratch");
            }
            if let Some(p) = a.validation {
                c.arg("--validation").arg(p);
            }
            if let Some(p) = a.holdout_run {
                c.arg("--holdout-run").arg(p);
            }
            if let Some(p) = a.config {
                c.arg("--config").arg(p);
            }
            ensure!(
                c.status()
                    .context("install optional uno-rt Python package and select --python")?
                    .success(),
                "IM training failed; see Python error above"
            );
            cnn(&a.out)?;
            Ok(())
        }
        ImCommand::Inspect { model } => {
            let m = cnn(&model)?;
            let mut info = m.cnn_info().unwrap().clone();
            if let Some(meta) = info["metadata"].as_object_mut() {
                if let Some(h) = meta.remove("seen_peptide_hashes") {
                    meta.insert(
                        "seen_peptide_count".into(),
                        json!(h.as_array().map_or(0, Vec::len)),
                    );
                }
            }
            println!("{}", serde_json::to_string_pretty(&info)?);
            Ok(())
        }
        ImCommand::Predict { model, data, out } => {
            let m = cnn(&model)?;
            let rows = read(&data, false)?;
            write(&out, &rows, &predictions(&m, &rows)?, false)
        }
        ImCommand::Evaluate(a) => {
            let m = cnn(&a.model)?;
            let rows = read(&a.data, true)?;
            let meta = &m.cnn_info().unwrap()["metadata"];
            let hashes: HashSet<&str> = meta["seen_peptide_hashes"]
                .as_array()
                .context("missing model lineage")?
                .iter()
                .filter_map(|v| v.as_str())
                .collect();
            for r in &rows {
                ensure!(
                    !hashes
                        .contains(format!("{:x}", Sha256::digest(r.sequence.as_bytes())).as_str()),
                    "evaluation peptide {} was used in training/model selection",
                    r.sequence
                );
            }
            let mut p = predictions(&m, &rows)?;
            let mut scales = BTreeMap::new();
            if let Some(path) = a.calibration_data {
                let cal = read(&path, true)?;
                let eval_seqs: HashSet<_> = rows.iter().map(|r| &r.sequence).collect();
                ensure!(
                    cal.iter().all(|r| !eval_seqs.contains(&r.sequence)),
                    "calibration/evaluation peptide overlap"
                );
                let cp = predictions(&m, &cal)?;
                let mut ratios: BTreeMap<u8, Vec<f64>> = BTreeMap::new();
                for (r, p) in cal.iter().zip(cp) {
                    ratios.entry(r.charge).or_default().push(r.im / p);
                }
                for (z, mut v) in ratios {
                    ensure!(
                        v.len() >= 2,
                        "need at least two calibration observations per charge"
                    );
                    scales.insert(z, median(&mut v));
                }
                for (r, p) in rows.iter().zip(&mut p) {
                    *p *= scales
                        .get(&r.charge)
                        .context("missing charge calibration")?;
                }
            } else if meta["im_units"] != a.units.name() {
                bail!("IM label units differ; supply separate --calibration-data");
            }
            let mut report =
                crate::rt_cli::metrics(&rows.iter().map(|r| r.im).collect::<Vec<_>>(), &p);
            report["units"] = json!(a.units.name());
            report["calibration_scales"] = json!(scales);
            if let Some(path) = a.out_predictions {
                write(&path, &rows, &p, true)?;
            }
            println!("{}", serde_json::to_string_pretty(&report)?);
            Ok(())
        }
    }
}
