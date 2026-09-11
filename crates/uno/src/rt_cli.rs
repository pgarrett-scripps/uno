//! Public CNN workflow. Only training launches optional Python; inspection and
//! independent evaluation use exactly the Rust predictor used during search.
use anyhow::{bail, Context, Result};
use clap::{Args, Subcommand, ValueEnum};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::collections::{BTreeMap, HashSet};
use std::path::{Path, PathBuf};
use uno_core::RtModel;

#[derive(Debug, Subcommand)]
pub enum Command {
    /// Train, inspect, or evaluate a retention-time CNN.
    Rt {
        #[command(subcommand)]
        command: RtCommand,
    },
}

#[derive(Debug, Subcommand)]
pub enum RtCommand {
    /// Adapt the bundled CNN, or explicitly train from scratch (requires uno-rt Python package).
    Train(TrainArgs),
    /// Inspect a CNN's architecture, normalization, and training provenance (no Python).
    Inspect { model: PathBuf },
    /// Evaluate independent, headered sequence/rt data using Rust inference (no Python).
    Evaluate(EvaluateArgs),
}

#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum Units {
    Minutes,
    Seconds,
    Irt,
    Hi,
}
impl Units {
    fn name(self) -> &'static str {
        match self {
            Self::Minutes => "minutes",
            Self::Seconds => "seconds",
            Self::Irt => "irt",
            Self::Hi => "hi",
        }
    }
}

#[derive(Debug, Args)]
pub struct TrainArgs {
    /// Headered TSV: sequence, rt, optional run_id. One separation method per model.
    #[arg(long)]
    data: PathBuf,
    /// New model file (existing files are never overwritten).
    #[arg(long)]
    out: PathBuf,
    /// Explicit label units; multiple raw-time runs must first be normalized.
    #[arg(long, value_enum)]
    units: Units,
    /// Descriptive column/separation name, recorded with the model.
    #[arg(long)]
    method: String,
    /// Starting CNN; architecture and normalization are read from the file.
    #[arg(long, conflicts_with = "from_scratch")]
    init: Option<PathBuf>,
    /// Initialize a fresh CNN instead of adapting the bundled model.
    #[arg(long)]
    from_scratch: bool,
    /// Separate validation TSV; peptides must not overlap training or parent lineage.
    #[arg(long, conflicts_with = "holdout_run")]
    validation: Option<PathBuf>,
    /// Hold out this run_id; its peptides are also removed from other training runs.
    #[arg(long)]
    holdout_run: Option<String>,
    #[arg(long, default_value_t = 0)]
    seed: u32,
    #[arg(long, default_value_t = 40)]
    epochs: u32,
    /// Advanced JSON: lr, batch_size, patience, val_fraction, device, threads, architecture.
    #[arg(long)]
    config: Option<PathBuf>,
    /// Python interpreter with the optional uno-rt package installed.
    #[arg(long, default_value = "python3")]
    python: PathBuf,
}

#[derive(Debug, Args)]
pub struct EvaluateArgs {
    #[arg(long)]
    model: PathBuf,
    /// Independent TSV with sequence/rt headers, unused in training or model selection.
    #[arg(long)]
    data: PathBuf,
    #[arg(long, value_enum)]
    units: Units,
    /// Separate observations for affine calibration; cannot overlap evaluation peptides.
    #[arg(long)]
    calibration_data: Option<PathBuf>,
    /// Optional starting CNN to compare; requires separate calibration observations.
    #[arg(long, requires = "calibration_data")]
    baseline: Option<PathBuf>,
    /// Save sequence/observed_rt/predicted_rt for inspection or parity checks.
    #[arg(long)]
    out_predictions: Option<PathBuf>,
}

pub fn run(command: Command) -> Result<()> {
    let Command::Rt { command } = command;
    match command {
        RtCommand::Train(args) => train(args),
        RtCommand::Inspect { model } => {
            let model = cnn(&model)?;
            let mut info = model.cnn_info().unwrap().clone();
            if let Some(metadata) = info["metadata"].as_object_mut() {
                if let Some(hashes) = metadata.remove("seen_peptide_hashes") {
                    metadata.insert(
                        "seen_peptide_count".into(),
                        json!(hashes.as_array().map_or(0, Vec::len)),
                    );
                }
            }
            println!("{}", serde_json::to_string_pretty(&info)?);
            Ok(())
        }
        RtCommand::Evaluate(args) => evaluate(args),
    }
}

fn train(args: TrainArgs) -> Result<()> {
    // The installed Python package may be from a different release. Always
    // adapt this binary's actual default when no starting model was selected.
    let bundled = if args.init.is_none() && !args.from_scratch {
        use std::io::Write;
        let mut file =
            tempfile::NamedTempFile::new().context("preparing bundled CNN for training")?;
        file.write_all(uno_core::rt::DEFAULT_RT_CNN)?;
        Some(file)
    } else {
        None
    };
    // Argument arrays avoid shell interpolation, including for paths/method names.
    let mut command = std::process::Command::new(&args.python);
    command
        .args(["-m", "uno_rt", "train", "--data"])
        .arg(&args.data)
        .arg("--out")
        .arg(&args.out)
        .arg("--units")
        .arg(args.units.name())
        .arg("--method")
        .arg(&args.method)
        .arg("--seed")
        .arg(args.seed.to_string())
        .arg("--epochs")
        .arg(args.epochs.to_string());
    if let Some(p) = args
        .init
        .as_deref()
        .or_else(|| bundled.as_ref().map(|f| f.path()))
    {
        command.arg("--init").arg(p);
    }
    if args.from_scratch {
        command.arg("--from-scratch");
    }
    if let Some(p) = args.validation {
        command.arg("--validation").arg(p);
    }
    if let Some(run) = args.holdout_run {
        command.arg("--holdout-run").arg(run);
    }
    if let Some(p) = args.config {
        command.arg("--config").arg(p);
    }
    let status = command.status().with_context(|| format!(
        "starting {}; install the optional trainer with `python -m pip install ./crates/uno-rt-train/python` and select that interpreter with --python",
        args.python.display()))?;
    if !status.success() {
        bail!("RT training failed ({status}); see the Python error above. The optional trainer is installed with `python -m pip install ./crates/uno-rt-train/python`; --python selects its interpreter");
    }
    // Verify the exported artifact with the production Rust loader as well.
    cnn(&args.out).context("training output failed Rust compatibility validation")?;
    Ok(())
}

fn cnn(path: &Path) -> Result<RtModel> {
    let model = RtModel::load(path)?;
    if model.cnn_info().is_none() {
        bail!("uno rt supports URTC CNNs only; legacy loaders remain available for search");
    }
    Ok(model)
}

#[derive(Debug)]
struct Observations {
    sequences: Vec<String>,
    rts: Vec<f64>,
}

fn read_data(path: &Path, units: Units) -> Result<Observations> {
    let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').from_path(path)?;
    let headers = reader.headers()?.clone();
    if headers.iter().collect::<HashSet<_>>().len() != headers.len() {
        bail!("duplicate TSV headers");
    }
    let seq_col = headers
        .iter()
        .position(|h| h == "sequence")
        .context("TSV requires sequence header")?;
    let rt_col = headers
        .iter()
        .position(|h| h == "rt")
        .context("TSV requires rt header")?;
    let run_col = headers.iter().position(|h| h == "run_id");
    let mut groups = BTreeMap::<(String, String), Vec<f64>>::new();
    let mut runs = HashSet::new();
    for (i, record) in reader.records().enumerate() {
        let record = record?;
        let sequence = record[seq_col].trim();
        if sequence.is_empty()
            || sequence.len() > 50
            || !sequence
                .bytes()
                .all(|c| b"ACDEFGHIKLMNPQRSTVWY".contains(&c))
        {
            bail!(
                "{}:{}: require 1–50 canonical amino acids; variable modifications are unsupported",
                path.display(),
                i + 2
            );
        }
        let rt: f64 = record[rt_col]
            .trim()
            .parse()
            .with_context(|| format!("invalid RT at row {}", i + 2))?;
        if !rt.is_finite() || (matches!(units, Units::Minutes | Units::Seconds) && rt < 0.) {
            bail!("invalid RT at row {}", i + 2);
        }
        let run = run_col.map_or("", |c| record[c].trim());
        if run_col.is_some() && run.is_empty() {
            bail!("run_id must be populated at row {}", i + 2);
        }
        runs.insert(run.to_owned());
        groups
            .entry((sequence.to_owned(), run.to_owned()))
            .or_default()
            .push(rt);
    }
    if matches!(units, Units::Minutes | Units::Seconds) && runs.len() > 1 {
        bail!("normalize multiple runs to shared irt/hi units before pooling");
    }
    if groups.len() < 2 {
        bail!("need at least two observations");
    }
    let mut sequences = Vec::new();
    let mut rts = Vec::new();
    for ((sequence, _), mut values) in groups {
        values.sort_by(f64::total_cmp);
        let n = values.len();
        let rt = if n % 2 == 0 {
            values[n / 2 - 1] / 2. + values[n / 2] / 2.
        } else {
            values[n / 2]
        };
        sequences.push(sequence);
        rts.push(rt);
    }
    Ok(Observations { sequences, rts })
}

fn check_unseen(model: &RtModel, data: &Observations) -> Result<()> {
    let meta = &model.cnn_info().unwrap()["metadata"];
    if let Some(seen) = meta["seen_peptide_hashes"].as_array() {
        let hashes: HashSet<&str> = seen.iter().filter_map(Value::as_str).collect();
        if data
            .sequences
            .iter()
            .any(|s| hashes.contains(format!("{:x}", Sha256::digest(s.as_bytes())).as_str()))
        {
            bail!("evaluation peptides overlap this model's recorded training/validation lineage; provide independent data");
        }
    } else {
        log::warn!(
            "legacy model has no peptide provenance; evaluation independence cannot be checked"
        );
    }
    Ok(())
}

fn affine(predictions: &[f64], y: &[f64]) -> Result<(f64, f64)> {
    let n = y.len() as f64;
    let px = predictions.iter().sum::<f64>() / n;
    let py = y.iter().sum::<f64>() / n;
    let variance = predictions.iter().map(|p| (p - px).powi(2)).sum::<f64>();
    if variance <= 1e-12 {
        bail!("calibration predictions have no variance; add representative peptides");
    }
    let slope = predictions
        .iter()
        .zip(y)
        .map(|(p, y)| (p - px) * (y - py))
        .sum::<f64>()
        / variance;
    Ok((slope, py - slope * px))
}

fn predictions(
    model: &RtModel,
    data: &Observations,
    units: Units,
    calibration: Option<&Observations>,
) -> Result<Vec<f64>> {
    let (slope, intercept) = if let Some(calibration) = calibration {
        affine(
            &model.predict_many(&calibration.sequences),
            &calibration.rts,
        )?
    } else {
        let model_units = model.cnn_info().unwrap()["metadata"]["rt_units"]
            .as_str()
            .unwrap_or("unknown");
        match (model_units, units.name()) {
            (a, b) if a == b => (1., 0.),
            ("seconds", "minutes") => (1. / 60., 0.),
            ("minutes", "seconds") => (60., 0.),
            _ => bail!("model units are {model_units}, data units are {}; provide separate --calibration-data to compare these scales", units.name()),
        }
    };
    let result: Vec<f64> = model
        .predict_many(&data.sequences)
        .into_iter()
        .map(|p| slope * p + intercept)
        .collect();
    if result.iter().any(|p| !p.is_finite()) {
        bail!("nonfinite model predictions");
    }
    Ok(result)
}

fn metrics(y: &[f64], predictions: &[f64]) -> Value {
    let n = y.len() as f64;
    let mean = y.iter().sum::<f64>() / n;
    let mut errors: Vec<f64> = y
        .iter()
        .zip(predictions)
        .map(|(a, b)| (a - b).abs())
        .collect();
    let mse = errors.iter().map(|e| e * e).sum::<f64>() / n;
    let variance = y.iter().map(|y| (y - mean).powi(2)).sum::<f64>() / n;
    errors.sort_by(f64::total_cmp);
    let index = 0.95 * (errors.len() - 1) as f64;
    let low = index.floor() as usize;
    let high = index.ceil() as usize;
    json!({"observations": y.len(), "mae": errors.iter().sum::<f64>() / n, "rmse": mse.sqrt(),
        "r2": if variance > 0. { Some(1. - mse / variance) } else { None },
        "p95_absolute_error": errors[low] + (errors[high] - errors[low]) * index.fract()})
}

fn evaluate(args: EvaluateArgs) -> Result<()> {
    let data = read_data(&args.data, args.units)?;
    let calibration = args
        .calibration_data
        .as_ref()
        .map(|p| read_data(p, args.units))
        .transpose()?;
    if let Some(calibration) = &calibration {
        let held: HashSet<_> = data.sequences.iter().collect();
        if calibration.sequences.iter().any(|s| held.contains(s)) {
            bail!("calibration and evaluation peptides overlap");
        }
    }
    let model = cnn(&args.model)?;
    check_unseen(&model, &data)?;
    let predicted = predictions(&model, &data, args.units, calibration.as_ref())?;
    let mut report = json!({"units": args.units.name(), "calibrated": calibration.is_some(),
        "model": metrics(&data.rts, &predicted),
        "independence": "Checked against recorded lineage only; legacy/pretrained corpora may be unknown.",
        "interpretation": "RT error alone does not establish search quality; validate identifications and FDR before adoption."});
    if let Some(path) = &args.baseline {
        let baseline = cnn(path)?;
        check_unseen(&baseline, &data)?;
        report["baseline"] = metrics(
            &data.rts,
            &predictions(&baseline, &data, args.units, calibration.as_ref())?,
        );
    }
    if let Some(path) = args.out_predictions {
        let file = std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(path)?;
        let mut writer = csv::WriterBuilder::new().delimiter(b'\t').from_writer(file);
        writer.write_record(["sequence", "observed_rt", "predicted_rt"])?;
        for ((s, y), p) in data.sequences.iter().zip(&data.rts).zip(&predicted) {
            writer.write_record([s.as_str(), &y.to_string(), &p.to_string()])?;
        }
        writer.flush()?;
    }
    println!("{}", serde_json::to_string_pretty(&report)?);
    Ok(())
}
