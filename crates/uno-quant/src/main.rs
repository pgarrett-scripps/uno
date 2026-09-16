use anyhow::Result;
use clap::Parser;
use std::path::PathBuf;
use uno_quant::{
    input::{self, ReferenceOptions},
    Config,
};
#[derive(Parser)]
#[command(
    version,
    about = "Experimental Rust DirectMS1Quant-style contrast analysis on native reference evidence"
)]
struct Args {
    #[arg(long,num_args=1..,required_unless_present="samples",conflicts_with="samples")]
    s1: Vec<PathBuf>,
    #[arg(long,num_args=1..,required_unless_present="samples",conflicts_with="samples")]
    s2: Vec<PathBuf>,
    /// Uno manifest: sample, condition (S1/S2), evidence, proteins. Experimental adaptation.
    #[arg(long, conflicts_with_all=["all_proteins","all_pfms","allowed_proteins","allowed_peptides"])]
    samples: Option<PathBuf>,
    /// Reject features supporting more than one admitted precursor identity (symmetric target/decoy).
    #[arg(long, requires = "samples")]
    exclusive_features: bool,
    /// Require this many observed runs in EACH condition before reference statistics.
    #[arg(long, requires = "samples", default_value_t = 0)]
    min_per_condition: usize,
    #[arg(long, default_value_t = 0.01)]
    max_protein_qvalue: f64,
    #[arg(long)]
    out: PathBuf,
    #[arg(long, default_value_t = 0)]
    min_samples: usize,
    #[arg(long, default_value_t = 7)]
    min_length: usize,
    #[arg(long, default_value_t = 2)]
    normalization: u8,
    #[arg(long, default_value_t = 2.)]
    fold_change: f64,
    #[arg(long)]
    absolute_fold_change: bool,
    #[arg(long, default_value_t = 0.05)]
    qvalue: f64,
    #[arg(long, default_value = "DECOY_")]
    decoy_prefix: String,
    #[arg(long, default_value_t = 42)]
    seed: u32,
    #[arg(long)]
    all_proteins: bool,
    #[arg(long)]
    all_pfms: bool,
    #[arg(long)]
    allowed_proteins: Option<PathBuf>,
    #[arg(long)]
    allowed_peptides: Option<PathBuf>,
}
fn main() -> Result<()> {
    let a = Args::parse();
    let cfg = Config {
        min_samples: a.min_samples,
        min_length: a.min_length,
        normalization: a.normalization,
        fold_change: a.fold_change,
        absolute_fold_change: a.absolute_fold_change,
        qvalue: a.qvalue,
        decoy_prefix: a.decoy_prefix,
        seed: a.seed,
    };
    let options = ReferenceOptions {
        all_proteins: a.all_proteins,
        all_pfms: a.all_pfms,
        allowed_proteins: a.allowed_proteins,
        allowed_peptides: a.allowed_peptides,
    };
    let (data, hashes, mode, admission) = if let Some(path) = &a.samples {
        let policy = uno_quant::uno_adapter::AdmissionPolicy {
            exclusive_features: a.exclusive_features,
            min_per_condition: a.min_per_condition,
        };
        let (data, hashes, audit) =
            uno_quant::uno_adapter::load_with_policy(path, &cfg, a.max_protein_qvalue, &policy)?;
        (data, hashes, "uno_ms1_adaptation", Some(audit))
    } else {
        let (data, hashes) = input::reference(&a.s1, &a.s2, &cfg, &options)?;
        (data, hashes, "native_directms1quant", None)
    };
    let analysis = uno_quant::analyze(&data, &cfg)?;
    uno_quant::output::write_with_admission(&a.out, &analysis, &hashes, mode, admission.as_ref())?;
    eprintln!(
        "{} eligible precursors; {} target proteins; {} significant (experimental, unvalidated)",
        analysis.eligible_precursors,
        analysis.proteins.len(),
        analysis
            .proteins
            .iter()
            .filter(|p| p.bh_pass && p.fc_pass)
            .count()
    );
    Ok(())
}
