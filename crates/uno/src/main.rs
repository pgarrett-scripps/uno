use anyhow::{Context, Result};
use clap::Parser;
use std::path::PathBuf;
use std::time::Instant;
use std::collections::HashMap;
use uno_core::im::ImModel;
use uno_core::{
    match_features, match_features_topk, score,
    DedupConfig, DigestConfig, FeatureMatrix, MatchConfig, PeptideDb,
    ProteinConfig, RtModel, ScoreConfig,
};

/// uno — MS1-only protein identification + quantitation over a koth aligned matrix.
#[derive(Parser, Debug)]
#[command(name = "uno", version)]
struct Cli {
    /// koth `intensity_matrix.tsv`.
    #[arg(long)]
    intensity_matrix: PathBuf,

    /// koth `qvalue_matrix.tsv` (optional; enables per-sample q filtering).
    #[arg(long)]
    qvalue_matrix: Option<PathBuf>,

    /// Protein FASTA database.
    #[arg(long)]
    fasta: PathBuf,

    /// Distilled RT model (`uno-rt-train` output). Without it, RT is not used.
    #[arg(long)]
    rt_model: Option<PathBuf>,

    /// Distilled IM/CCS model (`uno-rt-train --im-labels` output). Used only when
    /// the data has ion mobility.
    #[arg(long)]
    im_model: Option<PathBuf>,

    /// TOML config (digest + matching + scoring). Defaults are used if omitted.
    #[arg(long)]
    config: Option<PathBuf>,

    /// Precomputed RT library TSV (`peptide<TAB>rt`), e.g. from DeepLC. When given,
    /// these replace the distilled model's predictions as the abstract RT (still
    /// calibrated on-data). DeepLC-quality RT is what makes RT discriminate.
    #[arg(long)]
    rt_library: Option<PathBuf>,

    /// Write the unique candidate peptide sequences (one per line) to this file
    /// and exit. Use to build a DeepLC RT library, then re-run with --rt-library.
    #[arg(long)]
    dump_candidate_peptides: Option<PathBuf>,

    /// Write protein-level identifications to TSV.
    #[arg(long)]
    out_proteins: Option<PathBuf>,
}

#[derive(Debug, Default, serde::Deserialize)]
#[serde(default)]
struct Config {
    dedup: DedupConfig,
    digest: DigestConfig,
    matching: MatchConfig,
    scoring: ScoreConfig,
    protein: ProteinConfig,
}

fn main() -> Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    let cli = Cli::parse();

    let config: Config = match &cli.config {
        Some(p) => {
            let s = std::fs::read_to_string(p).with_context(|| format!("reading {}", p.display()))?;
            toml::from_str(&s).with_context(|| "parsing config")?
        }
        None => Config::default(),
    };

    let t = Instant::now();
    let mut matrix = FeatureMatrix::load(&cli.intensity_matrix, cli.qvalue_matrix.as_ref())?;
    log::info!(
        "loaded {} aligned features across {} samples ({:.1}s)",
        matrix.len(),
        matrix.n_samples(),
        t.elapsed().as_secs_f64()
    );
    if config.dedup.enabled {
        let before = matrix.len();
        let removed = matrix.dedup(&config.dedup);
        log::info!(
            "dedup: removed {} near-duplicate features ({} -> {})",
            removed,
            before,
            matrix.len()
        );
    }

    let t = Instant::now();
    let db = PeptideDb::build(&cli.fasta, &config.digest)?;
    log::info!(
        "built peptide db: {} targets + {} decoys ({:.1}s)",
        db.n_targets,
        db.n_decoys,
        t.elapsed().as_secs_f64()
    );

    // The candidate-peptide dump needs the full (unreduced) PFM set so the RT
    // library covers every match; the scoring path only ever consumes the top-k
    // per feature, so we stream that reduction during generation and never
    // materialise the full ~41M-PFM vector (~1.3GB on run10). `reduce_top_k`
    // inside `score` is idempotent on this already-reduced set, so output is
    // byte-identical to generating-then-reducing.
    if let Some(path) = &cli.dump_candidate_peptides {
        let t = Instant::now();
        let pfms = match_features(&matrix, &db, &config.matching);
        log::info!(
            "generated {} PFMs at {} ppm, offsets {:?} ({:.1}s)",
            pfms.len(),
            config.matching.ppm_tolerance,
            config.matching.isotope_offsets,
            t.elapsed().as_secs_f64()
        );
        use std::io::Write;
        let mut seen = std::collections::HashSet::new();
        let mut w = std::io::BufWriter::new(std::fs::File::create(path)?);
        // Include decoys: they need RT predictions too, else they escape the RT
        // penalty and the target-decoy FDR collapses.
        for p in &pfms {
            let seq = db.peptides[p.peptide_idx as usize].sequence.as_ref();
            if seen.insert(seq.to_string()) {
                writeln!(w, "{}", seq)?;
            }
        }
        log::info!("dumped {} unique candidate peptides to {}", seen.len(), path.display());
        return Ok(());
    }

    let t = Instant::now();
    let pfms = match_features_topk(&matrix, &db, &config.matching, config.scoring.top_k);
    let n_decoy = pfms.iter().filter(|p| p.decoy).count();
    log::info!(
        "generated {} reduced PFMs ({} target, {} decoy) at {} ppm, offsets {:?}, top-{} ({:.1}s)",
        pfms.len(),
        pfms.len() - n_decoy,
        n_decoy,
        config.matching.ppm_tolerance,
        config.matching.isotope_offsets,
        config.scoring.top_k,
        t.elapsed().as_secs_f64()
    );

    // Optional precomputed RT library (DeepLC). peptide<TAB>rt per line.
    let rt_library: Option<HashMap<String, f64>> = match &cli.rt_library {
        Some(p) => {
            let s = std::fs::read_to_string(p).with_context(|| format!("reading {}", p.display()))?;
            let mut map = HashMap::new();
            for line in s.lines() {
                let mut it = line.split('\t');
                if let (Some(seq), Some(rt)) = (it.next(), it.next()) {
                    if let Ok(rt) = rt.parse::<f64>() {
                        map.insert(seq.to_string(), rt);
                    }
                }
            }
            log::info!("loaded RT library with {} peptides from {}", map.len(), p.display());
            Some(map)
        }
        None => None,
    };

    let rt_model = match &cli.rt_model {
        Some(p) => {
            let m = RtModel::load(p)?;
            log::info!("loaded RT model from {}", p.display());
            Some(m)
        }
        None => {
            if rt_library.is_none() {
                log::warn!("no --rt-model or --rt-library given; scoring without retention time");
            }
            None
        }
    };

    let im_model = match &cli.im_model {
        Some(p) => {
            let m = ImModel::load(p)?;
            log::info!("loaded IM model from {}", p.display());
            Some(m)
        }
        None => None,
    };

    let t = Instant::now();
    let result = score(
        &matrix,
        &db,
        &pfms,
        rt_model.as_ref(),
        im_model.as_ref(),
        rt_library.as_ref(),
        &config.scoring,
    );
    log::info!("scored in {:.1}s", t.elapsed().as_secs_f64());

    let proteins = result.proteins;
    let n_prot_1pct = proteins
        .iter()
        .filter(|p| !p.decoy && p.q_value <= 0.01)
        .count();
    let n_prot_decoy_1pct = proteins
        .iter()
        .filter(|p| p.decoy && p.q_value <= 0.01)
        .count();
    log::info!(
        "proteins: {} target at 1% FDR ({} decoy proteins below 1%)",
        n_prot_1pct,
        n_prot_decoy_1pct
    );

    if let Some(path) = &cli.out_proteins {
        use std::io::Write;
        let mut w = std::io::BufWriter::new(std::fs::File::create(path)?);
        writeln!(
            w,
            "accession\tdecoy\tobserved_peptides\ttheoretical_peptides\tcoverage\tbinomial_score\tmean_rt_diff\tq_value"
        )?;
        for p in &proteins {
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{:.4}\t{:.3}\t{:.3}\t{:.5}",
                p.accession,
                p.decoy as u8,
                p.observed_peptides,
                p.theoretical_peptides,
                p.coverage,
                p.binomial_score,
                p.mean_rt_diff,
                p.q_value
            )?;
        }
        log::info!("wrote {} proteins to {}", proteins.len(), path.display());
    }

    Ok(())
}
