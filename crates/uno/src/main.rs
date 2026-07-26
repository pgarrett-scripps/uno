use anyhow::{bail, Context, Result};
use clap::Parser;
use std::path::PathBuf;
use std::time::Instant;
use std::collections::HashMap;
use uno_core::im::ImModel;
use uno_core::{
    match_features, match_features_topk, score,
    DigestConfig, FeatureMatrix, MatchConfig, PeptideDb,
    ProteinConfig, RtModel, ScoreConfig,
};

mod koth_ingest;

/// uno — MS1-only protein identification over a single-run biosaur2 feature file.
#[derive(Parser, Debug)]
#[command(name = "uno", version)]
struct Cli {
    /// Single-run biosaur2 feature file (`*.features.tsv`) — the same single-run
    /// input ms1searchpy consumes. Aligned/consensus matrices are NOT accepted:
    /// alignment drops the per-feature quality columns that separate target from
    /// decoy. Mutually exclusive with `--raw`; exactly one is required.
    #[arg(long, required_unless_present_any = ["raw", "matrix"], conflicts_with_all = ["raw", "matrix"])]
    features: Option<PathBuf>,

    /// Raw run to find features from in-process, with no intermediate biosaur2
    /// pass: mzML (`.mzML`/`.mzML.gz`), Bruker `.d`, or Thermo `.raw` (the last
    /// needs uno built with `--features thermo`). koth's feature finder runs with
    /// the bundled biosaur2-matching config unless the uno config has a `[koth]`
    /// table. Mutually exclusive with `--features`.
    #[arg(long)]
    raw: Option<PathBuf>,

    /// koth-aligned consensus matrix (`intensity_matrix.tsv`): one search over all
    /// runs at once, with each consensus feature carrying its per-run intensity
    /// vector for the multi-sample quant rollup. Mutually exclusive with
    /// `--features`/`--raw`.
    #[arg(long, conflicts_with_all = ["features", "raw"])]
    matrix: Option<PathBuf>,

    /// Drop consensus features detected in fewer than this many runs (only with
    /// `--matrix`). Required floor: alignment drops the per-feature isotope quality
    /// uno uses to reject junk, so without a reproducibility gate the random-match
    /// background swamps the search. Default 2.
    #[arg(long, default_value_t = 2)]
    matrix_min_runs: u32,

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

    /// Write a per-peptide intensity table for the multi-sample quant rollup (only
    /// meaningful with `--matrix`). One row per matched peptide (best-mass feature),
    /// with its protein(s), decoy flag, and per-run intensities. Targets AND decoys
    /// are emitted so the downstream DE test has its decoy null. Columns:
    /// peptide, decoy, proteins, charge, mass, rt, <run1..runN>.
    #[arg(long)]
    out_peptides: Option<PathBuf>,

    /// Write protein-level differential-expression + abundance to TSV (needs
    /// `--matrix`, `--design`, `--contrast`). Co-direction binomial + BH FDR over
    /// the ID-confident proteins; peptide FDR is not used.
    #[arg(long)]
    out_quant: Option<PathBuf>,

    /// Design file for quant: lines `<run-name-substring>\t<condition>`. A matrix
    /// run column joins a condition when the run field is a substring of its name.
    #[arg(long)]
    design: Option<PathBuf>,

    /// Two conditions to contrast for `--out-quant`, e.g. `A,E` (B vs A = second
    /// vs first, so a positive log2FC means higher in the second).
    #[arg(long)]
    contrast: Option<String>,

    /// Dilution-series mode for `--out-quant`: ordered `cond:level,...`
    /// (e.g. `A:3.3,B:10,C:33,D:67,E:90`). Scores each protein by how consistently
    /// its peptides track the level (Spearman), stronger than a single contrast.
    /// Takes precedence over `--contrast`.
    #[arg(long)]
    dilution: Option<String>,

    /// After the feature competition, re-route VARYING features to the protein
    /// whose across-sample profile they best correlate with (profile-guided
    /// reassignment). Uses the quant pattern to fix mass+RT ties the ID stage
    /// couldn't — recovers minority-species features mis-assigned to abundant
    /// backgrounds. `--out-quant` only.
    #[arg(long)]
    profile_reassign: bool,
}

#[derive(Debug, Default, serde::Deserialize)]
#[serde(default)]
struct Config {
    digest: DigestConfig,
    matching: MatchConfig,
    scoring: ScoreConfig,
    protein: ProteinConfig,
    /// Optional override for the built-in koth feature finder (`--raw` only).
    /// Absent → the bundled biosaur2-matching profile is used.
    koth: Option<koth_ff::config::KothConfig>,
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
    let matrix = match (&cli.features, &cli.raw, &cli.matrix) {
        (Some(features), None, None) => FeatureMatrix::load(features)?,
        (None, Some(raw), None) => {
            let koth_cfg = match &config.koth {
                Some(c) => c.clone(),
                None => koth_ingest::default_config()?,
            };
            koth_ingest::features_from_raw(raw, &koth_cfg)
                .with_context(|| format!("finding features in {}", raw.display()))?
        }
        (None, None, Some(mp)) => {
            let m = FeatureMatrix::load_matrix(mp, cli.matrix_min_runs)
                .with_context(|| format!("reading aligned matrix {}", mp.display()))?;
            log::info!(
                "aligned matrix: {} runs, {} consensus features (min_runs={})",
                m.n_samples(),
                m.len(),
                cli.matrix_min_runs
            );
            m
        }
        // clap's required_unless_present/conflicts_with guarantees exactly one.
        _ => bail!("provide exactly one of --features, --raw, or --matrix"),
    };
    log::info!(
        "loaded {} features from single run ({:.1}s)",
        matrix.len(),
        t.elapsed().as_secs_f64()
    );

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

    // Per-peptide intensity table for the multi-sample quant rollup. Uses the
    // RT+mass-gated feature uno actually trusts for each peptide (`result.peptides`),
    // NOT the best-mass feature — a mass-coincident contaminant would corrupt the
    // per-run intensity profile and flatten real fold changes. Targets + decoys
    // (the decoys are the DE test's null).
    if let Some(path) = &cli.out_peptides {
        use std::io::Write;
        let mut w = std::io::BufWriter::new(std::fs::File::create(path)?);
        write!(w, "peptide\tdecoy\tproteins\tcharge\tmass\trt")?;
        for rn in &matrix.run_names {
            write!(w, "\t{rn}")?;
        }
        writeln!(w)?;
        for hit in &result.peptides {
            let pep = &db.peptides[hit.peptide_idx as usize];
            let f = &matrix.features[hit.feature_idx as usize];
            let prots = pep
                .proteins
                .iter()
                .map(|p| p.as_ref())
                .collect::<Vec<_>>()
                .join(";");
            write!(
                w,
                "{}\t{}\t{}\t{}\t{:.5}\t{:.4}",
                pep.sequence.as_ref(),
                hit.decoy as u8,
                prots,
                f.charge,
                f.mass,
                f.rt
            )?;
            for v in &f.intensities {
                write!(w, "\t{:.5e}", v)?;
            }
            writeln!(w)?;
        }
        log::info!(
            "wrote {} RT-gated peptide-intensity rows to {}",
            result.peptides.len(),
            path.display()
        );
    }

    // Protein-level differential expression + abundance (multi-sample quant).
    if let Some(path) = &cli.out_quant {
        use std::io::Write;
        let design_path = cli
            .design
            .as_ref()
            .context("--out-quant requires --design")?;
        let design = std::fs::read_to_string(design_path)
            .with_context(|| format!("reading {}", design_path.display()))?;
        // ID-confident target proteins (q <= 1%) define the tested universe.
        let id_confident: std::collections::HashSet<String> = proteins
            .iter()
            .filter(|p| !p.decoy && p.q_value <= 0.01)
            .map(|p| p.accession.clone())
            .collect();
        // Feature-level competition (parsimony) for quant: one feature quantifies
        // for only its best-scoring protein, so a minority peptide can't measure a
        // mass-coincident feature owned by an abundant background protein.
        let prot_score: std::collections::HashMap<(bool, &str), f64> = proteins
            .iter()
            .map(|p| ((p.decoy, p.accession.as_str()), p.binomial_score))
            .collect();
        let mut quant_peptides =
            uno_core::quant::resolve_feature_competition(&result.peptides, &db, &prot_score);
        log::info!(
            "quant feature competition: {} -> {} peptide-features after 1-feature-1-protein claim",
            result.peptides.len(),
            quant_peptides.len()
        );
        // Optional profile-guided reassignment: re-route varying features to the
        // protein whose across-sample profile they best track (mass-constrained,
        // decoy-calibrated). Candidates = mono PFMs of scored proteins.
        if cli.profile_reassign {
            let pfms = match_features(&matrix, &db, &config.matching);
            let candidates: Vec<(u32, u32, bool)> = pfms
                .iter()
                .filter(|p| p.isotope_offset == 0)
                .filter(|p| {
                    db.peptides[p.peptide_idx as usize]
                        .proteins
                        .iter()
                        .any(|pr| prot_score.contains_key(&(p.decoy, pr.as_ref())))
                })
                .map(|p| (p.feature_idx, p.peptide_idx, p.decoy))
                .collect();
            let before = quant_peptides.len();
            quant_peptides =
                uno_core::quant::profile_reassign(&matrix, &db, &quant_peptides, &candidates, 0.5, 0.25);
            log::info!(
                "profile-reassign: {} candidates, {} -> {} peptide-features",
                candidates.len(),
                before,
                quant_peptides.len()
            );
        }
        let (quant, col_a, col_b) = if let Some(spec) = &cli.dilution {
            let levels = uno_core::quant::resolve_levels(&matrix.run_names, &design, spec)?;
            log::info!("quant dilution: {} levels, {} runs", levels.len(), levels.iter().map(|(r, _)| r.len()).sum::<usize>());
            let q = uno_core::quant::monotonic(&matrix, &db, &quant_peptides, &id_confident, &levels, 2);
            (q, "low".to_string(), "high".to_string())
        } else {
            let contrast = cli.contrast.as_ref().context("--out-quant requires --contrast or --dilution")?;
            let (ca, cb) = contrast.split_once(',').context("--contrast must be `condA,condB`")?;
            let (runs_a, runs_b) = uno_core::quant::resolve_contrast(&matrix.run_names, &design, ca, cb)?;
            log::info!("quant contrast {ca}(n={}) vs {cb}(n={})", runs_a.len(), runs_b.len());
            let q = uno_core::quant::differential(&matrix, &db, &quant_peptides, &id_confident, &runs_a, &runs_b, 2);
            (q, ca.to_string(), cb.to_string())
        };
        let stat_col = if cli.dilution.is_some() { "median_rho" } else { "log2fc" };
        let mut w = std::io::BufWriter::new(std::fs::File::create(path)?);
        writeln!(
            w,
            "accession\tdecoy\tn_peptides\tn_significant\tk\t{stat_col}\tabundance_{col_a}\tabundance_{col_b}\tp_value\tq_value"
        )?;
        for q in &quant {
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.4e}\t{:.4e}\t{:.3e}\t{:.4e}",
                q.accession,
                q.decoy as u8,
                q.n_peptides,
                q.n_significant,
                q.k,
                q.log2_fc,
                q.abundance_a,
                q.abundance_b,
                q.p_value,
                q.q_value
            )?;
        }
        let de = quant.iter().filter(|q| !q.decoy && q.q_value <= 0.05).count();
        log::info!(
            "wrote {} proteins to {} ({} differentially expressed at 5% FDR)",
            quant.len(),
            path.display(),
            de
        );
    }

    Ok(())
}
