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
mod rt_cli;
mod im_cli;
mod quant_export;

#[cfg(test)]
mod rt_command_tests {
    use super::*;

    #[test]
    fn rt_commands_do_not_require_search_inputs() {
        for argv in [
            vec!["uno", "rt", "inspect", "custom.urtc"],
            vec!["uno", "rt", "train", "--data", "labels.tsv", "--out", "custom.urtc", "--units", "minutes", "--method", "C18"],
            vec!["uno", "rt", "evaluate", "--model", "custom.urtc", "--data", "heldout.tsv", "--units", "minutes"],
        ] { assert!(Cli::try_parse_from(argv).is_ok()); }
        assert!(Cli::try_parse_from(["uno", "--features", "run.tsv"]).is_err());
        assert!(Cli::try_parse_from(["uno", "--fasta", "db.fasta"]).is_err());
    }

    #[test]
    fn im_defaults_to_bundled_with_explicit_disable_and_override() {
        let base = ["uno", "--features", "run.tsv", "--fasta", "db.fasta"];
        let cli = Cli::try_parse_from(base).unwrap();
        assert!(load_models(&cli).unwrap().im.unwrap().cnn_info().is_some());
        let cli = Cli::try_parse_from(base.into_iter().chain(["--no-im"])).unwrap();
        assert!(load_models(&cli).unwrap().im.is_none());
        assert!(Cli::try_parse_from(base.into_iter().chain(["--no-im", "--im-model", "custom.uimc"])).is_err());
        let cli = Cli::try_parse_from(base.into_iter().chain(["--im-model", "missing.uimc"])).unwrap();
        assert!(load_models(&cli).is_err());
    }

    #[test]
    fn long_and_batch_inputs_are_explicit_and_exclusive() {
        for args in [vec!["uno","--matrix-long","bundle","--fasta","db"],vec!["uno","--matrix-long","bundle","--fasta","db","--matrix-batch-out","out"]] {assert!(Cli::try_parse_from(args).is_ok());}
        assert!(Cli::try_parse_from(["uno","--matrix-long","bundle","--features","f","--fasta","db"]).is_err());
        assert!(Cli::try_parse_from(["uno","--matrix-long","bundle","--matrix-batch-out","out","--out-proteins","p","--fasta","db"]).is_err());
    }

    #[test]
    fn joint_mode_requires_long_input_and_cannot_be_a_batch_search() {
        assert!(Cli::try_parse_from(["uno", "--matrix-long", "b", "--matrix-joint", "--fasta", "db"]).is_ok());
        assert!(Cli::try_parse_from(["uno", "--matrix", "b", "--matrix-joint", "--fasta", "db"]).is_err());
        assert!(Cli::try_parse_from(["uno", "--matrix-long", "b", "--matrix-joint", "--matrix-batch-out", "o", "--fasta", "db"]).is_err());
    }

    #[test]
    fn rt_sources_are_mutually_exclusive_and_default_needs_no_flag() {
        assert!(Cli::try_parse_from(["uno", "--features", "run.tsv", "--fasta", "db.fasta"]).is_ok());
        assert!(Cli::try_parse_from(["uno", "--features", "run.tsv", "--fasta", "db.fasta", "--rt-model", "a", "--rt-library", "b"]).is_err());
        assert!(Cli::try_parse_from(["uno", "rt", "train", "--data", "d", "--out", "o", "--units", "irt", "--method", "m", "--init", "a", "--from-scratch"]).is_err());
    }
}

/// uno — MS1-only protein identification over a single-run biosaur2 feature file.
#[derive(Parser, Debug)]
#[command(name = "uno", version, args_conflicts_with_subcommands = true, subcommand_negates_reqs = true)]
struct Cli {
    #[command(subcommand)]
    command: Option<rt_cli::Command>,
    /// Single-run feature table. Use --matrix for wide consensus input or
    /// --matrix-long for an evidence bundle. Select exactly one input route.
    #[arg(long, required_unless_present_any = ["raw", "matrix", "matrix_long"], conflicts_with_all = ["raw", "matrix", "matrix_long"])]
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

    /// Versioned KOTH long bundle manifest (or its directory).
    #[arg(long, conflicts_with_all = ["features", "raw", "matrix"])]
    matrix_long: Option<PathBuf>,
    /// Experimental: one joint search using all original aligned matrix members.
    #[arg(long, requires="matrix_long", conflicts_with_all=["features","raw","matrix","matrix_batch_out","matrix_no_feature_quality","matrix_legacy_values","matrix_prefer_isotope_support","dump_candidate_peptides"])]
    matrix_joint: bool,
    /// Search every original run in one process with a shared database/prediction cache.
    #[arg(long, requires="matrix_long", conflicts_with_all=["out_proteins","out_peptides","out_lfq","out_abundance","dump_candidate_peptides","matrix_no_feature_quality","matrix_legacy_values"])]
    matrix_batch_out: Option<PathBuf>,
    /// Diagnostic: disable original feature-quality descriptors on long input.
    #[arg(long, requires = "matrix_long")]
    matrix_no_feature_quality: bool,
    /// Diagnostic: reproduce legacy wide TSV coordinate/intensity rounding.
    #[arg(long, requires = "matrix_long")]
    matrix_legacy_values: bool,
    /// Diagnostic: use the best original member with >=3 isotopes when available.
    #[arg(long, requires="matrix_long", conflicts_with_all=["matrix_no_feature_quality","matrix_batch_out"])]
    matrix_prefer_isotope_support: bool,


    /// Minimum original-run support: defaults to 1 for --matrix-joint, 2 for
    /// representative matrix search. MBR-filled cells never count as support.
    #[arg(long)]
    matrix_min_runs: Option<u32>,

    /// Protein FASTA database.
    #[arg(long, required = true)]
    fasta: Option<PathBuf>,

    /// Custom URTC CNN. Defaults to the bundled CNN; legacy models are deprecated.
    #[arg(long, conflicts_with = "rt_library")]
    rt_model: Option<PathBuf>,

    /// Custom IM model. Defaults to the bundled sequence CNN; IM is a soft
    /// rescoring feature predicted after RT shortlisting, never a hard cutoff.
    #[arg(long, conflicts_with = "no_im")]
    im_model: Option<PathBuf>,

    /// Disable IM prediction and scoring, retaining the RT-only workflow.
    #[arg(long)]
    no_im: bool,

    /// TOML config (digest + matching + scoring). Defaults are used if omitted.
    #[arg(long)]
    config: Option<PathBuf>,

    /// External predictions (`peptide<TAB>rt`), calibrated to the run at search time.
    /// Must cover candidate target AND decoy sequences; replaces the bundled CNN.
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

    /// Export single-run target/decoy MS1 evidence for experimental uno-quant.
    /// Preserves feature IDs and mobility; does not claim peptide-level q-values.
    #[arg(long, conflicts_with_all = ["matrix", "matrix_long", "matrix_batch_out"])]
    out_quant_evidence: Option<PathBuf>,

    /// Write the precursor × sample LFQ matrix to TSV (`--matrix` only), Sage-lfq
    /// style: one row per precursor (peptide + charge) with its intensity in each
    /// run, and the protein(s) it maps to — NO rollup. Uses the corrected
    /// peptide→feature assignment (competition + optional profile reassignment).
    /// Downstream tools roll up to protein and do the stats.
    #[arg(long)]
    out_lfq: Option<PathBuf>,

    /// Write the protein-group × sample abundance matrix to TSV (`--matrix` only):
    /// one row per confident protein group, its indistinguishable members, and its
    /// summed intensity in each run. Optional protein-level rollup on top of the
    /// precursor `--out-lfq`. Normalization, differential testing, and
    /// interpretation are the job of downstream stats tools.
    #[arg(long)]
    out_abundance: Option<PathBuf>,

    /// Before building the abundance matrix, re-route VARYING features to the
    /// protein whose across-sample profile they best correlate with (profile-guided
    /// reassignment). Uses the quant pattern to fix mass+RT ties the ID stage
    /// couldn't — recovers minority-species features mis-assigned to abundant
    /// backgrounds. Design-free (uses across-sample variation, not labels).
    #[arg(long)]
    profile_reassign: bool,

    /// Seed for every result-affecting random component (decoy shuffle, RANSAC
    /// sampling, CV fold assignment, booster subsampling). Output is a pure
    /// function of (input, seed). Default 0 = the shipped streams; other values
    /// give an independent draw, for measuring seed-to-seed spread. Overrides
    /// `digest.decoy_seed` and `scoring.seed` in the config file. Explicit
    /// `scoring.rt_seed`, `scoring.cv_seed`, and `scoring.classifier_seed` overrides
    /// still take precedence for their individual stages.
    #[arg(long)]
    seed: Option<u64>,
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

impl Config {
    fn apply_cli_seed(&mut self, seed: Option<u64>) {
        if let Some(seed) = seed {
            self.digest.decoy_seed = seed;
            self.scoring.seed = seed;
        }
    }
}

fn main() -> Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    let cli = Cli::parse();
    if let Some(command) = cli.command {
        return rt_cli::run(command);
    }

    let mut config: Config = match &cli.config {
        Some(p) => {
            let s = std::fs::read_to_string(p).with_context(|| format!("reading {}", p.display()))?;
            toml::from_str(&s).with_context(|| "parsing config")?
        }
        None => Config::default(),
    };
    config.apply_cli_seed(cli.seed);

    if let Some(out)=&cli.matrix_batch_out { return search_batch(&cli,&config,out); }

    let matrix_min_runs = cli.matrix_min_runs.unwrap_or(if cli.matrix_joint { 1 } else { 2 });
    let t = Instant::now();
    let mut matrix = if let Some(path) = &cli.matrix_long {
        let bundle = uno_core::long_input::Bundle::open(path)?;
        if cli.matrix_joint { bundle.load_joint(matrix_min_runs)? }
        else { bundle.load_consensus_with_representative(matrix_min_runs, !cli.matrix_no_feature_quality, cli.matrix_legacy_values, cli.matrix_prefer_isotope_support)? }
    } else { match (&cli.features, &cli.raw, &cli.matrix) {
        (Some(features), None, None) => FeatureMatrix::load(features)?,
        (None, Some(raw), None) => {
            let koth_cfg = match &config.koth {
                Some(c) => c.clone(),
                None => koth_ingest::default_config_for(raw)?,
            };
            // A Bruker run needs BOTH halves: koth's timsTOF profile and uno's
            // feature-quality gate. Without the gate the search is not merely
            // shallower, it reports nothing at all (114,659 features -> background
            // p 0.554 -> 0 proteins; gated -> 17,587 -> p 0.206 -> ~300 proteins),
            // so defaulting it is the difference between `--raw sample.d` working
            // and silently returning an empty result. Only applied when BOTH knobs
            // are untouched: set either one explicitly to take full control (and
            // note that means `min_isotopes = 0` alone cannot express "gate off on
            // Bruker" — pair it with a nonzero `min_combined_score`, or vice versa).
            if koth_ingest::is_bruker(raw) {
                log::info!("raw input is a Bruker .d: using the bundled timsTOF koth profile");
                if config.matching.min_isotopes == 0 && config.matching.min_combined_score == 0.0 {
                    config.matching.min_isotopes = 4;
                    config.matching.min_combined_score = 0.85;
                    log::info!(
                        "Bruker defaults: matching.min_isotopes=4, matching.min_combined_score=0.85 \
                         (timsTOF mass accuracy is ~3 ppm, so the search must be restricted to \
                         well-characterised features; override in [matching])"
                    );
                }
            }
            koth_ingest::features_from_raw(raw, &koth_cfg)
                .with_context(|| format!("finding features in {}", raw.display()))?
        }
        (None, None, Some(mp)) => {
            let m = FeatureMatrix::load_matrix(mp, matrix_min_runs)
                .with_context(|| format!("reading aligned matrix {}", mp.display()))?;
            log::info!(
                "aligned matrix: {} runs, {} consensus features (min_runs={})",
                m.n_samples(),
                m.len(),
                matrix_min_runs
            );
            m
        }
        // clap's required_unless_present/conflicts_with guarantees exactly one.
        _ => bail!("provide exactly one of --features, --raw, or --matrix"),
    } };
    uno_core::joint::restrict_observations(&mut matrix, &config.matching);
    log::info!(
        "loaded {} features from input ({:.1}s)",
        matrix.len(),
        t.elapsed().as_secs_f64()
    );

    let t = Instant::now();
    let db = PeptideDb::build(cli.fasta.as_ref().context("--fasta is required for search")?, &config.digest)?;
    log::info!(
        "built peptide db: {} targets + {} decoys ({:.1}s)",
        db.n_targets,
        db.n_decoys,
        t.elapsed().as_secs_f64()
    );

    let models=load_models(&cli)?;
    search_loaded(&cli,&config,&matrix,&db,&models,None)
}

struct Models {
    rt: Option<RtModel>,
    im: Option<ImModel>,
    library: Option<HashMap<String, f64>>,
}
fn load_models(cli: &Cli) -> Result<Models> {
    if cli.dump_candidate_peptides.is_some() {
        return Ok(Models {
            rt: None,
            im: None,
            library: None,
        });
    }
    // Optional precomputed RT library (DeepLC). peptide<TAB>rt per line.
    let rt_library: Option<HashMap<String, f64>> = match &cli.rt_library {
        Some(p) => {
            let s =
                std::fs::read_to_string(p).with_context(|| format!("reading {}", p.display()))?;
            let mut map = HashMap::new();
            for line in s.lines() {
                let mut it = line.split('\t');
                if let (Some(seq), Some(rt)) = (it.next(), it.next()) {
                    if let Ok(rt) = rt.parse::<f64>() {
                        map.insert(seq.to_string(), rt);
                    }
                }
            }
            log::info!(
                "loaded RT library with {} peptides from {}",
                map.len(),
                p.display()
            );
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
        None if rt_library.is_none() => {
            log::info!("using bundled RT CNN");
            Some(RtModel::bundled()?)
        }
        None => None,
    };

    let im_model = match &cli.im_model {
        Some(p) => {
            let m = ImModel::load(p)?;
            log::info!("loaded IM model from {}", p.display());
            Some(m)
        }
        None if cli.no_im => None,
        None => {
            log::info!("using bundled IM CNN (soft scoring after RT shortlisting)");
            Some(ImModel::bundled()?)
        }
    };

    Ok(Models {
        rt: rt_model,
        im: im_model,
        library: rt_library,
    })
}

fn search_loaded(cli:&Cli,config:&Config,matrix:&FeatureMatrix,db:&PeptideDb,models:&Models,dir:Option<&std::path::Path>)->Result<()> {
    let out_proteins=dir.map(|p|p.join("proteins.tsv")).or_else(||cli.out_proteins.clone());
    let out_peptides=dir.map(|p|p.join("peptides.tsv")).or_else(||cli.out_peptides.clone());
    let out_lfq=dir.map(|p|p.join("lfq.tsv")).or_else(||cli.out_lfq.clone());
    let out_abundance=dir.map(|p|p.join("abundance.tsv")).or_else(||cli.out_abundance.clone());
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

    let t = Instant::now();
    let result = score(
        &matrix,
        &db,
        &pfms,
        models.rt.as_ref(),
        models.im.as_ref(),
        models.library.as_ref(),
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

    if let Some(path) = &out_proteins {
        use std::io::Write;
        // Group indistinguishable proteins so the output shows the leading protein
        // AND its indistinguishable twins, instead of silently reporting one.
        let confident_all: std::collections::HashSet<String> =
            proteins.iter().filter(|p| !p.decoy).map(|p| p.accession.clone()).collect();
        let (groups, _) = uno_core::group::build_groups(&db, &result.peptides, &confident_all);
        let members = uno_core::group::members_by_leading(&groups);
        let mut w = std::io::BufWriter::new(std::fs::File::create(path)?);
        writeln!(
            w,
            "accession\tdecoy\tn_group_members\tgroup_members\tobserved_peptides\ttheoretical_peptides\tcoverage\tbinomial_score\tmean_rt_diff\tq_value"
        )?;
        for p in &proteins {
            let mem = members.get(&(p.decoy, p.accession.clone()));
            let n_mem = mem.map(|m| m.len()).unwrap_or(1);
            let mem_str = mem.map(|m| m.join(";")).unwrap_or_else(|| p.accession.clone());
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.3}\t{:.3}\t{:.5}",
                p.accession,
                p.decoy as u8,
                n_mem,
                mem_str,
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
    if let Some(path) = &out_peptides {
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

    if let Some(path) = &cli.out_quant_evidence {
        quant_export::write(path, &matrix, &db, &result.peptides)?;
    }

    // Multi-sample quant. Correct the peptide→feature assignment (feature
    // competition + optional profile reassignment), then emit the precursor LFQ
    // matrix and/or the protein-group rollup. Both are raw intensities: no design,
    // no differential testing — that's downstream stats' job.
    if out_lfq.is_some() || out_abundance.is_some() {
        use std::io::Write;
        // Feature-level competition (parsimony): one feature quantifies for only
        // its best-scoring protein, so a minority peptide can't measure a
        // mass-coincident feature owned by an abundant background protein.
        let prot_score: std::collections::HashMap<(bool, &str), f64> = proteins
            .iter()
            .map(|p| ((p.decoy, p.accession.as_str()), p.binomial_score))
            .collect();
        let mut quant_peptides =
            uno_core::quant::resolve_feature_competition(&result.peptides, &db, &prot_score);
        log::info!(
            "feature competition: {} -> {} peptide-features after 1-feature-1-protein claim",
            result.peptides.len(),
            quant_peptides.len()
        );
        // Optional profile-guided reassignment (design-free: uses across-sample
        // variation, not labels). Candidates = mono PFMs of scored proteins.
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
            log::info!("profile-reassign: {} candidates, {} -> {} peptide-features", candidates.len(), before, quant_peptides.len());
        }

        // Precursor × sample LFQ matrix (Sage-lfq style): one row per assigned
        // precursor (peptide + charge), no rollup.
        if let Some(path) = &out_lfq {
            let mut w = std::io::BufWriter::new(std::fs::File::create(path)?);
            write!(w, "peptide\tcharge\tdecoy\tproteins\tmass\trt")?;
            for rn in &matrix.run_names {
                write!(w, "\t{rn}")?;
            }
            writeln!(w)?;
            for h in &quant_peptides {
                let pep = &db.peptides[h.peptide_idx as usize];
                let f = &matrix.features[h.feature_idx as usize];
                let prots = pep.proteins.iter().map(|p| p.as_ref()).collect::<Vec<_>>().join(";");
                write!(
                    w,
                    "{}\t{}\t{}\t{}\t{:.5}\t{:.4}",
                    pep.sequence.as_ref(), f.charge, h.decoy as u8, prots, f.mass, f.rt
                )?;
                for v in &f.intensities {
                    write!(w, "\t{:.5e}", v)?;
                }
                writeln!(w)?;
            }
            log::info!("wrote {} precursor LFQ rows to {}", quant_peptides.len(), path.display());
        }

        // Optional protein-group rollup.
        if let Some(path) = &out_abundance {
            let id_confident: std::collections::HashSet<String> = proteins
                .iter()
                .filter(|p| !p.decoy && p.q_value <= 0.01)
                .map(|p| p.accession.clone())
                .collect();
            let (groups, pep_group) =
                uno_core::group::build_groups(&db, &quant_peptides, &id_confident);
            let members = uno_core::group::members_by_leading(&groups);
            let n_runs = matrix.run_names.len();
            let ab = uno_core::quant::group_abundance(&matrix, &quant_peptides, &pep_group, n_runs);
            let mut w = std::io::BufWriter::new(std::fs::File::create(path)?);
            write!(w, "accession\tn_group_members\tgroup_members\tn_peptides")?;
            for rn in &matrix.run_names {
                write!(w, "\t{rn}")?;
            }
            writeln!(w)?;
            let mut rows: Vec<(&String, &(Vec<f64>, std::collections::HashSet<u32>))> = ab.iter().collect();
            rows.sort_by(|a, b| a.0.cmp(b.0));
            for (leading, (ints, peps)) in rows {
                let mem = members.get(&(false, leading.to_string()));
                let n_mem = mem.map(|m| m.len()).unwrap_or(1);
                let mem_str = mem.map(|m| m.join(";")).unwrap_or_else(|| leading.clone());
                write!(w, "{}\t{}\t{}\t{}", leading, n_mem, mem_str, peps.len())?;
                for v in ints {
                    write!(w, "\t{:.5e}", v)?;
                }
                writeln!(w)?;
            }
            log::info!("wrote {} protein groups to {}", ab.len(), path.display());
        }
    }

    Ok(())
}

#[cfg(test)]
mod seed_config_tests {
    use super::*;

    #[test]
    fn cli_seed_keeps_legacy_joint_behavior() {
        let cli = Cli::try_parse_from([
            "uno", "--features", "run.tsv", "--fasta", "db.fasta", "--seed", "9",
        ]).unwrap();
        let mut cfg: Config = toml::from_str("[digest]\ndecoy_seed = 2\n[scoring]\nseed = 3\n").unwrap();
        cfg.apply_cli_seed(cli.seed);
        assert_eq!(cfg.digest.decoy_seed, 9);
        assert_eq!(cfg.scoring.seed, 9);
        assert_eq!(cfg.scoring.effective_rt_seed(), 9);
        assert_eq!(cfg.scoring.effective_cv_seed(), 9);
        assert_eq!(cfg.scoring.effective_classifier_seed(), 9);
    }

    #[test]
    fn explicit_stage_seeds_survive_cli_override() {
        let mut cfg: Config = toml::from_str(
            "[digest]\ndecoy_seed = 2\n[scoring]\nseed = 3\nrt_seed = 0\ncv_seed = 5\nclassifier_seed = 6\n",
        ).unwrap();
        cfg.apply_cli_seed(Some(9));
        assert_eq!(cfg.digest.decoy_seed, 9);
        assert_eq!(cfg.scoring.seed, 9);
        assert_eq!(cfg.scoring.effective_rt_seed(), 0);
        assert_eq!(cfg.scoring.effective_cv_seed(), 5);
        assert_eq!(cfg.scoring.effective_classifier_seed(), 6);
    }

    #[test]
    fn absent_cli_seed_preserves_independent_config_seeds() {
        let mut cfg: Config = toml::from_str("[digest]\ndecoy_seed = 2\n[scoring]\nseed = 3\ncv_seed = 5\n").unwrap();
        cfg.apply_cli_seed(None);
        assert_eq!(cfg.digest.decoy_seed, 2);
        assert_eq!(cfg.scoring.effective_rt_seed(), 3);
        assert_eq!(cfg.scoring.effective_cv_seed(), 5);
        assert_eq!(cfg.scoring.effective_classifier_seed(), 3);
    }
}

fn search_batch(cli: &Cli, config: &Config, out: &std::path::Path) -> Result<()> {
    use std::io::Write;
    let bundle = uno_core::long_input::Bundle::open(cli.matrix_long.as_ref().unwrap())?;
    std::fs::create_dir_all(out)?;
    anyhow::ensure!(
        !out.join("batch_manifest.json").exists(),
        "batch output already complete; select a new output directory"
    );
    let db = PeptideDb::build(cli.fasta.as_ref().unwrap(), &config.digest)?;
    let models = load_models(cli)?;
    let mut runs = Vec::new();
    let mut union: std::collections::BTreeMap<String, Vec<usize>> =
        std::collections::BTreeMap::new();
    uno_core::rt_cache::with_memory_cache(|| -> Result<()> {
        for r in &bundle.runs {
            let started = Instant::now();
            let matrix = bundle.load_run(r.run_id)?;
            let dir = out.join(format!("run{:03}", r.run_id));
            std::fs::create_dir_all(&dir)?;
            search_loaded(cli, config, &matrix, &db, &models, Some(&dir))?;
            let mut reader = csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .from_path(dir.join("proteins.tsv"))?;
            let header = reader.headers()?.clone();
            let col = |name| header.iter().position(|h| h == name).unwrap();
            for row in reader.records() {
                let row = row?;
                if &row[col("decoy")] == "0" && row[col("q_value")].parse::<f64>()? <= 0.01 {
                    for acc in row[col("group_members")].split(';') {
                        let ids = union.entry(acc.to_owned()).or_default();
                        if !ids.contains(&r.run_id) {
                            ids.push(r.run_id);
                        }
                    }
                }
            }
            runs.push(serde_json::json!({"run_id":r.run_id,"name":r.name,"seconds":started.elapsed().as_secs_f64(),"n_features":matrix.len(),"directory":dir}));
        }
        Ok(())
    })?;
    let mut writer =
        std::io::BufWriter::new(std::fs::File::create(out.join("accession_union.tsv"))?);
    writeln!(writer, "accession\trun_ids")?;
    for (acc, ids) in &union {
        writeln!(
            writer,
            "{}\t{}",
            acc,
            ids.iter()
                .map(|i| i.to_string())
                .collect::<Vec<_>>()
                .join(";")
        )?;
    }
    writer.flush()?;
    std::fs::write(
        out.join("batch_manifest.json"),
        serde_json::to_vec_pretty(
            &serde_json::json!({"runs":runs,"union_accessions":union.len(),"interpretation":"Union of per-run calls at q <= 0.01; no calibrated cohort q-value is claimed."}),
        )?,
    )?;
    Ok(())
}
