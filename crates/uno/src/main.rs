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

    // Multi-sample quant. Correct the peptide→feature assignment (feature
    // competition + optional profile reassignment), then emit the precursor LFQ
    // matrix and/or the protein-group rollup. Both are raw intensities: no design,
    // no differential testing — that's downstream stats' job.
    if cli.out_lfq.is_some() || cli.out_abundance.is_some() {
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
        if let Some(path) = &cli.out_lfq {
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
        if let Some(path) = &cli.out_abundance {
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
