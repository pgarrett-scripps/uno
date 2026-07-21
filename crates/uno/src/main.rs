use anyhow::{Context, Result};
use clap::Parser;
use std::path::PathBuf;
use std::time::Instant;
use std::collections::HashMap;
use uno_core::im::ImModel;
use uno_core::{
    hill_search, infer_proteins, load_hills, match_features, match_features_topk, peptide_quant,
    protein_quant, score,
    DedupConfig, DigestConfig, FeatureMatrix, HillConfig, HillIndex, MatchConfig, PeptideDb,
    PeptideMatch, ProteinConfig, QuantRow, RtModel, ScoreConfig,
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

    /// Write all scored feature identifications (one row per feature) to TSV.
    #[arg(long)]
    out_peptides: Option<PathBuf>,

    /// Write protein-level identifications to TSV.
    #[arg(long)]
    out_proteins: Option<PathBuf>,

    /// Write the "matched peptides for identified proteins" set to TSV — the
    /// distinct peptides that map to a confident (q<=1%) protein and have a
    /// match within 3σ mass and 3σ RT. This is ms1searchpy's headline
    /// "Matched peptides for top-scored proteins" metric (not peptide-FDR
    /// controlled), reported for a like-for-like comparison.
    #[arg(long)]
    out_matched_peptides: Option<PathBuf>,

    /// Write protein quant matrix (protein x sample) to TSV.
    #[arg(long)]
    out_protein_quant: Option<PathBuf>,

    /// Write peptide quant matrix (feature ID x sample) to TSV.
    #[arg(long)]
    out_peptide_quant: Option<PathBuf>,

    /// FDR threshold for quant rollup output.
    #[arg(long, default_value_t = 0.01)]
    quant_fdr: f64,

    /// koth `hills.tsv` for the peptide-centric second pass. After the feature
    /// pass fits the on-data RT calibration, every peptide (target + decoy) is
    /// scored against the raw hills at its predicted m/z + RT; target/decoy
    /// competition -> q-values; survivors feed the protein rollup. Recovers
    /// low-abundance peptides that never assembled into a feature.
    #[arg(long)]
    hills: Option<PathBuf>,
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
    let ids = &result.ids;
    let n_target_1pct = ids
        .iter()
        .filter(|i| !i.decoy && i.q_value <= 0.01)
        .count();
    log::info!(
        "scored: {} target IDs at 1% FDR ({:.1}s)",
        n_target_1pct,
        t.elapsed().as_secs_f64()
    );

    if let Some(path) = &cli.out_peptides {
        use std::io::Write;
        let mut w = std::io::BufWriter::new(std::fs::File::create(path)?);
        writeln!(
            w,
            "feature_idx\tpeptide\tdecoy\tscore\tq_value\tmass_error_ppm\tisotope_offset\trt_diff\tim_diff\tcharge\trt\tcombined_score\tn_contributing_runs"
        )?;
        for id in ids {
            let f = &matrix.features[id.feature_idx as usize];
            let pep = &db.peptides[id.peptide_idx as usize];
            writeln!(
                w,
                "{}\t{}\t{}\t{:.4}\t{:.4}\t{:.3}\t{}\t{:.3}\t{:.4}\t{}\t{:.3}\t{:.4}\t{}",
                id.feature_idx as usize,
                pep.sequence,
                id.decoy as u8,
                id.score,
                id.q_value,
                id.mass_error_ppm,
                id.isotope_offset,
                id.rt_diff,
                id.im_diff,
                f.charge,
                f.rt,
                f.combined_score,
                f.n_contributing_runs
            )?;
        }
        log::info!("wrote {} feature identifications to {}", ids.len(), path.display());
    }

    // Phase 6: protein-level aggregation + FDR. When the mass-anchor mass+RT
    // binomial identifier ran (single-run / DirectMS1 mode), it already produced
    // the confident protein set — use it directly; otherwise fall back to the
    // aligned-matrix protein inference.
    let proteins = match result.proteins.clone() {
        Some(p) => {
            log::info!("using mass-anchor protein identifications ({} entries)", p.len());
            p
        }
        None => infer_proteins(&db, &result.ids, &result.peptide_matches, &config.protein),
    };
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

    // Matched peptides for identified proteins — ms1searchpy's "Matched peptides
    // for top-scored proteins". Distinct target peptides that (a) map to a
    // confident (q<=1%) protein and (b) have ANY feature within 3σ mass and 3σ
    // RT (ms1searchpy credits a peptide via any PFM in its window, not only its
    // best match). This is a lenient, non-peptide-FDR count — the number
    // ms1searchpy prints as its headline peptide figure — reported here so the
    // two tools can be compared on the same metric.
    let matched_peptides: Vec<usize> = {
        use std::collections::HashSet;
        let confident: HashSet<&str> = proteins
            .iter()
            .filter(|p| !p.decoy && p.q_value <= 0.01)
            .map(|p| p.accession.as_str())
            .collect();
        let maps_to_confident = |peptide_idx: usize| -> bool {
            db.peptides[peptide_idx]
                .proteins
                .iter()
                .any(|a| confident.contains(a.as_ref()))
        };
        // Robust RT sigma from confident-protein target matches. The representative
        // PFM per peptide is picked by MIN COMBINED (mass+RT) deviation -- the SAME
        // criterion mass-anchor's own scoring uses ("min-deviation any-feature", see
        // score.rs's `run_deviation`) to decide which feature actually counts for a
        // peptide. `result.peptide_matches` instead picks by mass alone (deliberately,
        // to avoid biasing ITS other consumer, `infer_proteins`, in the non-mass-anchor
        // path -- left untouched here). Mass-only selection often finds the WRONG
        // feature among several similarly-massed candidates for a given peptide,
        // understating RT's real accuracy; using `result.pfm_dev` (every candidate
        // PFM, not just one representative) to pick by combined deviation instead
        // reports the RT residual of the SAME evidence mass-anchor actually credited.
        // half-normal median relation median(|x|) = 0.6745 σ.
        let rt_sig_for_gate = result.rt_sigma.max(1e-6);
        let mut best_dev: std::collections::HashMap<usize, (f64, f64)> =
            std::collections::HashMap::new();
        for &(peptide_idx, decoy, abs_mass_z, abs_rt_diff) in &result.pfm_dev {
            if decoy || abs_mass_z > 3.0 || !abs_rt_diff.is_finite() {
                continue;
            }
            if !maps_to_confident(peptide_idx) {
                continue;
            }
            let dev = abs_mass_z * abs_mass_z + (abs_rt_diff / rt_sig_for_gate).powi(2);
            let e = best_dev.entry(peptide_idx).or_insert((f64::INFINITY, f64::INFINITY));
            if dev < e.0 {
                *e = (dev, abs_rt_diff);
            }
        }
        let mut resid: Vec<f64> = best_dev.values().map(|&(_, ard)| ard).collect();
        let rt_gate = if resid.is_empty() {
            f64::INFINITY
        } else {
            resid.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let med = resid[resid.len() / 2];
            3.0 * (med / 0.6745).max(1e-6)
        };
        // Any-feature rule: over EVERY reduced PFM, keep each distinct target
        // peptide sequence that has at least one feature inside both windows and
        // maps to a confident protein.
        let mut seen_seq: HashSet<&str> = HashSet::new();
        let mut out = Vec::new();
        for &(peptide_idx, decoy, abs_mass_z, abs_rt_diff) in &result.pfm_dev {
            if decoy || abs_mass_z > 3.0 || abs_rt_diff > rt_gate {
                continue;
            }
            if maps_to_confident(peptide_idx)
                && seen_seq.insert(db.peptides[peptide_idx].sequence.as_ref())
            {
                out.push(peptide_idx);
            }
        }
        log::info!(
            "matched peptides for identified proteins: {} (RT gate ±{:.2} min = 3σ, mass ≤3σ)",
            out.len(),
            rt_gate
        );
        out
    };

    if let Some(path) = &cli.out_matched_peptides {
        use std::io::Write;
        let mut w = std::io::BufWriter::new(std::fs::File::create(path)?);
        writeln!(w, "peptide\tproteins")?;
        for &pi in &matched_peptides {
            let pep = &db.peptides[pi];
            let prots: Vec<&str> = pep.proteins.iter().map(|a| a.as_ref()).collect();
            writeln!(w, "{}\t{}", pep.sequence, prots.join(";"))?;
        }
        log::info!("wrote {} matched peptides to {}", matched_peptides.len(), path.display());
    }

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

    // Optional peptide-centric hill second pass. Reuses the feature pass's
    // on-data RT calibration (result.rt_cal): predict each peptide's RT with
    // uno's own CNN (or the RT library), score it against the raw hills, run
    // target/decoy competition, then re-roll proteins on the survivors.
    if let Some(hpath) = &cli.hills {
        match result.rt_cal.as_ref() {
            None => log::warn!("--hills given but no RT calibration available; skipping hill pass"),
            Some(cal) => {
                // Abstract RT per peptide: RT library if present, else the CNN.
                let abstract_rt: Vec<Option<f64>> = if let Some(lib) = rt_library.as_ref() {
                    db.peptides
                        .iter()
                        .map(|p| lib.get(p.sequence.as_ref()).copied())
                        .collect()
                } else if let Some(m) = rt_model.as_ref() {
                    let seqs: Vec<&str> = db.peptides.iter().map(|p| p.sequence.as_ref()).collect();
                    m.predict_many(&seqs).into_iter().map(Some).collect()
                } else {
                    vec![None; db.peptides.len()]
                };
                let pred_rt: Vec<Option<f64>> =
                    abstract_rt.iter().map(|a| a.map(|x| cal.apply(x))).collect();

                let t = Instant::now();
                let hills = load_hills(hpath)?;
                let index = HillIndex::new(hills);
                let mut hcfg = HillConfig::default();
                // Quick experiment knobs (env overrides) for sweeping the scorer.
                let envf = |k: &str, d: f64| std::env::var(k).ok().and_then(|s| s.parse().ok()).unwrap_or(d);
                hcfg.min_isotopes = envf("HILL_MIN_ISO", hcfg.min_isotopes as f64) as u32;
                hcfg.rt_window_sigma = envf("HILL_RT_WIN", hcfg.rt_window_sigma);
                hcfg.min_coelution = envf("HILL_MIN_COEL", hcfg.min_coelution);
                hcfg.w_envelope = envf("HILL_W_ENV", hcfg.w_envelope);
                hcfg.w_coelution = envf("HILL_W_COEL", hcfg.w_coelution);
                log::info!(
                    "hill cfg: min_iso={} rt_win={}σ min_coel={} w_env={} w_coel={} rt_sigma={:.3}min",
                    hcfg.min_isotopes, hcfg.rt_window_sigma, hcfg.min_coelution,
                    hcfg.w_envelope, hcfg.w_coelution, result.rt_sigma
                );
                let hits = hill_search(&db, &index, &pred_rt, result.rt_sigma, &hcfg);
                let mut pep_t = std::collections::HashSet::new();
                let mut pep_d = std::collections::HashSet::new();
                for h in &hits {
                    if h.decoy { pep_d.insert(h.peptide_idx as usize); } else { pep_t.insert(h.peptide_idx as usize); }
                }
                log::info!(
                    "hill pass: {} hills; {} clusters over {} target / {} decoy peptides ({:.1}s)",
                    index.len(),
                    hits.len(),
                    pep_t.len(),
                    pep_d.len(),
                    t.elapsed().as_secs_f64()
                );

                // Diagnostic: dump per-hit features + labels so we can check whether
                // target vs decoy hill-hits are separable at all (Percolator viability).
                if let Ok(path) = std::env::var("UNO_DUMP_HILL_FEATS") {
                    use std::io::Write;
                    let conf: std::collections::HashSet<&str> = proteins
                        .iter()
                        .filter(|p| !p.decoy && p.q_value <= 0.01)
                        .map(|p| p.accession.as_str())
                        .collect();
                    let mut w = std::io::BufWriter::new(std::fs::File::create(&path)?);
                    writeln!(w, "decoy\tconfident\tmass_z\trt_diff\tenvelope\tcoelution\tn_iso\tscore")?;
                    for h in &hits {
                        let c = !h.decoy
                            && db.peptides[h.peptide_idx as usize].proteins.iter().any(|a| conf.contains(a.as_ref()));
                        writeln!(w, "{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{}\t{:.4}",
                            h.decoy as u8, c as u8, h.abs_mass_z, h.abs_rt_diff,
                            h.envelope, h.coelution, h.n_isotopes, h.score)?;
                    }
                    log::info!("dumped {} hill-hit feature rows to {}", hits.len(), path);
                }

                // Hill evidence as PeptideMatches, EXACTLY like the feature pass:
                // real RT residual + mass z + isotope count. No peptide-level FDR —
                // the per-protein binomial + razor parsimony (infer_proteins) is what
                // separates target from decoy, by aggregating many peptides.
                let hill_pms: Vec<PeptideMatch> = hits
                    .iter()
                    .map(|h| PeptideMatch {
                        peptide_idx: h.peptide_idx as usize,
                        decoy: h.decoy,
                        abs_rt_diff: h.abs_rt_diff,
                        abs_im_diff: f64::INFINITY,
                        abs_mass_z: h.abs_mass_z,
                        score: h.score,
                        feat_quality: h.n_isotopes,
                    })
                    .collect();

                let is_yeast_prot = |a: &str| a.contains("YEAST");
                let has_yeast = db.peptides.iter().any(|p| p.proteins.iter().any(|a| is_yeast_prot(a)));
                let report = |label: &str, prots: &[uno_core::ProteinId]| {
                    let t = prots.iter().filter(|p| !p.decoy && p.q_value <= 0.01).count();
                    let y = prots
                        .iter()
                        .filter(|p| !p.decoy && p.q_value <= 0.01 && is_yeast_prot(&p.accession))
                        .count();
                    if has_yeast {
                        log::info!(
                            "{}: {} target proteins @1% ({} yeast entrapment -> {:.2}% FDP)",
                            label, t, y, 100.0 * y as f64 / t.max(1) as f64
                        );
                    } else {
                        log::info!("{}: {} target proteins @1% FDR", label, t);
                    }
                };

                let rt_sig = result.rt_sigma.max(0.02);

                // Baselines through infer_proteins (the weak rollup).
                let feat_only = infer_proteins(&db, &[], &result.peptide_matches, &config.protein);
                report("feature-only (infer_proteins)", &feat_only);
                let hill_infer = infer_proteins(&db, &[], &hill_pms, &config.protein);
                report("hill-only (infer_proteins)", &hill_infer);

                // Hill-only through the MASS-ANCHOR scorer (the strong one, 626 on
                // features). Deviation = mass_z² + rt_z² ONLY. Envelope + coelution
                // do NOT go in the protein deviation: a decoy that coincides with a
                // real peptide's cluster inherits that cluster's good env+coel, so
                // they don't separate target from decoy — folding them in cost 112
                // proteins and pushed entrapment FDP over 1%. They earn their keep at
                // the per-peptide selection/gate stage (rejecting non-peptide junk),
                // which already happened in hill_search. Opt-in override kept for
                // experiments via HILL_ENV_TOL / HILL_COEL_TOL.
                let env_tol = std::env::var("HILL_ENV_TOL").ok().and_then(|s| s.parse::<f64>().ok());
                let coel_tol = std::env::var("HILL_COEL_TOL").ok().and_then(|s| s.parse::<f64>().ok());
                let hill_dev: Vec<(usize, bool, f64, usize)> = hits
                    .iter()
                    .map(|h| {
                        let rt_z = h.abs_rt_diff / rt_sig;
                        let mut dev = h.abs_mass_z * h.abs_mass_z + rt_z * rt_z;
                        if let Some(t) = env_tol {
                            dev += ((1.0 - h.envelope) / t.max(1e-3)).powi(2);
                        }
                        if let Some(t) = coel_tol {
                            dev += ((1.0 - h.coelution) / t.max(1e-3)).powi(2);
                        }
                        (h.peptide_idx as usize, h.decoy, dev, h.mono_hill)
                    })
                    .collect();
                let hill_ma = uno_core::score_by_deviation(&db, &hill_dev, 0.01);
                report("hill-only (mass-anchor)", &hill_ma);

                // PERCOLATOR-style rescore: semi-supervised CV-LDA over the full hill
                // feature set (mass_z, rt_diff, envelope, coelution, n_iso), positives
                // = mass-anchor-confident peptides, negatives = decoys, folds by
                // peptide. Diagnostic showed confident-vs-decoy AUC up to 0.75, so
                // there's a multivariate boundary to learn. The learned discriminant
                // becomes the deviation into the same rollup.
                let conf: std::collections::HashSet<&str> = proteins
                    .iter()
                    .filter(|p| !p.decoy && p.q_value <= 0.01)
                    .map(|p| p.accession.as_str())
                    .collect();
                let feats: Vec<Vec<f32>> = hits
                    .iter()
                    .map(|h| vec![
                        h.abs_mass_z as f32,
                        h.abs_rt_diff as f32,
                        h.envelope as f32,
                        h.coelution as f32,
                        h.n_isotopes as f32,
                    ])
                    .collect();
                let labels: Vec<Option<f64>> = hits
                    .iter()
                    .map(|h| {
                        if h.decoy {
                            Some(0.0)
                        } else if db.peptides[h.peptide_idx as usize].proteins.iter().any(|a| conf.contains(a.as_ref())) {
                            Some(1.0)
                        } else {
                            None
                        }
                    })
                    .collect();
                let pep_idx: Vec<usize> = hits.iter().map(|h| h.peptide_idx as usize).collect();
                let disc = uno_core::cv_lda_rescore(&feats, &pep_idx, &labels, 3);
                // Map discriminant (higher = better) to a deviation (lower = better)
                // landing in the rollup's 1..25 cut range: standardize, then
                // dev = 10 - 3·z, clamped ≥ 0.
                let n = disc.len().max(1) as f64;
                let mean = disc.iter().sum::<f64>() / n;
                let sd = (disc.iter().map(|d| (d - mean).powi(2)).sum::<f64>() / n).sqrt().max(1e-9);
                let perc_dev: Vec<(usize, bool, f64, usize)> = hits
                    .iter()
                    .zip(&disc)
                    .map(|(h, &d)| {
                        let z = (d - mean) / sd;
                        ((h.peptide_idx as usize), h.decoy, (10.0 - 3.0 * z).max(0.0), h.mono_hill)
                    })
                    .collect();
                let hill_perc = uno_core::score_by_deviation(&db, &perc_dev, 0.01);
                report("hill-only (percolator LDA)", &hill_perc);
            }
        }
    }

    // Phase 7: quant rollup onto koth's per-sample intensities.
    let write_quant = |path: &PathBuf, rows: &[QuantRow]| -> Result<()> {
        use std::io::Write;
        let mut w = std::io::BufWriter::new(std::fs::File::create(path)?);
        write!(w, "key\tq_value")?;
        for name in &matrix.run_names {
            write!(w, "\t{}", name)?;
        }
        writeln!(w)?;
        for r in rows {
            write!(w, "{}\t{:.5}", r.key, r.q_value)?;
            for x in &r.intensities {
                write!(w, "\t{:.6e}", x)?;
            }
            writeln!(w)?;
        }
        Ok(())
    };

    if let Some(path) = &cli.out_peptide_quant {
        let rows = peptide_quant(&matrix, &db, ids, cli.quant_fdr);
        write_quant(path, &rows)?;
        log::info!("wrote {} peptides to {} (q<={})", rows.len(), path.display(), cli.quant_fdr);
    }
    if let Some(path) = &cli.out_protein_quant {
        let protein_q: HashMap<String, f64> = proteins
            .iter()
            .filter(|p| !p.decoy)
            .map(|p| (p.accession.clone(), p.q_value))
            .collect();
        let rows = protein_quant(&matrix, &db, ids, &protein_q, cli.quant_fdr);
        write_quant(path, &rows)?;
        log::info!("wrote {} proteins to {} (q<={})", rows.len(), path.display(), cli.quant_fdr);
    }

    Ok(())
}
