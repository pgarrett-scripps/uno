//! The *effective* configuration report.
//!
//! A search is only reproducible if the reader can see the settings the engine
//! actually used, which is not the same as the settings the user passed. Three
//! things make the two differ:
//!
//!   1. **Tuned shipped defaults.** Several values were chosen by maximizing
//!      entrapment depth on development data (run10, deephela) and are not
//!      user-facing config fields at all — the classifier's `min_child_weight`,
//!      the included band of the protein-score ensemble, the deviation cut
//!      ceiling. A reader cannot verify or vary what the engine never prints.
//!   2. **Overridden fields.** `scoring.max_depth` and `scoring.shrinkage` are
//!      clamped before the classifier sees them, so a user who sets them gets
//!      something else.
//!   3. **Inactive fields.** Some fields are still accepted by the config format
//!      but are no longer read by any code path. Setting them does nothing.
//!
//! [`render`] prints all three, each marked, so the block is a complete and
//! honest record of the run. It also returns a warning for every *user-set*
//! field that is overridden or ignored, which is the half of the problem that a
//! silent config cannot solve.
//!
//! This is a report, not a config file: `(unset)` values and `(inactive)` lines
//! are not meant to be fed back in.

use crate::database::DigestConfig;
use crate::matching::MatchConfig;
use crate::protein::ProteinConfig;
use crate::score::ScoreConfig;

/// Why a line reads the way it does.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mark {
    /// A plain default, or exactly what the user asked for.
    Plain,
    /// A shipped default chosen on development data, not a user input.
    Tuned,
    /// Computed at runtime from other settings or from the input.
    Derived,
    /// Read from the environment in this process, not from the config.
    Env,
    /// The engine uses a different value than this one.
    Overridden,
    /// Accepted by the config format but read by no code path.
    Inactive,
}

impl Mark {
    fn tag(self) -> &'static str {
        match self {
            Mark::Plain => "",
            Mark::Tuned => "(tuned)",
            Mark::Derived => "(derived)",
            Mark::Env => "(env)",
            Mark::Overridden => "(overridden)",
            Mark::Inactive => "(inactive)",
        }
    }
}

/// What the CLI knows and the config structs do not: which input route ran,
/// which models were loaded, and which file the settings came from.
#[derive(Debug, Clone, Default)]
pub struct RunContext {
    pub uno_version: String,
    /// `--features`, `--raw`, `--matrix`, `--matrix-long`, `--matrix-long --matrix-joint`.
    pub input_route: String,
    pub input_path: String,
    pub fasta_path: String,
    /// `--config PATH`, or None when the shipped defaults were used whole.
    pub config_path: Option<String>,
    /// Effective `--matrix-min-runs` (defaulted by route), where it applies.
    pub matrix_min_runs: Option<u32>,
    /// "bundled CNN", a custom `.urtc` path, or an external prediction library.
    pub rt_source: String,
    /// "bundled CNN", a custom path, or "disabled (--no-im)".
    pub im_source: String,
    /// `--seed N` was given, so it replaced `digest.decoy_seed`/`scoring.seed`.
    pub cli_seed: Option<u64>,
    /// The Bruker `.d` route auto-set the feature-quality gate.
    pub bruker_gate_applied: bool,
    /// Fields the user actually set, as `"section.field"`. Used to warn only
    /// about settings that were asked for and then ignored.
    pub user_set: Vec<String>,
}

/// The rendered block plus a warning per user-set field the engine does not honour.
#[derive(Debug, Clone)]
pub struct Report {
    pub text: String,
    pub warnings: Vec<String>,
}

/// Read an environment knob, as the engine reads it.
type EnvFn<'a> = &'a dyn Fn(&str) -> Option<String>;

/// Print a float the way the engine's own literal was written. Several shipped
/// constants are `f32`, and widening one to `f64` before formatting turns 0.05
/// into 0.05000000074505806 — true to the bits, useless in a report. Any value
/// that round-trips through `f32` is formatted as an `f32`.
fn num(v: f64) -> String {
    if (v as f32) as f64 == v {
        format!("{}", v as f32)
    } else {
        format!("{v}")
    }
}

fn env_f64(env: EnvFn, key: &str, default: f64) -> (f64, bool) {
    match env(key).and_then(|s| s.parse().ok()) {
        Some(v) => (v, true),
        None => (default, false),
    }
}

fn env_usize(env: EnvFn, key: &str, default: usize) -> (usize, bool) {
    match env(key).and_then(|s| s.parse().ok()) {
        Some(v) => (v, true),
        None => (default, false),
    }
}

/// Lines accumulate here so every section renders the same way: a `key = value`
/// padded out to a fixed column, then the mark, then the note.
#[derive(Default)]
struct Lines(Vec<String>);

impl Lines {
    fn section(&mut self, name: &str) {
        self.0.push(String::new());
        self.0.push(format!("[{name}]"));
    }

    fn put(&mut self, key: &str, value: impl std::fmt::Display, mark: Mark, note: &str) {
        let assign = format!("{key} = {value}");
        if mark == Mark::Plain && note.is_empty() {
            self.0.push(assign);
            return;
        }
        let tag = mark.tag();
        let comment = if note.is_empty() {
            tag.to_string()
        } else if tag.is_empty() {
            note.to_string()
        } else {
            format!("{tag} {note}")
        };
        self.0.push(format!("{assign:<44} # {comment}"));
    }

    fn plain(&mut self, key: &str, value: impl std::fmt::Display) {
        self.put(key, value, Mark::Plain, "");
    }
}

fn quote(s: &str) -> String {
    format!("\"{}\"", s.replace('\\', "\\\\").replace('"', "\\\""))
}

fn opt<T: std::fmt::Display>(v: &Option<T>) -> String {
    match v {
        Some(v) => v.to_string(),
        None => "(unset)".to_string(),
    }
}

fn i8_list(v: &[i8]) -> String {
    let items: Vec<String> = v.iter().map(|x| x.to_string()).collect();
    format!("[{}]", items.join(", "))
}

/// Inline table with keys sorted, so the report is byte-identical run to run
/// (a `HashMap`'s own order is not).
fn mod_table(m: &std::collections::HashMap<String, f32>) -> String {
    let mut keys: Vec<&String> = m.keys().collect();
    keys.sort();
    let items: Vec<String> = keys.iter().map(|k| format!("{} = {}", quote(k), m[*k])).collect();
    if items.is_empty() { "{}".to_string() } else { format!("{{ {} }}", items.join(", ")) }
}

fn var_mod_table(m: &std::collections::HashMap<String, Vec<f32>>) -> String {
    let mut keys: Vec<&String> = m.keys().collect();
    keys.sort();
    let items: Vec<String> = keys
        .iter()
        .map(|k| {
            let vals: Vec<String> = m[*k].iter().map(|v| v.to_string()).collect();
            format!("{} = [{}]", quote(k), vals.join(", "))
        })
        .collect();
    if items.is_empty() { "{}".to_string() } else { format!("{{ {} }}", items.join(", ")) }
}

/// Render the effective configuration, reading tuned knobs from the process
/// environment exactly as the engine does.
pub fn render(
    digest: &DigestConfig,
    matching: &MatchConfig,
    scoring: &ScoreConfig,
    protein: &ProteinConfig,
    run: &RunContext,
) -> Report {
    render_with(digest, matching, scoring, protein, run, &|k| std::env::var(k).ok())
}

/// [`render`] with the environment injected, so a test can pin the output.
pub fn render_with(
    digest: &DigestConfig,
    matching: &MatchConfig,
    scoring: &ScoreConfig,
    protein: &ProteinConfig,
    run: &RunContext,
    env: EnvFn,
) -> Report {
    let mut l = Lines::default();
    let mut warnings: Vec<String> = Vec::new();
    let was_set = |field: &str| run.user_set.iter().any(|f| f == field);
    let mut warn_if_set = |field: &str, why: &str, warnings: &mut Vec<String>| {
        if was_set(field) {
            warnings.push(format!("config sets `{field}` but {why}"));
        }
    };

    l.0.push("# uno effective configuration".to_string());
    l.0.push("# Every setting this search will actually use. Marks:".to_string());
    l.0.push("#   (tuned)      shipped default chosen on development data, not a user input".to_string());
    l.0.push("#   (derived)    computed at runtime from another setting or from the input".to_string());
    l.0.push("#   (env)        taken from an environment variable in this process".to_string());
    l.0.push("#   (overridden) the engine uses a different value than the one shown".to_string());
    l.0.push("#   (inactive)   accepted by the config format but read by no code path".to_string());
    l.0.push("# Unmarked lines are plain defaults or values you passed. This is a report,".to_string());
    l.0.push("# not a config file: (unset) and (inactive) lines are not meant to be fed back.".to_string());

    // ---- [run] --------------------------------------------------------------
    l.section("run");
    l.plain("uno_version", quote(&run.uno_version));
    l.plain("input_route", quote(&run.input_route));
    l.plain("input", quote(&run.input_path));
    l.plain("fasta", quote(&run.fasta_path));
    match &run.config_path {
        Some(p) => l.plain("config_file", quote(p)),
        None => l.put("config_file", "(unset)", Mark::Derived, "shipped defaults, no config file"),
    }
    if let Some(n) = run.matrix_min_runs {
        l.put("matrix_min_runs", n, Mark::Derived, "1 for --matrix-joint, else 2, unless passed");
    }
    l.plain("rt_source", quote(&run.rt_source));
    l.plain("im_source", quote(&run.im_source));
    match run.cli_seed {
        Some(s) => l.put("cli_seed", s, Mark::Plain, "--seed replaced digest.decoy_seed and scoring.seed"),
        None => l.put("cli_seed", "(unset)", Mark::Plain, "shipped seed streams"),
    }

    // ---- [digest] -----------------------------------------------------------
    l.section("digest");
    l.plain("missed_cleavages", digest.missed_cleavages);
    l.plain("min_len", digest.min_len);
    l.plain("max_len", digest.max_len);
    l.plain("cleave_at", quote(&digest.cleave_at));
    l.plain("restrict", quote(&digest.restrict));
    l.plain("c_terminal", digest.c_terminal);
    l.plain("peptide_min_mass", digest.peptide_min_mass);
    l.plain("peptide_max_mass", digest.peptide_max_mass);
    l.plain("static_mods", mod_table(&digest.static_mods));
    l.plain("variable_mods", var_mod_table(&digest.variable_mods));
    if digest.variable_mods.is_empty() {
        l.put("max_variable_mods", digest.max_variable_mods, Mark::Inactive, "no variable_mods are configured");
    } else {
        l.plain("max_variable_mods", digest.max_variable_mods);
    }
    l.plain("decoy_tag", quote(&digest.decoy_tag));
    l.plain("generate_decoys", digest.generate_decoys);

    let protein_level_decoys = digest.decoy_method != "shuffle";
    if protein_level_decoys {
        l.put("decoy_method", quote(&digest.decoy_method), Mark::Plain,
            "protein-level generator; overrides decoy_shuffle");
        l.put("decoy_shuffle", digest.decoy_shuffle, Mark::Overridden,
            "decoy_method is not \"shuffle\"");
        warn_if_set("digest.decoy_shuffle", "decoy_method is not \"shuffle\", so it is overridden", &mut warnings);
    } else {
        l.plain("decoy_method", quote(&digest.decoy_method));
        l.plain("decoy_shuffle", digest.decoy_shuffle);
    }
    let markov = matches!(digest.decoy_method.as_str(), "markov" | "debruijn");
    if markov {
        l.plain("decoy_order", digest.decoy_order);
    } else {
        l.put("decoy_order", digest.decoy_order, Mark::Inactive,
            "only \"markov\" and \"debruijn\" read it");
        warn_if_set("digest.decoy_order", "decoy_method is neither \"markov\" nor \"debruijn\", so it is ignored", &mut warnings);
    }
    if digest.decoy_shuffle && !protein_level_decoys {
        l.put("decoy_mass_shift", digest.decoy_mass_shift, Mark::Overridden,
            "ignored while decoy_shuffle = true (shuffled decoys keep their own mass)");
        warn_if_set("digest.decoy_mass_shift", "decoy_shuffle = true, so it is ignored", &mut warnings);
    } else {
        l.plain("decoy_mass_shift", digest.decoy_mass_shift);
    }
    l.plain("decoy_seed", digest.decoy_seed);

    // ---- [matching] ---------------------------------------------------------
    l.section("matching");
    l.plain("ppm_tolerance", matching.ppm_tolerance);
    l.plain("isotope_offsets", i8_list(&matching.isotope_offsets));
    l.plain("min_charge", matching.min_charge);
    if run.bruker_gate_applied {
        l.put("min_isotopes", matching.min_isotopes, Mark::Derived,
            "Bruker .d default: both gates were untouched, so the timsTOF gate was applied");
        l.put("min_combined_score", matching.min_combined_score, Mark::Derived,
            "Bruker .d default; set either gate explicitly to take full control");
    } else {
        l.plain("min_isotopes", matching.min_isotopes);
        l.plain("min_combined_score", matching.min_combined_score);
    }

    // ---- [scoring] ----------------------------------------------------------
    l.section("scoring");
    l.put("top_k", scoring.top_k, Mark::Tuned,
        "50 chosen on development-data entrapment depth (25 -> ~half the IDs)");
    l.plain("n_folds", scoring.n_folds);
    l.plain("seed", scoring.seed);
    l.put("rt_seed", opt(&scoring.rt_seed), Mark::Derived,
        &format!("effective = {}", scoring.effective_rt_seed()));
    l.put("cv_seed", opt(&scoring.cv_seed), Mark::Derived,
        &format!("effective = {}", scoring.effective_cv_seed()));
    l.put("classifier_seed", opt(&scoring.classifier_seed), Mark::Derived,
        &format!("effective = {}", scoring.effective_classifier_seed()));
    l.plain("mass_anchor", scoring.mass_anchor);
    l.plain("mass_anchor_ppm", scoring.mass_anchor_ppm);
    l.plain("rt_calibrated", scoring.rt_calibrated);
    l.put("ml_rescore", scoring.ml_rescore, Mark::Tuned,
        "on by default; the cumulative-band protein ensemble below depends on it");
    l.plain("matrix_joint_agreement", scoring.matrix_joint_agreement);
    l.plain("matrix_joint_cell_quality", scoring.matrix_joint_cell_quality);
    l.plain("matrix_joint_cell_counts", scoring.matrix_joint_cell_counts);
    l.plain("matrix_joint_cell_q", scoring.matrix_joint_cell_q);
    l.plain("matrix_joint_alternative_groups", scoring.matrix_joint_alternative_groups);
    l.plain("matrix_recurrence", scoring.matrix_recurrence);
    l.plain("diagnostic_peptides_path", opt(&scoring.diagnostic_peptides_path.as_deref().map(quote)));

    // The classifier hyperparameters the config exposes are not the ones it gets.
    let (ml_depth, depth_from_env) = env_usize(env, "UNO_ML_DEPTH", scoring.max_depth.max(7) as usize);
    let (ml_iters, iters_from_env) = env_usize(env, "UNO_ML_ITERS", 300);
    let (ml_lr, lr_from_env) = env_f64(env, "UNO_ML_LR", scoring.shrinkage.min(0.05) as f64);
    l.put("iterations", scoring.iterations, Mark::Inactive,
        &format!("the classifier always trains {ml_iters} rounds; see [tuned]"));
    warn_if_set("scoring.iterations", &format!("the classifier always trains {ml_iters} rounds"), &mut warnings);
    l.put("max_depth", scoring.max_depth, Mark::Overridden,
        &format!("effective depth = max(max_depth, 7) = {ml_depth}"));
    if was_set("scoring.max_depth") && ml_depth != scoring.max_depth as usize {
        warnings.push(format!(
            "config sets `scoring.max_depth` = {} but the classifier trains at depth {ml_depth}",
            scoring.max_depth
        ));
    }
    l.put("shrinkage", scoring.shrinkage, Mark::Overridden,
        &format!("effective learning rate = min(shrinkage, 0.05) = {}", num(ml_lr)));
    if was_set("scoring.shrinkage") && (ml_lr - scoring.shrinkage as f64).abs() > f64::EPSILON {
        warnings.push(format!(
            "config sets `scoring.shrinkage` = {} but the classifier trains at learning rate {}",
            scoring.shrinkage, num(ml_lr)
        ));
    }

    // Fields the config still accepts and no code path reads. Each one is a
    // reproducibility trap: a user who sets it gets the shipped behaviour.
    for (field, value, why) in [
        ("calibration_fdr", scoring.calibration_fdr.to_string(), "the bootstrap selects its own provisional set"),
        ("bootstrap_rounds", scoring.bootstrap_rounds.to_string(), "the calibration round count is fixed in the engine"),
        ("max_train_rows", scoring.max_train_rows.to_string(), "training rows are capped by the 3x negative balance instead"),
        ("seed_min_runs_frac", scoring.seed_min_runs_frac.to_string(), "the RANSAC seed does not read it"),
        ("discriminant", quote(&scoring.discriminant), "the deviation scorer is unconditional"),
        ("ondata_rt", scoring.ondata_rt.to_string(), "no code path reads it"),
        ("ondata_rt_lambda", scoring.ondata_rt_lambda.to_string(), "no code path reads it"),
        ("ondata_rt_finetune", scoring.ondata_rt_finetune.to_string(), "no code path reads it"),
    ] {
        l.put(field, &value, Mark::Inactive, why);
        warn_if_set(&format!("scoring.{field}"), why, &mut warnings);
    }

    // ---- [protein] ----------------------------------------------------------
    l.section("protein");
    l.plain("rt_window", protein.rt_window);
    l.plain("mass_z_window", protein.mass_z_window);
    l.plain("fdr", protein.fdr);
    l.plain("ml_rescore", protein.ml_rescore);
    for (field, value, why) in [
        ("im_window", num(protein.im_window), "no code path reads it"),
        ("optimize_thresholds", protein.optimize_thresholds.to_string(), "the threshold scan is unconditional"),
        ("max_feature_isotopes_gate", protein.max_feature_isotopes_gate.to_string(), "no code path reads it"),
    ] {
        l.put(field, &value, Mark::Inactive, why);
        warn_if_set(&format!("protein.{field}"), why, &mut warnings);
    }

    // ---- [tuned] ------------------------------------------------------------
    // Values with no config field at all. The manuscript names three of these as
    // tuned defaults, so they have to be printable or the claim cannot be checked.
    l.section("tuned");
    l.0.push("# No config field exposes these. Each was chosen by maximizing entrapment".to_string());
    l.0.push("# depth on development data (run10, deephela) and is overridable only by its".to_string());
    l.0.push("# environment variable, shown in the note.".to_string());

    let lightgbm = env("UNO_ML_LIGHTGBM").is_some();
    l.put("ml_backend", quote(if lightgbm { "lightgbm" } else { "sequoia-boost" }),
        if lightgbm { Mark::Env } else { Mark::Tuned }, "UNO_ML_LIGHTGBM selects the C LightGBM trainer");
    l.put("ml_iterations", ml_iters, if iters_from_env { Mark::Env } else { Mark::Tuned }, "UNO_ML_ITERS");
    l.put("ml_max_depth", ml_depth, if depth_from_env { Mark::Env } else { Mark::Tuned },
        "UNO_ML_DEPTH; 7 was the fixed-path best on development data");
    l.put("ml_learning_rate", num(ml_lr), if lr_from_env { Mark::Env } else { Mark::Tuned }, "UNO_ML_LR");
    let (trials, trials_env) = env_usize(env, "UNO_ML_TUNE", 0);
    l.put("ml_hyperparameter_trials", trials, if trials_env { Mark::Env } else { Mark::Tuned },
        "UNO_ML_TUNE; 0 = the fixed-parameter path");
    let rank = env("UNO_ML_RANK").map(|s| s != "0").unwrap_or(false);
    l.put("ml_rank_normalize", rank, if env("UNO_ML_RANK").is_some() { Mark::Env } else { Mark::Tuned }, "UNO_ML_RANK");
    if lightgbm {
        let leaves = env_usize(env, "UNO_ML_LEAVES", (1usize << ml_depth.min(8)).saturating_sub(1).max(2));
        l.put("lightgbm_num_leaves", leaves.0, if leaves.1 { Mark::Env } else { Mark::Derived },
            "UNO_ML_LEAVES; derived as 2^min(depth, 8) - 1");
    } else {
        let mcw = env_f64(env, "UNO_SEQ_MCW", 80.0);
        l.put("sequoia_min_child_weight", num(mcw.0), if mcw.1 { Mark::Env } else { Mark::Tuned },
            "UNO_SEQ_MCW; the manuscript's tuned min_child_weight");
        for (key, var, default) in [
            ("sequoia_lambda", "UNO_SEQ_LAMBDA", 1.0),
            ("sequoia_alpha", "UNO_SEQ_ALPHA", 0.0),
            ("sequoia_gamma", "UNO_SEQ_GAMMA", 0.0),
            ("sequoia_subsample", "UNO_SEQ_SUB", 0.8),
            ("sequoia_colsample_bytree", "UNO_SEQ_COL", 0.8),
            ("sequoia_max_bin", "UNO_SEQ_BINS", 256.0),
        ] {
            let (v, from_env) = env_f64(env, var, default);
            l.put(key, num(v), if from_env { Mark::Env } else { Mark::Tuned }, var);
        }
    }

    let (ens_steps, steps_env) = env_usize(env, "UNO_ENS_STEPS", 40);
    let (ens_lo_raw, lo_env) = env_usize(env, "UNO_ENS_LO", 1);
    let (ens_hi_raw, hi_env) = env_usize(env, "UNO_ENS_HI", 14);
    let ens_lo = ens_lo_raw.max(1);
    let ens_hi = ens_hi_raw.min(ens_steps);
    l.put("protein_ensemble_steps", ens_steps, if steps_env { Mark::Env } else { Mark::Tuned }, "UNO_ENS_STEPS");
    l.put("protein_ensemble_band_lo", ens_lo, if lo_env { Mark::Env } else { Mark::Tuned }, "UNO_ENS_LO");
    l.put("protein_ensemble_band_hi", ens_hi, if hi_env { Mark::Env } else { Mark::Tuned }, "UNO_ENS_HI");
    l.put("protein_ensemble_bands_used", (ens_hi + 1).saturating_sub(ens_lo).max(1), Mark::Derived,
        &format!("the {} tightest of {ens_steps} cumulative inclusion bands", (ens_hi + 1).saturating_sub(ens_lo).max(1)));

    let (cut_max, cut_max_env) = env_f64(env, "UNO_CUT_MAX", 8.0);
    let (cut_step, cut_step_env) = env_f64(env, "UNO_CUT_STEP", 0.5);
    l.put("deviation_cut_max", num(cut_max), if cut_max_env { Mark::Env } else { Mark::Tuned }, "UNO_CUT_MAX");
    l.put("deviation_cut_step", num(cut_step), if cut_step_env { Mark::Env } else { Mark::Tuned },
        "UNO_CUT_STEP; the fine step over 1..12");

    // Any other UNO_* knob set in this process. Everything above is already
    // reported with its effective value; these are the diagnostic switches, and a
    // run with one set is not a default run.
    let reported = [
        "UNO_ML_LIGHTGBM", "UNO_ML_ITERS", "UNO_ML_DEPTH", "UNO_ML_LR", "UNO_ML_TUNE",
        "UNO_ML_RANK", "UNO_ML_LEAVES", "UNO_SEQ_MCW", "UNO_SEQ_LAMBDA", "UNO_SEQ_ALPHA",
        "UNO_SEQ_GAMMA", "UNO_SEQ_SUB", "UNO_SEQ_COL", "UNO_SEQ_BINS", "UNO_ENS_STEPS",
        "UNO_ENS_LO", "UNO_ENS_HI", "UNO_CUT_MAX", "UNO_CUT_STEP",
    ];
    let mut others: Vec<String> = Vec::new();
    for key in OTHER_ENV_KNOBS {
        if reported.contains(key) {
            continue;
        }
        if let Some(v) = env(key) {
            others.push(format!("{key}={v}"));
        }
    }
    if !others.is_empty() {
        l.section("environment_overrides");
        l.0.push("# Diagnostic switches set in this process. This is not a default run.".to_string());
        for (i, o) in others.iter().enumerate() {
            l.put(&format!("override_{}", i + 1), quote(o), Mark::Env, "");
            warnings.push(format!("environment override in effect: {o}"));
        }
    }

    if !warnings.is_empty() {
        l.section("warnings");
        for (i, w) in warnings.iter().enumerate() {
            l.put(&format!("warning_{}", i + 1), quote(w), Mark::Plain, "");
        }
    }

    l.0.push("# end of effective configuration".to_string());
    Report { text: l.0.join("\n"), warnings }
}

/// Every remaining `UNO_*` knob the engine reads. Listed so a run with a
/// diagnostic switch set cannot be mistaken for a default run; keep in step with
/// `std::env::var` calls in this crate (`env_knobs_are_all_listed` checks it).
pub const OTHER_ENV_KNOBS: &[&str] = &[
    "UNO_CONTINUOUS",
    "UNO_CS_GATE",
    "UNO_CS_GATE_ML",
    "UNO_DUMP_ML",
    "UNO_DUMP_PEPSET",
    "UNO_DYNAMIC_P",
    "UNO_ENSEMBLE",
    "UNO_ENS_BY",
    "UNO_IM_CACHE",
    "UNO_LOG_HITRATE",
    "UNO_LOG_PEPCOUNT",
    "UNO_ML_CORR",
    "UNO_ML_DECOYCUT",
    "UNO_MULTI_FEAT",
    "UNO_NO_FEATURE_COMPETE",
    "UNO_PROT_PERMUTE",
    "UNO_RT_CACHE",
    "UNO_RT_CACHE_DIR",
    "UNO_RT_NOISYGAUS",
    "UNO_RT_SIG",
    "UNO_RT_SIG_GCORE",
    "UNO_RT_STAGES",
    "UNO_SMOOTH",
    "UNO_STAGE_MASS",
];

#[cfg(test)]
mod tests {
    use super::*;

    fn ctx() -> RunContext {
        RunContext {
            uno_version: "0.0.0-test".to_string(),
            input_route: "--features".to_string(),
            input_path: "run.features.tsv".to_string(),
            fasta_path: "db.fasta".to_string(),
            config_path: None,
            matrix_min_runs: None,
            rt_source: "bundled CNN".to_string(),
            im_source: "bundled CNN".to_string(),
            cli_seed: None,
            bruker_gate_applied: false,
            user_set: Vec::new(),
        }
    }

    /// Render with a chosen environment, so the test never depends on the
    /// developer's shell (and parallel tests cannot perturb each other).
    fn render_env(run: &RunContext, pairs: &[(&str, &str)]) -> Report {
        let owned: Vec<(String, String)> =
            pairs.iter().map(|(k, v)| (k.to_string(), v.to_string())).collect();
        let lookup = move |k: &str| {
            owned.iter().find(|(key, _)| key == k).map(|(_, v)| v.clone())
        };
        render_with(
            &DigestConfig::default(),
            &MatchConfig::default(),
            &ScoreConfig::default(),
            &ProteinConfig::default(),
            run,
            &lookup,
        )
    }

    /// The point of the report: the tuned values the manuscript names, and the
    /// effective classifier hyperparameters, appear with their real values. This
    /// section is pinned whole — a change here changes what the paper can claim.
    #[test]
    fn tuned_section_is_pinned() {
        let text = render_env(&ctx(), &[]).text;
        let tuned = text
            .split("\n[tuned]\n")
            .nth(1)
            .expect("[tuned] section")
            .split("\n# end of effective configuration")
            .next()
            .unwrap()
            .to_string();
        let expected = "\
# No config field exposes these. Each was chosen by maximizing entrapment
# depth on development data (run10, deephela) and is overridable only by its
# environment variable, shown in the note.
ml_backend = \"sequoia-boost\"                 # (tuned) UNO_ML_LIGHTGBM selects the C LightGBM trainer
ml_iterations = 300                          # (tuned) UNO_ML_ITERS
ml_max_depth = 7                             # (tuned) UNO_ML_DEPTH; 7 was the fixed-path best on development data
ml_learning_rate = 0.05                      # (tuned) UNO_ML_LR
ml_hyperparameter_trials = 0                 # (tuned) UNO_ML_TUNE; 0 = the fixed-parameter path
ml_rank_normalize = false                    # (tuned) UNO_ML_RANK
sequoia_min_child_weight = 80                # (tuned) UNO_SEQ_MCW; the manuscript's tuned min_child_weight
sequoia_lambda = 1                           # (tuned) UNO_SEQ_LAMBDA
sequoia_alpha = 0                            # (tuned) UNO_SEQ_ALPHA
sequoia_gamma = 0                            # (tuned) UNO_SEQ_GAMMA
sequoia_subsample = 0.8                      # (tuned) UNO_SEQ_SUB
sequoia_colsample_bytree = 0.8               # (tuned) UNO_SEQ_COL
sequoia_max_bin = 256                        # (tuned) UNO_SEQ_BINS
protein_ensemble_steps = 40                  # (tuned) UNO_ENS_STEPS
protein_ensemble_band_lo = 1                 # (tuned) UNO_ENS_LO
protein_ensemble_band_hi = 14                # (tuned) UNO_ENS_HI
protein_ensemble_bands_used = 14             # (derived) the 14 tightest of 40 cumulative inclusion bands
deviation_cut_max = 8                        # (tuned) UNO_CUT_MAX
deviation_cut_step = 0.5                     # (tuned) UNO_CUT_STEP; the fine step over 1..12";
        assert_eq!(tuned.trim_end(), expected);
    }

    /// Every field of every config struct reaches the report. Serde gives the
    /// field list, so adding a field to a config and forgetting to print it
    /// fails here instead of silently shipping an incomplete record.
    #[test]
    fn every_config_field_is_printed() {
        let text = render_env(&ctx(), &[]).text;
        let sections: Vec<(&str, serde_json::Value)> = vec![
            ("digest", serde_json::to_value(DigestConfig::default()).unwrap()),
            ("matching", serde_json::to_value(MatchConfig::default()).unwrap()),
            ("scoring", serde_json::to_value(ScoreConfig::default()).unwrap()),
            ("protein", serde_json::to_value(ProteinConfig::default()).unwrap()),
        ];
        for (section, value) in sections {
            for field in value.as_object().expect("config serializes to a table").keys() {
                assert!(
                    text.contains(&format!("\n{field} = ")),
                    "[{section}] field `{field}` is missing from the effective-config report"
                );
            }
        }
    }

    /// A default run reports no warnings, so a warning always means something.
    #[test]
    fn a_default_run_warns_about_nothing() {
        let report = render_env(&ctx(), &[]);
        assert!(report.warnings.is_empty(), "unexpected warnings: {:?}", report.warnings);
        assert!(!report.text.contains("[warnings]"));
        assert!(!report.text.contains("[environment_overrides]"));
    }

    /// The reviewer's second point: a user who sets an ignored or clamped field
    /// must be told, not silently given something else.
    #[test]
    fn user_set_fields_that_are_ignored_produce_warnings() {
        let mut run = ctx();
        run.user_set = vec![
            "scoring.ondata_rt".to_string(),
            "scoring.iterations".to_string(),
            "digest.decoy_mass_shift".to_string(),
            "protein.optimize_thresholds".to_string(),
        ];
        let report = render_env(&run, &[]);
        assert_eq!(report.warnings.len(), 4, "warnings: {:?}", report.warnings);
        assert!(report.warnings.iter().any(|w| w.contains("`scoring.ondata_rt`")));
        assert!(report.warnings.iter().any(|w| w.contains("always trains 300 rounds")));
        assert!(report.warnings.iter().any(|w| w.contains("`digest.decoy_mass_shift`")));
        assert!(report.text.contains("[warnings]"));
    }

    /// Setting a field that IS honoured must stay silent.
    #[test]
    fn user_set_fields_that_are_honoured_do_not_warn() {
        let mut run = ctx();
        run.user_set = vec!["scoring.top_k".to_string(), "matching.ppm_tolerance".to_string()];
        assert!(render_env(&run, &[]).warnings.is_empty());
    }

    /// Clamped hyperparameters: the report shows the value the classifier gets,
    /// and warns when it is not the value that was asked for.
    #[test]
    fn clamped_classifier_hyperparameters_show_the_effective_value() {
        let text = render_env(&ctx(), &[]).text;
        assert!(text.contains("effective depth = max(max_depth, 7) = 7"), "{text}");
        assert!(text.contains("effective learning rate = min(shrinkage, 0.05) = 0.05"), "{text}");

        let mut run = ctx();
        run.user_set = vec!["scoring.max_depth".to_string(), "scoring.shrinkage".to_string()];
        let report = render_env(&run, &[]);
        assert_eq!(report.warnings.len(), 2, "warnings: {:?}", report.warnings);
        assert!(report.warnings.iter().any(|w| w.contains("trains at depth 7")));
        assert!(report.warnings.iter().any(|w| w.contains("learning rate 0.05")));
    }

    /// An environment override is reported as such, with its effective value, and
    /// a run carrying one is flagged as not a default run.
    #[test]
    fn environment_overrides_are_reported_and_warned() {
        let report = render_env(&ctx(), &[("UNO_ENS_HI", "10"), ("UNO_SEQ_MCW", "20"), ("UNO_RT_SIG", "0.7")]);
        let band_hi = format!("{:<44} # (env) UNO_ENS_HI", "protein_ensemble_band_hi = 10");
        assert!(report.text.contains(&band_hi), "{}", report.text);
        assert!(report.text.contains("sequoia_min_child_weight = 20"), "{}", report.text);
        assert!(report.text.contains("protein_ensemble_bands_used = 10"), "{}", report.text);
        assert!(report.text.contains("[environment_overrides]"), "{}", report.text);
        assert!(report.warnings.iter().any(|w| w.contains("UNO_RT_SIG=0.7")), "{:?}", report.warnings);
    }

    /// `UNO_ML_LIGHTGBM` swaps the whole backend, so the sequoia knobs stop
    /// applying and LightGBM's derived leaf count starts.
    #[test]
    fn the_lightgbm_backend_reports_its_own_knobs() {
        let text = render_env(&ctx(), &[("UNO_ML_LIGHTGBM", "1")]).text;
        assert!(text.contains("ml_backend = \"lightgbm\""), "{text}");
        assert!(text.contains("lightgbm_num_leaves = 127"), "{text}");
        assert!(!text.contains("sequoia_min_child_weight"), "{text}");
    }

    /// Every inactive field is marked, so the listing cannot imply that setting
    /// one would do something.
    #[test]
    fn inactive_fields_are_all_marked() {
        let text = render_env(&ctx(), &[]).text;
        for field in [
            "calibration_fdr", "bootstrap_rounds", "max_train_rows", "seed_min_runs_frac",
            "discriminant", "ondata_rt", "ondata_rt_lambda", "ondata_rt_finetune",
            "optimize_thresholds", "max_feature_isotopes_gate", "iterations", "im_window",
        ] {
            let line = text
                .lines()
                .find(|l| l.starts_with(&format!("{field} = ")))
                .unwrap_or_else(|| panic!("no line for {field}"));
            assert!(line.contains("(inactive)"), "{field} is not marked inactive: {line}");
        }
    }

    /// The report has to be byte-identical across runs or it is useless as a
    /// record; `static_mods` is a `HashMap`, whose own order is not stable.
    #[test]
    fn the_report_is_deterministic() {
        let mut digest = DigestConfig::default();
        digest.static_mods.insert("K".to_string(), 8.014_199);
        digest.static_mods.insert("M".to_string(), 15.994_9);
        digest.variable_mods.insert("M".to_string(), vec![15.994_9]);
        let run = ctx();
        let one = render_with(&digest, &MatchConfig::default(), &ScoreConfig::default(),
            &ProteinConfig::default(), &run, &|_| None);
        let two = render_with(&digest, &MatchConfig::default(), &ScoreConfig::default(),
            &ProteinConfig::default(), &run, &|_| None);
        assert_eq!(one.text, two.text);
        assert!(one.text.contains("static_mods = { \"C\" = 57.02146, \"K\" = 8.014199, \"M\" = 15.9949 }"),
            "{}", one.text);
        assert!(one.text.contains("variable_mods = { \"M\" = [15.9949] }"), "{}", one.text);
        // With variable mods configured, max_variable_mods is live again.
        assert!(one.text.lines().any(|l| l == "max_variable_mods = 2"), "{}", one.text);
    }

    /// The Bruker route sets the feature-quality gate itself. That is a derived
    /// value, not something the user passed, and the report must say so.
    #[test]
    fn the_bruker_gate_is_reported_as_derived() {
        let mut run = ctx();
        run.input_route = "--raw".to_string();
        run.bruker_gate_applied = true;
        let matching = MatchConfig { min_isotopes: 4, min_combined_score: 0.85, ..MatchConfig::default() };
        let text = render_with(&DigestConfig::default(), &matching, &ScoreConfig::default(),
            &ProteinConfig::default(), &run, &|_| None).text;
        let gate = format!("{:<44} # (derived)", "min_isotopes = 4");
        assert!(text.contains(&gate), "{text}");
        assert!(text.contains("min_combined_score = 0.85"), "{text}");
    }

    /// `OTHER_ENV_KNOBS` exists so a diagnostic switch cannot hide. It is a hand
    /// list, so check it against the `std::env::var` calls in this crate.
    #[test]
    fn env_knobs_are_all_listed() {
        let src_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src");
        let reported_here = [
            "UNO_ML_LIGHTGBM", "UNO_ML_ITERS", "UNO_ML_DEPTH", "UNO_ML_LR", "UNO_ML_TUNE",
            "UNO_ML_RANK", "UNO_ML_LEAVES", "UNO_SEQ_MCW", "UNO_SEQ_LAMBDA", "UNO_SEQ_ALPHA",
            "UNO_SEQ_GAMMA", "UNO_SEQ_SUB", "UNO_SEQ_COL", "UNO_SEQ_BINS", "UNO_ENS_STEPS",
            "UNO_ENS_LO", "UNO_ENS_HI", "UNO_CUT_MAX", "UNO_CUT_STEP",
            // Not a search setting: test-data strictness and the cache location.
            "UNO_REQUIRE_TEST_DATA",
        ];
        let mut missing: Vec<String> = Vec::new();
        for entry in std::fs::read_dir(&src_dir).expect("src dir") {
            let path = entry.expect("dir entry").path();
            if path.extension().and_then(|e| e.to_str()) != Some("rs") {
                continue;
            }
            let src = std::fs::read_to_string(&path).expect("read source");
            for (_, rest) in src.match_indices("std::env::var(\"") {
                let name: String = rest
                    .trim_start_matches("std::env::var(\"")
                    .chars()
                    .take_while(|c| *c != '"')
                    .collect();
                if !name.starts_with("UNO_") {
                    continue;
                }
                if OTHER_ENV_KNOBS.contains(&name.as_str()) || reported_here.contains(&name.as_str()) {
                    continue;
                }
                missing.push(format!("{name} (in {})", path.display()));
            }
        }
        missing.sort();
        missing.dedup();
        assert!(
            missing.is_empty(),
            "these UNO_* knobs are read but absent from the effective-config report: {missing:?}"
        );
    }
}
