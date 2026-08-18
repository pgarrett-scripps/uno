//! Built-in feature finding: run koth in-process on a raw file and hand uno a
//! `FeatureMatrix` directly, with no intermediate `*.features.tsv`.
//!
//! This is the same three-stage koth pipeline the standalone `koth_ff` binary
//! runs (hills -> features -> scoring), but we keep the scored features in
//! memory and convert them straight into uno's single-run `Feature` rows instead
//! of serialising to TSV and re-parsing. The field mapping mirrors koth's
//! biosaur2-compatible TSV writer exactly, so `--raw` and a `--features` file
//! produced by `koth_ff` on the same run yield byte-identical feature tables.
//!
//! Input routing (mzML / Bruker `.d` / Thermo `.raw`) is handled inside
//! `koth_ff::run_hills_streaming` by file extension. `.raw` needs uno built with
//! `--features thermo` (a .NET 8 runtime); `.d` needs koth's default `tdf`
//! feature (on by default); mzML always works.
//!
//! The default koth config is the biosaur2-matching profile (`min_scans = 2` +
//! `min_scan_overlap = 1`), bundled from `koth_default.toml`. Those two knobs
//! together bring koth's single-run feature depth to parity with biosaur2 for
//! uno's MS1 search; each alone underperforms. Override by adding a `[koth]`
//! table to the uno config TOML.

use anyhow::{Context, Result};
use std::path::Path;
use std::time::Instant;

use koth_ff::config::KothConfig;
use uno_core::{Feature, FeatureMatrix};

/// The biosaur2-matching koth config used when the uno config has no `[koth]`
/// table. See module docs for why these particular values.
pub const DEFAULT_KOTH_CONFIG: &str = include_str!("koth_default.toml");

/// Parse the bundled default koth config.
pub fn default_config() -> Result<KothConfig> {
    toml::from_str(DEFAULT_KOTH_CONFIG).context("parsing bundled default koth config")
}

/// Run koth's feature finder in-process on `path` (mzML / `.d` / `.raw`) and
/// return the single-run `FeatureMatrix` uno's search consumes.
pub fn features_from_raw(path: &Path, config: &KothConfig) -> Result<FeatureMatrix> {
    // uno generates its own decoys downstream; koth's spectrum-shuffle decoy
    // mode must stay off here or we'd be searching shuffled spectra.
    let mut cfg = config.clone();
    cfg.file.decoy_mode = false;

    // One call into koth's library seam: `run_pipeline` runs the exact staged
    // functions the `koth_ff` binary uses (hills -> features -> scoring), so
    // `--raw` stays byte-identical to a `--features` file that `koth_ff` produced
    // on the same run. Scoring on: uno's search reads the cosine/isotope/combined
    // scores. Going through this single entry point keeps uno off koth's internal
    // stage functions, which change shape more often than the pipeline seam.
    let t = Instant::now();
    // `emit_ms2: false` keeps the MS1 feature path byte-for-byte unchanged (koth
    // documents it as exactly reproducing the historic MS1-only `run_pipeline`);
    // uno is an MS1-only engine and never reads MS2 hills. The field was added
    // upstream after `7f2cece` and broke this build until it was named here.
    let out = koth_ff::run_pipeline(
        path,
        &cfg,
        &koth_ff::PipelineOptions { scoring: true, emit_ms2: false },
    )
        .with_context(|| format!("koth feature finding on {}", path.display()))?;
    log::info!(
        "koth: {} hills, {} scored features ({:.1}s)",
        out.hills.len(),
        out.features.len(),
        t.elapsed().as_secs_f64()
    );
    let scored = out.features;

    let run_name = path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("run")
        .to_string();

    let features = scored
        .iter()
        .filter_map(|sf| {
            let f = &sf.feature;
            // Charge-0 / mass-less features carry no usable precursor mass; koth's
            // TSV writer leaves massCalib blank for them, so drop them here too.
            let mass = sf.monoisotopic_neutral_mass()?;
            let im = f.im_apex();
            Some(Feature {
                mass,
                mz: sf.monoisotopic_mz(),
                charge: f.charge,
                rt: f.rt_apex(),
                im: if im != 0.0 { Some(im) } else { None },
                combined_score: sf.combined_score,
                // Single-run input: always 1 (kept for ML-featurizer stability).
                n_contributing_runs: 1,
                n_isotopes: Some(f.hills.len() as f64),
                n_scans: Some(f.n_scans_total() as f64),
                // biosaur2-only shape columns koth does not produce.
                isoerror: None,
                isoerror2: None,
                rt_width: None,
                cosine_score: Some(sf.cosine_score),
                isotope_score: Some(sf.isotope_score),
                ppm_error: Some(f.ppm_error),
                neutron_offset: Some(sf.neutron_offset as f64),
                // uno's TSV loader reads the `intensityApex` column (apex-scan
                // total), falling back to `intensitySum` only if apex is absent.
                // Match that here — `total_intensity()` (the sum) is the wrong
                // quantity and shifts the intensity-derived ML features.
                intensities: vec![f.total_intensity_at_apex()],
            })
        })
        .collect::<Vec<_>>();

    Ok(FeatureMatrix {
        run_names: vec![run_name],
        features,
    })
}
