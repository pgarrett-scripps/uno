//! Built-in feature finding: run koth in-process on a raw file and hand uno a
//! `FeatureMatrix` directly, with no intermediate `*.features.tsv`.
//!
//! This is the same three-stage koth pipeline the standalone `koth_ff` binary
//! runs (hills -> features -> scoring), but we keep the scored features in
//! memory and convert them straight into uno's single-run `Feature` rows instead
//! of serialising to TSV and re-parsing. The fields match koth's biosaur2-style
//! TSV schema, but the TSV writer sorts by summed intensity and rounds numeric
//! values. In-process rows retain pipeline order and full precision. Equal
//! feature membership therefore does not guarantee identical search outputs.
//!
//! Input routing (mzML / Bruker `.d` / Thermo `.raw`) is handled inside
//! `koth_ms::run_pipeline` by file extension. `.raw` needs uno's `thermo`
//! feature (on by default; koth-ms's pure-Rust reader, no .NET runtime); `.d`
//! uses koth-ms's `tdf` feature (always enabled here); mzML always works.
//!
//! The default koth config is the biosaur2-matching profile (`min_scans = 2` +
//! `min_scan_overlap = 1`), bundled from `koth_default.toml`. Those two knobs
//! together bring koth's single-run feature depth to parity with biosaur2 for
//! uno's MS1 search; each alone underperforms. Override by adding a `[koth]`
//! table to the uno config TOML.

use anyhow::{Context, Result};
use std::path::Path;
use std::time::Instant;

use koth_ms::config::KothConfig;
use uno_core::{Feature, FeatureMatrix};

/// The biosaur2-matching koth config used when the uno config has no `[koth]`
/// table. See module docs for why these particular values.
pub const DEFAULT_KOTH_CONFIG: &str = include_str!("koth_default.toml");

/// The Bruker `.d` profile. The biosaur2-matching default above is tuned for
/// Orbitrap and is wrong on timsTOF (singly-charged background, ~3 ppm mass
/// accuracy); see the file header for the specific deltas and their evidence.
pub const BRUKER_KOTH_CONFIG: &str = include_str!("koth_bruker.toml");

/// Parse the bundled default koth config.
pub fn default_config() -> Result<KothConfig> {
    let cfg: KothConfig =
        toml::from_str(DEFAULT_KOTH_CONFIG).context("parsing bundled default koth config")?;
    cfg.features.validate()?;
    Ok(cfg)
}

/// Parse the bundled Bruker koth config.
pub fn bruker_config() -> Result<KothConfig> {
    let cfg: KothConfig =
        toml::from_str(BRUKER_KOTH_CONFIG).context("parsing bundled Bruker koth config")?;
    cfg.features.validate()?;
    Ok(cfg)
}

/// Is this raw input a Bruker `.d` directory? koth routes by extension, so match
/// on the same thing rather than probing the filesystem.
pub fn is_bruker(path: &Path) -> bool {
    path.extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("d"))
}

/// The bundled config for `path`: the Bruker profile for `.d`, the
/// biosaur2-matching default for mzML and Thermo `.raw`.
pub fn default_config_for(path: &Path) -> Result<KothConfig> {
    if is_bruker(path) {
        bruker_config()
    } else {
        default_config()
    }
}

/// Run koth's feature finder in-process on `path` (mzML / `.d` / `.raw`) and
/// return the single-run `FeatureMatrix` uno's search consumes.
pub fn features_from_raw(path: &Path, config: &KothConfig) -> Result<FeatureMatrix> {
    // uno generates its own decoys downstream; koth's spectrum-shuffle decoy
    // mode must stay off here or we'd be searching shuffled spectra.
    let mut cfg = config.clone();
    cfg.file.decoy_mode = false;

    // One call into koth's library seam: `run_pipeline` runs the exact staged
    // functions the `koth_ff` binary uses (hills -> features -> scoring).
    // Serialization order and precision still differ from external TSV input.
    // Scoring on: uno's search reads the cosine/isotope/combined
    // scores. Going through this single entry point keeps uno off koth's internal
    // stage functions, which change shape more often than the pipeline seam.
    let t = Instant::now();
    // `emit_ms2: false` keeps the MS1 feature path byte-for-byte unchanged (koth
    // documents it as exactly reproducing the historic MS1-only `run_pipeline`);
    // uno is an MS1-only engine and never reads MS2 hills. The field was added
    // upstream after `7f2cece` and broke this build until it was named here.
    let out = koth_ms::run_pipeline(
        path,
        &cfg,
        &koth_ms::PipelineOptions {
            scoring: true,
            emit_ms2: false,
        },
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
            // koth 0.4.0 made ion polarity explicit; take it from the run config so
            // `--raw` matches koth_ff's own TSV writer rather than assuming positive.
            let mass = sf.monoisotopic_neutral_mass(cfg.file.polarity)?;
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
        matrix_evidence: None,
        joint_evidence: None,
        joint_cell_evidence: None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn koth070_profiles_and_lfq_config_are_valid() {
        let orbitrap = default_config().unwrap();
        let bruker = bruker_config().unwrap();
        for cfg in [&orbitrap, &bruker] {
            assert_eq!(cfg.features.isotope_evidence_ratio_sigma, 0.75);
            assert_eq!(cfg.features.isotope_evidence_cosine_shape, 2.0);
        }
        let standalone: KothConfig =
            toml::from_str(include_str!("../../../tests/fixtures/koth_uno.toml")).unwrap();
        assert_eq!(
            toml::to_string(&standalone).unwrap(),
            toml::to_string(&orbitrap).unwrap()
        );
        let align: koth_ms::config::AlignConfig =
            toml::from_str(include_str!("../../../tests/fixtures/koth_align.toml")).unwrap();
        align.lfq.consensus.validate().unwrap();
        assert_eq!(align.lfq.consensus.max_group_qvalue, 0.05);
        assert_eq!(align.lfq.rt_window_pct, 0.005);
        assert_eq!(align.lfq.im_tolerance, 0.015);
        assert!(align.output.export_long);
    }

    #[test]
    fn legacy_lfq_admission_flags_are_rejected() {
        for setting in [
            "min_group_size = 1",
            "allow_replicated_weak_seeds = true",
            "min_seed_combined_score = 0.75",
            "min_member_combined_score = 0.5",
        ] {
            let text = format!("[lfq.consensus]\n{setting}\n");
            assert!(toml::from_str::<koth_ms::config::AlignConfig>(&text).is_err());
        }
    }
}
