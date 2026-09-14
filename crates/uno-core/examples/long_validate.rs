//! Validate an exported bundle against its wide export and original source runs.
use anyhow::{ensure, Result};
use uno_core::{long_input::Bundle, FeatureMatrix};
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    ensure!(
        args.len() >= 2,
        "usage: long_validate BUNDLE WIDE [SOURCE_FEATURES ...]"
    );
    let bundle = Bundle::open(std::path::Path::new(&args[0]))?;
    let legacy = bundle.load_consensus(2, false, true)?;
    let wide = FeatureMatrix::load_matrix(&args[1], 2)?;
    ensure!(
        legacy.run_names == wide.run_names && legacy.features == wide.features,
        "long/wide consumed feature mismatch"
    );
    let restored = bundle.load_consensus(2, true, false)?;
    ensure!(
        restored.features.iter().all(|f| f.n_isotopes.is_some()),
        "missing original seed isotope count"
    );
    let mut run_counts = Vec::new();
    for (run, path) in args.iter().skip(2).enumerate() {
        let original = FeatureMatrix::load(path)?;
        let roundtrip = bundle.load_run(run)?;
        ensure!(
            original.features == roundtrip.features,
            "original feature round-trip mismatch in run {}",
            run
        );
        run_counts.push(original.len());
    }
    println!(
        "{}",
        serde_json::json!({"long_wide_exact":true,"consensus_features":legacy.len(),"original_quality_restored":true,"original_run_feature_counts":run_counts})
    );
    Ok(())
}
