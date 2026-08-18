//! uno-core: MS1-only protein identification.
//!
//! Consumes a SINGLE-RUN biosaur2 feature file (the same input ms1searchpy uses),
//! assigns peptide/protein identities, and controls FDR. Single-run input is
//! required: it carries the per-feature quality columns that separate target from
//! decoy, which cross-run alignment drops.

pub mod calibrate;
pub mod database;
pub mod decoy;
pub mod group;
pub mod im;
pub mod input;
pub mod mass;
pub mod matching;
pub mod protein;
pub mod quant;
pub mod rt;
pub mod rt_cache;
pub mod rt_cnn;
pub mod rt_mlp;
pub mod score;

pub use calibrate::Calibration;
pub use database::{DigestConfig, PeptideDb, TheoPeptide};
pub use im::{im_features, ImModel, IM_N_FEATURES};
pub use input::{Feature, FeatureMatrix};
pub use matching::{match_features, match_features_topk, MatchConfig, Pfm};
pub use protein::{ProteinConfig, ProteinId};
pub use rt::{rt_features, RtModel, RT_N_FEATURES};
pub use rt_cnn::RtCnn;
pub use rt_mlp::RtMlp;
pub use score::{score, ScoreConfig, ScoreResult};

/// Locate a large test input, portably and without lying about coverage.
///
/// The benchmark inputs these tests need are gitignored (multi-gigabyte), so a
/// fresh clone cannot run them. Three tests used to hard-code
/// `/home/patrick-garrett/...` and silently `return` when it was missing, which
/// reports as PASSED. That is the failure mode that let `uno` sit un-compilable at
/// HEAD behind a green `cargo check`.
///
/// This resolves paths relative to the repository instead of one developer's home
/// directory, prints a clearly marked SKIP when the input is absent, and honours
/// `UNO_REQUIRE_TEST_DATA=1` by failing instead of skipping. Set that in CI and on
/// the pre-benchmark checklist, where the data IS expected to be present.
#[cfg(test)]
pub(crate) fn test_input(rel: &str) -> Option<std::path::PathBuf> {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .expect("crates/uno-core -> repo root");
    let p = root.join(rel);
    if p.exists() {
        return Some(p);
    }
    if std::env::var("UNO_REQUIRE_TEST_DATA").as_deref() == Ok("1") {
        panic!("UNO_REQUIRE_TEST_DATA=1 but test input is missing: {}", p.display());
    }
    eprintln!("SKIP: test input not present: {} (set UNO_REQUIRE_TEST_DATA=1 to make this a failure)", p.display());
    None
}
