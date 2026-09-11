//! uno-core: MS1-only protein identification.
//!
//! Consumes a SINGLE-RUN biosaur2 feature file (the same input ms1searchpy uses),
//! assigns peptide/protein identities, and controls FDR. Single-run input is
//! required: it carries the per-feature quality columns that separate target from
//! decoy, which cross-run alignment drops.

pub mod calibrate;
mod broad_rt;
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

/// Mix a user-chosen seed into the fixed constants behind uno's stochastic
/// stages: the global-pool decoy shuffle, RANSAC pair sampling, CV fold
/// assignment, and the tree booster's row/column subsampling. Seed 0 is the
/// shipped default and returns 0, so every constant is left untouched and
/// historical output is reproduced byte for byte; any other seed perturbs all
/// of those stages together (`uno --seed N`). Used to measure uno's seed-to-seed
/// spread the same way ms1searchpy's run-to-run spread is measured.
pub fn seed_mix(seed: u64) -> u64 {
    if seed == 0 {
        return 0;
    }
    // splitmix64 finaliser: a small integer becomes a well-mixed 64-bit word.
    let mut z = seed.wrapping_add(0x9e37_79b9_7f4a_7c15);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

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
