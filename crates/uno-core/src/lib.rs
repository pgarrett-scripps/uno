//! uno-core: MS1-only protein identification.
//!
//! Consumes a SINGLE-RUN biosaur2 feature file (the same input ms1searchpy uses),
//! assigns peptide/protein identities, and controls FDR. Single-run input is
//! required: it carries the per-feature quality columns that separate target from
//! decoy, which cross-run alignment drops.

pub mod calibrate;
pub mod database;
pub mod decoy;
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
