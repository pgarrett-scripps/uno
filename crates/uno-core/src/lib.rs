//! uno-core: MS1-only protein identification + quantitation.
//!
//! Consumes the aligned LFQ matrix from `koth_align`, assigns peptide/protein
//! identities, controls FDR, and rolls koth's intensities up into peptide and
//! protein quant matrices. See `PLAN.md` at the workspace root.

pub mod calibrate;
pub mod database;
pub mod decoy;
pub mod im;
pub mod input;
pub mod mass;
pub mod matching;
pub mod protein;
pub mod rt;
pub mod rt_cnn;
pub mod rt_mlp;
pub mod score;

pub use calibrate::Calibration;
pub use database::{DigestConfig, PeptideDb, TheoPeptide};
pub use im::{im_features, ImModel, IM_N_FEATURES};
pub use input::{DedupConfig, Feature, FeatureMatrix};
pub use matching::{match_features, match_features_topk, MatchConfig, Pfm};
pub use protein::{ProteinConfig, ProteinId};
pub use rt::{rt_features, RtModel, RT_N_FEATURES};
pub use rt_cnn::RtCnn;
pub use rt_mlp::RtMlp;
pub use score::{score, ScoreConfig, ScoreResult};
