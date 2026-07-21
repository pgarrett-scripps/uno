//! uno-core: MS1-only protein identification + quantitation.
//!
//! Consumes the aligned LFQ matrix from `koth_align`, assigns peptide/protein
//! identities, controls FDR, and rolls koth's intensities up into peptide and
//! protein quant matrices. See `PLAN.md` at the workspace root.

pub mod calibrate;
pub mod database;
pub mod decoy;
pub mod hill;
pub mod im;
pub mod input;
pub mod mass;
pub mod matching;
pub mod protein;
pub mod quant;
pub mod rt;
pub mod rt_cnn;
pub mod rt_mlp;
pub mod score;

pub use calibrate::Calibration;
pub use database::{DigestConfig, PeptideDb, TheoPeptide};
pub use hill::{hill_search, load_hills, Hill, HillConfig, HillHit, HillIndex};
pub use im::{im_features, ImModel, IM_N_FEATURES};
pub use input::{DedupConfig, Feature, FeatureMatrix};
pub use matching::{match_features, match_features_topk, MatchConfig, Pfm};
pub use protein::{infer_proteins, ProteinConfig, ProteinId};
pub use quant::{peptide_quant, protein_quant, QuantRow};
pub use rt::{rt_features, RtModel, RT_N_FEATURES};
pub use rt_cnn::RtCnn;
pub use rt_mlp::RtMlp;
pub use score::{
    cv_lda_rescore, score, score_by_deviation, Identification, PeptideMatch, ScoreConfig,
    ScoreResult,
};
