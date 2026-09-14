//! Ion-mobility (1/K0) prediction: a gbdt model distilled from im2deep.
//!
//! Like the RT model, but charge-aware — 1/K0 depends strongly on charge state.
//! im2deep (`-i`) predicts 1/K0 directly in the same units as koth's `im`
//! column, so the model predicts 1/K0 and a per-dataset calibration maps the
//! im2deep reference scale to the instrument's. The sequence-specific deviation
//! from the mass/charge trend is what discriminates the correct peptide.

use crate::rt::rt_features;
use anyhow::Result;
use gbdt::decision_tree::Data;
use gbdt::gradient_boost::GBDT;
use std::path::Path;

/// CCS/IM features = RT composition features plus charge (the dominant factor).
pub const IM_N_FEATURES: usize = crate::rt::RT_N_FEATURES + 1;

pub fn im_features(seq: &str, charge: u8) -> Vec<f64> {
    let mut v = rt_features(seq);
    v.push(charge as f64);
    v
}

fn to_f32(v: &[f64]) -> Vec<f32> {
    v.iter().map(|&x| x as f32).collect()
}

/// Either backend. `gbdt` is the legacy single-threaded crate the shipped
/// `im_im2deep.gbdt` was trained with; `sequoia-boost` is uno's own pure-Rust
/// XGBoost reimplementation (already the default ML backend for scoring) and
/// trains on the ambient Rayon pool, which is ~20x faster here — the legacy
/// crate needed ~28 min on the 270k-row public CCS corpus.
pub enum ImModel {
    Gbdt(GBDT),
    Sequoia(Box<sequoia_boost::prelude::BoostedModel>),
}

impl ImModel {
    /// Load by content: a sequoia native model first, else the legacy gbdt.
    pub fn load(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        if let Ok(m) = sequoia_boost::prelude::BoostedModel::load_binary(path) {
            return Ok(ImModel::Sequoia(Box::new(m)));
        }
        let metadata: serde_json::Value = serde_json::from_slice(&std::fs::read(path)?)?;
        let width = metadata["conf"]["feature_size"].as_u64()
            .ok_or_else(|| anyhow::anyhow!("IM GBDT has no feature-size metadata"))? as usize;
        anyhow::ensure!(width == IM_N_FEATURES,
            "IM model expects {} features, current encoder supplies {}. Retrain with the current encoder; do not silently reinterpret an old schema.", width, IM_N_FEATURES);
        let model = GBDT::load_model(&path.to_string_lossy())
            .map_err(|e| anyhow::anyhow!("loading im model {}: {e}", path.display()))?;
        Ok(ImModel::Gbdt(model))
    }

    /// Predict each distinct sequence/charge once. Bound temporary feature
    /// matrices so a large PFM list cannot allocate several gigabytes of input.
    pub fn predict_many(&self, items: &[(&str, u8)]) -> Vec<f64> {
        let mut lookup = std::collections::HashMap::new();
        let mut unique = Vec::new();
        let mut indices = Vec::with_capacity(items.len());
        for &item in items {
            let idx = *lookup.entry(item).or_insert_with(|| {
                unique.push(item);
                unique.len() - 1
            });
            indices.push(idx);
        }
        let mut predictions = Vec::with_capacity(unique.len());
        for chunk in unique.chunks(16_384) {
            predictions.extend(self.predict_distinct(chunk));
        }
        indices.into_iter().map(|i| predictions[i]).collect()
    }

    fn predict_distinct(&self, items: &[(&str, u8)]) -> Vec<f64> {
        match self {
            ImModel::Gbdt(model) => {
                let data: Vec<Data> = items.iter()
                    .map(|&(s, z)| Data::new_test_data(to_f32(&im_features(s, z)), None))
                    .collect();
                model.predict(&data).into_iter().map(|v| v as f64).collect()
            }
            ImModel::Sequoia(model) => {
                use sequoia_boost::prelude::*;
                let mut flat = Vec::with_capacity(items.len() * IM_N_FEATURES);
                for &(s, z) in items { flat.extend_from_slice(&to_f32(&im_features(s, z))); }
                let matrix = DMatrix::from_dense(&flat, items.len(), IM_N_FEATURES)
                    .expect("validated IM feature matrix dimensions");
                model.predict(&matrix).expect("validated IM model dimensions")
                    .into_iter().map(|v| v as f64).collect()
            }
        }
    }

}

/// Read an im2deep `-i` output CSV: `modified_seq,charge,predicted IM`, where
/// modified_seq is `SEQUENCE/charge`. Returns (stripped_seq, charge, im).
pub fn read_im2deep_csv(path: impl AsRef<Path>) -> Result<(Vec<String>, Vec<u8>, Vec<f64>)> {
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .flexible(true)
        .from_path(path.as_ref())?;
    let (mut seqs, mut charges, mut ims) = (Vec::new(), Vec::new(), Vec::new());
    for rec in rdr.records() {
        let rec = rec?;
        let modseq = match rec.get(0) {
            Some(s) => s,
            None => continue,
        };
        let seq = modseq.split('/').next().unwrap_or("").to_string();
        let charge: u8 = match rec.get(1).and_then(|s| s.trim().parse().ok()) {
            Some(z) => z,
            None => continue,
        };
        let im: f64 = match rec.get(2).and_then(|s| s.trim().parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        if seq.is_empty() || !seq.bytes().all(|b| b.is_ascii_uppercase()) {
            continue;
        }
        seqs.push(seq);
        charges.push(charge);
        ims.push(im);
    }
    Ok((seqs, charges, ims))
}

#[cfg(test)]
mod schema_tests {
    use super::*;
    #[test]
    fn incompatible_legacy_schema_is_rejected_before_prediction() {
        let p = std::env::temp_dir().join(format!("uno-im-schema-{}.json", std::process::id()));
        std::fs::write(&p, r#"{"conf":{"feature_size":70}}"#).unwrap();
        let err = ImModel::load(&p).err().expect("old feature width must fail");
        std::fs::remove_file(p).unwrap();
        assert!(err.to_string().contains("expects 70 features, current encoder supplies 100"));
    }
    #[test]
    fn charge_is_last_and_sequence_descriptors_are_charge_independent() {
        let a = im_features("PEPTIDEK", 2);
        let b = im_features("PEPTIDEK", 3);
        assert_eq!(a.len(), IM_N_FEATURES);
        assert_eq!(&a[..IM_N_FEATURES-1], &b[..IM_N_FEATURES-1]);
        assert_eq!(a[IM_N_FEATURES-1], 2.0);
        assert_eq!(b[IM_N_FEATURES-1], 3.0);
    }
}
