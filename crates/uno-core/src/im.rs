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

pub struct ImModel {
    model: GBDT,
}

impl ImModel {
    pub fn load(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let model = GBDT::load_model(&path.to_string_lossy())
            .map_err(|e| anyhow::anyhow!("loading im model {}: {e}", path.display()))?;
        Ok(ImModel { model })
    }

    pub fn predict_many(&self, items: &[(&str, u8)]) -> Vec<f64> {
        let data: Vec<Data> = items
            .iter()
            .map(|&(s, z)| Data::new_test_data(to_f32(&im_features(s, z)), None))
            .collect();
        self.model.predict(&data).into_iter().map(|v| v as f64).collect()
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
