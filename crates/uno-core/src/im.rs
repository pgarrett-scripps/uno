//! Charge-aware IM prediction with the shared RT sequence CNN or legacy tree models.
use crate::rt::rt_features;
use anyhow::Result;
use gbdt::decision_tree::Data;
use gbdt::gradient_boost::GBDT;
use std::path::Path;

/// Portable sequence CNN, trained on public measured TIMS labels.
pub const DEFAULT_IM_CNN: &[u8] = include_bytes!("../../../models/im_cnn_public_20260914.uimc");

/// Tree-model inputs: RT composition descriptors plus charge.
pub const IM_N_FEATURES: usize = crate::rt::RT_N_FEATURES + 1;

pub fn im_features(seq: &str, charge: u8) -> Vec<f64> {
    let mut v = rt_features(seq);
    v.push(charge as f64);
    v
}

fn to_f32(v: &[f64]) -> Vec<f32> {
    v.iter().map(|&x| x as f32).collect()
}

/// Supported sequence CNN and legacy tree prediction backends.
enum ImBackend {
    Cnn(crate::rt_cnn::RtCnn),
    Gbdt(GBDT),
    Sequoia(Box<sequoia_boost::prelude::BoostedModel>),
}

pub struct ImModel {
    backend: ImBackend,
    fingerprint: String,
    cache: std::sync::Mutex<crate::im_cache::PredictionCache>,
}

impl ImModel {
    fn from_backend(backend: ImBackend, bytes: &[u8]) -> Self {
        use sha2::{Digest, Sha256};
        let mut hash = Sha256::new();
        // Bump this namespace if encoding/inference semantics change.
        hash.update(b"uno-im-prediction-v1\0");
        hash.update(bytes);
        let fingerprint = format!("{:x}", hash.finalize());
        let cache = std::sync::Mutex::new(crate::im_cache::PredictionCache::new(&fingerprint));
        Self { backend, fingerprint, cache }
    }

    pub fn fingerprint(&self) -> &str { &self.fingerprint }

    /// The default charge-aware model needs no model path or Python runtime.
    pub fn bundled() -> Result<Self> {
        Ok(Self::from_backend(ImBackend::Cnn(crate::rt_cnn::RtCnn::from_im_bytes(DEFAULT_IM_CNN)?), DEFAULT_IM_CNN))
    }

    /// Load by content: UIMC sequence CNN, sequoia tree, or legacy gbdt.
    pub fn load(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let bytes = std::fs::read(path)?;
        if bytes.starts_with(b"UIMC") {
            return Ok(Self::from_backend(ImBackend::Cnn(crate::rt_cnn::RtCnn::from_im_bytes(&bytes)?), &bytes));
        }
        if let Ok(m) = sequoia_boost::prelude::BoostedModel::load_binary(path) {
            anyhow::ensure!(m.n_features() == IM_N_FEATURES,
                "IM model expects {} features, current encoder supplies {}. Retrain with the current encoder; do not silently reinterpret an old schema.", m.n_features(), IM_N_FEATURES);
            return Ok(Self::from_backend(ImBackend::Sequoia(Box::new(m)), &bytes));
        }
        let metadata: serde_json::Value = serde_json::from_slice(&std::fs::read(path)?)?;
        let width = metadata["conf"]["feature_size"].as_u64()
            .ok_or_else(|| anyhow::anyhow!("IM GBDT has no feature-size metadata"))? as usize;
        anyhow::ensure!(width == IM_N_FEATURES,
            "IM model expects {} features, current encoder supplies {}. Retrain with the current encoder; do not silently reinterpret an old schema.", width, IM_N_FEATURES);
        let model = GBDT::load_model(&path.to_string_lossy())
            .map_err(|e| anyhow::anyhow!("loading im model {}: {e}", path.display()))?;
        Ok(Self::from_backend(ImBackend::Gbdt(model), &bytes))
    }

    pub fn cnn_info(&self) -> Option<&serde_json::Value> {
        match &self.backend { ImBackend::Cnn(m) => Some(m.info()), _ => None }
    }

    /// Raw predictions shared across calibration attempts and search runs.
    pub fn predict_many_cached(&self, items: &[(&str, u8)]) -> Vec<f64> {
        let t = std::time::Instant::now();
        let (values, misses) = self.cache.lock().unwrap().predict(items, |q| self.predict_many(q));
        log::info!("IM predictions: {} unique misses of {} requested ({:.3}s)",
            misses, items.len(), t.elapsed().as_secs_f64());
        values
    }

    /// Prime every target/decoy library sequence at the supplied observed charges.
    /// Bounded query batches avoid materializing the full Cartesian product.
    pub fn precompute_library(&self, seqs: &[&str], charges: &[u8]) {
        if seqs.is_empty() || charges.is_empty() { return; }
        let t = std::time::Instant::now();
        let mut cache = self.cache.lock().unwrap();
        // Disabled means the original uncached, on-demand path, including no disk IO.
        if !cache.enabled { return; }
        let mut seqs = seqs.to_vec();
        seqs.sort_unstable(); seqs.dedup();
        let mut charges = charges.to_vec();
        charges.retain(|&z| z > 0); charges.sort_unstable(); charges.dedup();
        let mut misses = 0;
        for &z in &charges {
            for chunk in seqs.chunks(crate::im_cache::BATCH_SIZE) {
                let queries: Vec<_> = chunk.iter().map(|&s| (s, z)).collect();
                misses += cache.predict(&queries, |q| self.predict_many(q)).1;
            }
        }
        log::info!("IM library cache: {} sequences x {} charges = {} pairs; {} unique misses ({:.3}s)",
            seqs.len(), charges.len(), seqs.len() * charges.len(), misses, t.elapsed().as_secs_f64());
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
        match &self.backend {
            ImBackend::Cnn(model) => model.predict_im_many(items),
            ImBackend::Gbdt(model) => {
                let data: Vec<Data> = items.iter()
                    .map(|&(s, z)| Data::new_test_data(to_f32(&im_features(s, z)), None))
                    .collect();
                model.predict(&data).into_iter().map(|v| v as f64).collect()
            }
            ImBackend::Sequoia(model) => {
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
    fn incompatible_sequoia_schema_is_rejected_at_load() {
        use sequoia_boost::prelude::*;
        let d = DMatrix::from_dense(&[0.0, 1.0, 2.0, 3.0], 4, 1).unwrap()
            .with_labels(&[0.0, 1.0, 2.0, 3.0]).unwrap();
        let params = TrainingParams::builder().objective("reg:squarederror").build().unwrap();
        let model = train(&params, &d, 1).unwrap();
        let p = std::env::temp_dir().join(format!("uno-im-schema-{}.sequoia", std::process::id()));
        model.save_binary(&p).unwrap();
        let err = ImModel::load(&p).err().expect("wrong feature width must fail");
        std::fs::remove_file(p).unwrap();
        assert!(err.to_string().contains("expects 1 features, current encoder supplies 100"));
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

/// Label-free per-charge calibration on the RT shortlist, treating targets and
/// decoys identically. Missing/invalid IM never triggers prediction or rejection.
/// The relative residual is a classifier input, never an acceptance threshold.
pub fn soft_residuals(model: &ImModel, items: &[(&str, u8, Option<f64>)]) -> Vec<Option<f64>> {
    let valid: Vec<usize> = items.iter().enumerate().filter_map(|(i, &(s, z, o))| {
        (z > 0 && !s.is_empty() && o.is_some_and(|v| v.is_finite() && v > 0.0)).then_some(i)
    }).collect();
    let queries: Vec<(&str, u8)> = valid.iter().map(|&i| (items[i].0, items[i].1)).collect();
    let predictions = model.predict_many_cached(&queries);
    let observed: Vec<(u8, f64, f64)> = valid.iter().zip(&predictions)
        .map(|(&i, &p)| (items[i].1, items[i].2.unwrap(), p)).collect();
    let scaled = calibrated_residuals(&observed);
    let mut result = vec![None; items.len()];
    for (&i, r) in valid.iter().zip(scaled) { result[i] = r; }
    let unique = queries.iter().collect::<std::collections::HashSet<_>>().len();
    log::info!("IM soft feature: {} unique sequence/charge predictions for {} RT-shortlisted candidates; {} calibrated, no IM rejection",
        unique, items.len(), result.iter().filter(|r| r.is_some()).count());
    result
}

fn calibrated_residuals(items: &[(u8, f64, f64)]) -> Vec<Option<f64>> {
    let mut ratios: std::collections::BTreeMap<u8, Vec<f64>> = Default::default();
    for &(z, o, p) in items {
        if z > 0 && o.is_finite() && o > 0.0 && p.is_finite() && p > 0.0 {
            ratios.entry(z).or_default().push(o / p);
        }
    }
    let scales: std::collections::BTreeMap<u8, f64> = ratios.into_iter().filter_map(|(z, mut v)| {
        // Do not infer a charge calibration from a handful of candidates.
        if v.len() < 20 { return None; }
        v.sort_by(f64::total_cmp);
        let n = v.len();
        let median = (v[(n - 1) / 2] + v[n / 2]) / 2.0;
        (median.is_finite() && median > 0.0).then_some((z, median))
    }).collect();
    items.iter().map(|&(z, o, p)| {
        if !(o.is_finite() && o > 0.0 && p.is_finite() && p > 0.0) { return None; }
        scales.get(&z).and_then(|s| {
            let residual = (o - p * s) / o;
            residual.is_finite().then_some(residual)
        })
    }).collect()
}

#[cfg(test)]
mod soft_tests {
    use super::*;
    #[test]
    fn full_library_cache_preserves_predictions_and_recalibrates_each_shortlist() {
        let dir = std::env::temp_dir().join(format!("uno-im-library-{}-{}",
            std::process::id(), std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos()));
        let mut model = ImModel::bundled().unwrap();
        let path = dir.join(format!("{}.tsv", model.fingerprint()));
        model.cache = std::sync::Mutex::new(crate::im_cache::PredictionCache::isolated(path.clone(), true));
        let queries = [("PEPTIDEK", 2), ("AAAAAAAR", 3), ("PEPTIDEK", 3), ("AAAAAAAR", 2)];
        let direct = model.predict_many(&queries);
        model.precompute_library(&["PEPTIDEK", "AAAAAAAR", "PEPTIDEK"], &[2, 3, 2, 0]);
        let bytes = std::fs::read(&path).unwrap();
        // A backend that cannot infer proves both memory and fresh-model disk hits.
        model.backend = ImBackend::Gbdt(GBDT::new(&gbdt::config::Config::new()));
        let cached = model.predict_many_cached(&queries);
        assert_eq!(cached.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            direct.iter().map(|v| v.to_bits()).collect::<Vec<_>>());
        let mut reloaded = ImModel::bundled().unwrap();
        reloaded.cache = std::sync::Mutex::new(crate::im_cache::PredictionCache::isolated(path.clone(), true));
        reloaded.backend = ImBackend::Gbdt(GBDT::new(&gbdt::config::Config::new()));
        assert_eq!(reloaded.predict_many_cached(&queries), direct);
        for scale in [2.0, 3.0] {
            let rows: Vec<_> = (0..30).flat_map(|_| queries.iter().zip(&direct)
                .map(|(&(s, z), &p)| (s, z, Some(p * scale)))).collect();
            let expected: Vec<_> = (0..30).flat_map(|_| queries.iter().zip(&direct)
                .map(|(&(_, z), &p)| (z, p * scale, p))).collect();
            assert_eq!(soft_residuals(&reloaded, &rows), calibrated_residuals(&expected));
        }
        assert_eq!(std::fs::read(&path).unwrap(), bytes);
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn model_content_controls_cache_identity() {
        let a = ImModel::bundled().unwrap();
        let same = ImModel::bundled().unwrap();
        assert_eq!(a.fingerprint(), same.fingerprint());
        let different = ImModel::from_backend(
            ImBackend::Gbdt(GBDT::new(&gbdt::config::Config::new())), b"different-weights");
        assert_ne!(a.fingerprint(), different.fingerprint());
    }

    #[test]
    fn bundled_cnn_is_portable_and_charge_aware() {
        let m = ImModel::bundled().unwrap();
        assert_eq!(m.cnn_info().unwrap()["format"], "UIMC");
        let p = m.predict_many(&[("PEPTIDEK", 2), ("PEPTIDEK", 3), ("PEPTIDEK", 2)]);
        assert!(p.iter().all(|v| v.is_finite() && *v > 0.0));
        assert_eq!(p[0], p[2]);
        assert!((p[0] - p[1]).abs() > 0.01);
    }
    #[test]
    fn missing_im_never_invokes_the_untrained_predictor() {
        let model = ImModel::from_backend(ImBackend::Gbdt(GBDT::new(&gbdt::config::Config::new())), b"untrained-test");
        let rows = [("PEPTIDEK", 2, None), ("PEPTIDER", 3, Some(f64::NAN)), ("AAAAK", 2, Some(0.0))];
        assert_eq!(soft_residuals(&model, &rows), vec![None, None, None]);
    }
    #[test]
    fn calibration_is_charge_specific_and_permutation_equivariant() {
        let mut items = vec![(2, 2.0, 1.0); 30];
        items.extend(vec![(3, 3.0, 1.0); 30]);
        items.extend([(2, 4.0, 1.0), (4, 1.0, 1.0), (2, f64::NAN, 1.0), (2, 2.0, -1.0)]);
        let r = calibrated_residuals(&items);
        assert_eq!(r[0], Some(0.0)); assert_eq!(r[30], Some(0.0));
        assert_eq!(r[60], Some(0.5)); assert_eq!(&r[61..], &[None, None, None]);
        items.reverse();
        let mut reverse = calibrated_residuals(&items); reverse.reverse();
        assert_eq!(r, reverse);
    }
}
