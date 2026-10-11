//! Minimal, hand-written FFI bindings to LightGBM's stable C API
//! (`LGBM_Dataset*` / `LGBM_Booster*`), linked against the vendored
//! `lib_lightgbm.so` (see `vendor/README.md`). No `bindgen`/`cmake` — the C API
//! surface used here (create dataset from a row-major matrix, set the label
//! field, train N boosting iterations, predict) has been stable for years, and
//! is validated against the real library in the test below.
//!
//! Why this exists: uno's pure-Rust `gbdt` crate PFM classifier repeatedly
//! failed to separate decoys from real-but-unconfirmed peptides where
//! ms1searchpy's LightGBM classifier does (see uno-core's
//! `ML_RESCORE_PORT.md`). This binds the real thing instead of continuing to
//! approximate it.

use anyhow::{bail, Result};
use std::ffi::{c_void, CStr, CString};
use std::os::raw::{c_char, c_int};

type DatasetHandle = *mut c_void;
type BoosterHandle = *mut c_void;

const C_API_DTYPE_FLOAT32: c_int = 0;
const C_API_PREDICT_NORMAL: c_int = 0;

#[allow(non_snake_case)]
extern "C" {
    fn LGBM_DatasetCreateFromMat(
        data: *const c_void,
        data_type: c_int,
        nrow: i32,
        ncol: i32,
        is_row_major: c_int,
        parameters: *const c_char,
        reference: DatasetHandle,
        out: *mut DatasetHandle,
    ) -> c_int;
    fn LGBM_DatasetSetField(
        handle: DatasetHandle,
        field_name: *const c_char,
        field_data: *const c_void,
        num_element: c_int,
        field_type: c_int,
    ) -> c_int;
    fn LGBM_DatasetFree(handle: DatasetHandle) -> c_int;
    fn LGBM_BoosterCreate(
        train_data: DatasetHandle,
        parameters: *const c_char,
        out: *mut BoosterHandle,
    ) -> c_int;
    fn LGBM_BoosterUpdateOneIter(handle: BoosterHandle, is_finished: *mut c_int) -> c_int;
    fn LGBM_BoosterFree(handle: BoosterHandle) -> c_int;
    #[allow(clippy::too_many_arguments)]
    fn LGBM_BoosterPredictForMat(
        handle: BoosterHandle,
        data: *const c_void,
        data_type: c_int,
        nrow: i32,
        ncol: i32,
        is_row_major: c_int,
        predict_type: c_int,
        start_iteration: c_int,
        num_iteration: c_int,
        parameter: *const c_char,
        out_len: *mut i64,
        out_result: *mut f64,
    ) -> c_int;
    fn LGBM_GetLastError() -> *const c_char;
}

fn last_error() -> String {
    unsafe {
        let p = LGBM_GetLastError();
        if p.is_null() {
            "unknown LightGBM error".to_string()
        } else {
            CStr::from_ptr(p).to_string_lossy().into_owned()
        }
    }
}

fn check(rc: c_int) -> Result<()> {
    if rc != 0 {
        bail!("LightGBM error: {}", last_error());
    }
    Ok(())
}

fn flatten(feats: &[Vec<f32>]) -> (Vec<f32>, usize, usize) {
    let ncol = feats.first().map_or(0, |r| r.len());
    let nrow = feats.len();
    let mut flat = Vec::with_capacity(nrow * ncol);
    for row in feats {
        debug_assert_eq!(row.len(), ncol, "ragged feature matrix");
        flat.extend_from_slice(row);
    }
    (flat, nrow, ncol)
}

/// A trained LightGBM binary classifier. Owns the C booster handle.
pub struct LgbmBooster {
    handle: BoosterHandle,
}

// The handle is an opaque pointer into LightGBM's own heap; LightGBM's C API
// has no thread affinity requirement on a booster used from a single owning
// thread at a time, which is how uno uses it (trained and predicted from
// within one CV fold's sequential code path).
unsafe impl Send for LgbmBooster {}

impl Drop for LgbmBooster {
    fn drop(&mut self) {
        unsafe {
            let _ = LGBM_BoosterFree(self.handle);
        }
    }
}

impl LgbmBooster {
    /// Train a binary classifier. `feats[i]` is row i's feature vector (all
    /// rows must share the same length); `labels[i]` is 0.0 or 1.0. `params` is
    /// a LightGBM parameter string, e.g. `"objective=binary num_leaves=31
    /// learning_rate=0.1 verbosity=-1"`. Runs exactly `num_iterations` boosting
    /// rounds (no early stopping -- the caller's CV loop controls that).
    pub fn train(
        feats: &[Vec<f32>],
        labels: &[f32],
        params: &str,
        num_iterations: usize,
    ) -> Result<Self> {
        if feats.is_empty() {
            bail!("LightGBM train: empty feature matrix");
        }
        let (flat, nrow, ncol) = flatten(feats);
        let cparams = CString::new(params)?;
        let mut dataset: DatasetHandle = std::ptr::null_mut();
        unsafe {
            check(LGBM_DatasetCreateFromMat(
                flat.as_ptr() as *const c_void,
                C_API_DTYPE_FLOAT32,
                nrow as i32,
                ncol as i32,
                1, // row-major
                cparams.as_ptr(),
                std::ptr::null_mut(),
                &mut dataset,
            ))?;
        }
        let label_name = CString::new("label")?;
        let set_label = unsafe {
            check(LGBM_DatasetSetField(
                dataset,
                label_name.as_ptr(),
                labels.as_ptr() as *const c_void,
                labels.len() as c_int,
                C_API_DTYPE_FLOAT32,
            ))
        };
        if let Err(e) = set_label {
            unsafe { LGBM_DatasetFree(dataset) };
            return Err(e);
        }
        let mut booster: BoosterHandle = std::ptr::null_mut();
        let created = unsafe { check(LGBM_BoosterCreate(dataset, cparams.as_ptr(), &mut booster)) };
        if let Err(e) = created {
            unsafe { LGBM_DatasetFree(dataset) };
            return Err(e);
        }
        let mut is_finished: c_int = 0;
        for _ in 0..num_iterations {
            let step = unsafe { check(LGBM_BoosterUpdateOneIter(booster, &mut is_finished)) };
            if let Err(e) = step {
                unsafe {
                    LGBM_BoosterFree(booster);
                    LGBM_DatasetFree(dataset);
                }
                return Err(e);
            }
            if is_finished != 0 {
                break;
            }
        }
        unsafe { LGBM_DatasetFree(dataset) };
        Ok(LgbmBooster { handle: booster })
    }

    /// Predict class-1 probabilities for each row.
    pub fn predict(&self, feats: &[Vec<f32>]) -> Result<Vec<f64>> {
        if feats.is_empty() {
            return Ok(Vec::new());
        }
        let (flat, nrow, ncol) = flatten(feats);
        let mut out_len: i64 = 0;
        let mut out = vec![0.0f64; nrow];
        let empty = CString::new("")?;
        unsafe {
            check(LGBM_BoosterPredictForMat(
                self.handle,
                flat.as_ptr() as *const c_void,
                C_API_DTYPE_FLOAT32,
                nrow as i32,
                ncol as i32,
                1,
                C_API_PREDICT_NORMAL,
                0,
                -1,
                empty.as_ptr(),
                &mut out_len,
                out.as_mut_ptr(),
            ))?;
        }
        debug_assert_eq!(out_len as usize, nrow);
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn separates_linearly_separable_classes() {
        // Deterministic synthetic data (no RNG crate needed): label = 1 iff
        // x+y > 0, plus an irrelevant third feature. A working classifier
        // should assign high scores to label==1 and low to label==0.
        let mut feats = Vec::new();
        let mut labels = Vec::new();
        for i in 0..400u32 {
            let x = ((i.wrapping_mul(2654435761)) % 1000) as f32 / 500.0 - 1.0;
            let y = ((i.wrapping_mul(40503).wrapping_add(7)) % 1000) as f32 / 500.0 - 1.0;
            let z = ((i.wrapping_mul(2246822519)) % 1000) as f32 / 500.0 - 1.0;
            let label = if x + y > 0.0 { 1.0 } else { 0.0 };
            feats.push(vec![x, y, z]);
            labels.push(label);
        }
        let booster = LgbmBooster::train(
            &feats,
            &labels,
            "objective=binary num_leaves=15 learning_rate=0.1 verbosity=-1 min_data_in_leaf=5",
            50,
        )
        .expect("training failed");
        let preds = booster.predict(&feats).expect("predict failed");
        let (mut pos_sum, mut pos_n, mut neg_sum, mut neg_n) = (0.0, 0usize, 0.0, 0usize);
        for (p, &l) in preds.iter().zip(&labels) {
            if l > 0.5 {
                pos_sum += p;
                pos_n += 1;
            } else {
                neg_sum += p;
                neg_n += 1;
            }
        }
        let pos_mean = pos_sum / pos_n as f64;
        let neg_mean = neg_sum / neg_n as f64;
        assert!(
            pos_mean > neg_mean + 0.3,
            "expected clear separation, got pos_mean={pos_mean} neg_mean={neg_mean}"
        );
    }
}
