//! On-disk cache of abstract RT predictions.
//!
//! The RT model's *abstract* RT is a pure function of `(model, peptide
//! sequence)` — it does not depend on the run being searched (the run-specific
//! part is the on-data calibration applied downstream). So predicting ~1.2M
//! unique peptides (~18s) is repeated in full on every search over the same
//! FASTA/digest/model. This module memoizes those predictions to a file keyed by
//! the model's content fingerprint, so repeat runs only predict sequences not
//! seen before.
//!
//! Cache values are the model's exact `f64` output, written with Rust's
//! round-trippable float formatting, so a cache hit yields a bit-identical
//! prediction to recomputing — search output is unchanged whether cached or not.
//!
//! Location: `$UNO_RT_CACHE_DIR` if set, else `$XDG_CACHE_HOME/uno/rt`, else
//! `~/.cache/uno/rt`. Set `UNO_RT_CACHE=0` to disable (always recompute, never
//! read or write).

use std::collections::HashMap;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

use crate::rt::RtModel;

/// Predict abstract RT for `seqs`, reusing an on-disk cache keyed by the model
/// fingerprint. Only sequences missing from the cache are run through the model;
/// newly predicted `(seq, rt)` pairs are appended to the cache file. The returned
/// vec is aligned to `seqs` and bit-identical to `model.predict_many(seqs)`.
pub fn predict_many_cached(model: &RtModel, seqs: &[&str]) -> Vec<f64> {
    let cache_path = match cache_file(model.fingerprint()) {
        Some(p) => p,
        // Caching disabled or no home dir: fall straight through to the model.
        None => return model.predict_many(seqs),
    };

    BATCH_CACHE.with(|state| {
        let mut state = state.borrow_mut();
        if let Some(models) = state.as_mut() {
            let cached = models
                .entry(model.fingerprint())
                .or_insert_with(|| load_cache(&cache_path));
            predict_into(model, seqs, &cache_path, cached)
        } else {
            predict_into(model, seqs, &cache_path, &mut load_cache(&cache_path))
        }
    })
}

thread_local! {
    static BATCH_CACHE: std::cell::RefCell<Option<HashMap<u64,HashMap<String,f64>>>> = const { std::cell::RefCell::new(None) };
}

/// Retain prediction maps across sequential per-run searches on this thread.
/// The guard restores the previous scope even on early return or panic.
pub fn with_memory_cache<T>(work: impl FnOnce() -> T) -> T {
    struct Restore(Option<HashMap<u64, HashMap<String, f64>>>);
    impl Drop for Restore {
        fn drop(&mut self) {
            BATCH_CACHE.with(|s| *s.borrow_mut() = self.0.take());
        }
    }
    let previous = BATCH_CACHE.with(|s| s.replace(Some(HashMap::new())));
    let _guard = Restore(previous);
    work()
}

fn predict_into(
    model: &RtModel,
    seqs: &[&str],
    path: &std::path::Path,
    cached: &mut HashMap<String, f64>,
) -> Vec<f64> {
    let mut miss_seqs = Vec::new();
    let mut seen = std::collections::HashSet::new();
    for &s in seqs {
        if !cached.contains_key(s) && seen.insert(s) {
            miss_seqs.push(s);
        }
    }
    log::info!(
        "RT cache {}: {} unique misses of {} requested ({} rows cached)",
        path.display(),
        miss_seqs.len(),
        seqs.len(),
        cached.len()
    );
    if !miss_seqs.is_empty() {
        let predictions = model.predict_many(&miss_seqs);
        if let Err(e) = append_cache(path, &miss_seqs, &predictions) {
            log::warn!("RT cache write {}: {}", path.display(), e);
        }
        for (s, p) in miss_seqs.into_iter().zip(predictions) {
            cached.insert(s.to_owned(), p);
        }
    }
    seqs.iter().map(|s| cached[*s]).collect()
}

/// The cache file for a model fingerprint, or `None` if caching is disabled or no
/// cache directory can be determined.
fn cache_file(fingerprint: u64) -> Option<PathBuf> {
    if std::env::var("UNO_RT_CACHE").as_deref() == Ok("0") {
        return None;
    }
    let dir = cache_dir()?;
    Some(dir.join(format!("{fingerprint:016x}.tsv")))
}

fn cache_dir() -> Option<PathBuf> {
    if let Ok(d) = std::env::var("UNO_RT_CACHE_DIR") {
        return Some(PathBuf::from(d));
    }
    if let Ok(d) = std::env::var("XDG_CACHE_HOME") {
        if !d.is_empty() {
            return Some(PathBuf::from(d).join("uno").join("rt"));
        }
    }
    std::env::var("HOME")
        .ok()
        .map(|h| PathBuf::from(h).join(".cache").join("uno").join("rt"))
}

/// Load `seq<TAB>rt` rows into a map. A missing/corrupt file yields an empty map
/// (the run simply recomputes and rewrites).
fn load_cache(path: &PathBuf) -> HashMap<String, f64> {
    let mut map = HashMap::new();
    let f = match std::fs::File::open(path) {
        Ok(f) => f,
        Err(_) => return map,
    };
    for line in BufReader::new(f).lines().map_while(Result::ok) {
        let mut it = line.split('\t');
        if let (Some(seq), Some(rt)) = (it.next(), it.next()) {
            if let Ok(rt) = rt.parse::<f64>() {
                map.insert(seq.to_string(), rt);
            }
        }
    }
    map
}

/// Append newly predicted `(seq, rt)` pairs to the cache, creating the directory
/// and file on first use. Append-only, so concurrent searches can both add
/// entries without clobbering each other's rows.
fn append_cache(path: &std::path::Path, seqs: &[&str], preds: &[f64]) -> std::io::Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let f = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)?;
    let mut w = BufWriter::new(f);
    for (s, p) in seqs.iter().zip(preds) {
        // Rust's default f64 formatting is the shortest round-trippable decimal,
        // so parsing it back yields the identical bit pattern.
        writeln!(w, "{s}\t{p}")?;
    }
    w.flush()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn warm_predictions_are_bit_exact_and_append_only_misses() {
        let dir = std::env::temp_dir().join(format!(
            "uno-cache-test-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let path = dir.join("model.tsv");
        let model = RtModel::bundled().unwrap();
        let seqs = ["PEPTIDEK", "AAAAAAAR", "PEPTIDEK"];
        let expected = model.predict_many(&seqs);
        let mut map = HashMap::new();
        let first = predict_into(&model, &seqs, &path, &mut map);
        let before = std::fs::read(&path).unwrap();
        let second = predict_into(&model, &seqs, &path, &mut map);
        assert_eq!(first, expected);
        assert_eq!(second, expected);
        assert_eq!(std::fs::read(&path).unwrap(), before);
        let disk = load_cache(&path);
        assert_eq!(disk.len(), 2);
        assert_eq!(disk["PEPTIDEK"].to_bits(), expected[0].to_bits());
        std::fs::remove_dir_all(dir).unwrap();
    }
    #[test]
    fn batch_cache_scope_restores_previous_state() {
        with_memory_cache(|| {
            BATCH_CACHE.with(|s| {
                s.borrow_mut().as_mut().unwrap().insert(42, HashMap::new());
            });
            with_memory_cache(|| {
                assert!(BATCH_CACHE.with(|s| s.borrow().as_ref().unwrap().is_empty()))
            });
            assert!(BATCH_CACHE.with(|s| s.borrow().as_ref().unwrap().contains_key(&42)));
        });
        assert!(BATCH_CACHE.with(|s| s.borrow().is_none()));
    }
}
