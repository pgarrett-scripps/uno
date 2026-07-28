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

    let cached = load_cache(&cache_path);

    // Unique misses, in first-seen order (stable, independent of dedup impl).
    let mut miss_seqs: Vec<&str> = Vec::new();
    let mut miss_seen: HashMap<&str, ()> = HashMap::new();
    for &s in seqs {
        if !cached.contains_key(s) && miss_seen.insert(s, ()).is_none() {
            miss_seqs.push(s);
        }
    }

    let hits = seqs.len() - seqs.iter().filter(|s| !cached.contains_key(**s)).count();
    log::info!(
        "RT cache {}: {} unique misses of {} requested ({} rows already cached)",
        cache_path.display(),
        miss_seqs.len(),
        seqs.len(),
        cached.len(),
    );

    let mut cached = cached;
    if !miss_seqs.is_empty() {
        let preds = model.predict_many(&miss_seqs);
        if let Err(e) = append_cache(&cache_path, &miss_seqs, &preds) {
            log::warn!("RT cache: failed to write {}: {e}", cache_path.display());
        }
        for (s, p) in miss_seqs.iter().zip(preds) {
            cached.insert((*s).to_string(), p);
        }
    }
    let _ = hits;

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
fn append_cache(path: &PathBuf, seqs: &[&str], preds: &[f64]) -> std::io::Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let f = std::fs::OpenOptions::new().create(true).append(true).open(path)?;
    let mut w = BufWriter::new(f);
    for (s, p) in seqs.iter().zip(preds) {
        // Rust's default f64 formatting is the shortest round-trippable decimal,
        // so parsing it back yields the identical bit pattern.
        writeln!(w, "{s}\t{p}")?;
    }
    w.flush()
}
