//! Persistent, model-specific raw IM predictions. Calibration is never cached.
//!
//! Keys include the model bytes, cache/encoder version, sequence, and charge.
//! A loaded model retains its map across calibration attempts and batch runs.
//! File locks protect concurrent readers/writers; duplicate computations by
//! separate cold processes are harmless. Failed disk writes retain memory hits.

use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

type Predictions = HashMap<u8, HashMap<String, f64>>;
pub(crate) const BATCH_SIZE: usize = 16_384;

pub(crate) struct PredictionCache {
    pub enabled: bool,
    path: Option<PathBuf>,
    loaded: bool,
    values: Predictions,
}

impl PredictionCache {
    #[cfg(test)]
    pub fn isolated(path: PathBuf, enabled: bool) -> Self {
        Self {
            enabled,
            path: Some(path),
            loaded: false,
            values: HashMap::new(),
        }
    }

    pub fn new(fingerprint: &str) -> Self {
        Self {
            enabled: std::env::var("UNO_IM_CACHE").as_deref() != Ok("0"),
            path: cache_dir().map(|d| d.join(format!("{fingerprint}.tsv"))),
            loaded: false,
            values: HashMap::new(),
        }
    }

    /// Return predictions in input order, plus the number of unique misses.
    pub fn predict(
        &mut self,
        items: &[(&str, u8)],
        infer: impl Fn(&[(&str, u8)]) -> Vec<f64>,
    ) -> (Vec<f64>, usize) {
        if items.is_empty() {
            return (Vec::new(), 0);
        }
        if !self.enabled {
            return (infer(items), items.iter().collect::<HashSet<_>>().len());
        }
        if !self.loaded {
            if let Some(path) = &self.path {
                match load(path) {
                    Ok(values) => self.values = values,
                    Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
                    Err(e) => log::warn!("IM cache read {}: {e}", path.display()),
                }
            }
            self.loaded = true;
        }
        let mut seen = HashSet::new();
        let missing: Vec<_> = items
            .iter()
            .copied()
            .filter(|&(s, z)| {
                !self.values.get(&z).is_some_and(|m| m.contains_key(s)) && seen.insert((s, z))
            })
            .collect();
        for chunk in missing.chunks(BATCH_SIZE) {
            let predictions = infer(chunk);
            assert_eq!(predictions.len(), chunk.len(), "IM predictor output length");
            if let Some(path) = &self.path {
                if let Err(e) = append(path, chunk, &predictions) {
                    log::warn!(
                        "IM cache write {}: {e}; retaining in-memory predictions",
                        path.display()
                    );
                }
            }
            for (&(s, z), p) in chunk.iter().zip(predictions) {
                self.values.entry(z).or_default().insert(s.to_owned(), p);
            }
        }
        (
            items.iter().map(|&(s, z)| self.values[&z][s]).collect(),
            missing.len(),
        )
    }
}

fn cache_dir() -> Option<PathBuf> {
    if let Ok(d) = std::env::var("UNO_IM_CACHE_DIR") {
        return Some(PathBuf::from(d));
    }
    if let Ok(d) = std::env::var("XDG_CACHE_HOME") {
        if !d.is_empty() {
            return Some(PathBuf::from(d).join("uno/im"));
        }
    }
    std::env::var_os("HOME").map(|h| PathBuf::from(h).join(".cache/uno/im"))
}

fn load(path: &Path) -> std::io::Result<Predictions> {
    let file = std::fs::File::open(path)?;
    file.lock_shared()?;
    let mut map: Predictions = HashMap::new();
    for line in BufReader::new(&file).lines() {
        let line = line?;
        let fields: Vec<_> = line.split('\t').collect();
        // Fixed-width bits and a terminal marker reject incomplete appends.
        if fields.len() != 4
            || fields[3] != "IM1"
            || fields[2].len() != 16
            || fields[0].is_empty()
            || !fields[0].bytes().all(|b| b.is_ascii_uppercase())
        {
            continue;
        }
        if let (Ok(z), Ok(bits)) = (fields[1].parse::<u8>(), u64::from_str_radix(fields[2], 16)) {
            let p = f64::from_bits(bits);
            if z > 0 && p.is_finite() {
                map.entry(z).or_default().insert(fields[0].to_owned(), p);
            }
        }
    }
    Ok(map)
}

fn append(path: &Path, items: &[(&str, u8)], values: &[f64]) -> std::io::Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let file = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)?;
    file.lock()?;
    let mut out = BufWriter::new(&file);
    // Delimit any partial record left by a killed writer before adding rows.
    writeln!(out)?;
    for (&(s, z), p) in items.iter().zip(values) {
        writeln!(out, "{s}\t{z}\t{:016x}\tIM1", p.to_bits())?;
    }
    out.flush()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn directory() -> PathBuf {
        let path = std::env::temp_dir().join(format!(
            "uno-im-cache-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&path).unwrap();
        path
    }

    fn cache(path: PathBuf) -> PredictionCache {
        PredictionCache {
            enabled: true,
            path: Some(path),
            loaded: false,
            values: HashMap::new(),
        }
    }

    #[test]
    fn memory_and_disk_hits_preserve_bits_order_charge_and_only_infer_misses() {
        let dir = directory();
        let path = dir.join("model.tsv");
        let rows = [("PEPTIDEK", 2), ("PEPTIDEK", 3), ("PEPTIDEK", 2)];
        let mut first = cache(path.clone());
        let (pred, misses) = first.predict(&rows, |queries| {
            assert_eq!(queries, &rows[..2]);
            vec![f64::from_bits(0x3ff0000000000001), -0.0]
        });
        assert_eq!(misses, 2);
        let bytes = std::fs::read(&path).unwrap();
        let no_inference = |_: &[(&str, u8)]| -> Vec<f64> { panic!("cache hit must not infer") };
        let (warm, n) = first.predict(&rows, no_inference);
        assert_eq!(n, 0);
        let mut second = cache(path.clone());
        let (disk, n) = second.predict(&rows, no_inference);
        assert_eq!(n, 0);
        for (a, b) in warm.iter().chain(&disk).zip(pred.iter().cycle()) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
        assert_eq!(std::fs::read(&path).unwrap(), bytes);
        let (_, n) = second.predict(&[("PEPTIDEK", 4), rows[0]], |q| {
            assert_eq!(q, &[("PEPTIDEK", 4)]);
            vec![1.5]
        });
        assert_eq!(n, 1);
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn disabled_and_empty_queries_never_touch_disk() {
        let dir = directory();
        let path = dir.join("unused.tsv");
        let mut c = cache(path.clone());
        c.predict(&[], |_| panic!("empty request"));
        c.enabled = false;
        for _ in 0..2 {
            let (p, n) = c.predict(&[("PEPTIDEK", 2)], |_| vec![1.25]);
            assert_eq!((p, n), (vec![1.25], 1));
        }
        assert!(!path.exists());
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn partial_rows_and_concurrent_appends_are_safe() {
        let dir = directory();
        let path = dir.join("model.tsv");
        std::fs::write(&path, "BAD\t2\t3ff0000000000000\tIM").unwrap();
        std::thread::scope(|scope| {
            for z in 1..=8 {
                let path = &path;
                scope.spawn(move || append(path, &[("PEPTIDEK", z)], &[z as f64]).unwrap());
            }
        });
        let disk = load(&path).unwrap();
        assert_eq!(disk.values().map(HashMap::len).sum::<usize>(), 8);
        for z in 1..=8 {
            assert_eq!(disk[&z]["PEPTIDEK"], z as f64);
        }
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn unwritable_cache_still_retains_memory_predictions() {
        let dir = directory();
        let mut c = cache(dir.clone()); // A directory cannot be opened for append.
        c.predict(&[("PEPTIDEK", 2)], |_| vec![1.2]);
        let (p, n) = c.predict(&[("PEPTIDEK", 2)], |_| panic!("memory hit"));
        assert_eq!((p, n), (vec![1.2], 0));
        std::fs::remove_dir_all(dir).unwrap();
    }
}
