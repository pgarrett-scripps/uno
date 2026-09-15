//! Model-specific raw IM predictions; calibration is never cached.
//!
//! Immutable Zstd-compressed Parquet shards are published only after the writer
//! closes. Readers ignore temporary files, and concurrent writers never replace
//! each other's shards. Failed disk operations retain in-memory predictions.

use parquet::basic::Compression;
use parquet::data_type::{ByteArray, ByteArrayType, DoubleType, Int32Type};
use parquet::file::properties::WriterProperties;
use parquet::file::reader::{FileReader, SerializedFileReader};
use parquet::file::writer::SerializedFileWriter;
use parquet::record::RowAccessor;
use parquet::schema::parser::parse_message_type;
use parquet::schema::types::{ColumnPath, TypePtr};
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::sync::Arc;

type Predictions = HashMap<u8, HashMap<String, f64>>;
pub(crate) const BATCH_SIZE: usize = 16_384;
const FORMAT: &str = "parquet-v2";

pub(crate) struct PredictionCache {
    pub enabled: bool,
    path: PathBuf,
    private_root: Option<PathBuf>,
    loaded: bool,
    values: Predictions,
}

impl PredictionCache {
    #[cfg(test)]
    pub fn isolated(path: PathBuf, enabled: bool) -> Self {
        Self {
            enabled,
            path,
            private_root: None,
            loaded: false,
            values: HashMap::new(),
        }
    }

    pub fn new(fingerprint: &str) -> Self {
        let configured = std::env::var_os("UNO_IM_CACHE_DIR").filter(|s| !s.is_empty());
        let private_root = configured.is_none().then(default_root);
        let root = configured
            .map(PathBuf::from)
            .unwrap_or_else(|| private_root.clone().unwrap());
        Self {
            enabled: std::env::var("UNO_IM_CACHE").as_deref() != Ok("0"),
            path: root.join(FORMAT).join(fingerprint),
            private_root,
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
            // No storage access until there is an enabled, nonempty query.
            match self.prepare_root().and_then(|()| load(&self.path)) {
                Ok(values) => self.values = values,
                Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
                Err(e) => log::warn!("IM cache read {}: {e}", self.path.display()),
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
            if let Err(e) = self
                .prepare_root()
                .and_then(|()| append(&self.path, chunk, &predictions))
            {
                log::warn!(
                    "IM cache write {}: {e}; retaining in-memory predictions",
                    self.path.display()
                );
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

    fn prepare_root(&self) -> std::io::Result<()> {
        match &self.private_root {
            Some(root) => prepare_private_root(root),
            None => Ok(()),
        }
    }
}

fn default_root() -> PathBuf {
    // Unix temp directories are shared; keep each user's cache private. On
    // Windows the OS-provided temp directory is already scoped to the user.
    #[cfg(unix)]
    let name = format!("uno-im-{}", rustix::process::geteuid().as_raw());
    #[cfg(not(unix))]
    let name = "uno-im";
    std::env::temp_dir().join(name)
}

fn prepare_private_root(root: &Path) -> std::io::Result<()> {
    let mut builder = std::fs::DirBuilder::new();
    #[cfg(unix)]
    {
        use std::os::unix::fs::DirBuilderExt;
        builder.mode(0o700);
    }
    match builder.create(root) {
        Ok(()) => {}
        Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => {}
        Err(e) => return Err(e),
    }
    let metadata = std::fs::symlink_metadata(root)?;
    let mut safe = metadata.is_dir() && !metadata.file_type().is_symlink();
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        safe &=
            metadata.uid() == rustix::process::geteuid().as_raw() && metadata.mode() & 0o077 == 0;
    }
    if !safe {
        return Err(std::io::Error::new(
            std::io::ErrorKind::PermissionDenied,
            "IM temporary cache root must be a private directory owned by this user",
        ));
    }
    Ok(())
}

fn schema() -> TypePtr {
    Arc::new(
        parse_message_type(
            "message uno_im_prediction_cache_v2 {
        REQUIRED BINARY sequence (UTF8);
        REQUIRED INT32 charge;
        REQUIRED DOUBLE prediction;
    }",
        )
        .expect("static IM cache schema"),
    )
}

fn invalid_data(error: impl std::fmt::Display) -> std::io::Error {
    std::io::Error::new(std::io::ErrorKind::InvalidData, error.to_string())
}

fn load(path: &Path) -> std::io::Result<Predictions> {
    let mut shards = Vec::new();
    for entry in std::fs::read_dir(path)? {
        let entry = entry?;
        if entry.file_type()?.is_file() && entry.path().extension().is_some_and(|e| e == "parquet")
        {
            shards.push(entry.path());
        }
    }
    shards.sort();
    let mut map: Predictions = HashMap::new();
    for shard in shards {
        // A corrupt/evicted shard is a cache miss; other complete shards survive.
        match load_shard(&shard) {
            Ok(rows) => {
                for (charge, values) in rows {
                    map.entry(charge).or_default().extend(values);
                }
            }
            Err(e) => log::warn!("IM cache ignoring shard {}: {e}", shard.display()),
        }
    }
    Ok(map)
}

fn load_shard(path: &Path) -> std::io::Result<Predictions> {
    let reader = SerializedFileReader::new(std::fs::File::open(path)?).map_err(invalid_data)?;
    if reader
        .metadata()
        .file_metadata()
        .schema_descr()
        .root_schema()
        != schema().as_ref()
    {
        return Err(invalid_data("unsupported IM cache schema"));
    }
    let mut map: Predictions = HashMap::new();
    for row in reader.get_row_iter(None).map_err(invalid_data)? {
        let row = row.map_err(invalid_data)?;
        let seq = row.get_string(0).map_err(invalid_data)?;
        let charge = row.get_int(1).map_err(invalid_data)?;
        let prediction = row.get_double(2).map_err(invalid_data)?;
        if (1..=255).contains(&charge)
            && !seq.is_empty()
            && seq.bytes().all(|b| b.is_ascii_uppercase())
            && prediction.is_finite()
        {
            map.entry(charge as u8)
                .or_default()
                .insert(seq.clone(), prediction);
        }
    }
    Ok(map)
}

fn append(path: &Path, items: &[(&str, u8)], values: &[f64]) -> std::io::Result<()> {
    std::fs::create_dir_all(path)?;
    // Same-directory temporary files keep publication on the same filesystem.
    // Normal errors remove the temporary file; a killed process may leave a
    // .tmp file, which readers never consider a completed Parquet shard.
    let mut temp = tempfile::Builder::new()
        .prefix(".pending-")
        .suffix(".tmp")
        .tempfile_in(path)?;
    let properties = WriterProperties::builder()
        .set_compression(Compression::ZSTD(Default::default()))
        // Avoid floating-point dictionary equality merging +0 and -0.
        .set_column_dictionary_enabled(ColumnPath::from("prediction"), false)
        .build();
    let mut writer = SerializedFileWriter::new(temp.as_file_mut(), schema(), Arc::new(properties))
        .map_err(invalid_data)?;
    let mut group = writer.next_row_group().map_err(invalid_data)?;
    let sequences: Vec<ByteArray> = items.iter().map(|&(s, _)| ByteArray::from(s)).collect();
    let charges: Vec<i32> = items.iter().map(|&(_, z)| i32::from(z)).collect();
    let mut col = group
        .next_column()
        .map_err(invalid_data)?
        .expect("sequence column");
    col.typed::<ByteArrayType>()
        .write_batch(&sequences, None, None)
        .map_err(invalid_data)?;
    col.close().map_err(invalid_data)?;
    let mut col = group
        .next_column()
        .map_err(invalid_data)?
        .expect("charge column");
    col.typed::<Int32Type>()
        .write_batch(&charges, None, None)
        .map_err(invalid_data)?;
    col.close().map_err(invalid_data)?;
    let mut col = group
        .next_column()
        .map_err(invalid_data)?
        .expect("prediction column");
    col.typed::<DoubleType>()
        .write_batch(values, None, None)
        .map_err(invalid_data)?;
    col.close().map_err(invalid_data)?;
    group.close().map_err(invalid_data)?;
    writer.close().map_err(invalid_data)?;
    temp.as_file().sync_all()?;
    let published = path.join(format!(
        "part-{}.parquet",
        temp.path()
            .file_stem()
            .unwrap()
            .to_string_lossy()
            .trim_start_matches(".pending-")
    ));
    temp.persist_noclobber(published).map_err(|e| e.error)?;
    Ok(())
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    fn cache(path: PathBuf) -> PredictionCache {
        PredictionCache::isolated(path, true)
    }

    pub(crate) fn disk_snapshot(path: &Path) -> Vec<(PathBuf, Vec<u8>)> {
        let mut files: Vec<_> = std::fs::read_dir(path)
            .unwrap()
            .map(|e| {
                let p = e.unwrap().path();
                let bytes = std::fs::read(&p).unwrap();
                (p, bytes)
            })
            .collect();
        files.sort();
        files
    }

    #[test]
    fn memory_and_disk_hits_preserve_bits_order_charge_and_only_infer_misses() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("model");
        let rows = [
            ("PEPTIDEK", 2),
            ("PEPTIDEK", 3),
            ("PEPTIDEK", 2),
            ("AAAAAAAR", 2),
        ];
        let mut first = cache(path.clone());
        let (pred, misses) = first.predict(&rows, |queries| {
            assert_eq!(queries, &[rows[0], rows[1], rows[3]]);
            vec![f64::from_bits(0x3ff0000000000001), -0.0, 0.0]
        });
        assert_eq!(misses, 3);
        let files = disk_snapshot(&path);
        assert_eq!(files.len(), 1);
        assert_eq!(files[0].0.extension().unwrap(), "parquet");
        assert_eq!(&files[0].1[..4], b"PAR1");
        let reader = SerializedFileReader::new(std::fs::File::open(&files[0].0).unwrap()).unwrap();
        for column in reader.metadata().row_group(0).columns() {
            assert!(matches!(column.compression(), Compression::ZSTD(_)));
        }
        let no_inference = |_: &[(&str, u8)]| -> Vec<f64> { panic!("cache hit must not infer") };
        let (warm, n) = first.predict(&rows, no_inference);
        assert_eq!(n, 0);
        let mut second = cache(path.clone());
        let (disk, n) = second.predict(&rows, no_inference);
        assert_eq!(n, 0);
        for (a, b) in warm.iter().chain(&disk).zip(pred.iter().cycle()) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
        assert_eq!(disk_snapshot(&path), files);
        let (_, n) = second.predict(&[("PEPTIDEK", 4), rows[0]], |q| {
            assert_eq!(q, &[("PEPTIDEK", 4)]);
            vec![1.5]
        });
        assert_eq!(n, 1);
        assert_eq!(disk_snapshot(&path).len(), 2);
        assert!(disk_snapshot(&path).contains(&files[0]));
    }

    #[test]
    fn disabled_and_empty_queries_never_touch_disk() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("unused");
        let mut c = cache(path.clone());
        c.private_root = Some(dir.path().join("private-unused"));
        c.predict(&[], |_| panic!("empty request"));
        c.enabled = false;
        for _ in 0..2 {
            let (p, n) = c.predict(&[("PEPTIDEK", 2)], |_| vec![1.25]);
            assert_eq!((p, n), (vec![1.25], 1));
        }
        assert!(!path.exists());
        assert!(!c.private_root.unwrap().exists());
    }

    #[test]
    fn partial_and_corrupt_shards_do_not_hide_concurrent_writes() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path();
        std::fs::write(path.join(".pending-interrupted.tmp"), b"PAR1incomplete").unwrap();
        std::fs::write(path.join("corrupt.parquet"), b"PAR1incomplete").unwrap();
        std::thread::scope(|scope| {
            for z in 1..=8 {
                scope.spawn(move || append(path, &[("PEPTIDEK", z)], &[z as f64]).unwrap());
            }
        });
        let disk = load(path).unwrap();
        assert_eq!(disk.values().map(HashMap::len).sum::<usize>(), 8);
        for z in 1..=8 {
            assert_eq!(disk[&z]["PEPTIDEK"], z as f64);
        }
        assert_eq!(std::fs::read_dir(path).unwrap().count(), 10);
    }

    #[test]
    fn unwritable_cache_still_retains_memory_predictions() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("not-a-directory");
        std::fs::write(&path, b"keep me").unwrap();
        let mut c = cache(path.clone());
        c.predict(&[("PEPTIDEK", 2)], |_| vec![1.2]);
        let (p, n) = c.predict(&[("PEPTIDEK", 2)], |_| panic!("memory hit"));
        assert_eq!((p, n), (vec![1.2], 0));
        assert_eq!(std::fs::read(path).unwrap(), b"keep me");
    }

    #[test]
    fn evicted_cache_is_recreated_from_new_misses() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("evictable");
        let mut c = cache(path.clone());
        c.predict(&[("PEPTIDEK", 2)], |_| vec![1.2]);
        std::fs::remove_dir_all(&path).unwrap();
        let (p, n) = c.predict(&[("PEPTIDEK", 2), ("PEPTIDEK", 3)], |q| {
            assert_eq!(q, &[("PEPTIDEK", 3)]);
            vec![1.8]
        });
        assert_eq!((p, n), (vec![1.2, 1.8], 1));
        assert_eq!(load(&path).unwrap()[&3]["PEPTIDEK"], 1.8);
    }

    #[cfg(unix)]
    #[test]
    fn default_temp_root_is_private_and_rejects_symlinks_or_shared_directories() {
        use std::os::unix::fs::{symlink, PermissionsExt};
        assert!(default_root().starts_with(std::env::temp_dir()));
        let dir = tempfile::tempdir().unwrap();
        let private = dir.path().join("private");
        prepare_private_root(&private).unwrap();
        assert_eq!(
            std::fs::metadata(&private).unwrap().permissions().mode() & 0o777,
            0o700
        );
        let link = dir.path().join("link");
        symlink(&private, &link).unwrap();
        assert!(prepare_private_root(&link).is_err());
        std::fs::set_permissions(&private, std::fs::Permissions::from_mode(0o755)).unwrap();
        assert!(prepare_private_root(&private).is_err());
    }
}
