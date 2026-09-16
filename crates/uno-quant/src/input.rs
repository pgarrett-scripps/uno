use crate::{Config, Dataset, InputPeptide};
use anyhow::{bail, ensure, Context, Result};
use serde::Serialize;
use sha2::{Digest, Sha256};
use std::{
    collections::{BTreeMap, BTreeSet},
    path::{Path, PathBuf},
};

#[derive(Clone, Debug, Serialize)]
pub struct InputHash {
    pub path: PathBuf,
    pub sha256: String,
}
#[derive(Default)]
pub struct Reader {
    pub hashes: Vec<InputHash>,
}
pub type Row = BTreeMap<String, String>;
impl Reader {
    pub fn table(&mut self, path: &Path, required: &[&str]) -> Result<Vec<Row>> {
        let bytes = std::fs::read(path).with_context(|| format!("reading {}", path.display()))?;
        self.hashes.push(InputHash {
            path: path.to_path_buf(),
            sha256: format!("{:x}", Sha256::digest(&bytes)),
        });
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(bytes.as_slice());
        let headers = reader.headers()?.clone();
        ensure!(
            headers.iter().collect::<BTreeSet<_>>().len() == headers.len(),
            "duplicate columns in {}",
            path.display()
        );
        for name in required {
            ensure!(
                headers.iter().any(|h| h == *name),
                "missing column {name} in {}",
                path.display()
            );
        }
        reader
            .records()
            .map(|r| {
                Ok(headers
                    .iter()
                    .map(str::to_owned)
                    .zip(r?.iter().map(str::to_owned))
                    .collect())
            })
            .collect()
    }
    pub fn lines(&mut self, path: &Path) -> Result<BTreeSet<String>> {
        let bytes = std::fs::read(path)?;
        self.hashes.push(InputHash {
            path: path.to_path_buf(),
            sha256: format!("{:x}", Sha256::digest(&bytes)),
        });
        Ok(std::str::from_utf8(&bytes)?
            .lines()
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .map(str::to_owned)
            .collect())
    }
}
fn companion(path: &Path, suffix: &str) -> Result<PathBuf> {
    let name = path.to_str().context("non-UTF8 reference path")?;
    let base = name
        .strip_suffix("_proteins_full.tsv")
        .context("reference inputs must end in _proteins_full.tsv")?;
    Ok(PathBuf::from(format!("{base}{suffix}")))
}
fn float(row: &Row, key: &str) -> Result<f64> {
    row[key]
        .parse()
        .with_context(|| format!("invalid {key}: {}", row[key]))
}
fn py_float(value: f64) -> String {
    if value.is_nan() {
        "nan".to_string()
    } else if value.fract() == 0. {
        format!("{value:.1}")
    } else {
        value.to_string()
    }
}
#[derive(Default)]
pub struct ReferenceOptions {
    pub all_proteins: bool,
    pub all_pfms: bool,
    pub allowed_proteins: Option<PathBuf>,
    pub allowed_peptides: Option<PathBuf>,
}
pub fn reference(
    s1: &[PathBuf],
    s2: &[PathBuf],
    cfg: &Config,
    options: &ReferenceOptions,
) -> Result<(Dataset, Vec<InputHash>)> {
    ensure!(
        s1.len() >= 2 && s2.len() >= 2,
        "at least two input files per condition required"
    );
    let paths: Vec<_> = s1.iter().chain(s2).collect();
    let canonical: Result<BTreeSet<_>> = paths
        .iter()
        .map(|p| Ok(std::fs::canonicalize(p)?))
        .collect();
    ensure!(
        canonical?.len() == paths.len(),
        "a reference file occurs in multiple samples or conditions"
    );
    let mut reader = Reader::default();
    let mut allowed_proteins = BTreeSet::new();
    let mut allowed_peptides = BTreeSet::new();
    for path in &paths {
        if options.allowed_proteins.is_none() {
            let p = if options.all_proteins {
                (*path).clone()
            } else {
                companion(path, "_proteins.tsv")?
            };
            for row in reader.table(&p, &["dbname"])? {
                let name = row["dbname"].clone();
                if !options.all_proteins {
                    allowed_proteins.insert(format!("{}{name}", cfg.decoy_prefix));
                }
                allowed_proteins.insert(name);
            }
        }
        if options.allowed_peptides.is_none() {
            for row in reader.table(
                &companion(path, "_PFMs_ML.tsv")?,
                &["seqs", "qpreds", "plen"],
            )? {
                if (options.all_pfms || float(&row, "qpreds")? <= 10.)
                    && float(&row, "plen")? >= cfg.min_length as f64
                {
                    allowed_peptides.insert(row["seqs"].clone());
                }
            }
        }
    }
    if let Some(p) = &options.allowed_proteins {
        allowed_proteins = reader.lines(p)?;
        for p in allowed_proteins.clone() {
            allowed_proteins.insert(format!("{}{p}", cfg.decoy_prefix));
        }
    }
    if let Some(p) = &options.allowed_peptides {
        allowed_peptides = reader.lines(p)?;
    }
    let mut tables = Vec::new();
    let mut expanded = BTreeSet::new();
    for path in &paths {
        let rows = reader.table(
            &companion(path, "_PFMs.tsv")?,
            &[
                "sequence",
                "proteins",
                "charge",
                "ion_mobility",
                "Intensity",
                "RT",
            ],
        )?;
        for row in &rows {
            if allowed_peptides.contains(&row["sequence"])
                && row["proteins"]
                    .split(';')
                    .any(|p| allowed_proteins.contains(p))
            {
                expanded.extend(row["proteins"].split(';').map(str::to_owned));
            }
        }
        tables.push(rows);
    }
    let mut merged: BTreeMap<String, InputPeptide> = BTreeMap::new();
    for (sample, table) in tables.iter().enumerate() {
        let mut best: BTreeMap<String, (f64, InputPeptide)> = BTreeMap::new();
        for row in table {
            let seq = &row["sequence"];
            if !allowed_peptides.contains(seq) {
                continue;
            }
            let proteins: Vec<_> = row["proteins"]
                .split(';')
                .filter(|p| expanded.contains(*p))
                .map(str::to_owned)
                .collect();
            if proteins.is_empty() {
                continue;
            }
            let intensity = float(row, "Intensity")?;
            ensure!(intensity.is_finite() && intensity>0.,"nonpositive/nonfinite reference intensity for {seq}; supply missing rows rather than zeros");
            let charge: u8 = row["charge"].parse().context("invalid precursor charge")?;
            ensure!(charge > 0, "charge must be positive");
            // pandas.read_table represents absent native mobility as NaN; its
            // string "nan" is part of the reference precursor key.
            let im = if row["ion_mobility"].trim().is_empty() {
                f64::NAN
            } else {
                float(row, "ion_mobility")?
            };
            ensure!(!im.is_infinite(), "invalid mobility");
            let key = format!("{seq}{charge}{}", py_float(im));
            let entry = InputPeptide {
                key: key.clone(),
                sequence: seq.clone(),
                proteins,
                values: vec![None; paths.len()],
            };
            if let Some((v, old)) = best.get(&key) {
                ensure!(
                    old.sequence == *seq,
                    "ambiguous reference composite precursor identity"
                );
                if intensity <= *v {
                    continue;
                }
            }
            best.insert(key, (intensity, entry));
        }
        for (key, (intensity, row)) in best {
            let p = merged.entry(key).or_insert(row);
            p.values[sample] = Some(intensity);
        }
    }
    let samples = paths
        .iter()
        .map(|p| p.file_name().unwrap().to_string_lossy().to_string())
        .collect();
    if merged.is_empty() {
        bail!("no eligible peptide-feature matches");
    }
    Ok((
        Dataset {
            samples,
            n_s1: s1.len(),
            rows: merged.into_values().collect(),
        },
        reader.hashes,
    ))
}
