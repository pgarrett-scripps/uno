//! Uno-specific admission policy; deliberately distinct from native-reference mode.
use crate::{
    input::{InputHash, Reader},
    Config, Dataset, InputPeptide,
};
use anyhow::{ensure, Context, Result};
use std::{
    collections::{BTreeMap, BTreeSet},
    path::{Path, PathBuf},
};
fn relative(base: &Path, value: &str) -> PathBuf {
    let p = Path::new(value);
    if p.is_absolute() {
        p.into()
    } else {
        base.join(p)
    }
}
pub fn load(
    manifest: &Path,
    cfg: &Config,
    max_protein_qvalue: f64,
) -> Result<(Dataset, Vec<InputHash>)> {
    ensure!(
        max_protein_qvalue.is_finite() && (0. ..=1.).contains(&max_protein_qvalue),
        "invalid protein identification qvalue threshold"
    );
    let mut reader = Reader::default();
    let design = reader.table(manifest, &["sample", "condition", "evidence", "proteins"])?;
    let base = manifest.parent().unwrap_or(Path::new("."));
    ensure!(!design.is_empty(), "empty sample manifest");
    let mut seen = BTreeSet::new();
    let mut files = BTreeSet::new();
    let mut samples = Vec::new();
    // Contrast direction is explicit: condition must be S1 or S2, fold change is S2/S1.
    for condition in ["S1", "S2"] {
        for r in &design {
            ensure!(
                r["condition"] == "S1" || r["condition"] == "S2",
                "condition must be S1 or S2"
            );
            if r["condition"] != condition {
                continue;
            }
            ensure!(
                !r["sample"].is_empty() && seen.insert(r["sample"].clone()),
                "duplicate/empty sample identifier"
            );
            let path = relative(base, &r["evidence"]);
            ensure!(
                files.insert(std::fs::canonicalize(&path)?),
                "an evidence file is assigned to more than one sample"
            );
            samples.push((
                r["sample"].clone(),
                condition,
                path,
                relative(base, &r["proteins"]),
            ));
        }
    }
    let n_s1 = samples.iter().filter(|s| s.1 == "S1").count();
    ensure!(
        n_s1 >= 2 && samples.len() - n_s1 >= 2,
        "at least two runs per condition required"
    );
    let mut admitted = BTreeSet::new();
    for (_, _, _, path) in &samples {
        for row in reader.table(path, &["accession", "decoy", "q_value"])? {
            ensure!(
                row["decoy"] == "0" || row["decoy"] == "1",
                "invalid protein decoy flag"
            );
            let q: f64 = row["q_value"].parse()?;
            ensure!(
                q.is_finite() && (0. ..=1.).contains(&q),
                "invalid protein ID q-value"
            );
            if row["decoy"] == "0" && q <= max_protein_qvalue {
                let accession = &row["accession"];
                ensure!(
                    !accession.starts_with(&cfg.decoy_prefix),
                    "target accession uses decoy prefix"
                );
                admitted.insert(accession.clone());
                admitted.insert(format!("{}{accession}", cfg.decoy_prefix));
            }
        }
    }
    ensure!(!admitted.is_empty(), "no admitted target proteins");
    let mut merged: BTreeMap<String, InputPeptide> = BTreeMap::new();
    for (sample, (_, _, path, _)) in samples.iter().enumerate() {
        let table = reader.table(
            path,
            &[
                "schema_version",
                "run",
                "feature_id",
                "peptide_id",
                "sequence",
                "decoy",
                "proteins",
                "charge",
                "ion_mobility",
                "observed_mass",
                "rt_minutes",
                "intensity",
                "identity_evidence",
            ],
        )?;
        let mut run_names = BTreeSet::new();
        let mut observations = BTreeSet::new();
        let mut best: BTreeMap<String, (f64, InputPeptide)> = BTreeMap::new();
        for row in table {
            ensure!(
                row["schema_version"] == "1"
                    && row["identity_evidence"] == "uno_ms1_rt_mass_assignment",
                "unsupported evidence schema or identity type"
            );
            run_names.insert(row["run"].clone());
            ensure!(
                run_names.len() == 1,
                "evidence contains multiple runs; matrix input is unsupported by this adapter"
            );
            ensure!(
                observations.insert((
                    row["feature_id"].clone(),
                    row["peptide_id"].clone(),
                    row["decoy"].clone()
                )),
                "duplicate observation identifier"
            );
            ensure!(
                row["decoy"] == "0" || row["decoy"] == "1",
                "invalid peptide decoy flag"
            );
            let decoy = row["decoy"] == "1";
            // Uno strips decoy tags from accessions and uses (decoy, accession)
            // as the identity. The authoritative flag must define the null namespace.
            let mut proteins = Vec::new();
            for accession in row["proteins"].split(';') {
                ensure!(!accession.is_empty(), "empty protein accession");
                let name = if decoy {
                    format!(
                        "{}{}",
                        cfg.decoy_prefix,
                        accession
                            .strip_prefix(&cfg.decoy_prefix)
                            .unwrap_or(accession)
                    )
                } else {
                    ensure!(
                        !accession.starts_with(&cfg.decoy_prefix),
                        "target accession collides with the quantitative decoy namespace"
                    );
                    accession.to_owned()
                };
                if admitted.contains(&name) {
                    proteins.push(name);
                }
            }
            let sequence = &row["sequence"];
            if proteins.is_empty() || sequence.len() < cfg.min_length {
                continue;
            }
            let charge: u8 = row["charge"].parse().context("invalid charge")?;
            ensure!(charge > 0, "zero charge");
            let im: Option<f64> = if row["ion_mobility"].is_empty() {
                None
            } else {
                Some(row["ion_mobility"].parse()?)
            };
            ensure!(im.is_none_or(|x| x.is_finite()), "invalid mobility");
            let intensity: f64 = row["intensity"].parse()?;
            ensure!(
                intensity.is_finite() && intensity >= 0.,
                "invalid intensity"
            );
            if intensity == 0. {
                continue;
            }
            let key = format!(
                "{}|{}|{}|{}",
                sequence,
                charge,
                im.map(|v| v.to_string()).unwrap_or_default(),
                u8::from(decoy)
            );
            let p = InputPeptide {
                key: key.clone(),
                sequence: sequence.clone(),
                proteins,
                values: vec![None; samples.len()],
            };
            if best.get(&key).is_none_or(|(old, _)| intensity > *old) {
                best.insert(key, (intensity, p));
            }
        }
        for (key, (value, row)) in best {
            let dest = merged.entry(key).or_insert(row);
            dest.values[sample] = Some(value);
        }
    }
    Ok((
        Dataset {
            samples: samples.into_iter().map(|s| s.0).collect(),
            n_s1,
            rows: merged.into_values().collect(),
        },
        reader.hashes,
    ))
}
