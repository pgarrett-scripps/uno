use crate::{input::InputHash, Analysis};
use anyhow::{Context, Result};
use serde_json::json;
use sha2::{Digest, Sha256};
use std::path::Path;
fn writer(path: &Path) -> Result<csv::Writer<std::fs::File>> {
    Ok(csv::WriterBuilder::new().delimiter(b'\t').from_path(path)?)
}
fn tf(v: bool) -> String {
    if v { "True" } else { "False" }.into()
}
pub fn write(out: &Path, a: &Analysis, inputs: &[InputHash], mode: &str) -> Result<()> {
    write_with_admission(out, a, inputs, mode, None)
}
pub fn write_with_admission(
    out: &Path,
    a: &Analysis,
    inputs: &[InputHash],
    mode: &str,
    admission: Option<&crate::uno_adapter::AdmissionAudit>,
) -> Result<()> {
    // Exclusive directory creation prevents accidental replacement of a previous analysis.
    std::fs::create_dir(out)
        .with_context(|| format!("create a NEW output directory {}", out.display()))?;
    let mut p = writer(&out.join("peptides.tsv"))?;
    let mut headers = vec![
        "proteins",
        "origseq",
        "peptide",
        "decoy",
        "nummissing",
        "nonmissing_S1",
        "nonmissing_S2",
        "intensity_median",
        "iq",
        "FC_raw",
        "FC",
        "FC_corrected",
        "p-value",
        "FC_gr_mean",
        "bp",
        "sign",
        "up",
        "down",
    ]
    .into_iter()
    .map(str::to_owned)
    .collect::<Vec<_>>();
    headers.extend(a.samples.clone());
    p.write_record(headers)?;
    for x in &a.peptides {
        let mut row = vec![
            x.protein.clone(),
            x.sequence.clone(),
            x.peptide.clone(),
            tf(x.decoy),
            x.missing.to_string(),
            x.observed_s1.to_string(),
            x.observed_s2.to_string(),
            x.intensity.to_string(),
            x.bin.to_string(),
            x.raw_fc.to_string(),
            x.fc.to_string(),
            x.corrected_fc.to_string(),
            x.pvalue.to_string(),
            x.running_fc.to_string(),
            x.better.to_string(),
            tf(x.significant),
            tf(x.up),
            tf(x.down),
        ];
        row.extend(x.intensities.iter().map(ToString::to_string));
        p.write_record(row)?;
    }
    p.flush()?;
    let headers = [
        "score",
        "p-value",
        "dbname",
        "log2FoldChange(S2/S1)",
        "differentially expressed peptides",
        "identified peptides",
        "log2FoldChange(S2/S1) no normalization",
        "log2FoldChange(S2/S1) using all peptides",
        "log2FoldChange(S2/S1) using all peptides and no normalization",
        "protein_quant_group",
        "FC_pass",
        "FDR_pass",
        "BH_pass",
    ];
    let mut full = writer(&out.join("proteins.tsv"))?;
    full.write_record(headers)?;
    let mut significant = writer(&out.join("significant.tsv"))?;
    significant.write_record(headers)?;
    for x in &a.proteins {
        let row = [
            x.score.to_string(),
            x.pvalue.to_string(),
            x.name.clone(),
            x.fc.to_string(),
            x.n_changed.to_string(),
            x.n_peptides.to_string(),
            x.raw_fc.to_string(),
            x.all_fc.to_string(),
            x.all_raw_fc.to_string(),
            x.group.to_string(),
            tf(x.fc_pass),
            x.fdr_pass.map(tf).unwrap_or_default(),
            tf(x.bh_pass),
        ];
        full.write_record(&row)?;
        if x.bh_pass && x.fc_pass {
            significant.write_record(&row)?;
        }
    }
    full.flush()?;
    significant.flush()?;
    let mut files = std::collections::BTreeMap::new();
    for name in ["peptides.tsv", "proteins.tsv", "significant.tsv"] {
        files.insert(
            name,
            format!("{:x}", Sha256::digest(std::fs::read(out.join(name))?)),
        );
    }
    let manifest = json!({"schema_version":1,"admission_audit":admission,"executable_sha256":format!("{:x}",Sha256::digest(std::fs::read(std::env::current_exe()?)?)),"peptide_sample_columns":"normalized intensities with reference minimum-intensity filling; not raw observed measurements","status":"experimental_not_biologically_validated","mode":mode,"reference_commit":"9df4a941a901aac657b1c9ab5229cb49671e3c7d","reference_sha256":"cb15b7c43b276eb45e77da2b8cf7cd97cd20765f95264202283212fde85cfc90","crate_version":env!("CARGO_PKG_VERSION"),"arguments":std::env::args().collect::<Vec<_>>(),"config":a.config,"samples":a.samples,"condition_1_count":a.n_s1,"inputs":inputs,"outputs":files,"normalization_factors":a.normalization_factors,"calibration_mean":a.calibration_mean,"calibration_std":a.calibration_std,"fold_change_threshold":a.fold_change_threshold,"decoy_probability":a.decoy_probability,"decoy_sequences":a.n_decoy_sequences,"input_rows":a.input_rows,"eligible_precursors":a.eligible_precursors,"protein_rows":a.proteins.len(),"significant_proteins":a.proteins.iter().filter(|p|p.bh_pass&&p.fc_pass).count(),"semantics":"FDR_pass and BH_pass reproduce reference decisions, not newly validated q-values. Scalar NumPy 1.26 quicksort tie ordering; reference CPU SIMD sort must be disabled. Group numbers are opaque. Uno mode uses protein-admitted RT/mass assignments, not reference PFM-ML admission."});
    std::fs::write(
        out.join("manifest.json"),
        serde_json::to_string_pretty(&manifest)?,
    )?;
    Ok(())
}
