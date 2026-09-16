//! Opt-in, single-run evidence for uno-quant. No peptide q-value is invented.
use anyhow::{ensure, Result};
use std::path::Path;
use uno_core::score::PeptideHit;
use uno_core::{FeatureMatrix, PeptideDb};

pub fn write(
    path: &Path,
    matrix: &FeatureMatrix,
    db: &PeptideDb,
    hits: &[PeptideHit],
) -> Result<()> {
    ensure!(
        matrix.run_names.len() == 1,
        "quantitative evidence export requires a single run"
    );
    let file = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)?;
    let mut out = csv::WriterBuilder::new().delimiter(b'\t').from_writer(file);
    out.write_record([
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
    ])?;
    for h in hits {
        let p = &db.peptides[h.peptide_idx as usize];
        let f = &matrix.features[h.feature_idx as usize];
        out.write_record([
            "1".to_owned(),
            matrix.run_names[0].clone(),
            h.feature_idx.to_string(),
            h.peptide_idx.to_string(),
            p.sequence.to_string(),
            (h.decoy as u8).to_string(),
            p.proteins
                .iter()
                .map(|p| p.as_ref())
                .collect::<Vec<_>>()
                .join(";"),
            f.charge.to_string(),
            f.im.map(|v| v.to_string()).unwrap_or_default(),
            f.mass.to_string(),
            f.rt.to_string(),
            f.intensities[0].to_string(),
            "uno_ms1_rt_mass_assignment".into(),
        ])?;
    }
    out.flush()?;
    Ok(())
}
