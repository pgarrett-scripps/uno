//! Diagnostic compiled against the frozen UNO source; does not change the engine.
#[path = "../src/koth_ingest.rs"]
mod koth_ingest;
use std::collections::BTreeMap;
use std::path::Path;
use uno_core::{Feature, FeatureMatrix};

fn key(f: &Feature) -> String {
    format!("{:.6}|{:.6}|{}|{:.6}|{}|{:.6}|{}|{}|{:.6}|{:.6}|{:.6}|{}|{:.5e}",
        f.mass, f.mz, f.charge, f.rt, f.im.map(|x| format!("{x:.6}")).unwrap_or_default(),
        f.combined_score, f.n_isotopes.unwrap_or(0.), f.n_scans.unwrap_or(0.),
        f.cosine_score.unwrap_or(0.), f.isotope_score.unwrap_or(0.),
        f.ppm_error.unwrap_or(0.), f.neutron_offset.unwrap_or(0.), f.intensities[0])
}

fn multiset(m: &FeatureMatrix) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for f in &m.features { *counts.entry(key(f)).or_default() += 1; }
    counts
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    anyhow::ensure!(args.len() == 5, "route_audit raw features config output-dir");
    let config = toml::from_str(&std::fs::read_to_string(&args[3])?)?;
    let raw = koth_ingest::features_from_raw(Path::new(&args[1]), &config)?;
    let file = FeatureMatrix::load(&args[2])?;
    let a = multiset(&raw); let b = multiset(&file);
    let common: usize = a.iter().map(|(k,v)| v.min(b.get(k).unwrap_or(&0))).sum();
    let equal_positions = raw.features.iter().zip(&file.features).filter(|(a,b)| key(a)==key(b)).count();
    let out = Path::new(&args[4]); std::fs::create_dir_all(out)?;
    let report = serde_json::json!({"raw_features":raw.len(), "tsv_features":file.len(),
        "matched_after_tsv_precision":common, "same_row_position":equal_positions,
        "comparison":"all UNO-consumed feature fields, rounded to KOTH TSV precision"});
    std::fs::write(out.join("feature_comparison.json"), serde_json::to_string_pretty(&report)?)?;
    let mut w = csv::WriterBuilder::new().delimiter(b'\t').from_path(out.join("raw_full_precision.tsv"))?;
    w.write_record(["massCalib","mz","charge","rtApex","im","combined_score","nIsotopes",
        "nScans","cosine_score","isotope_score","ppm_error","neutron_offset","intensityApex"])?;
    for f in &raw.features {
        w.write_record([f.mass.to_string(),f.mz.to_string(),f.charge.to_string(),f.rt.to_string(),
            f.im.map(|x|x.to_string()).unwrap_or_default(),f.combined_score.to_string(),
            f.n_isotopes.unwrap().to_string(),f.n_scans.unwrap().to_string(),
            f.cosine_score.unwrap().to_string(),f.isotope_score.unwrap().to_string(),
            f.ppm_error.unwrap().to_string(),f.neutron_offset.unwrap().to_string(),f.intensities[0].to_string()])?;
    }
    w.flush()?;
    println!("{report}");
    Ok(())
}
