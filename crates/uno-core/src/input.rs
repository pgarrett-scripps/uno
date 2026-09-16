//! Read single-run biosaur2/KOTH feature TSVs or aligned KOTH intensity matrices.
//!
//! Single-run tables retain feature-finding quality descriptors. Aligned
//! matrices retain cross-run contributions and intensities, but omit several
//! individual-feature descriptors. Optional matrix evidence columns supply
//! validated summaries of KOTH's LFQ extraction diagnostics. Raw-file feature
//! finding is handled separately by `uno::koth_ingest`.
//!
//! Columns are located BY NAME (order-independent). Required: `massCalib`, `mz`,
//! `charge`, `rtApex`. Everything else is optional and read if present, with
//! biosaur2 name aliases (`nIsotopes`↔`n_isotopes`, `nScans`↔`n_scans`).

use anyhow::{anyhow, Context, Result};
use std::path::Path;

/// Optional, label-independent summaries derived from KOTH's LFQ details.
/// All columns must precede `n_contributing_runs`; partial schemas are rejected.
/// These describe extraction confidence, not peptide identification confidence.
pub const MATRIX_EVIDENCE_COLUMNS: [&str; 9] = [
    "matrix_confident_fraction",
    "matrix_confident_direct_fraction",
    "matrix_confident_mbr_fraction",
    "matrix_median_log_confidence",
    "matrix_median_isotope_pattern",
    "matrix_median_coelution",
    "matrix_median_isotopes",
    "matrix_median_abs_ppm",
    "matrix_median_abs_rt",
];

/// One LC-MS feature from a single run: a deisotoped, charge-resolved precursor
/// mass at an apex RT, that we try to assign a peptide to.
#[derive(Debug, Clone, PartialEq)]
pub struct Feature {
    /// Monoisotopic neutral mass (Da).
    pub mass: f64,
    /// Monoisotopic m/z.
    pub mz: f64,
    pub charge: u8,
    /// Apex retention time (minutes).
    pub rt: f64,
    /// Ion mobility (1/K0); `None` when the data has no IM dimension.
    pub im: Option<f64>,
    /// Feature-quality score (isotope x coelution), in [0, 1]. Defaults to 1.0
    /// when the input has no `combined_score` column.
    pub combined_score: f64,
    /// Number of original runs contributing to a consensus feature; always 1
    /// for single-run input. Does not count MBR-filled intensity cells.
    pub n_contributing_runs: u32,
    /// Per-feature shape descriptors: number of isotopic peaks and number of
    /// scans the feature spans (biosaur2 `nIsotopes`/`nScans`). `None` if absent.
    pub n_isotopes: Option<f64>,
    pub n_scans: Option<f64>,
    /// biosaur2 isotope-pattern fit error (deviation of the observed isotope
    /// envelope from theory) and RT peak width (minutes). Real peptides fit
    /// theory tightly; random mass coincidences (and shuffled decoys) do not.
    pub isoerror: Option<f64>,
    pub isoerror2: Option<f64>,
    pub rt_width: Option<f64>,
    /// Per-feature PFM-quality columns. These are the target/decoy discriminators
    /// the PFM ML rescorer relies on (the same signal ms1searchpy's LightGBM PFM
    /// filter uses to drop the random-match rate). `None` if the column is absent.
    pub cosine_score: Option<f64>,
    pub isotope_score: Option<f64>,
    pub ppm_error: Option<f64>,
    pub neutron_offset: Option<f64>,
    /// One apex intensity for single-run input, or per-run integrated LFQ
    /// intensities for matrix input. The estimators are not interchangeable.
    pub intensities: Vec<f64>,
}

/// A single-run or aligned multi-run feature set.
#[derive(Debug, Clone)]
pub struct FeatureMatrix {
    /// Run intensity column names in vector order.
    pub run_names: Vec<String>,
    pub features: Vec<Feature>,
    /// Optional LFQ evidence, one vector per retained consensus feature.
    /// Absent for ordinary matrices and single-run inputs: their classifier
    /// vectors keep the exact original dimensions and column order.
    pub matrix_evidence: Option<Vec<[f32; 9]>>,
    /// Original, aligned observations per consensus row for one joint search.
    /// Quantitative cells remain on `features`; observations never duplicate them.
    pub joint_evidence: Option<Vec<Vec<crate::joint::Observation>>>,
    /// Extraction summaries per joint row, independent of candidate identity.
    pub joint_cell_evidence: Option<Vec<[f32; 12]>>,
}

impl FeatureMatrix {
    pub fn len(&self) -> usize {
        self.features.len()
    }
    pub fn is_empty(&self) -> bool {
        self.features.is_empty()
    }
    /// Number of runs; 1 for single-run input.
    pub fn n_samples(&self) -> usize {
        self.run_names.len()
    }

    /// Load a single-run biosaur2 feature file (`*.features.tsv`).
    pub fn load(features_path: impl AsRef<Path>) -> Result<Self> {
        let (intensity_name, rows) =
            read_features(features_path.as_ref()).with_context(|| "reading feature file")?;
        let features = rows
            .into_iter()
            .map(MatrixRow::into_feature)
            .collect::<Vec<_>>();
        Ok(FeatureMatrix {
            run_names: vec![intensity_name],
            features,
            matrix_evidence: None,
            joint_evidence: None,
            joint_cell_evidence: None,
        })
    }

    /// Load a koth-aligned consensus matrix (`intensity_matrix.tsv`): one row per
    /// consensus feature, with the leading columns `massCalib, mz, charge, rtApex,
    /// im, combined_score, seed_run, n_contributing_runs` followed by one
    /// intensity column per run. Unlike single-run input, each feature carries the
    /// FULL per-run intensity vector (`intensities`) — that's what the multi-sample
    /// quant rollup consumes, and the ID featurizer already reads `n_detected`
    /// (runs with signal) and `sum_int` from it. The per-feature single-run quality
    /// columns (cosine/isotope/ppm/nScans) are absent here and read as `None`
    /// (inert in the featurizer); `combined_score` and `n_contributing_runs`
    /// survive alignment and are used.
    ///
    /// `min_runs` drops consensus features detected in fewer than that many runs.
    /// This filter is REQUIRED for aligned-matrix search, not cosmetic: alignment
    /// discards the per-feature `nIsotopes` column that uno's single-run path uses
    /// to reject junk, so without a reproducibility floor the random-match
    /// background `p` climbs to ~0.83 and the binomial finds nothing. A real
    /// peptide is detected across replicates; single-run features are dominated by
    /// noise. `min_runs = 2` is a safe floor (won't drop a condition-specific
    /// peptide that has ≥2 replicates).
    pub fn load_matrix(path: impl AsRef<Path>, min_runs: u32) -> Result<Self> {
        let path = path.as_ref();
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .flexible(false)
            .from_path(path)
            .with_context(|| format!("opening {}", path.display()))?;
        let header = rdr.headers()?.clone();
        let col = |name: &str| -> Result<usize> {
            header
                .iter()
                .position(|h| h == name)
                .ok_or_else(|| anyhow!("aligned matrix missing column `{name}`"))
        };
        let (c_mass, c_mz, c_charge, c_rt) =
            (col("massCalib")?, col("mz")?, col("charge")?, col("rtApex")?);
        let c_im = header.iter().position(|h| h == "im");
        let c_cs = header.iter().position(|h| h == "combined_score");
        let c_ncr = col("n_contributing_runs")?;
        let evidence_cols: Vec<Option<usize>> = MATRIX_EVIDENCE_COLUMNS.iter()
            .map(|name| header.iter().position(|h| h == *name)).collect();
        let has_evidence = evidence_cols.iter().any(Option::is_some);
        if has_evidence && evidence_cols.iter().any(|c| c.is_none_or(|c| c >= c_ncr)) {
            return Err(anyhow!("matrix evidence requires all nine named columns before n_contributing_runs"));
        }
        let mut matrix_evidence = has_evidence.then(Vec::new);
        // Run intensity columns are everything after `n_contributing_runs`.
        let first_run = c_ncr + 1;
        if first_run >= header.len() {
            return Err(anyhow!("aligned matrix has no run intensity columns"));
        }
        let run_names: Vec<String> = header.iter().skip(first_run).map(|s| s.to_string()).collect();

        let parse = |rec: &csv::StringRecord, c: usize, what: &str, line: usize| -> Result<f64> {
            rec.get(c)
                .unwrap_or("")
                .parse::<f64>()
                .map_err(|_| anyhow!("line {line}: cannot parse {what}"))
        };
        let mut features = Vec::new();
        for (n, rec) in rdr.records().enumerate() {
            let line = n + 2;
            let rec = rec.with_context(|| format!("reading record at line {line}"))?;
            let mut evidence = [0.0f32; 9];
            if has_evidence {
                for (i, c) in evidence_cols.iter().enumerate() {
                    let value = parse(&rec, c.unwrap(), MATRIX_EVIDENCE_COLUMNS[i], line)?;
                    if !value.is_finite() || value < 0.0 || value > f32::MAX as f64
                        || (matches!(i, 0 | 1 | 2 | 4 | 5) && value > 1.0) {
                        return Err(anyhow!("line {line}: invalid matrix evidence {}", MATRIX_EVIDENCE_COLUMNS[i]));
                    }
                    evidence[i] = value as f32;
                }
            }
            let intensities: Vec<f64> = (first_run..header.len())
                .map(|c| rec.get(c).unwrap_or("0").trim().parse::<f64>().unwrap_or(0.0))
                .collect();
            // Reproducibility floor uses koth's `n_contributing_runs` — the number
            // of runs that DETECTED the feature, NOT the count of non-zero cells.
            // With MBR, nearly every cell is filled, so counting non-zero cells
            // would pass ~all junk (defeating the floor); the detected count still
            // discriminates. Falls back to the non-zero count if the column is
            // absent (e.g. an older matrix without it).
            let n_detected = rec
                .get(c_ncr)
                .and_then(|s| s.trim().parse::<u32>().ok())
                .unwrap_or_else(|| intensities.iter().filter(|&&x| x > 0.0).count() as u32);
            if n_detected < min_runs {
                continue;
            }
            if let Some(values) = &mut matrix_evidence {
                values.push(evidence);
            }
            let im = c_im
                .and_then(|c| rec.get(c))
                .map(str::trim)
                .filter(|s| !s.is_empty())
                .and_then(|s| s.parse::<f64>().ok());
            let combined_score = c_cs
                .and_then(|c| rec.get(c))
                .and_then(|s| s.trim().parse::<f64>().ok())
                .unwrap_or(1.0);
            features.push(Feature {
                mass: parse(&rec, c_mass, "massCalib", line)?,
                mz: parse(&rec, c_mz, "mz", line)?,
                charge: parse(&rec, c_charge, "charge", line)? as u8,
                rt: parse(&rec, c_rt, "rtApex", line)?,
                im,
                combined_score,
                n_contributing_runs: n_detected.max(1),
                n_isotopes: None,
                n_scans: None,
                isoerror: None,
                isoerror2: None,
                rt_width: None,
                cosine_score: None,
                isotope_score: None,
                ppm_error: None,
                neutron_offset: None,
                intensities,
            });
        }
        Ok(FeatureMatrix { run_names, features, matrix_evidence, joint_evidence: None, joint_cell_evidence: None })
    }
}

/// A parsed feature row.
struct MatrixRow {
    mass: f64,
    mz: f64,
    charge: u8,
    rt: f64,
    im: Option<f64>,
    combined_score: f64,
    n_isotopes: Option<f64>,
    n_scans: Option<f64>,
    isoerror: Option<f64>,
    isoerror2: Option<f64>,
    rt_width: Option<f64>,
    cosine_score: Option<f64>,
    isotope_score: Option<f64>,
    ppm_error: Option<f64>,
    neutron_offset: Option<f64>,
    intensity: f64,
}

impl MatrixRow {
    fn into_feature(self) -> Feature {
        Feature {
            mass: self.mass,
            mz: self.mz,
            charge: self.charge,
            rt: self.rt,
            im: self.im,
            combined_score: self.combined_score,
            n_contributing_runs: 1,
            n_isotopes: self.n_isotopes,
            n_scans: self.n_scans,
            isoerror: self.isoerror,
            isoerror2: self.isoerror2,
            rt_width: self.rt_width,
            cosine_score: self.cosine_score,
            isotope_score: self.isotope_score,
            ppm_error: self.ppm_error,
            neutron_offset: self.neutron_offset,
            intensities: vec![self.intensity],
        }
    }
}

/// Column indices within a single-run feature-file header, located by name.
struct Layout {
    mass: usize,
    mz: usize,
    charge: usize,
    rt: usize,
    im: Option<usize>,
    combined_score: Option<usize>,
    n_isotopes: Option<usize>,
    n_scans: Option<usize>,
    isoerror: Option<usize>,
    isoerror2: Option<usize>,
    rt_width: Option<usize>,
    cosine_score: Option<usize>,
    isotope_score: Option<usize>,
    ppm_error: Option<usize>,
    neutron_offset: Option<usize>,
    intensity: Option<usize>,
    intensity_name: String,
}

impl Layout {
    fn from_header(header: &csv::StringRecord) -> Result<Self> {
        let find = |name: &str| -> Result<usize> {
            header
                .iter()
                .position(|h| h == name)
                .ok_or_else(|| anyhow!("missing required column `{name}` in feature-file header"))
        };
        // Find the first present of a set of alias names.
        let find_any = |names: &[&str]| -> Option<usize> {
            names.iter().find_map(|&n| header.iter().position(|h| h == n))
        };

        // Intensity: prefer apex, fall back to sum.
        let (intensity, intensity_name) = match find_any(&["intensityApex", "intensitySum", "intensity"]) {
            Some(i) => (Some(i), header.get(i).unwrap_or("intensity").to_string()),
            None => (None, "intensity".to_string()),
        };

        Ok(Layout {
            mass: find("massCalib")?,
            mz: find("mz")?,
            charge: find("charge")?,
            rt: find("rtApex")?,
            im: find_any(&["im", "ion_mobility", "FAIMS"]),
            combined_score: find_any(&["combined_score"]),
            n_isotopes: find_any(&["n_isotopes", "nIsotopes"]),
            n_scans: find_any(&["n_scans", "nScans"]),
            isoerror: find_any(&["isoerror"]),
            isoerror2: find_any(&["isoerror2"]),
            rt_width: find_any(&["rt_width"]),
            cosine_score: find_any(&["cosine_score"]),
            isotope_score: find_any(&["isotope_score"]),
            ppm_error: find_any(&["ppm_error"]),
            neutron_offset: find_any(&["neutron_offset"]),
            intensity,
            intensity_name,
        })
    }
}

/// Read a single-run feature TSV. Returns (intensity_column_name, rows).
fn read_features(path: &Path) -> Result<(String, Vec<MatrixRow>)> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .flexible(false)
        .from_path(path)
        .with_context(|| format!("opening {}", path.display()))?;

    let header = rdr.headers()?.clone();
    let layout = Layout::from_header(&header)?;

    let parse_f64 = |rec: &csv::StringRecord, col: usize, what: &str, line: usize| -> Result<f64> {
        let s = rec.get(col).unwrap_or("");
        s.parse::<f64>()
            .map_err(|_| anyhow!("line {line}: cannot parse {what} `{s}` as number"))
    };

    let mut rows = Vec::new();
    for (n, rec) in rdr.records().enumerate() {
        let line = n + 2; // 1-based, +1 for header
        let rec = rec.with_context(|| format!("reading record at line {line}"))?;

        let opt_f64 = |col: Option<usize>| -> Option<f64> {
            col.and_then(|c| rec.get(c))
                .map(str::trim)
                .filter(|s| !s.is_empty())
                .and_then(|s| s.parse::<f64>().ok())
        };

        let charge = rec
            .get(layout.charge)
            .unwrap_or("")
            .parse::<u8>()
            .map_err(|_| anyhow!("line {line}: cannot parse charge"))?;

        rows.push(MatrixRow {
            mass: parse_f64(&rec, layout.mass, "massCalib", line)?,
            mz: parse_f64(&rec, layout.mz, "mz", line)?,
            charge,
            rt: parse_f64(&rec, layout.rt, "rtApex", line)?,
            im: opt_f64(layout.im),
            // No combined_score column → treat as 1.0 (neutral quality).
            combined_score: opt_f64(layout.combined_score).unwrap_or(1.0),
            n_isotopes: opt_f64(layout.n_isotopes),
            n_scans: opt_f64(layout.n_scans),
            isoerror: opt_f64(layout.isoerror),
            isoerror2: opt_f64(layout.isoerror2),
            rt_width: opt_f64(layout.rt_width),
            cosine_score: opt_f64(layout.cosine_score),
            isotope_score: opt_f64(layout.isotope_score),
            ppm_error: opt_f64(layout.ppm_error),
            neutron_offset: opt_f64(layout.neutron_offset),
            intensity: opt_f64(layout.intensity).unwrap_or(0.0),
        });
    }

    Ok((layout.intensity_name, rows))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn matrix_fixture(evidence: bool, invalid: bool) -> std::path::PathBuf {
        let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
        let path = std::env::temp_dir().join(format!("uno-matrix-{}-{stamp}.tsv", std::process::id()));
        let extra = if evidence { format!("{}\t", MATRIX_EVIDENCE_COLUMNS.join("\t")) } else { String::new() };
        let values = if evidence {
            if invalid { "NaN\t0\t0\t0\t0\t0\t0\t0\t0\t" } else { "0.5\t0.25\t0.25\t2\t0.8\t0.9\t3\t1.2\t0.1\t" }
        } else { "" };
        std::fs::write(&path, format!(
            "massCalib\tmz\tcharge\trtApex\t{extra}n_contributing_runs\tA\tB\n\
             1000\t501\t2\t10\t{values}1\t100\t0\n\
             1100\t551\t2\t12\t{values}2\t200\t300\n"
        )).unwrap();
        path
    }

    #[test]
    fn ordinary_matrix_keeps_original_classifier_schema() {
        let path = matrix_fixture(false, false);
        let m = FeatureMatrix::load_matrix(&path, 2).unwrap();
        std::fs::remove_file(path).unwrap();
        assert_eq!(m.run_names, ["A", "B"]);
        assert_eq!(m.len(), 1);
        assert_eq!(m.features[0].intensities, [200.0, 300.0]);
        assert!(m.matrix_evidence.is_none());
    }

    #[test]
    fn matrix_evidence_tracks_filtered_feature_rows() {
        let path = matrix_fixture(true, false);
        let m = FeatureMatrix::load_matrix(&path, 2).unwrap();
        std::fs::remove_file(path).unwrap();
        assert_eq!(m.run_names, ["A", "B"]);
        assert_eq!(m.features[0].mass, 1100.0);
        let values = m.matrix_evidence.unwrap();
        assert_eq!(values.len(), 1);
        assert_eq!(values[0], [0.5, 0.25, 0.25, 2.0, 0.8, 0.9, 3.0, 1.2, 0.1]);
    }

    #[test]
    fn matrix_evidence_rejects_nonfinite_values_even_in_filtered_rows() {
        let path = matrix_fixture(true, true);
        let error = FeatureMatrix::load_matrix(&path, 2).unwrap_err();
        std::fs::remove_file(path).unwrap();
        assert!(error.to_string().contains("invalid matrix evidence"));
    }

    #[test]
    fn matrix_evidence_rejects_partial_schema() {
        let path = matrix_fixture(true, false);
        let text = std::fs::read_to_string(&path).unwrap().replace("matrix_median_abs_rt", "unexpected");
        std::fs::write(&path, text).unwrap();
        let error = FeatureMatrix::load_matrix(&path, 2).unwrap_err();
        std::fs::remove_file(path).unwrap();
        assert!(error.to_string().contains("all nine named columns"));
    }

    #[test]
    fn loads_single_run_features() {
        let Some(path) = crate::test_input("tests/data/09744.features.tsv") else {
            return;
        };
        let m = FeatureMatrix::load(&path).expect("load features");
        assert_eq!(m.n_samples(), 1, "single-run input");
        assert!(m.len() > 1000, "expected many features, got {}", m.len());

        let f = &m.features[0];
        assert!(f.charge >= 1);
        assert!(f.mass > 0.0);
        assert_eq!(f.intensities.len(), 1);
        // The quality columns that separate target from decoy must be present on
        // single-run input (the reason we do not use aligned matrices).
        assert!(
            f.cosine_score.is_some() || f.isotope_score.is_some(),
            "single-run features must carry PFM-quality columns"
        );
    }
}
