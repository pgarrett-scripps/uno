//! Input contract: read the aligned LFQ matrix produced by `koth_align`.
//!
//! `uno` does not read raw data or find features. Its input is koth's
//! consensus/aligned output. Three sibling files share the same row order and
//! the same leading metadata columns:
//!   - `intensity_matrix.tsv` — metadata + one integrated-intensity column per sample
//!   - `qvalue_matrix.tsv`    — metadata + one target-decoy q-value column per sample
//!   - `consensus_features.tsv` — metadata only (+ `n_runs_detected`)
//!
//! Metadata columns (koth): massCalib, mz, charge, rtApex, im, combined_score,
//! seed_run, n_contributing_runs[, n_runs_detected]. `im` is empty for data
//! without ion mobility (e.g. Orbitrap).

use anyhow::{anyhow, bail, Context, Result};
use std::path::Path;

/// Columns koth writes as feature metadata. Everything else in a matrix file is
/// treated as a per-sample value column, in file order.
const META_COLUMNS: &[&str] = &[
    "massCalib",
    "mz",
    "charge",
    "rtApex",
    "im",
    "combined_score",
    "seed_run",
    "n_contributing_runs",
    "n_runs_detected",
    // Optional per-feature shape descriptors (e.g. from biosaur2). Present for
    // single-run searches; absent on koth aligned matrices. Used by the learned
    // discriminant to downweight junk features (which contribute random
    // target+decoy matches) relative to well-formed ones.
    "n_isotopes",
    "n_scans",
    // Optional biosaur2 discriminators (single-run): isotope-pattern fit error and
    // RT peak width. These separate real peptides from random mass-coincidence
    // (incl. shuffled decoys) where n_isotopes/n_scans cannot — the signal the PFM
    // ML rescorer needs. Absent on koth aligned matrices.
    "isoerror",
    "isoerror2",
    "rt_width",
];

/// One consensus LC-MS feature: a deisotoped, charge-resolved, RT-aligned
/// precursor mass shared across runs. This is what we try to assign a peptide to.
#[derive(Debug, Clone)]
pub struct Feature {
    /// Monoisotopic neutral mass (Da), in koth reference-run space.
    pub mass: f64,
    /// Monoisotopic m/z.
    pub mz: f64,
    pub charge: u8,
    /// Apex retention time (minutes), reference-run space.
    pub rt: f64,
    /// Ion mobility (1/K0); `None` when the data has no IM dimension.
    pub im: Option<f64>,
    /// koth feature-quality score (isotope x coelution), in [0, 1].
    pub combined_score: f64,
    /// Number of runs contributing a detection to this consensus feature.
    pub n_contributing_runs: u32,
    /// Optional per-feature shape descriptors (number of isotopic peaks and
    /// number of scans the feature spans). `None` when the matrix lacks them.
    pub n_isotopes: Option<f64>,
    pub n_scans: Option<f64>,
    /// Optional biosaur2 isotope-pattern fit error (deviation of the observed
    /// isotope envelope from theory) and RT peak width (minutes). Real peptides
    /// fit theory tightly; random mass coincidences (and shuffled decoys) do not.
    pub isoerror: Option<f64>,
    pub isoerror2: Option<f64>,
    pub rt_width: Option<f64>,
    /// Per-sample integrated intensities (aligned to [`FeatureMatrix::run_names`]).
    /// `0.0` means not detected in that run.
    pub intensities: Vec<f64>,
    /// Per-sample koth q-values, if a qvalue matrix was joined. Same order as
    /// `intensities`. `None` when no qvalue file was supplied.
    pub q_values: Option<Vec<f64>>,
}

/// The full aligned feature set: features plus the sample (run) names that index
/// every per-sample vector.
#[derive(Debug, Clone)]
pub struct FeatureMatrix {
    pub run_names: Vec<String>,
    pub features: Vec<Feature>,
}

/// Tolerances for collapsing near-duplicate consensus features (the same
/// analyte split across the cross-run alignment).
#[derive(Debug, Clone, serde::Deserialize)]
#[serde(default)]
pub struct DedupConfig {
    pub enabled: bool,
    pub ppm: f64,
    pub rt_window: f64,
    pub im_tol: f64,
}

impl Default for DedupConfig {
    fn default() -> Self {
        DedupConfig {
            enabled: false,
            ppm: 10.0,
            rt_window: 0.15,
            im_tol: 0.01,
        }
    }
}

impl FeatureMatrix {
    pub fn len(&self) -> usize {
        self.features.len()
    }

    /// Merge near-duplicate consensus features (same charge; mass/RT/IM within
    /// tolerance). Keeps the highest-`combined_score` member as representative,
    /// taking the per-sample MAX intensity (and MIN q-value) across the cluster
    /// so a split analyte's evidence is consolidated, not double-counted.
    /// Returns the number of features removed.
    pub fn dedup(&mut self, cfg: &DedupConfig) -> usize {
        let n = self.features.len();
        // index sorted by (charge, mass)
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| {
            let fa = &self.features[a];
            let fb = &self.features[b];
            fa.charge
                .cmp(&fb.charge)
                .then(fa.mass.partial_cmp(&fb.mass).unwrap())
        });

        let mut removed = vec![false; n];
        let mut i = 0;
        while i < n {
            let a = order[i];
            if removed[a] {
                i += 1;
                continue;
            }
            // gather a cluster of following features within tolerance of `a`
            let (ma, za, rta, ima) = (
                self.features[a].mass,
                self.features[a].charge,
                self.features[a].rt,
                self.features[a].im,
            );
            let mut members = vec![a];
            let mut j = i + 1;
            while j < n {
                let b = order[j];
                if self.features[b].charge != za {
                    break;
                }
                if (self.features[b].mass - ma) / ma * 1e6 > cfg.ppm {
                    break;
                }
                if !removed[b]
                    && (self.features[b].rt - rta).abs() <= cfg.rt_window
                    && match (ima, self.features[b].im) {
                        (Some(x), Some(y)) => (x - y).abs() <= cfg.im_tol,
                        _ => true,
                    }
                {
                    members.push(b);
                    removed[b] = true;
                }
                j += 1;
            }
            if members.len() > 1 {
                // representative = highest combined_score
                let rep = *members
                    .iter()
                    .max_by(|&&x, &&y| {
                        self.features[x]
                            .combined_score
                            .partial_cmp(&self.features[y].combined_score)
                            .unwrap()
                    })
                    .unwrap();
                let n_samples = self.features[rep].intensities.len();
                let mut merged_int = vec![0.0f64; n_samples];
                let mut merged_q = self.features[rep].q_values.clone();
                let mut max_contrib = 0u32;
                for &m in &members {
                    max_contrib = max_contrib.max(self.features[m].n_contributing_runs);
                    for s in 0..n_samples {
                        merged_int[s] = merged_int[s].max(self.features[m].intensities[s]);
                    }
                    if let (Some(mq), Some(q)) = (merged_q.as_mut(), self.features[m].q_values.as_ref()) {
                        for s in 0..n_samples {
                            mq[s] = mq[s].min(q[s]);
                        }
                    }
                }
                self.features[rep].intensities = merged_int;
                self.features[rep].q_values = merged_q;
                self.features[rep].n_contributing_runs = max_contrib;
                removed[rep] = false; // keep the representative
                for &m in &members {
                    if m != rep {
                        removed[m] = true;
                    }
                }
            }
            i += 1;
        }

        let kept: Vec<Feature> = self
            .features
            .drain(..)
            .enumerate()
            .filter_map(|(idx, f)| if removed[idx] { None } else { Some(f) })
            .collect();
        let n_removed = n - kept.len();
        self.features = kept;
        n_removed
    }
    pub fn is_empty(&self) -> bool {
        self.features.is_empty()
    }
    pub fn n_samples(&self) -> usize {
        self.run_names.len()
    }

    /// Load the aligned matrix from koth's `intensity_matrix.tsv`. If
    /// `qvalue_path` is given, per-sample q-values are joined by row index
    /// (the files must share row order and sample columns, which koth guarantees).
    pub fn load(
        intensity_path: impl AsRef<Path>,
        qvalue_path: Option<impl AsRef<Path>>,
    ) -> Result<Self> {
        let (run_names, rows) = read_matrix(intensity_path.as_ref())
            .with_context(|| "reading intensity matrix")?;

        let q_rows = match qvalue_path {
            Some(p) => {
                let (q_runs, q_rows) =
                    read_matrix(p.as_ref()).with_context(|| "reading qvalue matrix")?;
                if q_runs != run_names {
                    bail!(
                        "qvalue matrix sample columns ({}) do not match intensity matrix ({})",
                        q_runs.len(),
                        run_names.len()
                    );
                }
                if q_rows.len() != rows.len() {
                    bail!(
                        "qvalue matrix has {} rows but intensity matrix has {}",
                        q_rows.len(),
                        rows.len()
                    );
                }
                Some(q_rows)
            }
            None => None,
        };

        let features = rows
            .into_iter()
            .enumerate()
            .map(|(i, row)| {
                let q = q_rows.as_ref().map(|qr| qr[i].values.clone());
                row.into_feature(q)
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(FeatureMatrix {
            run_names,
            features,
        })
    }
}

/// A parsed matrix row: metadata fields plus the per-sample value vector.
struct MatrixRow {
    mass: f64,
    mz: f64,
    charge: u8,
    rt: f64,
    im: Option<f64>,
    combined_score: f64,
    n_contributing_runs: u32,
    n_isotopes: Option<f64>,
    n_scans: Option<f64>,
    isoerror: Option<f64>,
    isoerror2: Option<f64>,
    rt_width: Option<f64>,
    values: Vec<f64>,
}

impl MatrixRow {
    fn into_feature(self, q_values: Option<Vec<f64>>) -> Result<Feature> {
        Ok(Feature {
            mass: self.mass,
            mz: self.mz,
            charge: self.charge,
            rt: self.rt,
            im: self.im,
            combined_score: self.combined_score,
            n_contributing_runs: self.n_contributing_runs,
            n_isotopes: self.n_isotopes,
            n_scans: self.n_scans,
            isoerror: self.isoerror,
            isoerror2: self.isoerror2,
            rt_width: self.rt_width,
            intensities: self.values,
            q_values,
        })
    }
}

/// Index of each metadata column within a matrix header, plus the list of
/// sample-column indices (everything that is not metadata).
struct Layout {
    mass: usize,
    mz: usize,
    charge: usize,
    rt: usize,
    im: usize,
    combined_score: usize,
    n_contributing_runs: usize,
    n_isotopes: Option<usize>,
    n_scans: Option<usize>,
    isoerror: Option<usize>,
    isoerror2: Option<usize>,
    rt_width: Option<usize>,
    sample_cols: Vec<usize>,
    sample_names: Vec<String>,
}

impl Layout {
    fn from_header(header: &csv::StringRecord) -> Result<Self> {
        let find = |name: &str| -> Result<usize> {
            header
                .iter()
                .position(|h| h == name)
                .ok_or_else(|| anyhow!("missing required column `{name}` in koth matrix header"))
        };
        let find_opt = |name: &str| -> Option<usize> { header.iter().position(|h| h == name) };

        let mut sample_cols = Vec::new();
        let mut sample_names = Vec::new();
        for (i, h) in header.iter().enumerate() {
            if !META_COLUMNS.contains(&h) {
                sample_cols.push(i);
                sample_names.push(h.to_string());
            }
        }
        if sample_names.is_empty() {
            bail!("koth matrix has no per-sample value columns");
        }

        Ok(Layout {
            mass: find("massCalib")?,
            mz: find("mz")?,
            charge: find("charge")?,
            rt: find("rtApex")?,
            im: find("im")?,
            combined_score: find("combined_score")?,
            n_contributing_runs: find("n_contributing_runs")?,
            n_isotopes: find_opt("n_isotopes"),
            n_scans: find_opt("n_scans"),
            isoerror: find_opt("isoerror"),
            isoerror2: find_opt("isoerror2"),
            rt_width: find_opt("rt_width"),
            sample_cols,
            sample_names,
        })
    }
}

/// Read a koth matrix TSV. Returns (sample_names, rows).
fn read_matrix(path: &Path) -> Result<(Vec<String>, Vec<MatrixRow>)> {
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

        let im_raw = rec.get(layout.im).unwrap_or("").trim();
        let im = if im_raw.is_empty() {
            None
        } else {
            Some(
                im_raw
                    .parse::<f64>()
                    .map_err(|_| anyhow!("line {line}: cannot parse im `{im_raw}`"))?,
            )
        };

        let charge = rec
            .get(layout.charge)
            .unwrap_or("")
            .parse::<u8>()
            .map_err(|_| anyhow!("line {line}: cannot parse charge"))?;
        let n_contributing_runs = rec
            .get(layout.n_contributing_runs)
            .unwrap_or("")
            .parse::<u32>()
            .map_err(|_| anyhow!("line {line}: cannot parse n_contributing_runs"))?;
        // Optional shape descriptors: parse if the column exists and is non-empty.
        let opt_f64 = |col: Option<usize>| -> Option<f64> {
            col.and_then(|c| rec.get(c))
                .map(str::trim)
                .filter(|s| !s.is_empty())
                .and_then(|s| s.parse::<f64>().ok())
        };
        let n_isotopes = opt_f64(layout.n_isotopes);
        let n_scans = opt_f64(layout.n_scans);
        let isoerror = opt_f64(layout.isoerror);
        let isoerror2 = opt_f64(layout.isoerror2);
        let rt_width = opt_f64(layout.rt_width);

        let values = layout
            .sample_cols
            .iter()
            .map(|&c| {
                let s = rec.get(c).unwrap_or("").trim();
                if s.is_empty() {
                    Ok(0.0)
                } else {
                    s.parse::<f64>()
                        .map_err(|_| anyhow!("line {line}: cannot parse sample value `{s}`"))
                }
            })
            .collect::<Result<Vec<_>>>()?;

        rows.push(MatrixRow {
            mass: parse_f64(&rec, layout.mass, "massCalib", line)?,
            mz: parse_f64(&rec, layout.mz, "mz", line)?,
            charge,
            rt: parse_f64(&rec, layout.rt, "rtApex", line)?,
            im,
            combined_score: parse_f64(&rec, layout.combined_score, "combined_score", line)?,
            n_contributing_runs,
            n_isotopes,
            n_scans,
            isoerror,
            isoerror2,
            rt_width,
            values,
        });
    }

    Ok((layout.sample_names, rows))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn benchmark_dir() -> PathBuf {
        PathBuf::from("/home/patrick-garrett/Repos/koth_rust/benchmark/data/koth/align_output")
    }

    #[test]
    fn loads_real_koth_matrix() {
        let dir = benchmark_dir();
        let intensity = dir.join("intensity_matrix.tsv");
        if !intensity.exists() {
            eprintln!("skipping: benchmark matrix not present at {}", intensity.display());
            return;
        }
        let qvalue = dir.join("qvalue_matrix.tsv");
        let m = FeatureMatrix::load(&intensity, Some(&qvalue)).expect("load matrix");

        assert_eq!(m.n_samples(), 20, "20 sample columns expected");
        assert_eq!(m.len(), 55431, "row count");

        let f = &m.features[0];
        assert_eq!(f.charge, 2);
        assert!((f.mass - 798.373104).abs() < 1e-4);
        assert!((f.rt - 31.237699).abs() < 1e-4);
        assert!(f.im.is_none(), "Orbitrap data has no ion mobility");
        assert_eq!(f.intensities.len(), 20);
        assert!(f.intensities[0] > 0.0);
        let q = f.q_values.as_ref().expect("q-values joined");
        assert_eq!(q.len(), 20);
    }

    #[test]
    fn loads_without_qvalues() {
        let intensity = benchmark_dir().join("intensity_matrix.tsv");
        if !intensity.exists() {
            return;
        }
        let m = FeatureMatrix::load(&intensity, None::<&str>).expect("load");
        assert!(m.features[0].q_values.is_none());
    }
}
