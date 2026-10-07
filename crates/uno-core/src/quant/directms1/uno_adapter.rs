//! Uno-specific admission policy; deliberately distinct from native-reference mode.
//!
//! With `AdmissionPolicy::default()` this is the `reference` adaptation (the
//! former `uno-quant` binary) and must stay bit-for-bit identical to it. Each
//! further option is independent so it can be ablated alone:
//! - `feature_competition`: one feature quantifies for one peptide, the row
//!   whose best backing protein has the highest per-run binomial score
//!   (`crate::quant::feature_winners`, the rule the matrix path uses).
//! - `peptide_confidence`: keep only assignments the search's selected
//!   calibration scored (`--out-peptide-scores`), either those in the selected
//!   ensemble or those passing a per-run target-decoy q-value on `ml`.
//! - `intensity`: the feature apex (reference) or its integrated `intensitySum`
//!   looked up in the KOTH features table.
//! - `reduction`: which observation represents a precursor in a run.
use super::{
    input::{InputHash, Reader},
    Config, Dataset, InputPeptide,
};
use anyhow::{ensure, Context, Result};
use serde::Serialize;
use std::{
    collections::{BTreeMap, BTreeSet, HashMap},
    path::{Path, PathBuf},
    str::FromStr,
};
fn relative(base: &Path, value: &str) -> PathBuf {
    let p = Path::new(value);
    if p.is_absolute() {
        p.into()
    } else {
        base.join(p)
    }
}
fn is_default<T: Default + PartialEq>(v: &T) -> bool {
    *v == T::default()
}
/// Which search assignments count as confident peptide evidence.
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum PeptideConfidence {
    /// Every exported assignment (reference).
    #[default]
    None,
    /// Assignments inside the selected calibration's inclusion ensemble.
    Ensemble,
    /// Assignments whose per-run target-decoy q-value on `ml` passes `peptide_qvalue`.
    Qvalue,
}
/// Which feature intensity quantifies an observation.
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum IntensityKind {
    /// `intensityApex`, as exported in the evidence table (reference).
    #[default]
    Apex,
    /// `intensitySum` of the same feature.
    Integrated,
}
/// Which observation represents a precursor within one run.
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Reduction {
    /// The most intense observation (reference; DirectMS1Quant's rule).
    #[default]
    MaxIntensity,
    /// The observation with the smallest search deviation; ties go to intensity.
    MinDeviation,
}
macro_rules! parse_enum {
    ($t:ty { $($name:literal => $v:expr),+ $(,)? }) => {
        impl FromStr for $t {
            type Err = String;
            fn from_str(s: &str) -> Result<Self, String> {
                match s {
                    $($name => Ok($v),)+
                    _ => Err(format!("expected one of: {}", [$($name),+].join(", "))),
                }
            }
        }
    };
}
parse_enum!(PeptideConfidence { "none" => Self::None, "ensemble" => Self::Ensemble, "qvalue" => Self::Qvalue });
parse_enum!(IntensityKind { "apex" => Self::Apex, "integrated" => Self::Integrated });
parse_enum!(Reduction { "max-intensity" => Self::MaxIntensity, "min-deviation" => Self::MinDeviation });

#[derive(Clone, Debug, Default, PartialEq, Serialize)]
pub struct AdmissionPolicy {
    pub exclusive_features: bool,
    pub min_per_condition: usize,
    #[serde(skip_serializing_if = "is_default")]
    pub feature_competition: bool,
    #[serde(skip_serializing_if = "is_default")]
    pub peptide_confidence: PeptideConfidence,
    /// Threshold for `PeptideConfidence::Qvalue`.
    #[serde(skip_serializing_if = "is_default")]
    pub peptide_qvalue: f64,
    #[serde(skip_serializing_if = "is_default")]
    pub intensity: IntensityKind,
    #[serde(skip_serializing_if = "is_default")]
    pub reduction: Reduction,
}
impl AdmissionPolicy {
    fn needs_scores(&self) -> bool {
        self.peptide_confidence != PeptideConfidence::None
            || self.reduction == Reduction::MinDeviation
    }
}
#[derive(Clone, Debug, Default, Serialize)]
pub struct RunAudit {
    pub sample: String,
    pub candidate_assignments: usize,
    pub candidate_features: usize,
    pub ambiguous_features: usize,
    pub excluded_ambiguous_assignments: usize,
    pub retained_precursors: usize,
    #[serde(skip_serializing_if = "is_default")]
    pub removed_by_competition: usize,
    #[serde(skip_serializing_if = "is_default")]
    pub removed_by_confidence: usize,
    /// Precursors with more than one surviving observation, i.e. where the
    /// reduction rule had a choice to make. Always reported.
    pub precursors_with_multiple_observations: usize,
}
#[derive(Clone, Debug, Default, Serialize)]
pub struct AdmissionAudit {
    pub policy: AdmissionPolicy,
    pub runs: Vec<RunAudit>,
    pub precursors_before_condition_filter: usize,
    pub precursors_after_condition_filter: usize,
}
/// Per-assignment scores from `uno --out-peptide-scores`, keyed by
/// (sequence, decoy, feature index). `q` is the per-run target-decoy q-value on `ml`.
struct PeptideScore {
    in_ensemble: bool,
    deviation: f64,
    q: f64,
}
fn peptide_scores(
    reader: &mut Reader,
    path: &Path,
) -> Result<HashMap<(String, bool, u32), PeptideScore>> {
    let mut rows: Vec<((String, bool, u32), f64, f64, bool)> = Vec::new();
    reader.each(
        path,
        &[
            "peptide",
            "decoy",
            "feature_idx",
            "ml",
            "deviation",
            "in_ensemble",
        ],
        |r| {
            ensure!(
                matches!(&r["decoy"], "0" | "1") && matches!(&r["in_ensemble"], "0" | "1"),
                "invalid peptide-score flag"
            );
            let ml: f64 = r["ml"].parse().context("invalid ml")?;
            let deviation: f64 = r["deviation"].parse().context("invalid deviation")?;
            ensure!(
                ml.is_finite() && deviation.is_finite(),
                "non-finite peptide score"
            );
            rows.push((
                (
                    r["peptide"].to_owned(),
                    &r["decoy"] == "1",
                    r["feature_idx"].parse()?,
                ),
                ml,
                deviation,
                &r["in_ensemble"] == "1",
            ));
            Ok(())
        },
    )?;
    // q = decoys / targets at or above each score, made monotone from the bottom.
    let mut order: Vec<usize> = (0..rows.len()).collect();
    order.sort_by(|&a, &b| rows[b].1.total_cmp(&rows[a].1).then(a.cmp(&b)));
    let mut q = vec![0.; rows.len()];
    let (mut targets, mut decoys) = (0usize, 0usize);
    for &i in &order {
        if rows[i].0 .1 {
            decoys += 1;
        } else {
            targets += 1;
        }
        q[i] = decoys as f64 / targets.max(1) as f64;
    }
    let mut floor = f64::INFINITY;
    for &i in order.iter().rev() {
        floor = floor.min(q[i]);
        q[i] = floor;
    }
    let mut out = HashMap::with_capacity(rows.len());
    for (i, (key, _, deviation, in_ensemble)) in rows.into_iter().enumerate() {
        let score = PeptideScore {
            in_ensemble,
            deviation,
            q: q[i],
        };
        ensure!(
            out.insert(key, score).is_none(),
            "duplicate peptide-score assignment"
        );
    }
    Ok(out)
}
/// (intensityApex, intensitySum, massCalib) per KOTH feature, by data-row index.
fn feature_table(reader: &mut Reader, path: &Path) -> Result<Vec<(f64, f64, f64)>> {
    let mut out = Vec::new();
    reader.each(path, &["intensityApex", "intensitySum", "massCalib"], |r| {
        out.push((
            r["intensityApex"].parse()?,
            r["intensitySum"].parse()?,
            r["massCalib"].parse()?,
        ));
        Ok(())
    })?;
    Ok(out)
}
struct Sample {
    name: String,
    condition: &'static str,
    evidence: PathBuf,
    proteins: PathBuf,
    peptide_scores: Option<PathBuf>,
    features: Option<PathBuf>,
}
struct Candidate {
    feature: String,
    key: String,
    intensity: f64,
    deviation: f64,
    confident: bool,
    peptide: InputPeptide,
}
pub fn load(
    manifest: &Path,
    cfg: &Config,
    max_protein_qvalue: f64,
) -> Result<(Dataset, Vec<InputHash>)> {
    let (data, hashes, _) = load_with_policy(
        manifest,
        cfg,
        max_protein_qvalue,
        &AdmissionPolicy::default(),
    )?;
    Ok((data, hashes))
}
const EVIDENCE_COLUMNS: [&str; 13] = [
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
];
pub fn load_with_policy(
    manifest: &Path,
    cfg: &Config,
    max_protein_qvalue: f64,
    policy: &AdmissionPolicy,
) -> Result<(Dataset, Vec<InputHash>, AdmissionAudit)> {
    let mut audit = AdmissionAudit {
        policy: policy.clone(),
        ..AdmissionAudit::default()
    };
    ensure!(
        max_protein_qvalue.is_finite() && (0. ..=1.).contains(&max_protein_qvalue),
        "invalid protein identification qvalue threshold"
    );
    if policy.peptide_confidence == PeptideConfidence::Qvalue {
        ensure!(
            policy.peptide_qvalue.is_finite() && (0. ..=1.).contains(&policy.peptide_qvalue),
            "invalid peptide qvalue threshold"
        );
    }
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
                &r["condition"] == "S1" || &r["condition"] == "S2",
                "condition must be S1 or S2"
            );
            if &r["condition"] != condition {
                continue;
            }
            ensure!(
                !r["sample"].is_empty() && seen.insert(r["sample"].to_owned()),
                "duplicate/empty sample identifier"
            );
            let path = relative(base, &r["evidence"]);
            ensure!(
                files.insert(std::fs::canonicalize(&path)?),
                "an evidence file is assigned to more than one sample"
            );
            let optional = |column: &str| -> Result<Option<PathBuf>> {
                match r.get(column).filter(|v| !v.is_empty()) {
                    Some(v) => Ok(Some(relative(base, v))),
                    None => Ok(None),
                }
            };
            let sample = Sample {
                name: r["sample"].to_owned(),
                condition,
                evidence: path,
                proteins: relative(base, &r["proteins"]),
                peptide_scores: optional("peptide_scores")?,
                features: optional("features")?,
            };
            ensure!(
                !policy.needs_scores() || sample.peptide_scores.is_some(),
                "this policy needs a peptide_scores column (uno --out-peptide-scores)"
            );
            ensure!(
                policy.intensity == IntensityKind::Apex || sample.features.is_some(),
                "integrated intensity needs a features column (KOTH features.tsv)"
            );
            samples.push(sample);
        }
    }
    let n_s1 = samples.iter().filter(|s| s.condition == "S1").count();
    ensure!(
        n_s1 >= 2 && samples.len() - n_s1 >= 2,
        "at least two runs per condition required"
    );
    ensure!(
        policy.min_per_condition <= n_s1 && policy.min_per_condition <= samples.len() - n_s1,
        "minimum per-condition support exceeds replicate count"
    );
    let mut admitted = BTreeSet::new();
    // Per-run protein binomial scores keyed by (decoy, accession), as uno keys them.
    let mut protein_scores: Vec<HashMap<(bool, String), f64>> = Vec::new();
    for s in &samples {
        let mut scores = HashMap::new();
        let mut required = vec!["accession", "decoy", "q_value"];
        if policy.feature_competition {
            required.push("binomial_score");
        }
        reader.each(&s.proteins, &required, |row| {
            ensure!(
                &row["decoy"] == "0" || &row["decoy"] == "1",
                "invalid protein decoy flag"
            );
            let q: f64 = row["q_value"].parse()?;
            ensure!(
                q.is_finite() && (0. ..=1.).contains(&q),
                "invalid protein ID q-value"
            );
            if policy.feature_competition {
                let score: f64 = row["binomial_score"]
                    .parse()
                    .context("invalid binomial score")?;
                scores.insert((&row["decoy"] == "1", row["accession"].to_owned()), score);
            }
            if &row["decoy"] == "0" && q <= max_protein_qvalue {
                let accession = &row["accession"];
                ensure!(
                    !accession.starts_with(&cfg.decoy_prefix),
                    "target accession uses decoy prefix"
                );
                admitted.insert(accession.to_owned());
                admitted.insert(format!("{}{accession}", cfg.decoy_prefix));
            }
            Ok(())
        })?;
        protein_scores.push(scores);
    }
    ensure!(!admitted.is_empty(), "no admitted target proteins");
    let mut merged: BTreeMap<String, InputPeptide> = BTreeMap::new();
    for (sample, s) in samples.iter().enumerate() {
        let scores = match &s.peptide_scores {
            Some(path) if policy.needs_scores() => Some(peptide_scores(&mut reader, path)?),
            _ => None,
        };
        let features = match &s.features {
            Some(path) if policy.intensity == IntensityKind::Integrated => {
                Some(feature_table(&mut reader, path)?)
            }
            _ => None,
        };
        let mut run_names = BTreeSet::new();
        let mut observations = BTreeSet::new();
        // Every evidence row: its feature, its competition score, and the
        // candidate it yields if it passes protein/length/intensity admission.
        let mut rows: Vec<(u32, f64, Option<Candidate>)> = Vec::new();
        reader.each(&s.evidence, &EVIDENCE_COLUMNS, |row| {
            ensure!(
                &row["schema_version"] == "1"
                    && &row["identity_evidence"] == "uno_ms1_rt_mass_assignment",
                "unsupported evidence schema or identity type"
            );
            run_names.insert(row["run"].to_owned());
            ensure!(
                run_names.len() == 1,
                "evidence contains multiple runs; matrix input is unsupported by this adapter"
            );
            ensure!(
                observations.insert((
                    row["feature_id"].to_owned(),
                    row["peptide_id"].to_owned(),
                    row["decoy"].to_owned()
                )),
                "duplicate observation identifier"
            );
            ensure!(
                &row["decoy"] == "0" || &row["decoy"] == "1",
                "invalid peptide decoy flag"
            );
            let decoy = &row["decoy"] == "1";
            let feature_idx: u32 = row["feature_id"].parse().context("invalid feature_id")?;
            let competition = if policy.feature_competition {
                row["proteins"]
                    .split(';')
                    .filter_map(|a| protein_scores[sample].get(&(decoy, a.to_owned())).copied())
                    .fold(f64::NEG_INFINITY, f64::max)
            } else {
                0.
            };
            rows.push((feature_idx, competition, None));
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
                return Ok(());
            }
            let charge: u8 = row["charge"].parse().context("invalid charge")?;
            ensure!(charge > 0, "zero charge");
            let im: Option<f64> = if row["ion_mobility"].is_empty() {
                None
            } else {
                Some(row["ion_mobility"].parse()?)
            };
            ensure!(im.is_none_or(|x| x.is_finite()), "invalid mobility");
            let mut intensity: f64 = row["intensity"].parse()?;
            ensure!(
                intensity.is_finite() && intensity >= 0.,
                "invalid intensity"
            );
            if let Some(table) = &features {
                let &(apex, sum, mass) = table
                    .get(feature_idx as usize)
                    .context("feature_id outside the features table")?;
                let observed: f64 = row["observed_mass"].parse()?;
                ensure!(
                    (apex - intensity).abs() <= 1e-5 * apex.abs().max(1.)
                        && (mass - observed).abs() <= 1e-5,
                    "features table does not match evidence at feature {feature_idx}"
                );
                ensure!(sum.is_finite() && sum >= 0., "invalid integrated intensity");
                intensity = sum;
            }
            if intensity == 0. {
                return Ok(());
            }
            let score = scores
                .as_ref()
                .and_then(|s| s.get(&(sequence.to_owned(), decoy, feature_idx)));
            // Unscored assignments fail every confidence rule; with no score
            // file, deviation is unknown (NaN) and never consulted.
            let confident = match (policy.peptide_confidence, score) {
                (PeptideConfidence::None, _) => true,
                (_, None) => false,
                (PeptideConfidence::Ensemble, Some(p)) => p.in_ensemble,
                (PeptideConfidence::Qvalue, Some(p)) => p.q <= policy.peptide_qvalue,
            };
            let deviation = score.map_or(f64::NAN, |p| p.deviation);
            let key = format!(
                "{}|{}|{}|{}",
                sequence,
                charge,
                im.map(|v| v.to_string()).unwrap_or_default(),
                u8::from(decoy)
            );
            let peptide = InputPeptide {
                key: key.clone(),
                sequence: sequence.to_owned(),
                proteins,
                values: vec![None; samples.len()],
            };
            rows.last_mut().unwrap().2 = Some(Candidate {
                feature: row["feature_id"].to_owned(),
                key,
                intensity,
                deviation,
                confident,
                peptide,
            });
            Ok(())
        })?;
        let winners: Option<BTreeSet<usize>> = policy.feature_competition.then(|| {
            let (features, scores): (Vec<u32>, Vec<f64>) =
                rows.iter().map(|(f, s, _)| (*f, *s)).unzip();
            crate::quant::feature_winners(&features, &scores)
                .into_iter()
                .collect()
        });
        let (mut removed_by_competition, mut removed_by_confidence) = (0, 0);
        let mut candidates = Vec::new();
        let mut feature_keys: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
        for (i, (_, _, candidate)) in rows.into_iter().enumerate() {
            let Some(c) = candidate else { continue };
            if winners.as_ref().is_some_and(|w| !w.contains(&i)) {
                removed_by_competition += 1;
                continue;
            }
            if !c.confident {
                removed_by_confidence += 1;
                continue;
            }
            feature_keys
                .entry(c.feature.clone())
                .or_default()
                .insert(c.key.clone());
            candidates.push(c);
        }
        let mut run_audit = RunAudit {
            sample: s.name.clone(),
            candidate_assignments: candidates.len(),
            candidate_features: feature_keys.len(),
            ambiguous_features: feature_keys.values().filter(|keys| keys.len() > 1).count(),
            removed_by_competition,
            removed_by_confidence,
            ..RunAudit::default()
        };
        // Symmetric target/decoy exclusion, before choosing the best remaining
        // feature for a precursor. This permits an unambiguous lower-intensity
        // observation to replace an excluded ambiguous maximum.
        let mut best: BTreeMap<String, Candidate> = BTreeMap::new();
        let mut seen_twice = BTreeSet::new();
        for c in candidates {
            if policy.exclusive_features && feature_keys[&c.feature].len() > 1 {
                run_audit.excluded_ambiguous_assignments += 1;
                continue;
            }
            let replace = match best.get(&c.key) {
                None => true,
                Some(old) => {
                    seen_twice.insert(c.key.clone());
                    match policy.reduction {
                        Reduction::MaxIntensity => c.intensity > old.intensity,
                        Reduction::MinDeviation => {
                            c.deviation < old.deviation
                                || (c.deviation == old.deviation && c.intensity > old.intensity)
                        }
                    }
                }
            };
            if replace {
                best.insert(c.key.clone(), c);
            }
        }
        run_audit.precursors_with_multiple_observations = seen_twice.len();
        run_audit.retained_precursors = best.len();
        audit.runs.push(run_audit);
        for (key, c) in best {
            let dest = merged.entry(key).or_insert(c.peptide);
            dest.values[sample] = Some(c.intensity);
        }
    }
    audit.precursors_before_condition_filter = merged.len();
    merged.retain(|_, row| {
        row.values[..n_s1].iter().flatten().count() >= policy.min_per_condition
            && row.values[n_s1..].iter().flatten().count() >= policy.min_per_condition
    });
    audit.precursors_after_condition_filter = merged.len();
    Ok((
        Dataset {
            samples: samples.into_iter().map(|s| s.name).collect(),
            n_s1,
            rows: merged.into_values().collect(),
        },
        reader.hashes,
        audit,
    ))
}
