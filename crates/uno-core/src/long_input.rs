//! KOTH long-format bundle ingestion. Quantities are attached once per cell;
//! original observations remain separately available for per-run batch search.
use crate::input::{Feature, FeatureMatrix};
use anyhow::{ensure, Context, Result};
use serde::Deserialize;
use sha2::{Digest, Sha256};
use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::Read,
    path::{Path, PathBuf},
};

#[derive(Deserialize)]
struct Table {
    path: String,
    sha256: String,
}
#[derive(Deserialize)]
pub struct Run {
    pub run_id: usize,
    pub name: String,
    pub n_observations: usize,
    #[serde(default)]
    sha256: Option<String>,
}
#[derive(Deserialize)]
pub struct Bundle {
    #[serde(default)]
    config: Option<serde_json::Value>,
    schema_version: u32,
    n_consensus: usize,
    reference_run_id: usize,
    pub runs: Vec<Run>,
    observations: Table,
    cells: Table,
    #[serde(skip)]
    root: PathBuf,
}

pub fn file_sha(path: &Path) -> Result<String> {
    let mut f = File::open(path)?;
    let mut h = Sha256::new();
    let mut buf = [0; 65536];
    loop {
        let n = f.read(&mut buf)?;
        if n == 0 {
            break;
        }
        h.update(&buf[..n]);
    }
    Ok(format!("{:x}", h.finalize()))
}
impl Bundle {
    pub fn open(path: &Path) -> Result<Self> {
        let path = if path.is_dir() {
            path.join("matrix_manifest.json")
        } else {
            path.to_owned()
        };
        let mut b: Self = serde_json::from_reader(File::open(&path)?)?;
        b.root = path
            .parent()
            .context("manifest requires parent directory")?
            .to_owned();
        ensure!(
            (1..=3).contains(&b.schema_version),
            "unsupported long-matrix schema version {} (supported: 1, 2, 3)",
            b.schema_version
        );
        ensure!(
            !b.runs.is_empty() && b.reference_run_id < b.runs.len(),
            "invalid run/reference list"
        );
        let mut names = HashSet::new();
        for (i, r) in b.runs.iter().enumerate() {
            ensure!(
                r.run_id == i && names.insert(&r.name),
                "duplicate or invalid run ID/name"
            );
        }
        for t in [&b.observations, &b.cells] {
            ensure!(
                matches!(
                    Path::new(&t.path).components().next(),
                    Some(std::path::Component::Normal(_))
                ) && Path::new(&t.path).components().count() == 1,
                "table path must be a filename"
            );
            ensure!(
                file_sha(&b.root.join(&t.path))? == t.sha256,
                "bundle table checksum mismatch: {}",
                t.path
            );
        }
        Ok(b)
    }
    fn observations(&self, mut visit: impl FnMut(Observation) -> Result<()>) -> Result<()> {
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(self.root.join(&self.observations.path))?;
        let mut counts = vec![0; self.runs.len()];
        let mut previous = None;
        for row in rdr.deserialize::<Observation>() {
            let o = row?;
            ensure!(o.run_id < counts.len(), "unknown observation run");
            ensure!(
                o.observation_id == counts[o.run_id],
                "duplicate, missing, or unordered observation ID"
            );
            ensure!(
                previous.is_none_or(|r| o.run_id >= r),
                "observations must be ordered by run"
            );
            previous = Some(o.run_id);
            counts[o.run_id] += 1;
            o.validate()?;
            if let Some(c) = o.consensus_id {
                ensure!(c < self.n_consensus, "invalid observation consensus");
            }
            visit(o)?;
        }
        ensure!(
            counts
                .iter()
                .zip(&self.runs)
                .all(|(n, r)| *n == r.n_observations),
            "incomplete observations table"
        );
        Ok(())
    }
    /// Original source order and measurements reproduce the single-file route.
    pub fn load_run(&self, run: usize) -> Result<FeatureMatrix> {
        ensure!(run < self.runs.len(), "unknown run");
        let mut features = Vec::new();
        self.observations(|o| {
            if o.run_id == run {
                features.push(o.feature());
            }
            Ok(())
        })?;
        Ok(FeatureMatrix {
            run_names: vec![self.runs[run].name.clone()],
            features,
            matrix_evidence: None,
            joint_evidence: None,
            joint_cell_evidence: None,
        })
    }
    /// One matrix search retaining every original member (including alternatives).
    /// KOTH-aligned coordinates supply a common calibration frame; no independent
    /// run searches or protein lists are involved. Missing/MBR cells supply no IDs.
    pub fn load_joint(&self, min_runs: u32) -> Result<FeatureMatrix> {
        // Reuse complete cell/membership validation, keeping consensus index order.
        let mut matrix = self.load_consensus(0, true, false)?;
        let mut groups = vec![Vec::new(); self.n_consensus];
        // Schema 1 used a seed-quality admission rule. New bundles have already
        // passed KOTH's cross-run group-confidence gate; do not silently restore
        // the removed seed floor. Original feature anchor checks still apply.
        let seed_floor = if self.schema_version == 1 {
            self.config
                .as_ref()
                .and_then(|c| c.pointer("/lfq/consensus/min_seed_combined_score"))
                .and_then(serde_json::Value::as_f64)
                .unwrap_or(0.75)
        } else {
            0.0
        };
        ensure!(
            seed_floor.is_finite() && (0.0..=1.0).contains(&seed_floor),
            "invalid joint bootstrap quality floor"
        );
        let mut source_ids = HashMap::new();
        let run_ids: Vec<usize> = self
            .runs
            .iter()
            .map(|r| {
                r.sha256
                    .as_ref()
                    .map_or(r.run_id, |hash| *source_ids.entry(hash).or_insert(r.run_id))
            })
            .collect();
        self.observations(|o| {
            if let Some(c) = o.consensus_id {
                let mut f = o.feature();
                f.mass = o
                    .aligned_mass
                    .context("joint input requires aligned mass for every member")?;
                f.rt = o
                    .aligned_rt
                    .context("joint input requires aligned RT for every member")?;
                f.im = o.aligned_im;
                ensure!(
                    f.charge == matrix.features[c].charge,
                    "member/consensus charge mismatch"
                );
                f.mz = f.mass / f.charge as f64 + 1.007_276_466_621;
                groups[c].push(crate::joint::Observation {
                    bootstrap_eligible: matrix.features[c].combined_score >= seed_floor,
                    run_id: run_ids[o.run_id],
                    observation_id: o.observation_id,
                    feature: f,
                });
            }
            Ok(())
        })?;
        // Cell table was fully validated by load_consensus above. Canonical source
        // IDs prevent duplicate input files from inflating extraction support.
        let mut cell_summaries = vec![crate::cell_evidence::Summary::default(); self.n_consensus];
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(self.root.join(&self.cells.path))?;
        for row in rdr.deserialize::<Cell>() {
            let e = row?;
            if e.extraction_kind == "target" && run_ids[e.run_id] == e.run_id {
                cell_summaries[e.consensus_id].add(e.lfq_is_mbr, e.lfq_intensity, e.lfq_q_value);
            }
        }
        let n_sources = run_ids.iter().copied().collect::<HashSet<_>>().len();
        let mut retained_cells = Vec::new();
        let mut retained = Vec::new();
        let mut retained_groups = Vec::new();
        for ((cell, mut group), summary) in
            matrix.features.into_iter().zip(groups).zip(cell_summaries)
        {
            group.sort_by_key(|o| (o.run_id, o.observation_id));
            // The same source table submitted twice must not become two pieces of
            // evidence; retain the first source's aligned observation deterministically.
            group.dedup_by_key(|o| (o.run_id, o.observation_id));
            let n_runs = group.iter().map(|o| o.run_id).collect::<HashSet<_>>().len();
            if n_runs < min_runs.max(1) as usize {
                continue;
            }
            // A real member supplies calibration coordinates AND its own shape.
            // Three-isotope members qualify for UNO's existing anchor criterion.
            let representative = group
                .iter()
                .max_by(|a, b| {
                    let a_iso = a.feature.n_isotopes.is_some_and(|n| n >= 3.0);
                    let b_iso = b.feature.n_isotopes.is_some_and(|n| n >= 3.0);
                    a_iso
                        .cmp(&b_iso)
                        .then(
                            a.feature
                                .combined_score
                                .total_cmp(&b.feature.combined_score),
                        )
                        .then_with(|| {
                            (b.run_id, b.observation_id).cmp(&(a.run_id, a.observation_id))
                        })
                })
                .context("joint group has no original members")?;
            let mut f = representative.feature.clone();
            f.intensities = cell.intensities; // exactly one integrated LFQ value per cell
            f.n_contributing_runs = n_runs as u32;
            retained.push(f);
            retained_groups.push(group);
            retained_cells.push(summary.values(n_sources));
        }
        log::info!(
            "joint matrix: {} groups, {} original observations, {} runs, min_runs={}",
            retained.len(),
            retained_groups.iter().map(Vec::len).sum::<usize>(),
            self.runs.len(),
            min_runs.max(1)
        );
        matrix.features = retained;
        matrix.joint_evidence = Some(retained_groups);
        matrix.joint_cell_evidence = Some(retained_cells);
        Ok(matrix)
    }
    pub fn load_consensus(
        &self,
        min_runs: u32,
        quality: bool,
        legacy_values: bool,
    ) -> Result<FeatureMatrix> {
        self.load_consensus_with_representative(min_runs, quality, legacy_values, false)
    }
    /// Optional diagnostic: choose the best actual member that clears UNO's
    /// existing three-isotope anchor criterion, falling back to the seed.
    pub fn load_consensus_with_representative(
        &self,
        min_runs: u32,
        quality: bool,
        legacy_values: bool,
        prefer_isotopes: bool,
    ) -> Result<FeatureMatrix> {
        ensure!(
            !prefer_isotopes || quality,
            "isotope representative requires original quality"
        );
        let nr = self.runs.len();
        let mut features: Vec<Option<Feature>> = vec![None; self.n_consensus];
        let mut seeds = vec![None; self.n_consensus];
        let mut cells_seen = HashSet::new();
        let mut primary: HashMap<(usize, usize), usize> = HashMap::new();
        let mut direct = vec![0usize; self.n_consensus];
        let mut group_evidence = vec![None; self.n_consensus];
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(self.root.join(&self.cells.path))?;
        let release_columns = [
            "group_score",
            "group_qvalue",
            "lfq_ownership_status",
            "lfq_owned_samples",
            "lfq_excluded_samples",
            "lfq_competing_consensus_id",
            "lfq_preceding_signal_fraction",
        ];
        let required = match self.schema_version {
            1 => 0,
            2 => 2,
            _ => release_columns.len(),
        };
        for field in &release_columns[..required] {
            ensure!(
                rdr.headers()?.iter().any(|h| h == *field),
                "schema {} requires column {field}",
                self.schema_version
            );
        }
        for row in rdr.deserialize::<Cell>() {
            let e = row?;
            let c = e.consensus_id;
            let r = e.run_id;
            ensure!(
                c < features.len() && r < nr && e.seed_run_id < nr,
                "invalid cell/seed identity"
            );
            self.validate_release_evidence(&e)?;
            if self.schema_version >= 2 {
                let evidence = (e.group_score.unwrap(), e.group_qvalue.unwrap());
                if let Some(previous) = group_evidence[c] {
                    ensure!(
                        previous == evidence,
                        "inconsistent group confidence across cells"
                    );
                } else {
                    group_evidence[c] = Some(evidence);
                }
            }
            ensure!(
                matches!(e.extraction_kind.as_str(), "target" | "decoy"),
                "invalid extraction kind"
            );
            ensure!(
                cells_seen.insert((c, r, e.extraction_kind == "decoy")),
                "duplicate LFQ cell"
            );
            ensure!(
                [e.mass, e.mz, e.rt, e.combined_score, e.lfq_intensity]
                    .iter()
                    .all(|x| x.is_finite())
                    && e.mass > 0.0
                    && e.mz > 0.0
                    && e.charge > 0
                    && e.lfq_intensity >= 0.0,
                "invalid cell coordinates/intensity"
            );
            ensure!(
                e.im.is_none_or(|x| x.is_finite()) && (0.0..=1.0).contains(&e.combined_score),
                "invalid cell quality/IM"
            );
            for value in [
                e.lfq_expected_rt,
                e.lfq_apex_rt,
                e.lfq_expected_mz,
                e.lfq_observed_mz,
                e.lfq_hybrid_score,
                e.lfq_isotope_pattern,
                e.lfq_coelution,
            ]
            .into_iter()
            .flatten()
            {
                ensure!(value.is_finite(), "nonfinite extraction measurement");
            }
            ensure!(
                e.lfq_q_value
                    .is_none_or(|q| q.is_finite() && (0.0..=1.0).contains(&q)),
                "invalid extraction q-value"
            );
            ensure!(
                e.n_contributing_runs > 0 && e.n_contributing_runs as usize <= nr,
                "invalid contributing-run count"
            );
            let seed = (e.seed_run_id, e.seed_observation_id);
            if let Some(f) = &features[c] {
                ensure!(
                    (
                        f.mass,
                        f.mz,
                        f.charge,
                        f.rt,
                        f.im,
                        f.combined_score,
                        f.n_contributing_runs
                    ) == (
                        e.mass,
                        e.mz,
                        e.charge,
                        e.rt,
                        e.im,
                        e.combined_score,
                        e.n_contributing_runs
                    ) && seeds[c] == Some(seed),
                    "inconsistent consensus metadata"
                );
            } else {
                features[c] = Some(Feature {
                    mass: e.mass,
                    mz: e.mz,
                    charge: e.charge,
                    rt: e.rt,
                    im: e.im,
                    combined_score: e.combined_score,
                    n_contributing_runs: e.n_contributing_runs,
                    n_isotopes: None,
                    n_scans: None,
                    isoerror: None,
                    isoerror2: None,
                    rt_width: None,
                    cosine_score: None,
                    isotope_score: None,
                    ppm_error: None,
                    neutron_offset: None,
                    intensities: vec![0.0; nr],
                });
                seeds[c] = Some(seed);
            }
            if e.extraction_kind == "decoy" {
                ensure!(
                    e.observation_id.is_none() && e.lfq_q_value.is_none(),
                    "synthetic decoy has original observation or target q-value"
                );
                continue;
            }
            ensure!(
                e.lfq_is_mbr == e.observation_id.is_none(),
                "direct/MBR observation mismatch"
            );
            if let Some(o) = e.observation_id {
                ensure!(
                    primary.insert((r, o), c).is_none(),
                    "primary observation reused across cells"
                );
                direct[c] += 1;
            }
            features[c].as_mut().unwrap().intensities[r] = e.lfq_intensity;
        }
        for c in 0..features.len() {
            ensure!(
                features[c].is_some() && (0..nr).all(|r| cells_seen.contains(&(c, r, false))),
                "missing consensus/target cell"
            );
            ensure!(
                direct[c] == features[c].as_ref().unwrap().n_contributing_runs as usize,
                "direct-cell/support mismatch"
            );
        }
        let seed_lookup: HashMap<_, _> = seeds
            .iter()
            .enumerate()
            .filter_map(|(i, s)| s.map(|s| (s, i)))
            .collect();
        ensure!(
            seed_lookup.len() == seeds.len(),
            "seed reused across consensus groups"
        );
        let mut seed_seen = HashSet::new();
        let mut primary_seen = HashSet::new();
        let mut member_runs = HashSet::new();
        let mut representatives: HashMap<usize, Observation> = HashMap::new();
        self.observations(|o| {
            let key = (o.run_id, o.observation_id);
            if let Some(c) = o.consensus_id {
                member_runs.insert((c, o.run_id));
            }
            ensure!(
                o.is_primary == primary.contains_key(&key),
                "primary observation flag/link mismatch"
            );
            ensure!(
                o.is_seed == seed_lookup.contains_key(&key),
                "seed observation flag/link mismatch"
            );
            if let Some(&c) = primary.get(&key) {
                ensure!(o.consensus_id == Some(c), "primary links wrong group");
                primary_seen.insert(key);
            }
            if let Some(&c) = seed_lookup.get(&key) {
                ensure!(
                    o.consensus_id == Some(c) && o.is_primary,
                    "seed is not a primary member"
                );
                let f = features[c].as_mut().unwrap();
                ensure!(
                    o.feature_charge == f.charge && o.feature_combined_score == f.combined_score,
                    "seed quality/charge mismatch"
                );
                if quality {
                    o.apply_quality(f);
                }
                seed_seen.insert(key);
            }
            if prefer_isotopes && o.feature_n_isotopes.is_some_and(|n| n >= 3.0) {
                if let Some(c) = o.consensus_id {
                    if representatives
                        .get(&c)
                        .is_none_or(|prior| o.feature_combined_score > prior.feature_combined_score)
                    {
                        representatives.insert(c, o);
                    }
                }
            }
            Ok(())
        })?;
        ensure!(
            primary_seen.len() == primary.len() && seed_seen.len() == seed_lookup.len(),
            "missing original observation link"
        );
        let mut support = vec![0usize; features.len()];
        for (c, _) in member_runs {
            support[c] += 1;
        }
        ensure!(support == direct, "member/direct run support mismatch");
        for (c, o) in representatives {
            let f = features[c].as_mut().unwrap();
            f.mass = o
                .aligned_mass
                .context("representative requires aligned mass")?;
            f.mz = (f.mass + f.charge as f64 * 1.007_276_466_621) / f.charge as f64;
            f.rt = o.aligned_rt.context("representative requires aligned RT")?;
            f.im = o.aligned_im;
            f.combined_score = o.feature_combined_score;
            o.apply_quality(f);
        }
        let mut retained = Vec::new();
        for f in features.into_iter().flatten() {
            if f.n_contributing_runs < min_runs {
                continue;
            }
            let mut f = f;
            if legacy_values {
                f.mass = format!("{:.6}", f.mass).parse()?;
                f.mz = format!("{:.6}", f.mz).parse()?;
                f.rt = format!("{:.6}", f.rt).parse()?;
                f.im = f.im.map(|x| format!("{x:.6}").parse().unwrap());
                f.combined_score = format!("{:.6}", f.combined_score).parse()?;
                for x in &mut f.intensities {
                    *x = format!("{x:.5e}").parse()?;
                }
            }
            retained.push(f);
        }
        Ok(FeatureMatrix {
            run_names: self.runs.iter().map(|r| r.name.clone()).collect(),
            features: retained,
            matrix_evidence: None,
            joint_evidence: None,
            joint_cell_evidence: None,
        })
    }

    fn validate_release_evidence(&self, e: &Cell) -> Result<()> {
        if self.schema_version >= 2 {
            let score = e.group_score.context("schema 2/3 requires group_score")?;
            let q = e.group_qvalue.context("schema 2/3 requires group_qvalue")?;
            let limit = self
                .config
                .as_ref()
                .and_then(|c| c.pointer("/lfq/consensus/max_group_qvalue"))
                .and_then(serde_json::Value::as_f64)
                .context("schema 2/3 requires max_group_qvalue in manifest config")?;
            ensure!(
                score.is_finite()
                    && q.is_finite()
                    && (0.0..=1.0).contains(&q)
                    && limit.is_finite()
                    && (0.0..=1.0).contains(&limit),
                "invalid group confidence"
            );
            ensure!(
                q <= limit && e.n_contributing_runs >= 2,
                "group fails declared confidence or original-run support requirement"
            );
        }
        if self.schema_version >= 3 {
            let status = e
                .lfq_ownership_status
                .as_deref()
                .context("schema 3 requires ownership status")?;
            ensure!(
                matches!(
                    status,
                    "exclusive" | "residual" | "shared_signal_only" | "ambiguous_residual"
                ),
                "invalid LFQ ownership status"
            );
            let owned = e
                .lfq_owned_samples
                .context("schema 3 requires owned-sample count")?;
            e.lfq_excluded_samples
                .context("schema 3 requires excluded-sample count")?;
            let fraction = e
                .lfq_preceding_signal_fraction
                .context("schema 3 requires preceding signal fraction")?;
            ensure!(
                fraction.is_finite() && (0.0..=1.0).contains(&fraction),
                "invalid preceding signal fraction"
            );
            ensure!(
                e.lfq_competing_consensus_id
                    .is_none_or(|c| c < self.n_consensus && c != e.consensus_id),
                "invalid competing consensus identity"
            );
            if matches!(status, "shared_signal_only" | "ambiguous_residual") {
                ensure!(
                    e.lfq_intensity == 0.0 && owned == 0,
                    "suppressed LFQ signal must have zero intensity and ownership"
                );
            }
        }
        Ok(())
    }
}

#[derive(Deserialize)]
struct Cell {
    consensus_id: usize,
    run_id: usize,
    extraction_kind: String,
    observation_id: Option<usize>,
    seed_run_id: usize,
    seed_observation_id: usize,
    mass: f64,
    mz: f64,
    charge: u8,
    rt: f64,
    im: Option<f64>,
    combined_score: f64,
    n_contributing_runs: u32,
    group_score: Option<f64>,
    group_qvalue: Option<f64>,
    lfq_intensity: f64,
    lfq_q_value: Option<f64>,
    lfq_is_mbr: bool,
    lfq_expected_rt: Option<f64>,
    lfq_apex_rt: Option<f64>,
    lfq_expected_mz: Option<f64>,
    lfq_observed_mz: Option<f64>,
    lfq_hybrid_score: Option<f64>,
    lfq_isotope_pattern: Option<f64>,
    lfq_coelution: Option<f64>,
    lfq_ownership_status: Option<String>,
    lfq_owned_samples: Option<usize>,
    lfq_excluded_samples: Option<usize>,
    lfq_competing_consensus_id: Option<usize>,
    lfq_preceding_signal_fraction: Option<f64>,
}
#[derive(Deserialize)]
struct Observation {
    run_id: usize,
    observation_id: usize,
    consensus_id: Option<usize>,
    is_primary: bool,
    is_seed: bool,
    feature_mass: Option<f64>,
    feature_mz: f64,
    feature_charge: u8,
    feature_rt: f64,
    feature_im: Option<f64>,
    feature_combined_score: f64,
    feature_n_isotopes: Option<f64>,
    feature_n_scans: Option<f64>,
    feature_cosine_score: Option<f64>,
    feature_isotope_score: Option<f64>,
    feature_ppm_error: Option<f64>,
    feature_neutron_offset: Option<f64>,
    feature_rt_width: Option<f64>,
    feature_isoerror: Option<f64>,
    feature_isoerror2: Option<f64>,
    feature_intensity_apex: f64,
    aligned_mass: Option<f64>,
    aligned_rt: Option<f64>,
    aligned_im: Option<f64>,
}
impl Observation {
    fn validate(&self) -> Result<()> {
        ensure!(
            [
                self.feature_mz,
                self.feature_rt,
                self.feature_combined_score,
                self.feature_intensity_apex
            ]
            .iter()
            .all(|x| x.is_finite()),
            "nonfinite original feature"
        );
        ensure!(
            self.feature_mz > 0.0
                && self.feature_intensity_apex >= 0.0
                && (0.0..=1.0).contains(&self.feature_combined_score),
            "invalid original feature"
        );
        ensure!(
            self.feature_charge == 0 || self.feature_mass.is_some_and(|x| x > 0.0),
            "missing original mass"
        );
        for x in [
            self.feature_mass,
            self.aligned_mass,
            self.aligned_rt,
            self.aligned_im,
            self.feature_im,
            self.feature_n_isotopes,
            self.feature_n_scans,
            self.feature_cosine_score,
            self.feature_isotope_score,
            self.feature_ppm_error,
            self.feature_neutron_offset,
            self.feature_rt_width,
            self.feature_isoerror,
            self.feature_isoerror2,
        ]
        .into_iter()
        .flatten()
        {
            ensure!(x.is_finite(), "nonfinite original descriptor");
        }
        for count in [self.feature_n_isotopes, self.feature_n_scans]
            .into_iter()
            .flatten()
        {
            ensure!(
                count >= 0.0 && count.fract() == 0.0,
                "invalid original feature count"
            );
        }
        for fraction in [self.feature_cosine_score, self.feature_isotope_score]
            .into_iter()
            .flatten()
        {
            ensure!(
                (0.0..=1.0).contains(&fraction),
                "invalid original quality fraction"
            );
        }
        ensure!(
            self.feature_rt_width.is_none_or(|x| x >= 0.0),
            "negative original peak width"
        );
        ensure!(
            self.feature_charge == 0 || self.aligned_mass.is_none_or(|x| x > 0.0),
            "invalid aligned mass"
        );
        Ok(())
    }
    fn apply_quality(&self, f: &mut Feature) {
        f.n_isotopes = self.feature_n_isotopes;
        f.n_scans = self.feature_n_scans;
        f.cosine_score = self.feature_cosine_score;
        f.isotope_score = self.feature_isotope_score;
        f.ppm_error = self.feature_ppm_error;
        f.neutron_offset = self.feature_neutron_offset;
        f.rt_width = self.feature_rt_width;
        f.isoerror = self.feature_isoerror;
        f.isoerror2 = self.feature_isoerror2;
    }
    fn feature(&self) -> Feature {
        let mut f = Feature {
            mass: self.feature_mass.unwrap_or(0.0),
            mz: self.feature_mz,
            charge: self.feature_charge,
            rt: self.feature_rt,
            im: self.feature_im,
            combined_score: self.feature_combined_score,
            n_contributing_runs: 1,
            n_isotopes: None,
            n_scans: None,
            cosine_score: None,
            isotope_score: None,
            ppm_error: None,
            neutron_offset: None,
            rt_width: None,
            isoerror: None,
            isoerror2: None,
            intensities: vec![self.feature_intensity_apex],
        };
        self.apply_quality(&mut f);
        f
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    struct Fixture(PathBuf);
    impl Drop for Fixture {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }
    impl Fixture {
        fn new() -> Self {
            let p = std::env::temp_dir().join(format!(
                "uno-long-{}-{}",
                std::process::id(),
                std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .unwrap()
                    .as_nanos()
            ));
            std::fs::create_dir(&p).unwrap();
            std::fs::write(p.join("observations.tsv"),concat!(
                "run_id\tobservation_id\tconsensus_id\tis_primary\tis_seed\tfeature_mass\tfeature_mz\tfeature_charge\tfeature_rt\tfeature_combined_score\tfeature_n_isotopes\tfeature_intensity_apex\n",
                "0\t0\t0\ttrue\ttrue\t1000.123456789\t501.069004861\t2\t10.123456789\t0.9\t3\t123456.789\n",
                "0\t1\t0\tfalse\tfalse\t1000.124\t501.069276466\t2\t10.2\t0.7\t2\t55\n",
                "1\t0\t0\ttrue\tfalse\t1000.125\t501.069776466\t2\t10.3\t0.8\t4\t65\n",
                "2\t0\t\tfalse\tfalse\t1100\t551.007276466\t2\t12\t0.4\t2\t75\n")).unwrap();
            std::fs::write(p.join("cells.tsv"),concat!(
                "consensus_id\trun_id\textraction_kind\tobservation_id\tseed_run_id\tseed_observation_id\tmass\tmz\tcharge\trt\tcombined_score\tn_contributing_runs\tlfq_intensity\tlfq_q_value\tlfq_is_mbr\n",
                "0\t0\ttarget\t0\t0\t0\t1000.123456789\t501.069004861\t2\t10.123456789\t0.9\t2\t123456.789\t0.001\tfalse\n",
                "0\t1\ttarget\t0\t0\t0\t1000.123456789\t501.069004861\t2\t10.123456789\t0.9\t2\t65\t0.001\tfalse\n",
                "0\t2\ttarget\t\t0\t0\t1000.123456789\t501.069004861\t2\t10.123456789\t0.9\t2\t5\t0.9\ttrue\n")).unwrap();
            Self(p)
        }
        fn open(&self) -> Bundle {
            let m = serde_json::json!({"schema_version":1,"n_consensus":1,"reference_run_id":0,
                "runs":[{"run_id":0,"name":"A","n_observations":2},{"run_id":1,"name":"B","n_observations":1},{"run_id":2,"name":"C","n_observations":1}],
                "observations":{"path":"observations.tsv","sha256":file_sha(&self.0.join("observations.tsv")).unwrap()},
                "cells":{"path":"cells.tsv","sha256":file_sha(&self.0.join("cells.tsv")).unwrap()}});
            std::fs::write(
                self.0.join("matrix_manifest.json"),
                serde_json::to_vec(&m).unwrap(),
            )
            .unwrap();
            Bundle::open(&self.0).unwrap()
        }
    }
    fn with_aligned_coordinates(f: &Fixture) {
        let p = f.0.join("observations.tsv");
        let text = std::fs::read_to_string(&p).unwrap();
        let mut out = String::new();
        for (i, line) in text.lines().enumerate() {
            if i == 0 {
                out.push_str(&format!("{line}\taligned_mass\taligned_rt\n"));
            } else {
                let cols: Vec<_> = line.split('\t').collect();
                out.push_str(&format!("{line}\t{}\t{}\n", cols[5], cols[8]));
            }
        }
        std::fs::write(p, out).unwrap();
    }

    fn release_fixture(version: u32) -> Fixture {
        let f = Fixture::new();
        with_aligned_coordinates(&f);
        let p = f.0.join("cells.tsv");
        let text = std::fs::read_to_string(&p).unwrap();
        let mut out = String::new();
        for (i, line) in text.lines().enumerate() {
            let columns = if i == 0 {
                "group_score\tgroup_qvalue"
            } else {
                "2.3\t0.03"
            };
            out.push_str(&format!("{line}\t{columns}"));
            if version >= 3 {
                out.push_str(if i == 0 {
                    "\tlfq_ownership_status\tlfq_owned_samples\tlfq_excluded_samples\tlfq_competing_consensus_id\tlfq_preceding_signal_fraction"
                } else { "\texclusive\t10\t0\t\t0" });
            }
            out.push('\n');
        }
        std::fs::write(p, out).unwrap();
        f
    }

    fn open_release(f: &Fixture, version: u32) -> Result<Bundle> {
        f.open(); // refresh table hashes after each fixture edit
        let path = f.0.join("matrix_manifest.json");
        let mut manifest: serde_json::Value =
            serde_json::from_slice(&std::fs::read(&path).unwrap()).unwrap();
        manifest["schema_version"] = version.into();
        manifest["config"] = serde_json::json!({"lfq":{"consensus":{"max_group_qvalue":0.05}}});
        std::fs::write(path, serde_json::to_vec(&manifest).unwrap()).unwrap();
        Bundle::open(&f.0)
    }

    #[test]
    fn release_schemas_preserve_originals_and_quantities() {
        for version in [2, 3] {
            let f = release_fixture(version);
            let old = f.open().load_joint(1).unwrap();
            let b = open_release(&f, version).unwrap();
            let new = b.load_joint(1).unwrap();
            assert_eq!(old.features, new.features);
            assert_eq!(old.joint_cell_evidence, new.joint_cell_evidence);
            assert_eq!(new.joint_evidence.unwrap()[0].len(), 3);
            assert_eq!(b.load_run(0).unwrap().features.len(), 2);
        }
        assert!(open_release(&release_fixture(3), 4).is_err());
    }

    #[test]
    fn release_group_confidence_replaces_only_the_removed_seed_floor() {
        let f = release_fixture(3);
        for name in ["cells.tsv", "observations.tsv"] {
            let p = f.0.join(name);
            let text = std::fs::read_to_string(&p)
                .unwrap()
                .replace("\t0.9\t", "\t0.6\t");
            std::fs::write(p, text).unwrap();
        }
        let legacy = f.open().load_joint(1).unwrap();
        assert!(legacy.joint_evidence.unwrap()[0]
            .iter()
            .all(|o| !o.bootstrap_eligible));
        let release = open_release(&f, 3).unwrap().load_joint(1).unwrap();
        assert!(release.joint_evidence.unwrap()[0]
            .iter()
            .all(|o| o.bootstrap_eligible));
        assert_eq!(legacy.features, release.features);
    }

    #[test]
    fn suppressed_lfq_retains_original_evidence_but_not_quantitative_support() {
        let f = release_fixture(3);
        let p = f.0.join("cells.tsv");
        let text = std::fs::read_to_string(&p)
            .unwrap()
            .replace("\t123456.789\t0.001\tfalse", "\t0\t0.001\tfalse")
            .replacen(
                "\texclusive\t10\t0\t\t0",
                "\tshared_signal_only\t0\t10\t\t1",
                1,
            );
        std::fs::write(p, text).unwrap();
        let m = open_release(&f, 3).unwrap().load_joint(1).unwrap();
        assert_eq!(m.features[0].intensities, vec![0.0, 65.0, 5.0]);
        assert_eq!(m.joint_evidence.unwrap()[0].len(), 3);
        assert_eq!(m.joint_cell_evidence.unwrap()[0][3], 1.0 / 3.0);
    }

    #[test]
    fn release_rejects_bad_group_confidence_and_ownership() {
        for (old, new) in [
            ("\t2.3\t0.03", "\t2.3\t0.06"), // beyond declared group threshold
            ("\t2.3\t0.03", "\tNaN\t0.03"),
            ("group_qvalue", "missing_group_qvalue"),
            ("lfq_competing_consensus_id", "missing_competitor_column"),
            (
                "\texclusive\t10\t0\t\t0",
                "\tambiguous_residual\t0\t10\t\t1",
            ), // nonzero suppressed intensity
            ("\texclusive\t10\t0\t\t0", "\texclusive\t10\t0\t1\t0"), // nonexistent competitor
            ("\texclusive\t10\t0\t\t0", "\texclusive\t10\t0\t\t1.1"),
        ] {
            let f = release_fixture(3);
            let p = f.0.join("cells.tsv");
            let text = std::fs::read_to_string(&p).unwrap().replace(old, new);
            std::fs::write(p, text).unwrap();
            assert!(
                open_release(&f, 3).unwrap().load_joint(1).is_err(),
                "accepted {new}"
            );
        }
        let f = release_fixture(3);
        let p = f.0.join("cells.tsv");
        let text = std::fs::read_to_string(&p)
            .unwrap()
            .replacen("\t2.3\t0.03", "\t2.4\t0.03", 1);
        std::fs::write(p, text).unwrap();
        assert!(open_release(&f, 3).unwrap().load_joint(1).is_err());
    }
    #[test]
    fn joint_loader_preserves_alternatives_and_quant_cells_without_mbr_evidence() {
        let f = Fixture::new();
        with_aligned_coordinates(&f);
        let m = f.open().load_joint(1).unwrap();
        assert_eq!(m.len(), 1);
        let g = &m.joint_evidence.as_ref().unwrap()[0];
        assert_eq!(g.len(), 3); // two originals in run A, one in B, no MBR observation
        assert_eq!(m.features[0].intensities, vec![123456.789, 65.0, 5.0]);
        assert_eq!(m.features[0].n_contributing_runs, 2);
        assert!(f.open().load_joint(3).unwrap().is_empty());
    }
    #[test]
    fn joint_loader_rejects_missing_aligned_coordinates() {
        let f = Fixture::new();
        assert!(f.open().load_joint(1).is_err());
    }
    #[test]
    fn joint_cells_use_extraction_q_without_adding_observations_or_decoy_cells() {
        let f = Fixture::new();
        with_aligned_coordinates(&f);
        let before = f.open().load_joint(1).unwrap();
        let values = before.joint_cell_evidence.as_ref().unwrap()[0];
        assert_eq!(values[0], 1.0);
        assert_eq!(values[1], 1.0 / 3.0);
        assert_eq!(values[2], 0.0);
        assert_eq!(values[3], 2.0 / 3.0);
        assert_eq!(&values[4..], &[0.001, 0.001, 1.0, 1.0, 0.9, 0.9, 1.0, 0.0]);
        // Synthetic extraction decoys have no target q and never enter summaries.
        let p = f.0.join("cells.tsv");
        let mut text = std::fs::read_to_string(&p).unwrap();
        text.push_str("0\t2\tdecoy\t\t0\t0\t1000.123456789\t501.069004861\t2\t10.123456789\t0.9\t2\t999999\t\ttrue\n");
        std::fs::write(p, text).unwrap();
        let after = f.open().load_joint(1).unwrap();
        assert_eq!(before.features, after.features);
        assert_eq!(before.joint_cell_evidence, after.joint_cell_evidence);
        assert_eq!(after.joint_evidence.as_ref().unwrap()[0].len(), 3);
        assert!(f
            .open()
            .load_joint(3)
            .unwrap()
            .joint_cell_evidence
            .unwrap()
            .is_empty());
    }
    #[test]
    fn duplicate_sources_do_not_inflate_cell_support() {
        let f = Fixture::new();
        with_aligned_coordinates(&f);
        let mut b = f.open();
        b.runs[0].sha256 = Some("same-source".into());
        b.runs[1].sha256 = Some("same-source".into());
        let m = b.load_joint(1).unwrap();
        let v = m.joint_cell_evidence.unwrap()[0];
        assert_eq!(v[0], 1.0);
        assert_eq!(v[1], 0.5);
        assert_eq!(v[3], 0.5);
    }
    #[test]
    fn weak_groups_remain_searchable_without_defining_bootstrap_calibration() {
        let f = Fixture::new();
        with_aligned_coordinates(&f);
        for name in ["observations.tsv", "cells.tsv"] {
            let p = f.0.join(name);
            let text = std::fs::read_to_string(&p)
                .unwrap()
                .replace("\t0.9\t", "\t0.7\t")
                .replace("\t0.8\t", "\t0.6\t");
            std::fs::write(p, text).unwrap();
        }
        let m = f.open().load_joint(1).unwrap();
        assert_eq!(m.len(), 1);
        assert!(m.joint_evidence.as_ref().unwrap()[0]
            .iter()
            .all(|o| !o.bootstrap_eligible));
    }
    #[test]
    fn joint_loader_rejects_member_charge_mismatch() {
        let f = Fixture::new();
        with_aligned_coordinates(&f);
        let p = f.0.join("observations.tsv");
        let text = std::fs::read_to_string(&p)
            .unwrap()
            .replace("501.069276466\t2", "501.069276466\t3");
        std::fs::write(p, text).unwrap();
        assert!(f.open().load_joint(1).is_err());
    }
    #[test]
    fn joint_matching_searches_member_windows_once_and_treats_decoys_equally() {
        let f = Fixture::new();
        with_aligned_coordinates(&f);
        let m = f.open().load_joint(1).unwrap();
        let db = crate::PeptideDb {
            n_targets: 2,
            n_decoys: 1,
            peptides: [(1000.123456789, false), (1000.125, false), (1000.125, true)]
                .iter()
                .map(|&(mass, decoy)| crate::TheoPeptide {
                    sequence: "PEPTIDE".into(),
                    mass,
                    decoy,
                    missed_cleavages: 0,
                    proteins: vec!["protein".into()],
                })
                .collect(),
        };
        let cfg = crate::MatchConfig {
            ppm_tolerance: 0.1,
            isotope_offsets: vec![0],
            ..Default::default()
        };
        let matches = crate::match_features_topk(&m, &db, &cfg, 50);
        assert_eq!(
            matches.iter().map(|p| p.peptide_idx).collect::<Vec<_>>(),
            vec![0, 1, 2]
        );
        let mut single = m.clone();
        single.joint_evidence = None;
        assert_eq!(crate::match_features_topk(&single, &db, &cfg, 50).len(), 1);
        let strict = crate::MatchConfig {
            min_isotopes: 4,
            ..cfg
        };
        let matches = crate::match_features_topk(&m, &db, &strict, 50);
        assert_eq!(
            matches.iter().map(|p| p.peptide_idx).collect::<Vec<_>>(),
            vec![1, 2]
        );
    }
    #[test]
    fn matching_quality_gates_also_restrict_joint_scoring_observations() {
        let f = Fixture::new();
        with_aligned_coordinates(&f);
        let mut m = f.open().load_joint(1).unwrap();
        let before = m.features.clone();
        crate::joint::restrict_observations(&mut m, &crate::MatchConfig::default());
        assert_eq!(before, m.features);
        let cells = m.joint_cell_evidence.clone();
        let intensities = m.features[0].intensities.clone();
        let cfg = crate::MatchConfig {
            min_isotopes: 4,
            ..Default::default()
        };
        crate::joint::restrict_observations(&mut m, &cfg);
        assert_eq!(m.features[0].n_isotopes, Some(4.0));
        assert_eq!(m.features[0].intensities, intensities);
        assert_eq!(m.joint_evidence.as_ref().unwrap()[0].len(), 1);
        assert_eq!(m.joint_cell_evidence, cells);
        let cfg = crate::MatchConfig {
            min_isotopes: 5,
            ..Default::default()
        };
        crate::joint::restrict_observations(&mut m, &cfg);
        assert!(m.is_empty());
        assert!(m.joint_cell_evidence.unwrap().is_empty());
    }

    #[test]
    fn original_quality_and_alternatives_survive_without_double_counting() {
        let f = Fixture::new();
        let b = f.open();
        let m = b.load_consensus(2, true, false).unwrap();
        assert_eq!(m.features[0].n_isotopes, Some(3.0));
        assert_eq!(m.features[0].intensities, vec![123456.789, 65.0, 5.0]);
        assert_eq!(m.features[0].n_contributing_runs, 2);
        let original = b.load_run(0).unwrap();
        assert_eq!(original.len(), 2);
        assert_eq!(original.features[1].n_isotopes, Some(2.0));
        assert_eq!(b.load_run(2).unwrap().features[0].mass, 1100.0); // filtered original remains available
        assert!(b.load_consensus(2, false, false).unwrap().features[0]
            .n_isotopes
            .is_none());
        let legacy = b.load_consensus(2, false, true).unwrap();
        assert_eq!(legacy.features[0].mass, 1000.123457);
        assert_eq!(legacy.features[0].intensities[0], 123457.0);
        assert!(b.load_consensus(3, true, false).unwrap().is_empty());
    }
    #[test]
    fn duplicate_or_missing_cells_are_rejected() {
        let f = Fixture::new();
        let path = f.0.join("cells.tsv");
        let text = std::fs::read_to_string(&path).unwrap();
        std::fs::write(&path, format!("{}{}\n", text, text.lines().nth(1).unwrap())).unwrap();
        assert!(f
            .open()
            .load_consensus(1, true, false)
            .unwrap_err()
            .to_string()
            .contains("duplicate"));
        std::fs::write(
            &path,
            text.lines().take(3).collect::<Vec<_>>().join("\n") + "\n",
        )
        .unwrap();
        assert!(f
            .open()
            .load_consensus(1, true, false)
            .unwrap_err()
            .to_string()
            .contains("missing"));
    }
    #[test]
    fn corrupt_hash_and_broken_member_links_are_rejected() {
        let f = Fixture::new();
        let _ = f.open();
        let path = f.0.join("observations.tsv");
        let text = std::fs::read_to_string(&path)
            .unwrap()
            .replace("0\t1\t0\tfalse", "0\t1\t0\ttrue");
        std::fs::write(path, text).unwrap();
        assert!(Bundle::open(&f.0)
            .err()
            .unwrap()
            .to_string()
            .contains("checksum"));
        assert!(f
            .open()
            .load_consensus(1, true, false)
            .unwrap_err()
            .to_string()
            .contains("primary"));
    }
    #[test]
    fn isotope_representative_uses_one_real_members_coordinates_and_quality() {
        let fixture = Fixture::new();
        let path = fixture.0.join("observations.tsv");
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(&path)
            .unwrap();
        let header = reader.headers().unwrap().clone();
        let rows: Vec<_> = reader.records().map(|r| r.unwrap()).collect();
        drop(reader);
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(&path)
            .unwrap();
        let mut names: Vec<_> = header.iter().map(str::to_owned).collect();
        names.extend([
            "aligned_mass".into(),
            "aligned_rt".into(),
            "aligned_im".into(),
        ]);
        writer.write_record(names).unwrap();
        let col = |n| header.iter().position(|h| h == n).unwrap();
        for (i, row) in rows.iter().enumerate() {
            let mut values: Vec<_> = row.iter().map(str::to_owned).collect();
            if i == 0 {
                values[col("feature_n_isotopes")] = "2".into();
            }
            values.extend([
                row[col("feature_mass")].to_owned(),
                row[col("feature_rt")].to_owned(),
                String::new(),
            ]);
            writer.write_record(values).unwrap();
        }
        writer.flush().unwrap();
        drop(writer);
        let bundle = fixture.open();
        let seed = bundle.load_consensus(2, true, false).unwrap();
        let actual = bundle
            .load_consensus_with_representative(2, true, false, true)
            .unwrap();
        assert_eq!(actual.features[0].mass, 1000.125);
        assert_eq!(actual.features[0].rt, 10.3);
        assert_eq!(actual.features[0].combined_score, 0.8);
        assert_eq!(actual.features[0].n_isotopes, Some(4.0));
        assert_eq!(actual.features[0].intensities, seed.features[0].intensities);
        assert_eq!(
            actual.features[0].n_contributing_runs,
            seed.features[0].n_contributing_runs
        );
    }

    #[test]
    fn negative_counts_are_rejected_before_row_filtering() {
        let fixture = Fixture::new();
        let path = fixture.0.join("observations.tsv");
        let text = std::fs::read_to_string(&path)
            .unwrap()
            .replace("0.9\t3\t123456.789", "0.9\t-1\t123456.789");
        std::fs::write(path, text).unwrap();
        assert!(fixture
            .open()
            .load_consensus(99, false, false)
            .unwrap_err()
            .to_string()
            .contains("feature count"));
    }
}
