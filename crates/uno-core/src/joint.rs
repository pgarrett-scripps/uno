//! Candidate-specific evidence from the original observations of an aligned row.
//! Measurements are features of one discriminant, NOT independent probabilities.
use crate::input::Feature;

#[derive(Debug, Clone)]
pub struct Observation {
    /// Group clears KOTH's original seed floor; weak replicated groups can be
    /// searched without defining the initial mass/RT calibration.
    pub bootstrap_eligible: bool,
    /// Duplicate source-table hashes share an evidence run ID at ingestion.
    pub run_id: usize,
    pub observation_id: usize,
    /// Original quality, aligned mass/RT/IM. Intensities contain original apex only.
    pub feature: Feature,
}

#[derive(Debug, Clone)]
pub struct Evidence {
    pub best: usize,
    pub mass_error: f64,
    pub rt_error: f64,
    pub deviation: f64,
    /// Candidate-specific agreement/quality summaries (see NAMES).
    pub values: [f32; 14],
}

pub struct Match {
    pub best: usize,
    pub mass_error: f64,
    pub rt_error: f64,
    pub deviation: f64,
}

pub const NAMES: [&str; 14] = [
    "direct_runs",
    "compatible_runs",
    "compatible_fraction",
    "second_deviation",
    "median_deviation",
    "median_mass_error",
    "mass_mad",
    "median_rt_error",
    "rt_mad",
    "median_quality",
    "median_isotopes",
    "median_scans",
    "best_quality",
    "has_second_run",
];

/// Apply user-specified matching quality gates to both matching AND scoring.
/// Otherwise a candidate found on an eligible member could be scored using a
/// different observation that the user's isotope/quality settings excluded.
pub fn restrict_observations(matrix: &mut crate::FeatureMatrix, cfg: &crate::MatchConfig) {
    let Some(groups) = matrix.joint_evidence.take() else {
        return;
    };
    let cells = matrix.joint_cell_evidence.take();
    let mut retained_cells = Vec::new();
    let mut features = Vec::new();
    let mut retained = Vec::new();
    for (i, (mut f, mut group)) in std::mem::take(&mut matrix.features)
        .into_iter()
        .zip(groups)
        .enumerate()
    {
        let before = group.len();
        group.retain(|o| crate::matching::passes(&o.feature, cfg));
        if group.is_empty() {
            continue;
        }
        if group.len() != before {
            let representative = group
                .iter()
                .max_by(|a, b| {
                    a.feature
                        .n_isotopes
                        .is_some_and(|n| n >= 3.0)
                        .cmp(&b.feature.n_isotopes.is_some_and(|n| n >= 3.0))
                        .then(
                            a.feature
                                .combined_score
                                .total_cmp(&b.feature.combined_score),
                        )
                        .then_with(|| {
                            (b.run_id, b.observation_id).cmp(&(a.run_id, a.observation_id))
                        })
                })
                .unwrap();
            let intensities = f.intensities;
            f = representative.feature.clone();
            f.intensities = intensities;
            f.n_contributing_runs = group
                .iter()
                .map(|o| o.run_id)
                .collect::<std::collections::HashSet<_>>()
                .len() as u32;
        }
        if let Some(cells) = &cells {
            retained_cells.push(cells[i]);
        }
        features.push(f);
        retained.push(group);
    }
    matrix.features = features;
    matrix.joint_evidence = Some(retained);
    matrix.joint_cell_evidence = cells.map(|_| retained_cells);
}

fn med(v: &mut [f64]) -> f64 {
    v.sort_unstable_by(f64::total_cmp);
    let n = v.len();
    if n.is_multiple_of(2) {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    } else {
        v[n / 2]
    }
}

/// Best *simultaneous* mass/RT match, at most one observation per source run.
/// A strong singleton survives; absent runs and MBR cells do not enter the model.
/// Alternatives can replace a poorer peak, but cannot create additional support.
pub fn evidence(
    observations: &[Observation],
    peptide_mass: f64,
    predicted_rt: f64,
    mass_shift: f64,
    mass_sigma: f64,
    rt_sigma: f64,
    mass_gate: f64,
) -> Option<Evidence> {
    let mut by_run: Vec<(usize, usize, f64, f64, f64)> = Vec::new();
    for (i, o) in observations.iter().enumerate() {
        let mass = (peptide_mass - o.feature.mass) / o.feature.mass * 1e6;
        let rt = o.feature.rt - predicted_rt;
        let mz = (mass - mass_shift) / mass_sigma;
        let dev = mz * mz + (rt / rt_sigma).powi(2);
        // Group observations are stored by run at load time. This remains correct
        // for arbitrary order in callers and tests, and small cohorts keep it cheap.
        if let Some(prior) = by_run.iter_mut().find(|p| p.0 == o.run_id) {
            let inside = (mass - mass_shift).abs() <= mass_gate;
            let prior_inside = (prior.3 - mass_shift).abs() <= mass_gate;
            let better_dev = dev < prior.2
                || (dev == prior.2
                    && o.feature.combined_score > observations[prior.1].feature.combined_score);
            if (inside && !prior_inside) || (inside == prior_inside && better_dev) {
                *prior = (o.run_id, i, dev, mass, rt);
            }
        } else {
            by_run.push((o.run_id, i, dev, mass, rt));
        }
    }
    by_run.sort_unstable_by(|a, b| {
        ((b.3 - mass_shift).abs() <= mass_gate)
            .cmp(&((a.3 - mass_shift).abs() <= mass_gate))
            .then(a.2.total_cmp(&b.2))
            .then_with(|| {
                observations[b.1]
                    .feature
                    .combined_score
                    .total_cmp(&observations[a.1].feature.combined_score)
            })
            .then(a.0.cmp(&b.0))
    });
    let best = by_run
        .iter()
        .find(|p| (p.3 - mass_shift).abs() <= mass_gate)?;
    let n = by_run.len();
    let compatible = by_run
        .iter()
        .filter(|p| p.2 <= 9.0 && (p.3 - mass_shift).abs() <= mass_gate)
        .count();
    let (mut mass, mut rt, mut devs, mut quality, mut isotopes, mut scans) = (
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    );
    for p in &by_run {
        let f = &observations[p.1].feature;
        mass.push(p.3 - mass_shift);
        rt.push(p.4);
        devs.push(p.2.min(100.0));
        quality.push(f.combined_score);
        isotopes.push(f.n_isotopes.unwrap_or(0.0));
        scans.push(f.n_scans.unwrap_or(0.0));
    }
    let mm = med(&mut mass);
    let rm = med(&mut rt);
    for x in &mut mass {
        *x = (*x - mm).abs();
    }
    for x in &mut rt {
        *x = (*x - rm).abs();
    }
    Some(Evidence {
        best: best.1,
        mass_error: best.3,
        rt_error: best.4,
        deviation: best.2,
        values: [
            (n as f64).ln_1p() as f32,
            (compatible as f64).ln_1p() as f32,
            compatible as f32 / n as f32,
            by_run.get(1).map_or(100.0, |p| p.2.min(100.0)) as f32,
            med(&mut devs) as f32,
            mm as f32,
            med(&mut mass) as f32,
            rm as f32,
            med(&mut rt) as f32,
            med(&mut quality) as f32,
            med(&mut isotopes).ln_1p() as f32,
            med(&mut scans).ln_1p() as f32,
            observations[best.1].feature.combined_score as f32,
            (n > 1) as u8 as f32,
        ],
    })
}

/// Joint deviation selection does not multiply replicate likelihoods. The ML
/// stage learns how agreement affects ranking using held-out peptide folds.
pub fn best_deviation(
    observations: &[Observation],
    peptide_mass: f64,
    predicted_rt: f64,
    mass_shift: f64,
    mass_sigma: f64,
    rt_sigma: f64,
    mass_gate: f64,
) -> Option<f64> {
    best_match(
        observations,
        peptide_mass,
        predicted_rt,
        mass_shift,
        mass_sigma,
        rt_sigma,
        mass_gate,
    )
    .map(|m| m.deviation)
}

pub fn best_match(
    observations: &[Observation],
    peptide_mass: f64,
    predicted_rt: f64,
    mass_shift: f64,
    mass_sigma: f64,
    rt_sigma: f64,
    mass_gate: f64,
) -> Option<Match> {
    observations
        .iter()
        .enumerate()
        .filter_map(|(i, o)| {
            let err = (peptide_mass - o.feature.mass) / o.feature.mass * 1e6;
            let rt = o.feature.rt - predicted_rt;
            ((err - mass_shift).abs() <= mass_gate).then(|| Match {
                best: i,
                mass_error: err,
                rt_error: rt,
                deviation: ((err - mass_shift) / mass_sigma).powi(2) + (rt / rt_sigma).powi(2),
            })
        })
        .min_by(|a, b| {
            a.deviation
                .total_cmp(&b.deviation)
                .then_with(|| {
                    observations[b.best]
                        .feature
                        .combined_score
                        .total_cmp(&observations[a.best].feature.combined_score)
                })
                .then_with(|| {
                    (
                        observations[a.best].run_id,
                        observations[a.best].observation_id,
                    )
                        .cmp(&(
                            observations[b.best].run_id,
                            observations[b.best].observation_id,
                        ))
                })
        })
}

/// One mass-compatible observed RT per candidate for cohort recalibration.
/// RT is taken from an actual member, never from an unrelated representative.
pub fn calibration_rt(
    observations: &[Observation],
    peptide_mass: f64,
    predicted_rt: f64,
    mass_shift: f64,
    mass_gate: f64,
) -> Option<f64> {
    observations
        .iter()
        .filter(|o| {
            ((peptide_mass - o.feature.mass) / o.feature.mass * 1e6 - mass_shift).abs() <= mass_gate
        })
        .min_by(|a, b| {
            (a.feature.rt - predicted_rt)
                .abs()
                .total_cmp(&(b.feature.rt - predicted_rt).abs())
                .then(a.feature.rt.total_cmp(&b.feature.rt))
        })
        .map(|o| o.feature.rt)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn obs(run: usize, mass: f64, rt: f64) -> Observation {
        Observation {
            bootstrap_eligible: true,
            run_id: run,
            observation_id: 0,
            feature: Feature {
                mass,
                mz: mass / 2.0 + 1.007276466621,
                charge: 2,
                rt,
                im: None,
                combined_score: 0.9,
                n_contributing_runs: 1,
                n_isotopes: Some(3.0),
                n_scans: Some(10.0),
                isoerror: None,
                isoerror2: None,
                rt_width: None,
                cosine_score: Some(0.9),
                isotope_score: Some(0.9),
                ppm_error: None,
                neutron_offset: None,
                intensities: vec![100.0],
            },
        }
    }
    fn ev(o: &[Observation]) -> Evidence {
        evidence(o, 1000.0, 10.0, 0.0, 1.0, 1.0, 3.0).unwrap()
    }
    #[test]
    fn alternatives_and_duplicate_sources_do_not_inflate_support() {
        let a = obs(0, 1000.0, 10.0);
        let b = obs(1, 1000.0005, 10.1);
        let base = ev(&[a.clone(), b.clone()]);
        let copies = ev(&[a.clone(), obs(0, 1000.0, 12.0), b, a]);
        assert_eq!(base.values, copies.values);
        assert_eq!(base.deviation, copies.deviation);
    }
    #[test]
    fn candidate_evidence_changes_when_an_additional_run_disagrees() {
        let good = ev(&[obs(0, 1000.0, 10.0), obs(1, 1000.0005, 10.1)]);
        let bad = ev(&[obs(0, 1000.0, 10.0), obs(1, 1000.006, 13.0)]);
        assert_eq!(good.deviation, bad.deviation); // strong first run retained
        assert!(good.values[2] > bad.values[2]);
        assert!(good.values[4] < bad.values[4]);
    }
    #[test]
    fn mass_and_rt_must_come_from_the_same_observation() {
        let o = [obs(0, 1000.0, 20.0), obs(0, 1000.01, 10.0)];
        assert!(best_deviation(&o, 1000.0, 10.0, 0.0, 1.0, 1.0, 3.0).unwrap() > 90.0);
        assert!(ev(&o).deviation > 90.0);
    }
    #[test]
    fn singleton_is_valid_and_order_does_not_change_summary() {
        assert_eq!(ev(&[obs(0, 1000.0, 10.0)]).values[13], 0.0);
        let a = obs(0, 1000.0, 10.0);
        let b = obs(1, 1000.0005, 10.1);
        assert_eq!(ev(&[a.clone(), b.clone()]).values, ev(&[b, a]).values);
    }
    #[test]
    fn calibration_uses_the_mass_compatible_original_observation() {
        let o = [obs(0, 1000.01, 10.0), obs(1, 1000.0, 10.2)];
        assert_eq!(calibration_rt(&o, 1000.0, 10.0, 0.0, 3.0), Some(10.2));
        assert_eq!(calibration_rt(&o[..1], 1000.0, 10.0, 0.0, 3.0), None);
    }
}
