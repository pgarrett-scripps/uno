//! Experimental label-blind, feature-balanced RT bootstrap. Parameters are frozen
//! in benchmark/results/stability_round2/rt_candidate/PROTOCOL.md.
use crate::Calibration;
use std::collections::{BTreeMap, BTreeSet};

const NX: usize = 32;
const NY: usize = 128;

#[derive(Clone, Copy, Debug)]
pub(crate) struct Evidence<'a> {
    pub position: usize,
    pub feature: u32,
    pub sequence: &'a str,
    pub x: f64,
    pub y: f64,
}

pub(crate) fn observation<'a>(position: usize, p: &crate::Pfm, f: &crate::Feature,
    peptide: &'a crate::TheoPeptide, rt: Option<f64>, shift: f64, gate: f64, has_iso: bool) -> Option<Evidence<'a>> {
    if p.isotope_offset != 0 || !p.mass_error_ppm.is_finite() || (p.mass_error_ppm-shift).abs() > gate
        || (has_iso && f.n_isotopes.unwrap_or(0.0) < 3.0) { return None; }
    rt.map(|x| Evidence { position,feature:p.feature_idx,sequence:peptide.sequence.as_ref(),x,y:f.rt })
}

pub(crate) struct BroadFit {
    pub calibration: Calibration,
    pub sigma_floor: f64,
    pub clean: Vec<usize>,
    #[cfg(test)]
    pub prior_only: bool,
    #[cfg(test)]
    pub histogram_weight: f64,
}

fn quantile(v: &mut [f64], q: f64) -> f64 {
    v.sort_by(f64::total_cmp);
    v[(q * (v.len() - 1) as f64).floor() as usize]
}

/// Route all later RT stages through this function. The fallback fitter is lazy:
/// an accepted broad map is never updated from a confident-protein list.
pub(crate) fn stage_map(frozen: Option<&Calibration>, fallback: impl FnOnce() -> Calibration) -> Calibration {
    frozen.cloned().unwrap_or_else(fallback)
}

pub(crate) fn fit(input: &[Evidence<'_>]) -> Option<BroadFit> {
    let mut rows: Vec<_> = input.iter().copied().filter(|r| r.x.is_finite() && r.y.is_finite()).collect();
    rows.sort_by(|a,b| a.feature.cmp(&b.feature).then(a.sequence.cmp(b.sequence))
        .then(a.x.total_cmp(&b.x)).then(a.y.total_cmp(&b.y)).then(a.position.cmp(&b.position)));
    rows.dedup_by(|a,b| a.feature == b.feature && a.sequence == b.sequence);
    let mut features = BTreeMap::new();
    let mut sequences = BTreeMap::new();
    for r in &rows {
        features.entry(r.feature).or_insert(r.y);
        sequences.entry(r.sequence).or_insert(r.x);
    }
    if features.len() < 500 || sequences.len() < 500 {
        log::info!("broad-RT fallback: insufficient evidence features={} sequences={}", features.len(), sequences.len());
        return None;
    }
    let mut xs: Vec<_> = sequences.values().copied().collect();
    let mut ys: Vec<_> = features.values().copied().collect();
    let x0 = quantile(&mut xs, 0.05);
    let x1 = quantile(&mut xs, 0.95);
    let y0 = quantile(&mut ys, 0.05);
    let y1 = quantile(&mut ys, 0.95);
    let dx = x1-x0;
    let dy = y1-y0;
    if !dx.is_finite() || !dy.is_finite() || dx <= 1e-9 || dy <= 1e-9 {
        log::info!("broad-RT fallback: degenerate quantile span x={} y={}", dx,dy);
        return None;
    }
    let mut counts = BTreeMap::<u32,usize>::new();
    for r in rows.iter().filter(|r| r.x >= x0 && r.x <= x1 && r.y >= y0 && r.y <= y1) {
        *counts.entry(r.feature).or_default() += 1;
    }
    // Canonical value order also keeps floating-point reductions independent of
    // feature row indices and input ordering.
    let mut weighted: Vec<_> = rows.iter().filter(|r| r.x >= x0 && r.x <= x1 && r.y >= y0 && r.y <= y1)
        .map(|r| (*r, 1.0 / counts[&r.feature] as f64)).collect();
    weighted.sort_by(|(a,wa),(b,wb)| a.x.total_cmp(&b.x).then(a.y.total_cmp(&b.y))
        .then(a.sequence.cmp(b.sequence)).then(wa.total_cmp(wb)));
    let mut hist = [[0.0; NY]; NX];
    for (r,w) in &weighted {
        let ix = ((((r.x-x0)/dx)*NX as f64) as usize).min(NX-1);
        let iy = ((((r.y-y0)/dy)*NY as f64) as usize).min(NY-1);
        hist[ix][iy] += w;
    }
    let hx: Vec<f64> = hist.iter().map(|row| row.iter().sum()).collect();
    let hy: Vec<f64> = (0..NY).map(|y| (0..NX).map(|x| hist[x][y]).sum()).collect();
    let total: f64 = hx.iter().sum();
    let strata: Vec<f64> = hx.chunks(4).map(|v| v.iter().sum()).collect();
    if strata.iter().any(|&w| w < 20.0) {
        log::info!("broad-RT fallback: insufficient x coverage strata={:?}", strata);
        return None;
    }
    let corridor = 0.5_f64.max(0.02*dy);
    let halfwidth = corridor/dy;
    #[derive(Clone)]
    struct Grid { slope: f64, intercept: f64, score: f64, excess: f64, supported: usize }
    let mut grid = Vec::new();
    for si in 0..=20 {
        for bi in -25..=25 {
            let slope = 0.5+si as f64*0.05;
            let intercept = bi as f64*0.01;
            let mut excess = [0.0;8];
            for x in 0..NX {
                let pred = slope*(x as f64+0.5)/NX as f64+intercept;
                for y in 0..NY {
                    if ((y as f64+0.5)/NY as f64-pred).abs() <= halfwidth {
                        excess[x/4] += hist[x][y]-hx[x]*hy[y]/total;
                    }
                }
            }
            let score = excess.iter().zip(&strata).map(|(e,n)| e/n).sum::<f64>()/8.0;
            let supported = excess.iter().zip(&strata).filter(|(e,n)| **e >= 0.01 * **n).count();
            grid.push(Grid { slope,intercept,score,excess:excess.iter().sum(),supported });
        }
    }
    grid.sort_by(|a,b| b.score.total_cmp(&a.score).then_with(|| {
        let distance = |g: &Grid| (g.slope-1.0).powi(2)+g.intercept.powi(2);
        distance(a).total_cmp(&distance(b))
    }).then(a.slope.total_cmp(&b.slope)).then(a.intercept.total_cmp(&b.intercept)));
    let best = &grid[0];
    let accepted = best.supported >= 6 && best.excess >= 20.0 && best.score > 0.0;
    let (slope, intercept, blend, sigma_floor) = if accepted {
        let blend = (best.excess/100.0).clamp(0.0,1.0);
        let ambiguity: Vec<_> = grid.iter().filter(|g| g.score >= 0.9*best.score).collect();
        let mut spread: f64 = 0.0;
        for x in [0.0,0.5,1.0] {
            let mut predictions: Vec<_> = ambiguity.iter().map(|g| (g.slope*x+g.intercept)*dy+y0).collect();
            let lo = quantile(&mut predictions,0.1);
            let hi = quantile(&mut predictions,0.9);
            spread = spread.max(hi-lo);
        }
        (1.0+blend*(best.slope-1.0),blend*best.intercept,blend,0.5_f64.max(corridor/2.0).max(spread/2.0).min(3.5))
    } else { (1.0,0.0,0.0,3.5) };
    let physical_slope = slope*dy/dx;
    let physical_intercept = y0+intercept*dy-physical_slope*x0;
    let calibration = Calibration::affine(physical_slope, physical_intercept, xs[0], xs[xs.len()-1]);
    let clean: BTreeSet<_> = rows.iter().filter(|r| (r.y-calibration.apply(r.x)).abs() <= corridor)
        .map(|r| r.position).collect();
    log::info!("broad-RT fit: mode={} features={} sequences={} interior_features={} weight={:.6} x0={:.6} x1={:.6} y0={:.6} y1={:.6} grid_slope={:.6} grid_intercept={:.6} score={:.6} excess={:.6} supported_strata={} blend={:.6} slope={:.6} intercept={:.6} sigma_floor={:.6} clean={} frozen=true",
        if accepted {"evidence"} else {"prior"}, features.len(), sequences.len(), counts.len(),total,x0,x1,y0,y1,
        best.slope,best.intercept,best.score,best.excess,best.supported,blend,physical_slope,physical_intercept,sigma_floor,clean.len());
    Some(BroadFit { calibration,sigma_floor,clean:clean.into_iter().collect(),
        #[cfg(test)] prior_only:!accepted,
        #[cfg(test)] histogram_weight:total })
}

#[cfg(test)]
mod tests {
    use super::*;
    type Row = (u32,String,f64,f64);
    fn input(rows: &[Row]) -> Vec<Evidence<'_>> {
        rows.iter().enumerate().map(|(position,(feature,sequence,x,y))|
            Evidence {position,feature:*feature,sequence,x:*x,y:*y}).collect()
    }
    fn signal() -> Vec<Row> {
        let mut rows=Vec::new();
        for f in 0..1600u32 {
            let x=(f as f64+0.5)/16.0;
            let y=2.0*x+10.0+((f%7) as f64-3.0)*0.025;
            rows.push((f,format!("true_{f}"),x,y));
            for j in 0..7u32 {
                let noise=((f.wrapping_mul(613)+j*239)%1600) as f64/16.0;
                rows.push((f,format!("noise_{f}_{j}"),noise,y));
            }
        }
        rows
    }
    #[test]
    fn independent_background_uses_global_prior() {
        // Every observed feature has the same predicted-RT distribution, so
        // the contingency table is its product-marginal null by construction.
        let mut rows=Vec::new();
        for f in 0..640u32 {
            for x in 0..32u32 {
                rows.push((f,format!("seq_{f}_{x}"),x as f64+0.5,(f/5) as f64+0.5));
            }
        }
        let result=fit(&input(&rows)).unwrap();
        assert!(result.prior_only);
        assert_eq!(result.sigma_floor,3.5);
        assert!(result.calibration.apply(25.0)>result.calibration.apply(5.0));
    }
    #[test]
    fn contaminated_signal_recovers_positive_affine_shape() {
        let rows=signal();
        let result=fit(&input(&rows)).unwrap();
        assert!(!result.prior_only,"12.5% coherent signal should exceed independent background");
        for x in [10.0,30.0,50.0,70.0,90.0] {
            assert!((result.calibration.apply(x)-(2.0*x+10.0)).abs()<4.0);
        }
        assert!(result.calibration.apply(90.0)>result.calibration.apply(10.0));
        assert!((0.5..=3.5).contains(&result.sigma_floor));
    }
    #[test]
    fn permutation_and_duplicate_pfms_do_not_change_fit() {
        let rows=signal();
        let a=fit(&input(&rows)).unwrap();
        let mut permuted=rows.clone();
        permuted.reverse();
        permuted.extend(rows.iter().step_by(3).cloned());
        // Renumber feature indices as a feature-file permutation would.
        for r in &mut permuted { r.0=1599-r.0; }
        let b=fit(&input(&permuted)).unwrap();
        for x in [0.0,10.0,50.0,90.0,100.0] {
            assert_eq!(a.calibration.apply(x).to_bits(),b.calibration.apply(x).to_bits());
        }
        assert_eq!(a.sigma_floor.to_bits(),b.sigma_floor.to_bits());
        assert_eq!(a.histogram_weight.to_bits(),b.histogram_weight.to_bits());
    }
    #[test]
    fn feature_mass_stays_one_despite_many_candidate_sequences() {
        let rows=signal();
        let a=fit(&input(&rows)).unwrap();
        let mut extra=rows.clone();
        // Add many different sequence identities with the same predictions to
        // every feature: this changes candidate multiplicity, not feature mass.
        for (f,s,x,y) in &rows {
            for copy in 0..3 { extra.push((*f,format!("{s}_copy{copy}"),*x,*y)); }
        }
        let b=fit(&input(&extra)).unwrap();
        assert!((a.histogram_weight-b.histogram_weight).abs()<1e-8);
        assert!((a.histogram_weight-a.histogram_weight.round()).abs()<1e-8);
    }
    #[test]
    fn insufficient_evidence_and_degenerate_span_fall_back() {
        assert!(fit(&input(&signal()[..100])).is_none());
        let rows: Vec<_>=(0..600).map(|f|(f,format!("s{f}"),1.0,f as f64)).collect();
        assert!(fit(&input(&rows)).is_none());
        // Sufficient identities and nonzero spans cannot substitute for support
        // across the gradient: two narrow islands leave interior strata empty.
        let islands: Vec<_>=(0..600).map(|f| (f,format!("s{f}"),
            if f < 300 { f as f64*1e-6 } else { 100.0+f as f64*1e-6 },f as f64)).collect();
        assert!(fit(&input(&islands)).is_none());
    }
    #[test]
    fn frozen_stage_never_invokes_confident_protein_fitter() {
        let fixed=Calibration::affine(2.0,10.0,0.0,100.0);
        for _ in 2..=4 {
            let result=stage_map(Some(&fixed),|| panic!("confidence-derived RT fitter must not run"));
            assert_eq!(result.apply(25.0),60.0);
        }
        let fallback=stage_map(None,|| Calibration::affine(3.0,1.0,0.0,100.0));
        assert_eq!(fallback.apply(25.0),76.0);
    }
    #[test]
    fn eligibility_ignores_labels_and_preserves_missing_prediction_positions() {
        let mut peptide=crate::TheoPeptide { sequence:"PEPTIDEK".into(),mass:900.0,decoy:false,
            missed_cleavages:0,proteins:vec!["P_HUMAN".into()] };
        let mut p=crate::Pfm {mass_error_ppm:0.2,feature_idx:7,peptide_idx:0,isotope_offset:0,decoy:false};
        let f=crate::Feature {mass:900.0,mz:451.0,charge:2,rt:42.0,im:None,combined_score:1.0,
            n_contributing_runs:1,n_isotopes:Some(3.0),n_scans:None,isoerror:None,isoerror2:None,
            rt_width:None,cosine_score:None,isotope_score:None,ppm_error:None,neutron_offset:None,intensities:vec![1.0]};
        let before=observation(29,&p,&f,&peptide,Some(21.0),0.0,1.0,true).unwrap();
        let original=(before.position,before.feature,before.sequence.to_string(),before.x,before.y);
        peptide.decoy=true;
        peptide.proteins=vec!["DECOY_P_YEAST".into()];
        p.decoy=true;
        let after=observation(29,&p,&f,&peptide,Some(21.0),0.0,1.0,true).unwrap();
        assert_eq!(original,(after.position,after.feature,after.sequence.to_string(),after.x,after.y));
        assert!(observation(28,&p,&f,&peptide,None,0.0,1.0,true).is_none());
        assert_eq!(after.position,29);
        p.isotope_offset=1;
        assert!(observation(29,&p,&f,&peptide,Some(21.0),0.0,1.0,true).is_none());
    }
}
