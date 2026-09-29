//! Numerical primitives used by the pinned DirectMS1Quant implementation.
use anyhow::{bail, ensure, Result};
use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};
use nalgebra::{storage::Owned, DMatrix, DVector, Dyn};
use statrs::distribution::{Binomial, ContinuousCDF, DiscreteCDF, StudentsT};

pub fn mean(x: &[f64]) -> f64 {
    x.iter().sum::<f64>() / x.len() as f64
}
pub fn quantile(x: &[f64], p: f64) -> f64 {
    let mut v = x.to_vec();
    v.sort_by(f64::total_cmp);
    if v.is_empty() {
        return f64::NAN;
    }
    let pos = (v.len() - 1) as f64 * p;
    let i = pos.floor() as usize;
    v[i] + (v[pos.ceil() as usize] - v[i]) * (pos - i as f64)
}
pub fn median(x: &[f64]) -> f64 {
    quantile(x, 0.5)
}
pub fn ttest(a: &[f64], b: &[f64]) -> f64 {
    if a.len() < 2 || b.len() < 2 {
        return f64::NAN;
    }
    let ma = mean(a);
    let mb = mean(b);
    let df = (a.len() + b.len() - 2) as f64;
    let ss = a.iter().map(|x| (x - ma).powi(2)).sum::<f64>()
        + b.iter().map(|x| (x - mb).powi(2)).sum::<f64>();
    let se = (ss / df * (1. / a.len() as f64 + 1. / b.len() as f64)).sqrt();
    if se == 0. {
        return if ma == mb { f64::NAN } else { 0. };
    }
    2. * StudentsT::new(0., 1., df)
        .unwrap()
        .sf(((ma - mb) / se).abs())
}
pub fn binomial_score(v: u64, n: u64, p: f64) -> f64 {
    if v <= 2 || n == 0 {
        return 0.;
    }
    let prob = Binomial::new(p, n).unwrap().sf(v - 1);
    if prob == 0. {
        20.
    } else {
        -prob.log10()
    }
}

struct Gaussian {
    x: Vec<f64>,
    y: Vec<f64>,
    params: DVector<f64>,
}
impl Gaussian {
    fn eval(&self, p: &DVector<f64>) -> DVector<f64> {
        DVector::from_iterator(
            self.x.len(),
            self.x
                .iter()
                .zip(&self.y)
                .map(|(x, y)| p[0] * (-(x - p[1]).powi(2) / (2. * p[2] * p[2])).exp() + p[3] - y),
        )
    }
}
impl LeastSquaresProblem<f64, Dyn, Dyn> for Gaussian {
    type ParameterStorage = Owned<f64, Dyn>;
    type ResidualStorage = Owned<f64, Dyn>;
    type JacobianStorage = Owned<f64, Dyn, Dyn>;
    fn set_params(&mut self, p: &DVector<f64>) {
        self.params.copy_from(p);
    }
    fn params(&self) -> DVector<f64> {
        self.params.clone()
    }
    fn residuals(&self) -> Option<DVector<f64>> {
        Some(self.eval(&self.params))
    }
    fn jacobian(&self) -> Option<DMatrix<f64>> {
        let base = self.eval(&self.params);
        let mut j = DMatrix::zeros(self.x.len(), 4);
        for c in 0..4 {
            let mut p = self.params.clone();
            let h = f64::EPSILON.sqrt() * if p[c] == 0. { 1. } else { p[c].abs() };
            p[c] += h;
            let diff = (self.eval(&p) - &base) / h;
            j.set_column(c, &diff);
        }
        Some(j)
    }
}
fn fit(values: &[f64], width: f64) -> Result<(f64, f64, bool)> {
    let lo = values.iter().copied().fold(f64::INFINITY, f64::min);
    let hi = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let count = ((hi - lo) / width).ceil() as usize;
    ensure!(
        (6..=100_000).contains(&count),
        "insufficient or excessive histogram bins for reference calibration"
    );
    // Match numpy.arange's floating-point step and histogram's inclusive final edge.
    let step = (lo + width) - lo;
    let edges: Vec<f64> = (0..count).map(|i| lo + i as f64 * step).collect();
    let mut y = vec![0.; count - 1];
    for &v in values {
        if v < edges[0] || v > edges[count - 1] {
            continue;
        }
        let bin = edges
            .partition_point(|x| *x <= v)
            .saturating_sub(1)
            .min(count - 2);
        y[bin] += 1.;
    }
    let problem = Gaussian {
        x: edges[..count - 1].iter().map(|x| x + width).collect(),
        y,
        params: DVector::from_vec(vec![1., median(values), 1., 1.]),
    };
    let (solution, report) = LevenbergMarquardt::new()
        .with_ftol(1.49012e-8)
        .with_xtol(1.49012e-8)
        .with_gtol(0.)
        .with_patience(200)
        .minimize(problem);
    ensure!(
        report.termination.was_successful(),
        "Gaussian calibration failed: {:?}",
        report.termination
    );
    let p = solution.params();
    ensure!(
        p.iter().all(|v| v.is_finite()) && p[2] != 0.,
        "nonfinite Gaussian fit"
    );
    // scipy's covariance is finite only for full-rank fits with residual degrees of freedom.
    let j = solution.jacobian().unwrap();
    let finite_cov = solution.x.len() > 4 && (j.transpose() * j).try_inverse().is_some();
    Ok((p[1], p[2].abs(), finite_cov))
}
pub fn calibrate(values: &[f64]) -> Result<(f64, f64)> {
    ensure!(
        !values.is_empty() && values.iter().all(|x| x.is_finite()),
        "no finite fold changes for calibration"
    );
    let first = fit(values, 0.05).and_then(|a| {
        fit(values, 0.1).map(|b| {
            if b.2 && b.0.abs() <= a.0.abs() / 10. {
                b
            } else {
                a
            }
        })
    });
    let x = first.or_else(|_| fit(values, 0.3))?;
    if !x.0.is_finite() || !x.1.is_finite() {
        bail!("invalid fold-change calibration");
    }
    Ok((x.0, x.1))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn exact_tail() {
        assert!((binomial_score(3, 3, 0.1) - 3.).abs() < 1e-12);
        assert_eq!(binomial_score(2, 8, 0.), 0.);
        assert_eq!(binomial_score(3, 8, 0.), 20.);
    }
    #[test]
    fn ttest_symmetry() {
        let a = [1., 2., 3.];
        let b = [3., 4., 5.];
        assert!((ttest(&a, &b) - 0.07048399691021993).abs() < 1e-12);
        assert_eq!(ttest(&a, &b), ttest(&b, &a));
    }
    #[test]
    fn median_even() {
        assert_eq!(median(&[4., 1., 2., 3.]), 2.5);
    }
}
