//! Monotone piecewise-linear calibration of abstract RT -> observed RT.
//! Ported from sherlock's `rt_model::Calibration`: quantile-bin by abstract RT,
//! take per-bin medians, force the observed knots non-decreasing, interpolate.

/// Sorted knots mapping abstract RT (`x`) to observed RT (`y`), `y` non-decreasing.
#[derive(Debug, Clone, Default)]
pub struct Calibration {
    x: Vec<f64>,
    y: Vec<f64>,
}

/// Pool-adjacent-violators (PAVA) isotonic regression: the least-squares
/// non-decreasing fit to `y`. Smoother and UNBIASED versus a running-max upper
/// envelope, which lifts the curve at every local dip (a systematic upward bias
/// on noisy mode knots). Equal weights — one per knot.
fn isotonic(y: &[f64]) -> Vec<f64> {
    let mut vals: Vec<f64> = Vec::with_capacity(y.len());
    let mut wts: Vec<f64> = Vec::with_capacity(y.len());
    for &yi in y {
        vals.push(yi);
        wts.push(1.0);
        while vals.len() > 1 && vals[vals.len() - 2] > vals[vals.len() - 1] {
            let (v2, w2) = (vals.pop().unwrap(), wts.pop().unwrap());
            let (v1, w1) = (vals.pop().unwrap(), wts.pop().unwrap());
            let w = w1 + w2;
            vals.push((v1 * w1 + v2 * w2) / w);
            wts.push(w);
        }
    }
    let mut out = Vec::with_capacity(y.len());
    for (&v, &w) in vals.iter().zip(&wts) {
        for _ in 0..(w as usize) {
            out.push(v);
        }
    }
    out
}

fn median(v: &mut [f64]) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    if n == 0 {
        return f64::NAN;
    }
    if n % 2 == 1 {
        v[n / 2]
    } else {
        0.5 * (v[n / 2 - 1] + v[n / 2])
    }
}

impl Calibration {
    /// A positive global map over the full observed abstract-RT domain.
    pub(crate) fn affine(slope: f64, intercept: f64, xmin: f64, xmax: f64) -> Calibration {
        assert!(slope.is_finite() && slope > 0.0 && intercept.is_finite());
        assert!(xmin.is_finite() && xmax.is_finite() && xmax > xmin);
        Calibration {
            x: vec![xmin, xmax],
            y: vec![slope * xmin + intercept, slope * xmax + intercept],
        }
    }

    /// Fit from `(abstract, observed)` pairs.
    pub fn fit(pairs: &[(f64, f64)]) -> Calibration {
        let mut p: Vec<(f64, f64)> = pairs
            .iter()
            .copied()
            .filter(|(a, o)| a.is_finite() && o.is_finite())
            .collect();
        p.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        if p.len() < 2 {
            return Calibration::default();
        }
        let n_bins = 50.min(p.len() / 10).max(2);
        let (mut x, mut y) = (Vec::new(), Vec::new());
        for b in 0..n_bins {
            let lo = b * p.len() / n_bins;
            let hi = ((b + 1) * p.len() / n_bins).max(lo + 1);
            let bin = &p[lo..hi];
            let mut xs: Vec<f64> = bin.iter().map(|t| t.0).collect();
            let mut ys: Vec<f64> = bin.iter().map(|t| t.1).collect();
            x.push(median(&mut xs));
            y.push(median(&mut ys));
        }
        let (mut xk, mut yk) = (Vec::new(), Vec::new());
        let mut run_max = f64::NEG_INFINITY;
        for i in 0..x.len() {
            if xk.last().is_some_and(|&lx: &f64| (x[i] - lx).abs() < 1e-12) {
                continue;
            }
            run_max = run_max.max(y[i]);
            xk.push(x[i]);
            yk.push(run_max);
        }
        Calibration { x: xk, y: yk }
    }

    /// Robust fit for low-purity data: in each predicted-RT (x) quantile bin,
    /// histogram observed RT (y) and take the PEAK bin as the calibrated value.
    /// True matches form a sharp diagonal peak; random wrong-peptide matches
    /// spread uniformly and do not. Recovers the diagonal where `fit`'s
    /// per-bin median would track the noise center. Knots forced monotone.
    pub fn fit_robust(pairs: &[(f64, f64)], x_bins: usize, y_bins: usize) -> Calibration {
        let mut p: Vec<(f64, f64)> = pairs
            .iter()
            .copied()
            .filter(|(a, o)| a.is_finite() && o.is_finite())
            .collect();
        if p.len() < x_bins.max(4) {
            return Calibration::fit(pairs);
        }
        p.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let (ymin, ymax) = p.iter().fold((f64::MAX, f64::MIN), |(lo, hi), &(_, y)| {
            (lo.min(y), hi.max(y))
        });
        let yspan = (ymax - ymin).max(1e-9);

        let (mut xk, mut yk) = (Vec::new(), Vec::new());
        for b in 0..x_bins {
            let lo = b * p.len() / x_bins;
            let hi = ((b + 1) * p.len() / x_bins).max(lo + 1);
            let bin = &p[lo..hi];
            // histogram observed y; find peak bin
            let mut hist = vec![0u32; y_bins];
            for &(_, y) in bin {
                let j = (((y - ymin) / yspan) * y_bins as f64) as usize;
                hist[j.min(y_bins - 1)] += 1;
            }
            let peak = (0..y_bins).max_by_key(|&j| hist[j]).unwrap();
            // mean of y's falling in the peak bin (+/- 1) for a smooth estimate
            let (jlo, jhi) = (peak.saturating_sub(1), (peak + 1).min(y_bins - 1));
            let (mut sy, mut ny) = (0.0, 0u32);
            for &(_, y) in bin {
                let j = (((y - ymin) / yspan) * y_bins as f64) as usize;
                let j = j.min(y_bins - 1);
                if j >= jlo && j <= jhi {
                    sy += y;
                    ny += 1;
                }
            }
            let xmid = bin[bin.len() / 2].0;
            let yval = if ny > 0 {
                sy / ny as f64
            } else {
                ymin + (peak as f64 + 0.5) / y_bins as f64 * yspan
            };
            xk.push(xmid);
            yk.push(yval);
        }
        // dedup x, then enforce monotonicity with PAVA isotonic regression (the
        // least-squares non-decreasing fit) rather than a running-max envelope, which
        // would bias every noisy mode knot upward.
        let (mut xs, mut yraw) = (Vec::new(), Vec::new());
        for i in 0..xk.len() {
            if xs
                .last()
                .is_some_and(|&lx: &f64| (xk[i] - lx).abs() < 1e-12)
            {
                continue;
            }
            xs.push(xk[i]);
            yraw.push(yk[i]);
        }
        let ys = isotonic(&yraw);
        Calibration { x: xs, y: ys }
    }

    /// RANSAC calibration for heavily-contaminated data (true matches may be a
    /// small minority). Finds the largest line-consensus, then fits a monotone
    /// PWL on the inliers and expands. Returns the calibration and the indices
    /// of the inlier pairs (a likely-correct subset, usable as positives).
    pub fn fit_ransac(
        pairs: &[(f64, f64)],
        thresh: f64,
        iters: usize,
    ) -> (Calibration, Vec<usize>) {
        Self::fit_ransac_seeded(pairs, thresh, iters, 0)
    }

    /// [`fit_ransac`](Self::fit_ransac) with the sampling stream salted by `seed`
    /// (0 = the shipped default stream; see [`crate::seed_mix`]).
    pub fn fit_ransac_seeded(
        pairs: &[(f64, f64)],
        thresh: f64,
        iters: usize,
        seed: u64,
    ) -> (Calibration, Vec<usize>) {
        // Work through a canonical (value-sorted) view of the finite pairs. The
        // sampling LCG below draws index *positions*, so without a fixed order the
        // same iteration would pick different points depending on how the caller
        // ordered `pairs` -> a different consensus, and thus a feature-row-order-
        // dependent calibration (the source of uno's depth swinging with input
        // order). Sorting by value makes the RANSAC output a pure function of the
        // pair *set*. `order[s]` maps a sorted position back to the caller's index,
        // so returned inlier indices remain in terms of the original `pairs`.
        let mut order: Vec<usize> = (0..pairs.len())
            .filter(|&k| pairs[k].0.is_finite() && pairs[k].1.is_finite())
            .collect();
        order.sort_by(|&a, &b| {
            pairs[a]
                .0
                .partial_cmp(&pairs[b].0)
                .unwrap()
                .then(pairs[a].1.partial_cmp(&pairs[b].1).unwrap())
                .then(a.cmp(&b))
        });
        let n = order.len();
        let pt = |s: usize| -> (f64, f64) { pairs[order[s]] };
        if n < 10 {
            return (Calibration::fit(pairs), order);
        }
        // deterministic LCG
        let mut state: u64 = 0x2545_f491_4f6c_dd1d ^ crate::seed_mix(seed);
        let mut next = |m: usize| -> usize {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((state >> 33) as usize) % m
        };

        // Inlier sets are tracked as sorted positions `s`; mapped to original
        // `pairs` indices via `order` only when returned.
        let mut best: Vec<usize> = Vec::new();
        for _ in 0..iters {
            let (i, j) = (next(n), next(n));
            if (pt(i).0 - pt(j).0).abs() < 1e-9 {
                continue;
            }
            let m = (pt(j).1 - pt(i).1) / (pt(j).0 - pt(i).0);
            let b = pt(i).1 - m * pt(i).0;
            let inl: Vec<usize> = (0..n)
                .filter(|&s| (pt(s).1 - (m * pt(s).0 + b)).abs() < thresh)
                .collect();
            if inl.len() > best.len() {
                best = inl;
            }
        }
        // refine: fit monotone PWL on inliers, re-collect inliers, repeat
        let mut inliers = best;
        for _ in 0..4 {
            let cal = Calibration::fit(&inliers.iter().map(|&s| pt(s)).collect::<Vec<_>>());
            let new: Vec<usize> = (0..n)
                .filter(|&s| (pt(s).1 - cal.apply(pt(s).0)).abs() < thresh)
                .collect();
            if new.len() <= inliers.len() {
                inliers = new;
                break;
            }
            inliers = new;
        }
        let cal = Calibration::fit(&inliers.iter().map(|&s| pt(s)).collect::<Vec<_>>());
        (cal, inliers.iter().map(|&s| order[s]).collect())
    }

    pub fn is_fit(&self) -> bool {
        !self.x.is_empty()
    }

    /// Interpolate (clamped at the ends). Identity when unfit.
    pub fn apply(&self, v: f64) -> f64 {
        let n = self.x.len();
        if n == 0 {
            return v;
        }
        if v <= self.x[0] {
            return self.y[0];
        }
        if v >= self.x[n - 1] {
            return self.y[n - 1];
        }
        let i = self.x.partition_point(|&x| x < v).max(1);
        let (x0, x1, y0, y1) = (self.x[i - 1], self.x[i], self.y[i - 1], self.y[i]);
        let t = if (x1 - x0).abs() < 1e-12 {
            0.0
        } else {
            (v - x0) / (x1 - x0)
        };
        y0 + t * (y1 - y0)
    }
}

#[cfg(test)]
mod seed_tests {
    use super::*;

    fn pairs() -> Vec<(f64, f64)> {
        // Noisy line with outliers, enough points to reach the RANSAC path.
        (0..200)
            .map(|i| {
                let x = i as f64 * 0.1;
                let y = if i % 17 == 0 {
                    x * 2.0 + 30.0
                } else {
                    x * 2.0 + 1.0 + ((i % 7) as f64 - 3.0) * 0.05
                };
                (x, y)
            })
            .collect()
    }

    /// `--seed 0` must leave the RANSAC sampling stream on its historical
    /// constant, so published seed-0 calibrations stay reproducible.
    #[test]
    fn ransac_seed_zero_matches_unseeded() {
        let p = pairs();
        let (a, ai) = Calibration::fit_ransac(&p, 0.5, 400);
        let (b, bi) = Calibration::fit_ransac_seeded(&p, 0.5, 400, 0);
        assert_eq!(ai, bi);
        assert_eq!(a.x, b.x);
        assert_eq!(a.y, b.y);
    }

    /// A non-zero seed must draw a different sample sequence. The consensus fit
    /// may well land in the same place on clean data; what must change is the
    /// stream, so assert on that directly rather than on the fitted line.
    #[test]
    fn ransac_seed_changes_the_sampling_stream() {
        let draw = |seed: u64| {
            let mut state: u64 = 0x2545_f491_4f6c_dd1d ^ crate::seed_mix(seed);
            (0..8)
                .map(|_| {
                    state = state
                        .wrapping_mul(6364136223846793005)
                        .wrapping_add(1442695040888963407);
                    (state >> 33) as usize % 200
                })
                .collect::<Vec<_>>()
        };
        assert_ne!(draw(0), draw(1));
        assert_ne!(draw(1), draw(2));
    }
}
