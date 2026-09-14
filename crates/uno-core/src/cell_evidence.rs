//! Optional extraction descriptors. These are learned covariates, not identity
//! probabilities or additional peptide observations. Zero-intensity cells are
//! unfilled; an unavailable q-value never counts as a confident extraction.
pub const NAMES: [&str; 12] = [
    "mbr_filled_count", "mbr_filled_fraction", "mbr_confident_count",
    "cell_confident_fraction", "direct_mean_q", "direct_worst_q",
    "direct_q_coverage", "direct_confident_fraction", "mbr_mean_q",
    "mbr_worst_q", "mbr_q_coverage", "mbr_confident_fraction",
];

/// Ablations keep column positions fixed. The first two fields require no q;
/// confident counts/fractions and q coverage remain with the q-derived block.
pub fn select(mut values: [f32; 12], counts: bool, quality: bool) -> [f32; 12] {
    if !counts { values[..2].fill(0.0); }
    if !quality { values[2..].fill(0.0); }
    values
}

#[derive(Clone, Default)]
struct Stratum { filled: usize, known: usize, good: usize, sum: f64, worst: f64 }
impl Stratum {
    fn add(&mut self, q: Option<f64>) {
        self.filled += 1;
        if let Some(q) = q {
            self.known += 1;
            self.good += usize::from(q <= 0.01);
            self.sum += q;
            self.worst = self.worst.max(q);
        }
    }
    fn values(&self) -> [f32; 4] {
        // Neutral conservative placeholders, distinguished by explicit coverage.
        [if self.known > 0 { (self.sum / self.known as f64) as f32 } else { 1.0 },
         if self.known > 0 { self.worst as f32 } else { 1.0 },
         self.known as f32 / self.filled.max(1) as f32,
         self.good as f32 / self.filled.max(1) as f32]
    }
}
#[derive(Clone, Default)]
pub(crate) struct Summary { direct: Stratum, mbr: Stratum }
impl Summary {
    pub fn add(&mut self, is_mbr: bool, intensity: f64, q: Option<f64>) {
        if intensity <= 0.0 { return; }
        if is_mbr { self.mbr.add(q); } else { self.direct.add(q); }
    }
    pub fn values(&self, n_sources: usize) -> [f32; 12] {
        let mut out = [0.0; 12];
        out[..4].copy_from_slice(&[self.mbr.filled as f32,
            self.mbr.filled as f32 / n_sources.max(1) as f32,
            self.mbr.good as f32,
            (self.direct.good + self.mbr.good) as f32 / n_sources.max(1) as f32]);
        out[4..8].copy_from_slice(&self.direct.values());
        out[8..12].copy_from_slice(&self.mbr.values());
        out
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn unfilled_and_unavailable_cells_do_not_become_confident() {
        let mut s = Summary::default();
        s.add(false, 10.0, Some(0.0));
        s.add(false, 20.0, Some(0.02));
        s.add(true, 5.0, Some(0.01));
        s.add(true, 5.0, None);
        s.add(true, 0.0, Some(0.0));
        let v = s.values(5);
        assert_eq!(v, [2.0, 0.4, 1.0, 0.4, 0.01, 0.02, 1.0, 0.5, 0.01, 0.01, 0.5, 0.5]);
    }
    #[test]
    fn ablations_separate_counts_from_q_derived_information() {
        let v = std::array::from_fn(|i| (i + 1) as f32);
        assert_eq!(select(v, true, true), v);
        assert_eq!(&select(v, true, false)[..2], &v[..2]);
        assert_eq!(&select(v, true, false)[2..], &[0.0; 10]);
        assert_eq!(&select(v, false, true)[..2], &[0.0; 2]);
        assert_eq!(&select(v, false, true)[2..], &v[2..]);
        assert_eq!(select(v, false, false), [0.0; 12]);
    }
    #[test]
    fn missing_q_is_distinct_from_zero_q() {
        let mut s = Summary::default(); s.add(true, 5.0, None);
        assert_eq!(&s.values(1)[8..], &[1.0, 1.0, 0.0, 0.0]);
        let mut s = Summary::default(); s.add(true, 5.0, Some(0.0));
        assert_eq!(&s.values(1)[8..], &[0.0, 0.0, 1.0, 1.0]);
    }
}
