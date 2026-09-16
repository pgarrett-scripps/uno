//! Scalar NumPy 1.26 argsort ordering, including equal-key partition behavior.
//! Adapted from NumPy's BSD-3-Clause npysort introsort (Charles R. Harris).
//! See THIRD_PARTY.md and licenses/NUMPY-BSD.txt.
pub fn argsort(values: &[f64], descending: bool) -> Vec<usize> {
    let mut order: Vec<usize> = (0..values.len()).collect();
    if descending {
        order.reverse();
    }
    if order.len() > 1 {
        let hi = order.len() - 1;
        let depth = 2 * (usize::BITS - 1 - order.len().leading_zeros()) as i32;
        sort(values, &mut order, 0, hi, depth);
    }
    if descending {
        order.reverse();
    }
    order
}
fn sort(v: &[f64], a: &mut [usize], mut lo: usize, mut hi: usize, mut depth: i32) {
    if depth < 0 {
        heap(v, &mut a[lo..=hi]);
        return;
    }
    while hi - lo > 15 {
        let mid = lo + (hi - lo) / 2;
        if v[a[mid]] < v[a[lo]] {
            a.swap(mid, lo);
        }
        if v[a[hi]] < v[a[mid]] {
            a.swap(hi, mid);
        }
        if v[a[mid]] < v[a[lo]] {
            a.swap(mid, lo);
        }
        let pivot = v[a[mid]];
        let mut i = lo;
        let mut j = hi - 1;
        a.swap(mid, j);
        loop {
            i += 1;
            while v[a[i]] < pivot {
                i += 1;
            }
            j -= 1;
            while pivot < v[a[j]] {
                j -= 1;
            }
            if i >= j {
                break;
            }
            a.swap(i, j);
        }
        a.swap(i, hi - 1);
        depth -= 1;
        // NumPy checks depth only when restoring the larger partition.
        if i - lo < hi - i {
            sort(v, a, i + 1, hi, depth);
            hi = i - 1;
        } else {
            sort(v, a, lo, i - 1, depth);
            lo = i + 1;
        }
    }
    for i in lo + 1..=hi {
        let x = a[i];
        let mut j = i;
        while j > lo && v[x] < v[a[j - 1]] {
            a[j] = a[j - 1];
            j -= 1;
        }
        a[j] = x;
    }
}
fn heap(v: &[f64], a: &mut [usize]) {
    fn sift(v: &[f64], a: &mut [usize], mut root: usize, end: usize) {
        loop {
            let mut child = 2 * root + 1;
            if child >= end {
                break;
            }
            if child + 1 < end && v[a[child]] < v[a[child + 1]] {
                child += 1;
            }
            if v[a[root]] >= v[a[child]] {
                break;
            }
            a.swap(root, child);
            root = child;
        }
    }
    for i in (0..a.len() / 2).rev() {
        sift(v, a, i, a.len());
    }
    for end in (1..a.len()).rev() {
        a.swap(0, end);
        sift(v, a, 0, end);
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn sorted_and_permutation() {
        let v: Vec<_> = (0..1000).map(|i| ((i * 71) % 19) as f64).collect();
        let a = argsort(&v, false);
        assert!(a.windows(2).all(|w| v[w[0]] <= v[w[1]]));
        let mut p = a;
        p.sort();
        assert_eq!(p, (0..1000).collect::<Vec<_>>());
    }
}
