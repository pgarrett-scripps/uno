//! CPython integer-seeded MT19937 and _randbelow used by random.choice.
//! This fixes resampling streams across Python and Rust, rather than substituting a different RNG.
pub struct PythonRandom {
    state: [u32; 624],
    index: usize,
}
impl PythonRandom {
    pub fn new(seed: u32) -> Self {
        let mut r = Self {
            state: [0; 624],
            index: 624,
        };
        r.state[0] = 19650218;
        for i in 1..624 {
            r.state[i] = 1812433253u32
                .wrapping_mul(r.state[i - 1] ^ (r.state[i - 1] >> 30))
                .wrapping_add(i as u32);
        }
        let mut i = 1;
        for _ in 0..624 {
            r.state[i] = (r.state[i]
                ^ (r.state[i - 1] ^ (r.state[i - 1] >> 30)).wrapping_mul(1664525))
            .wrapping_add(seed);
            i += 1;
            if i >= 624 {
                r.state[0] = r.state[623];
                i = 1;
            }
        }
        for _ in 0..623 {
            r.state[i] = (r.state[i]
                ^ (r.state[i - 1] ^ (r.state[i - 1] >> 30)).wrapping_mul(1566083941))
            .wrapping_sub(i as u32);
            i += 1;
            if i >= 624 {
                r.state[0] = r.state[623];
                i = 1;
            }
        }
        r.state[0] = 0x80000000;
        r
    }
    fn next(&mut self) -> u32 {
        if self.index == 624 {
            for i in 0..624 {
                let y = (self.state[i] & 0x80000000) | (self.state[(i + 1) % 624] & 0x7fffffff);
                self.state[i] = self.state[(i + 397) % 624]
                    ^ (y >> 1)
                    ^ if y & 1 != 0 { 0x9908b0df } else { 0 };
            }
            self.index = 0;
        }
        let mut y = self.state[self.index];
        self.index += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c5680;
        y ^= (y << 15) & 0xefc60000;
        y ^= y >> 18;
        y
    }
    pub fn choice(&mut self, n: usize) -> usize {
        assert!(n > 0 && n <= u32::MAX as usize);
        let bits = usize::BITS - n.leading_zeros();
        loop {
            let value = (self.next() >> (32 - bits)) as usize;
            if value < n {
                return value;
            }
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn python_stream() {
        let mut r = PythonRandom::new(42);
        assert_eq!(
            (0..10).map(|_| r.choice(1000)).collect::<Vec<_>>(),
            vec![654, 114, 25, 759, 281, 250, 228, 142, 754, 104]
        );
    }
}
