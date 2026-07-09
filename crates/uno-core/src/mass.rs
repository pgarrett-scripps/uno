//! High-precision (f64) monoisotopic masses.
//!
//! Sage stores residue masses as f32, which can accumulate ~0.1-0.5 ppm of
//! rounding over a 20-mer. For MS1-only search, precursor mass *is* the primary
//! identification signal, so we recompute neutral masses in f64 from the
//! stripped sequence plus Sage's modification deltas.

/// Monoisotopic mass of water (Da).
pub const H2O: f64 = 18.010_564_684;
/// Monoisotopic mass of a proton (Da).
pub const PROTON: f64 = 1.007_276_466_88;
/// C13 - C12 mass difference; the spacing between isotope peaks (Da).
pub const NEUTRON: f64 = 1.003_354_835;

/// Monoisotopic residue mass for an uppercase amino-acid byte, or `None` if the
/// character is not a standard residue.
pub fn residue(aa: u8) -> Option<f64> {
    Some(match aa {
        b'A' => 71.037_114_4,
        b'R' => 156.101_111_1,
        b'N' => 114.042_927_5,
        b'D' => 115.026_943_1,
        b'C' => 103.009_184_9,
        b'E' => 129.042_593_1,
        b'Q' => 128.058_577_5,
        b'G' => 57.021_463_7,
        b'H' => 137.058_911_9,
        b'I' => 113.084_064_0,
        b'L' => 113.084_064_0,
        b'K' => 128.094_963_0,
        b'M' => 131.040_484_9,
        b'F' => 147.068_413_9,
        b'P' => 97.052_763_9,
        b'S' => 87.032_028_4,
        b'T' => 101.047_678_5,
        b'W' => 186.079_312_9,
        b'Y' => 163.063_328_6,
        b'V' => 99.068_413_9,
        // Selenocysteine / pyrrolysine — rare; fall back to Sage if encountered.
        _ => return None,
    })
}

/// Neutral monoisotopic mass of a stripped sequence plus per-residue and
/// terminal modification deltas. Returns `None` if any residue is non-standard
/// (caller should fall back to Sage's f32 mass).
pub fn neutral_mass(
    sequence: &[u8],
    modifications: &[f32],
    nterm: Option<f32>,
    cterm: Option<f32>,
) -> Option<f64> {
    let mut m = H2O;
    for (i, &aa) in sequence.iter().enumerate() {
        m += residue(aa)?;
        if let Some(delta) = modifications.get(i) {
            m += *delta as f64;
        }
    }
    if let Some(n) = nterm {
        m += n as f64;
    }
    if let Some(c) = cterm {
        m += c as f64;
    }
    Some(m)
}
