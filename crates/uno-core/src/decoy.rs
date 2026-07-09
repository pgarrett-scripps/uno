//! Alternative decoy-protein generators (experimental).
//!
//! The default decoy strategy ([`crate::database::global_shuffle_decoy_db`]) works
//! at the *peptide* level: it shuffles each target peptide's interior. These
//! generators instead synthesize whole decoy *proteins* from the target proteome,
//! which are then digested by the same enzyme — so peptide count, length and mass
//! distributions emerge naturally from digestion rather than being imposed.
//!
//! The methods span a spectrum of how much real local sequence structure the
//! decoys inherit:
//!   - `RandomUniform`  — residues drawn uniformly from the 20 AAs (no structure)
//!   - `RandomFreq`     — residues drawn from the proteome AA frequency (order-0)
//!   - `Markov{order}`  — order-k residue Markov chain (local di/tri-peptide stats)
//!   - `DeBruijn{k}`    — random walk over the k-mer de Bruijn graph (recombines
//!                        real k-mers; dead-ends restart at a random observed node)
//!   - `Reverse`        — reverse each protein (classic; exact per-protein mass)
//!
//! Every generator is fully deterministic (seeded xorshift, no clock/RNG entropy)
//! so a search is byte-reproducible. Each target protein maps to exactly one decoy
//! protein of the same length, giving an approximately 1:1 peptide count.

use std::collections::HashMap;
use std::sync::Arc;

const AAS: &[u8; 20] = b"ARNDCQEGHILKMFPSTWYV";

/// Decoy generation method, parsed from the `decoy_method` config string.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DecoyMethod {
    /// Peptide-interior shuffle (handled in database.rs, not here).
    Shuffle,
    Reverse,
    RandomUniform,
    RandomFreq,
    Markov { order: usize },
    DeBruijn { k: usize },
    /// Cleavage-site-preserving de Bruijn: K/R/P stay at their exact positions
    /// (so the tryptic peptide lengths/count match the target exactly), and the
    /// remaining positions are filled by a back-off Markov walk (orders k..1..0)
    /// over real sequence — every residue drawn from real high-order context, no
    /// abrupt restarts. The most realistic decoy that still isn't a plain reversal.
    DeBruijnCleave { k: usize },
}

impl DecoyMethod {
    /// Parse `decoy_method` + `decoy_order` from config. Unknown → Shuffle.
    pub fn parse(name: &str, order: usize) -> DecoyMethod {
        match name.to_ascii_lowercase().as_str() {
            "reverse" => DecoyMethod::Reverse,
            "random_uniform" | "random" => DecoyMethod::RandomUniform,
            "random_freq" | "freq" => DecoyMethod::RandomFreq,
            "markov" => DecoyMethod::Markov {
                order: order.max(1),
            },
            "debruijn" => DecoyMethod::DeBruijn { k: order.max(2) },
            "debruijn_cleave" | "debruijn2" => DecoyMethod::DeBruijnCleave { k: order.max(2) },
            _ => DecoyMethod::Shuffle,
        }
    }
    pub fn is_protein_level(&self) -> bool {
        !matches!(self, DecoyMethod::Shuffle)
    }
}

#[inline]
fn xorshift(state: &mut u64) -> u64 {
    let mut x = *state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    *state = x;
    x
}

/// Seed a per-protein PRNG from a fixed constant + the accession bytes, so the
/// stream is deterministic but distinct per protein.
fn seed_for(acc: &str, salt: u64) -> u64 {
    let mut h: u64 = 0x9e37_79b9_7f4a_7c15 ^ salt;
    for b in acc.bytes() {
        h ^= b as u64;
        h = h.wrapping_mul(1099511628211);
    }
    h | 1
}

/// Sample an index in `[0, n)` weighted by `weights[0..n]` (sum > 0), else fall
/// back to a uniform pick. `r` is a fresh random u64.
#[inline]
fn weighted_pick(weights: &[u32], total: u64, r: u64) -> usize {
    if total == 0 {
        return (r % weights.len() as u64) as usize;
    }
    let mut t = (r % total) as i64;
    for (i, &w) in weights.iter().enumerate() {
        t -= w as i64;
        if t < 0 {
            return i;
        }
    }
    weights.len() - 1
}

#[inline]
fn aa_index(aa: u8) -> Option<usize> {
    AAS.iter().position(|&a| a == aa)
}

/// Generate decoy proteins for the given method as `(accession, sequence)` pairs.
/// `targets` are the real proteins (accession, sequence). Decoy accessions are
/// prefixed with `decoy_tag`.
pub fn generate(
    targets: &[(Arc<str>, String)],
    method: DecoyMethod,
    decoy_tag: &str,
) -> Vec<(String, String)> {
    match method {
        DecoyMethod::Shuffle => Vec::new(), // handled elsewhere
        DecoyMethod::Reverse => reverse_decoys(targets, decoy_tag),
        DecoyMethod::RandomUniform => random_decoys(targets, decoy_tag, None),
        DecoyMethod::RandomFreq => {
            let freq = global_aa_freq(targets);
            random_decoys(targets, decoy_tag, Some(&freq))
        }
        DecoyMethod::Markov { order } => markov_decoys(targets, decoy_tag, order, false),
        DecoyMethod::DeBruijn { k } => markov_decoys(targets, decoy_tag, k - 1, true),
        DecoyMethod::DeBruijnCleave { k } => debruijn_cleave_decoys(targets, decoy_tag, k),
    }
}

/// Trypsin-relevant residues: cut after K/R unless followed by P. Preserving the
/// positions of all three keeps the cleavage map (and thus every peptide length)
/// identical between target and decoy.
#[inline]
fn is_cleave(aa: u8) -> bool {
    aa == b'K' || aa == b'R' || aa == b'P'
}

/// Cleavage-site-preserving de Bruijn / back-off Markov decoys. K/R/P (and any
/// nonstandard residue) stay put; every other position is sampled from the
/// highest-order real context available (orders k down to 1, then global), with
/// K/R/P excluded from the fill alphabet so no cleavage site is created or
/// destroyed. Result: same digestion map, same peptide lengths, realistic interior.
fn debruijn_cleave_decoys(
    targets: &[(Arc<str>, String)],
    tag: &str,
    k: usize,
) -> Vec<(String, String)> {
    let k = k.max(1);
    // tables[o-1]: context of o residues -> next-residue counts (over real seq).
    let mut tables: Vec<HashMap<Vec<u8>, [u32; 20]>> = (0..k).map(|_| HashMap::new()).collect();
    let mut g0 = [0u32; 20]; // order-0 counts over NON-cleavage residues only.
    for (_, seq) in targets {
        let b = seq.as_bytes();
        for i in 0..b.len() {
            let ni = match aa_index(b[i]) {
                Some(x) => x,
                None => continue,
            };
            if !is_cleave(b[i]) {
                g0[ni] += 1;
            }
            for o in 1..=k {
                if i >= o {
                    let ctx = &b[i - o..i];
                    if ctx.iter().all(|&c| aa_index(c).is_some()) {
                        tables[o - 1].entry(ctx.to_vec()).or_insert([0u32; 20])[ni] += 1;
                    }
                }
            }
        }
    }
    let g0_total: u64 = g0.iter().map(|&w| w as u64).sum();

    // Zero the cleavage residues out of a count vector (fill alphabet = 17 AAs).
    let mask_cleave = |mut w: [u32; 20]| -> ([u32; 20], u64) {
        for &c in b"KRP" {
            if let Some(ix) = aa_index(c) {
                w[ix] = 0;
            }
        }
        let t: u64 = w.iter().map(|&x| x as u64).sum();
        (w, t)
    };

    targets
        .iter()
        .map(|(acc, seq)| {
            let b = seq.as_bytes();
            let n = b.len();
            let mut state = seed_for(acc, 0xDBC1);
            let mut out: Vec<u8> = Vec::with_capacity(n);
            for i in 0..n {
                if is_cleave(b[i]) || aa_index(b[i]).is_none() {
                    out.push(b[i]); // preserve cleavage residues & nonstandard chars
                    continue;
                }
                // Back off from the highest available order to the lowest.
                let mut chosen: Option<u8> = None;
                for o in (1..=k).rev() {
                    if out.len() < o {
                        continue;
                    }
                    let ctx = &out[out.len() - o..];
                    if let Some(counts) = tables[o - 1].get(ctx) {
                        let (w, total) = mask_cleave(*counts);
                        if total > 0 {
                            let r = xorshift(&mut state);
                            chosen = Some(AAS[weighted_pick(&w, total, r)]);
                            break;
                        }
                    }
                }
                let aa = chosen.unwrap_or_else(|| {
                    let r = xorshift(&mut state);
                    AAS[weighted_pick(&g0, g0_total.max(1), r)]
                });
                out.push(aa);
            }
            (format!("{tag}{acc}"), String::from_utf8(out).unwrap())
        })
        .collect()
}

/// Global amino-acid counts over the target proteome, indexed by AAS position.
fn global_aa_freq(targets: &[(Arc<str>, String)]) -> [u32; 20] {
    let mut f = [0u32; 20];
    for (_, seq) in targets {
        for &b in seq.as_bytes() {
            if let Some(i) = aa_index(b) {
                f[i] += 1;
            }
        }
    }
    f
}

fn reverse_decoys(targets: &[(Arc<str>, String)], tag: &str) -> Vec<(String, String)> {
    targets
        .iter()
        .map(|(acc, seq)| {
            let rev: String = seq.chars().rev().collect();
            (format!("{tag}{acc}"), rev)
        })
        .collect()
}

/// Random proteins, residues drawn either uniformly (`freq=None`) or from the
/// proteome AA frequency. Each decoy matches its target protein's length.
fn random_decoys(
    targets: &[(Arc<str>, String)],
    tag: &str,
    freq: Option<&[u32; 20]>,
) -> Vec<(String, String)> {
    let total: u64 = freq.map(|f| f.iter().map(|&w| w as u64).sum()).unwrap_or(0);
    targets
        .iter()
        .map(|(acc, seq)| {
            let mut state = seed_for(acc, 0xA5A5);
            let n = seq.len();
            let mut out = Vec::with_capacity(n);
            for _ in 0..n {
                let r = xorshift(&mut state);
                let idx = match freq {
                    Some(f) => weighted_pick(f, total, r),
                    None => (r % 20) as usize,
                };
                out.push(AAS[idx]);
            }
            (format!("{tag}{acc}"), String::from_utf8(out).unwrap())
        })
        .collect()
}

/// Order-k Markov / de Bruijn decoys. Builds an order-k transition table over the
/// proteome (context of k residues -> next-residue counts). Generation starts from
/// each protein's own first k residues and samples forward to matched length.
///
/// `debruijn=false`: on an unseen context, back off to the global AA frequency.
/// `debruijn=true` : on an unseen context (graph dead-end), restart from a random
///                   observed context — a true k-mer-graph walk.
fn markov_decoys(
    targets: &[(Arc<str>, String)],
    tag: &str,
    order: usize,
    debruijn: bool,
) -> Vec<(String, String)> {
    let order = order.max(1);
    // context (k bytes packed) -> [count per next AA]
    let mut table: HashMap<Vec<u8>, [u32; 20]> = HashMap::new();
    let mut contexts: Vec<Vec<u8>> = Vec::new();
    for (_, seq) in targets {
        let b = seq.as_bytes();
        if b.len() <= order {
            continue;
        }
        for w in b.windows(order + 1) {
            let ctx = &w[..order];
            let next = w[order];
            if ctx.iter().any(|&c| aa_index(c).is_none()) {
                continue;
            }
            if let Some(ni) = aa_index(next) {
                let e = table.entry(ctx.to_vec()).or_insert_with(|| {
                    contexts.push(ctx.to_vec());
                    [0u32; 20]
                });
                e[ni] += 1;
            }
        }
    }
    let gfreq = global_aa_freq(targets);
    let gtotal: u64 = gfreq.iter().map(|&w| w as u64).sum();

    targets
        .iter()
        .map(|(acc, seq)| {
            let b = seq.as_bytes();
            let n = b.len();
            let mut state = seed_for(acc, if debruijn { 0xDB12 } else { 0x3A7C });
            let mut out: Vec<u8> = Vec::with_capacity(n);
            // Seed with the protein's own first `order` residues (sanitized).
            let mut ctx: Vec<u8> = Vec::with_capacity(order);
            for &c in b.iter().take(order) {
                ctx.push(if aa_index(c).is_some() { c } else { AAS[0] });
                out.push(*ctx.last().unwrap());
            }
            while out.len() < n {
                let next = match table.get(&ctx) {
                    Some(counts) => {
                        let total: u64 = counts.iter().map(|&w| w as u64).sum();
                        let r = xorshift(&mut state);
                        AAS[weighted_pick(counts, total, r)]
                    }
                    None => {
                        if debruijn && !contexts.is_empty() {
                            // dead-end: jump to a random observed context
                            let r = xorshift(&mut state);
                            ctx = contexts[(r as usize) % contexts.len()].clone();
                            for &c in &ctx {
                                if out.len() < n {
                                    out.push(c);
                                }
                            }
                            continue;
                        }
                        let r = xorshift(&mut state);
                        AAS[weighted_pick(&gfreq, gtotal, r)]
                    }
                };
                out.push(next);
                // slide context window
                ctx.push(next);
                if ctx.len() > order {
                    ctx.remove(0);
                }
            }
            out.truncate(n);
            (format!("{tag}{acc}"), String::from_utf8(out).unwrap())
        })
        .collect()
}
