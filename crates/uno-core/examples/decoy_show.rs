//! Show example decoy peptides from each method, side by side with the target.
//! Run: cargo run --release -p uno-core --example decoy_show

use std::sync::Arc;
use uno_core::decoy::{generate, DecoyMethod};

/// Minimal tryptic digest: cut after K/R unless followed by P. mc=0, len 7..=50.
fn digest(seq: &str) -> Vec<String> {
    let b = seq.as_bytes();
    let mut peps = Vec::new();
    let mut start = 0usize;
    for i in 0..b.len() {
        let cut = (b[i] == b'K' || b[i] == b'R') && (i + 1 >= b.len() || b[i + 1] != b'P');
        if cut {
            let p = &seq[start..=i];
            if p.len() >= 7 && p.len() <= 50 {
                peps.push(p.to_string());
            }
            start = i + 1;
        }
    }
    peps
}

fn main() {
    // Load a chunk of real proteins so the Markov/de Bruijn models learn real
    // amino-acid statistics, but only DISPLAY one chosen protein's decoys.
    let path = "/home/patrick-garrett/Repos/koth_rust/benchmark/data/fasta/human_ecoli.fasta";
    let contents = std::fs::read_to_string(path).expect("fasta");
    let mut targets: Vec<(Arc<str>, String)> = Vec::new();
    let mut acc = String::new();
    let mut seq = String::new();
    for line in contents.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if !acc.is_empty() {
                targets.push((Arc::from(acc.as_str()), seq.clone()));
            }
            acc = rest.split_whitespace().next().unwrap_or("").to_string();
            seq.clear();
        } else {
            seq.push_str(line.trim());
        }
        if targets.len() >= 3000 {
            break;
        }
    }
    if !acc.is_empty() && targets.len() < 3000 {
        targets.push((Arc::from(acc.as_str()), seq.clone()));
    }

    // Choose a display protein with a handful of in-range tryptic peptides.
    let show_idx = targets
        .iter()
        .position(|(_, s)| {
            let n = digest(s).len();
            n >= 4 && n <= 8 && s.len() < 350
        })
        .unwrap_or(0);
    let (show_acc, show_seq) = &targets[show_idx];

    println!("=== TARGET protein: {} ({} aa) ===", show_acc, show_seq.len());
    println!("seq: {}\n", show_seq);
    let tpeps = digest(show_seq);
    println!("TARGET tryptic peptides ({}):", tpeps.len());
    for p in &tpeps {
        println!("    {}", p);
    }
    println!();

    let methods = [
        ("reverse", DecoyMethod::Reverse),
        ("random_uniform", DecoyMethod::RandomUniform),
        ("random_freq", DecoyMethod::RandomFreq),
        ("markov (order 1)", DecoyMethod::Markov { order: 1 }),
        ("markov (order 2)", DecoyMethod::Markov { order: 2 }),
        ("debruijn (k=4)", DecoyMethod::DeBruijn { k: 4 }),
        ("debruijn_cleave (k=4, keep K/R/P)", DecoyMethod::DeBruijnCleave { k: 4 }),
    ];

    for (name, m) in methods {
        let decoys = generate(&targets, m, "DECOY_");
        let (dacc, dseq) = &decoys[show_idx];
        println!("=== {} ===", name);
        println!("decoy protein {}:", dacc);
        println!("seq: {}", dseq);
        let dpeps = digest(dseq);
        println!("decoy tryptic peptides ({}):", dpeps.len());
        for p in &dpeps {
            println!("    {}", p);
        }
        println!();
    }

    // The default `shuffle` method works per-PEPTIDE (interior shuffle), not by
    // synthesizing a protein. Show it directly on the target's peptides.
    println!("=== shuffle (per-peptide interior shuffle; uno default) ===");
    println!("(keeps first+last residue, permutes the interior -> same residues, same mass)");
    let mut h: u64 = 0x9e37_79b9_7f4a_7c15;
    for p in &tpeps {
        let b = p.as_bytes();
        if b.len() <= 2 {
            println!("    {}  ->  {}", p, p);
            continue;
        }
        let mut interior: Vec<u8> = b[1..b.len() - 1].to_vec();
        for i in (1..interior.len()).rev() {
            h ^= h << 13;
            h ^= h >> 7;
            h ^= h << 17;
            interior.swap(i, (h % (i as u64 + 1)) as usize);
        }
        let mut d = vec![b[0]];
        d.extend_from_slice(&interior);
        d.push(b[b.len() - 1]);
        println!("    {}  ->  {}", p, String::from_utf8(d).unwrap());
    }
}
