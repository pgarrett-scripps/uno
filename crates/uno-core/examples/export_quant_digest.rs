//! Export the actual target digest for fixed quantitative protein grouping.
//! Usage: export_quant_digest FASTA DIGEST_CONFIG_JSON OUTPUT_TSV
use anyhow::{ensure, Result};
use std::io::Write;
use uno_core::{DigestConfig, PeptideDb};
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    ensure!(
        args.len() == 3,
        "usage: export_quant_digest FASTA DIGEST_CONFIG_JSON OUTPUT_TSV"
    );
    let mut config: DigestConfig = serde_json::from_str(&std::fs::read_to_string(&args[1])?)?;
    config.generate_decoys = false;
    let db = PeptideDb::build(&args[0], &config)?;
    let mut peptides: Vec<_> = db.peptides.iter().filter(|p| !p.decoy).collect();
    peptides.sort_by(|a, b| {
        a.sequence
            .cmp(&b.sequence)
            .then_with(|| a.mass.total_cmp(&b.mass))
    });
    let file = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&args[2])?;
    let mut out = std::io::BufWriter::new(file);
    writeln!(out, "sequence\tmass\tproteins")?;
    for p in &peptides {
        let mut proteins: Vec<_> = p.proteins.iter().map(|x| x.as_ref()).collect();
        proteins.sort();
        proteins.dedup();
        writeln!(out, "{}\t{}\t{}", p.sequence, p.mass, proteins.join(";"))?;
    }
    out.flush()?;
    eprintln!("{} target peptide entries", peptides.len());
    Ok(())
}
