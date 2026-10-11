//! Phase 2: in-silico peptide database.
//!
//! Digestion and modification expansion come from [`crate::digest`] (a port of
//! the Sage database builder); decoys are generated here. We keep only the lean subset an MS1 engine needs — stripped sequence, neutral
//! mass, target/decoy flag, source proteins — recompute the mass in f64 for
//! precision, and sort by mass for binary-search candidate lookup.

use crate::digest::{self, DigestParams, Fasta, Peptide};
use crate::mass;
use anyhow::{Context, Result};
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::sync::Arc;

/// Digestion + modification configuration. Mirrors Sage's database block.
#[derive(Debug, Clone, serde::Deserialize)]
#[serde(default)]
pub struct DigestConfig {
    pub missed_cleavages: u8,
    pub min_len: usize,
    pub max_len: usize,
    /// Residues to cleave at (trypsin = "KR").
    pub cleave_at: String,
    /// Residues that block cleavage when following a cleavage site ("P").
    pub restrict: String,
    pub c_terminal: bool,
    pub peptide_min_mass: f32,
    pub peptide_max_mass: f32,
    /// Fixed modifications, e.g. {"C": 57.02146}.
    pub static_mods: HashMap<String, f32>,
    /// Variable modifications, e.g. {"M": [15.9949]}.
    pub variable_mods: HashMap<String, Vec<f32>>,
    pub max_variable_mods: usize,
    pub decoy_tag: String,
    pub generate_decoys: bool,
    /// Mass (Da) added to every decoy peptide so decoys are NOT mass twins of
    /// their targets. A real feature then never competes against its own decoy;
    /// decoys match only by random mass coincidence, giving an honest MS1 null
    /// where retention time has full discriminating power. (koth uses 11 Da.)
    pub decoy_mass_shift: f64,
    /// Generate decoys by SHUFFLING the peptide interior (same mass as the
    /// target) instead of reversing + mass-shifting. Shuffled same-mass decoys
    /// model the random mass-coincidence background honestly (unlike mass-shifted
    /// decoys, which sample a different mass region and underestimate the false
    /// rate), while randomized order keeps them from being RT-twins (unlike plain
    /// reversal). When true, `decoy_mass_shift` is ignored.
    pub decoy_shuffle: bool,
    /// Protein-level decoy generator (experimental). When set to anything other
    /// than "shuffle", whole decoy *proteins* are synthesized from the target
    /// proteome and digested by the same enzyme, overriding `decoy_shuffle`.
    /// One of: "shuffle" (default), "reverse", "random_uniform", "random_freq",
    /// "markov", "debruijn". See [`crate::decoy`].
    pub decoy_method: String,
    /// Order k for "markov" / "debruijn" decoy methods (ignored otherwise).
    pub decoy_order: usize,
    /// Seed for the decoy shuffle (0 = the shipped default stream; see
    /// [`crate::seed_mix`]). Set with `uno --seed`.
    pub decoy_seed: u64,
}

/// Build the target + decoy database with ms1searchpy-style shuffled decoys.
///
/// The interior residues are pooled PER PROTEIN (not globally — see the grouping
/// below), shuffled, and dealt back out to that protein's decoys, keeping each
/// peptide's terminal residues. Unlike a per-peptide permutation, this generally
/// gives each decoy a different composition, and therefore a different mass, so
/// decoys sample the random mass-coincidence background instead of being
/// mass-twins of their targets. That is what keeps the per-protein binomial's
/// random rate p low and the FDR honest. Nothing here enforces non-collision with
/// a target: decoys are deduplicated against each other only.
fn global_shuffle_decoy_db(peptides_in: Vec<Peptide>, cfg: &DigestConfig) -> Vec<TheoPeptide> {
    // Per-residue static modification delta (e.g. C +57.021).
    let static_mod = |aa: u8| -> f64 {
        cfg.static_mods
            .get(&(aa as char).to_string())
            .copied()
            .unwrap_or(0.0) as f64
    };
    let seq_mass = |seq: &[u8]| -> f64 {
        let mut m = mass::H2O;
        for &aa in seq {
            m += mass::residue(aa).unwrap_or(0.0) + static_mod(aa);
        }
        m
    };

    // Targets only (drop the digest's reversed decoys). Group each target peptide
    // under its FIRST protein so the interior pool is built PER PROTEIN — this
    // preserves each organism's residue composition, so a decoy of a human protein
    // has human-like masses and coincides with features at the same rate an absent
    // human peptide does. A global pool would average organisms and make the
    // binomial p under-estimate the in-organism random rate.
    let mut targets: Vec<TheoPeptide> = Vec::with_capacity(peptides_in.len());
    // per-protein: accession -> Vec<(seq, mc)>; plus a stable protein order.
    let mut by_prot: HashMap<Arc<str>, Vec<(Vec<u8>, u8)>> = HashMap::new();
    let mut prot_order: Vec<Arc<str>> = Vec::new();
    for p in peptides_in.into_iter().filter(|p| !p.decoy) {
        let mass = mass::neutral_mass(&p.sequence, &p.modifications, p.nterm, p.cterm)
            .unwrap_or(p.monoisotopic as f64);
        let seq_vec = p.sequence.to_vec();
        let sequence: Arc<str> = Arc::from(std::str::from_utf8(&seq_vec).unwrap_or("").to_owned());
        targets.push(TheoPeptide {
            sequence,
            mass,
            decoy: false,
            missed_cleavages: p.missed_cleavages,
            proteins: p.proteins.clone(),
        });
        if let Some(first) = p.proteins.first() {
            by_prot
                .entry(first.clone())
                .or_insert_with(|| {
                    prot_order.push(first.clone());
                    Vec::new()
                })
                .push((seq_vec, p.missed_cleavages));
        }
    }
    prot_order.sort();

    // For each protein: pool its peptides' interior residues, deterministic
    // shuffle, deal back out to decoys (keeping terminal residues). Reuse a decoy
    // per unique target sequence so shared peptides map consistently.
    let mut shuf_map: HashMap<Vec<u8>, (Vec<u8>, u8)> = HashMap::new();
    let mut decoys: Vec<TheoPeptide> = Vec::with_capacity(targets.len());
    let mut decoy_prot: HashMap<Vec<u8>, Vec<Arc<str>>> = HashMap::new();
    // first, map each unique decoy sequence to the union of proteins it should
    // carry (mirror of the target peptide->proteins map).
    for t in &targets {
        decoy_prot
            .entry(t.sequence.as_bytes().to_vec())
            .or_insert_with(|| t.proteins.clone());
    }
    for acc in &prot_order {
        let peps = &by_prot[acc];
        let mut pool: Vec<u8> = Vec::new();
        for (seq, _) in peps {
            if seq.len() > 2 {
                pool.extend_from_slice(&seq[1..seq.len() - 1]);
            }
        }
        if pool.is_empty() {
            continue;
        }
        let mut h: u64 =
            0x9e37_79b9_7f4a_7c15 ^ (pool.len() as u64) ^ crate::seed_mix(cfg.decoy_seed);
        for b in acc.bytes() {
            h ^= b as u64;
            h = h.wrapping_mul(1099511628211);
        }
        for i in (1..pool.len()).rev() {
            h ^= h << 13;
            h ^= h >> 7;
            h ^= h << 17;
            pool.swap(i, (h % (i as u64 + 1)) as usize);
        }
        let mut cur = 0usize;
        for (seq, mc) in peps {
            if seq.len() <= 2 || shuf_map.contains_key(seq) {
                continue;
            }
            let mut d = Vec::with_capacity(seq.len());
            d.push(seq[0]);
            for _ in 1..seq.len() - 1 {
                d.push(pool[cur % pool.len()]);
                cur += 1;
            }
            d.push(seq[seq.len() - 1]);
            shuf_map.insert(seq.clone(), (d, *mc));
        }
    }
    // Deterministic order: HashMap iteration order is process-random, and decoy
    // sequence collisions resolve to whichever target is seen first — so iterate
    // shuf_map sorted by target sequence, or the decoy DB (and every downstream
    // result) varies run to run.
    let mut shuf_sorted: Vec<(&Vec<u8>, &(Vec<u8>, u8))> = shuf_map.iter().collect();
    shuf_sorted.sort_by(|a, b| a.0.cmp(b.0));
    let mut seen_decoy: HashSet<Vec<u8>> = HashSet::new();
    for (tseq, (dseq, mc)) in shuf_sorted {
        if !seen_decoy.insert(dseq.clone()) {
            continue;
        }
        let mass = seq_mass(dseq);
        let proteins = decoy_prot.get(tseq).cloned().unwrap_or_default();
        let sequence: Arc<str> = Arc::from(std::str::from_utf8(dseq).unwrap_or("").to_owned());
        decoys.push(TheoPeptide {
            sequence,
            mass,
            decoy: true,
            missed_cleavages: *mc,
            proteins,
        });
    }

    targets.extend(decoys);
    targets
}

/// Digest a FASTA text with `cfg`. For protein-level decoy methods, whole decoy
/// proteins are synthesized from the targets, appended tagged "DECOY_", and the
/// combined set is digested with decoy generation OFF: target peptides come out
/// decoy=false, DECOY_ peptides decoy=true, and decoy peptides equal to a target
/// digest are dropped. Returns (peptides, protein_level, decoy tag used).
fn digest_fasta(contents: &str, cfg: &DigestConfig) -> Result<(Vec<Peptide>, bool, String)> {
    use crate::decoy::DecoyMethod;
    let method = DecoyMethod::parse(&cfg.decoy_method, cfg.decoy_order);
    let protein_level = cfg.generate_decoys && method.is_protein_level();
    let (fasta, gen_decoys, decoy_tag) = if protein_level {
        let tfasta = Fasta::parse(contents, "DECOY_", false)?;
        let decoy_prots = crate::decoy::generate(&tfasta.targets, method, "DECOY_", cfg.decoy_seed);
        let mut combined = contents.to_string();
        if !combined.ends_with('\n') {
            combined.push('\n');
        }
        for (acc, seq) in &decoy_prots {
            combined.push('>');
            combined.push_str(acc);
            combined.push('\n');
            combined.push_str(seq);
            combined.push('\n');
        }
        log::info!(
            "protein-level decoys ({:?}): {} target proteins -> {} decoy proteins",
            method,
            tfasta.targets.len(),
            decoy_prots.len()
        );
        (
            Fasta::parse(&combined, "DECOY_", false)?,
            false,
            "DECOY_".to_string(),
        )
    } else {
        (
            Fasta::parse(contents, &cfg.decoy_tag, cfg.generate_decoys)?,
            cfg.generate_decoys,
            cfg.decoy_tag.clone(),
        )
    };
    let params = DigestParams {
        missed_cleavages: cfg.missed_cleavages,
        min_len: cfg.min_len,
        max_len: cfg.max_len,
        cleave_at: cfg.cleave_at.clone(),
        restrict: cfg.restrict.clone(),
        c_terminal: cfg.c_terminal,
        peptide_min_mass: cfg.peptide_min_mass,
        peptide_max_mass: cfg.peptide_max_mass,
        static_mods: cfg.static_mods.clone(),
        variable_mods: cfg.variable_mods.clone(),
        max_variable_mods: cfg.max_variable_mods,
        generate_decoys: gen_decoys,
    };
    Ok((digest::digest(&fasta, &params)?, protein_level, decoy_tag))
}

impl Default for DigestConfig {
    fn default() -> Self {
        let mut static_mods = HashMap::new();
        static_mods.insert("C".to_string(), 57.021_46);
        // DirectMS1/ms1searchpy convention: fully-tryptic (mc=0), NO variable mods.
        // MS1-only cannot localize mods, and oxidized/missed-cleavage variants just
        // inflate the search space and the random-match background. Opt in via config
        // if needed.
        let variable_mods = HashMap::new();
        DigestConfig {
            missed_cleavages: 0,
            min_len: 7,
            max_len: 50,
            cleave_at: "KR".to_string(),
            restrict: "P".to_string(),
            c_terminal: true,
            peptide_min_mass: 500.0,
            peptide_max_mass: 5000.0,
            static_mods,
            variable_mods,
            max_variable_mods: 2,
            decoy_tag: "rev_".to_string(),
            generate_decoys: true,
            decoy_mass_shift: 11.0,
            decoy_shuffle: true,
            decoy_method: "shuffle".to_string(),
            decoy_seed: 0,
            decoy_order: 2,
        }
    }
}

/// A theoretical peptide candidate. Modform identity is captured by
/// (`sequence`, `mass`): two variable-mod forms of one sequence differ in mass.
#[derive(Debug, Clone)]
pub struct TheoPeptide {
    /// Stripped amino-acid sequence (reverse-path decoys have their interior reversed).
    pub sequence: Arc<str>,
    /// Neutral monoisotopic mass (Da), recomputed in f64.
    pub mass: f64,
    pub decoy: bool,
    pub missed_cleavages: u8,
    /// Source protein accessions (shared `Arc`s; cheap to clone).
    pub proteins: Vec<Arc<str>>,
}

/// Mass-sorted peptide database with binary-search candidate lookup.
pub struct PeptideDb {
    /// Sorted ascending by `mass`.
    pub peptides: Vec<TheoPeptide>,
    pub n_targets: usize,
    pub n_decoys: usize,
}

impl PeptideDb {
    /// Build the database from a FASTA file using the given digest config.
    pub fn build(fasta_path: impl AsRef<Path>, cfg: &DigestConfig) -> Result<Self> {
        let fasta_path = fasta_path.as_ref();
        let contents = std::fs::read_to_string(fasta_path)
            .with_context(|| format!("reading fasta {}", fasta_path.display()))?;
        let (peptides_in, protein_level, decoy_tag) = digest_fasta(&contents, cfg)?;
        Ok(Self::assemble(peptides_in, cfg, protein_level, &decoy_tag))
    }

    /// Turn digested peptides into the sorted uno database (decoy handling + f64 masses).
    fn assemble(
        peptides_in: Vec<Peptide>,
        cfg: &DigestConfig,
        protein_level: bool,
        decoy_tag: &str,
    ) -> Self {
        let mut peptides: Vec<TheoPeptide> =
            if cfg.decoy_shuffle && cfg.generate_decoys && !protein_level {
                global_shuffle_decoy_db(peptides_in, cfg)
            } else {
                peptides_in
                    .into_iter()
                    .map(|p| {
                        let mut mass =
                            mass::neutral_mass(&p.sequence, &p.modifications, p.nterm, p.cterm)
                                .unwrap_or(p.monoisotopic as f64);
                        // Protein-level decoys carry real digestion masses; only the
                        // reverse+mass-shift path applies an artificial offset.
                        if p.decoy && !protein_level {
                            mass += cfg.decoy_mass_shift;
                        }
                        let sequence: Arc<str> =
                            Arc::from(std::str::from_utf8(&p.sequence).unwrap_or("").to_owned());
                        // Strip the decoy_tag prefix from decoy protein accessions so they
                        // read IDENTICALLY to their target counterpart's accession (only
                        // the `decoy` flag distinguishes them). Sage's own Fasta::parse
                        // keeps the tag baked into the accession string when reading an
                        // already-decoy-tagged FASTA (generate_decoys=false); without this,
                        // a decoy accession can never equal any target accession, so the
                        // (decoy, accession) pair-competition FDR step never fires for ANY
                        // protein in the run -- silently inflating the reported q-values
                        // (measured: 713 nominal @q<=1% vs 167 true on yeast entrapment).
                        // uno's own decoy_shuffle path already satisfies this by cloning
                        // the target's Arc<str> directly (see global_shuffle_decoy_db).
                        let proteins = if p.decoy {
                            p.proteins
                                .iter()
                                .map(|a| {
                                    a.strip_prefix(decoy_tag)
                                        .map(Arc::from)
                                        .unwrap_or_else(|| a.clone())
                                })
                                .collect()
                        } else {
                            p.proteins
                        };
                        TheoPeptide {
                            sequence,
                            mass,
                            decoy: p.decoy,
                            missed_cleavages: p.missed_cleavages,
                            proteins,
                        }
                    })
                    .collect()
            };

        peptides.sort_by(|a, b| a.mass.partial_cmp(&b.mass).unwrap());

        let n_decoys = peptides.iter().filter(|p| p.decoy).count();
        let n_targets = peptides.len() - n_decoys;
        PeptideDb {
            peptides,
            n_targets,
            n_decoys,
        }
    }

    pub fn len(&self) -> usize {
        self.peptides.len()
    }
    pub fn is_empty(&self) -> bool {
        self.peptides.is_empty()
    }

    /// Candidate peptides whose mass falls in `[lo, hi]` (inclusive), as a slice
    /// of the sorted vector.
    pub fn candidates(&self, lo: f64, hi: f64) -> &[TheoPeptide] {
        let start = self.peptides.partition_point(|p| p.mass < lo);
        let end = self.peptides.partition_point(|p| p.mass <= hi);
        &self.peptides[start..end]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builds_human_ecoli_db() {
        // The two-proteome entrapment FASTA the benchmark uses. Portable and
        // in-repo (gitignored), rather than a path into a sibling checkout.
        let Some(f) = crate::test_input("tests/data/entrap_human_yeast.fasta") else {
            return;
        };
        let cfg = DigestConfig::default();
        let db = PeptideDb::build(&f, &cfg).expect("build db");
        assert!(db.n_targets > 100_000, "targets: {}", db.n_targets);
        assert!(db.n_decoys > 0);

        // sorted ascending
        assert!(db.peptides.windows(2).all(|w| w[0].mass <= w[1].mass));

        // A 10 ppm window around a mid-range mass yields a small candidate set.
        // Probe an actual peptide mass rather than a round number: peptide masses
        // cluster by mass defect, so a nominal 1500.0 falls in the gap between
        // clusters and legitimately has no candidates within 10 ppm.
        let m = db.peptides[db.peptides.len() / 2].mass;
        let tol = m * 10e-6;
        let c = db.candidates(m - tol, m + tol);
        assert!(!c.is_empty());
        assert!(c.iter().all(|p| (p.mass - m).abs() <= tol + 1e-9));
        eprintln!(
            "db: {} targets, {} decoys; {} candidates within 10 ppm of {}",
            db.n_targets,
            db.n_decoys,
            c.len(),
            m
        );
    }
}
