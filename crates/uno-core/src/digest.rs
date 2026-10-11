//! In-silico digestion: FASTA parsing, enzymatic cleavage, and static/variable
//! modification expansion.
//!
//! This reproduces the peptide set that Sage's core library (`sage-core`, rev
//! 0639176) builds in `Parameters::digest`, so databases built before and after
//! dropping that dependency are identical: same sequences, modifications, f32
//! filter masses, protein mappings and order. Only the subset uno uses is kept
//! (no semi-enzymatic digestion, no fragment index).
//!
//! Behaviour worth knowing, all inherited from Sage:
//! - FASTA accessions are the first whitespace-delimited token of the header.
//! - No N-terminal methionine removal is performed.
//! - When decoys are generated, proteins whose accession contains the decoy tag
//!   are dropped and every target peptide gets a decoy with its interior residues
//!   reversed (termini kept). Decoys whose sequence equals any target digest are
//!   discarded. When decoys are not generated, peptides of proteins whose
//!   accession contains the tag are marked as decoys.
//! - Missed-cleavage windows are counted over the cleavage-site list, and a
//!   sequence seen twice in one protein keeps its first occurrence.
//! - Peptide mass filtering uses f32 monoisotopic masses.

use rayon::prelude::*;
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::sync::Arc;

/// f32 water mass used for the mass filter (Sage's constant).
const H2O_F32: f32 = 18.010565;

/// f32 residue masses indexed by `aa - b'A'`; 0.0 marks an invalid residue.
const MONOISOTOPIC_F32: [f32; 26] = [
    71.03711, 0.0, 103.00919, 115.02694, 129.04259, 147.0684, 57.02146, 137.05891, 113.08406, 0.0,
    128.09496, 113.08406, 131.0405, 114.04293, 237.14774, 97.05276, 128.05858, 156.1011, 87.03203,
    101.04768, 150.95363, 99.06841, 186.07932, 0.0, 163.06332, 0.0,
];

const VALID_AA: &[u8] = b"ACDEFGHIKLMNPQRSTVWYUO";

fn monoisotopic_f32(aa: u8) -> f32 {
    if aa.is_ascii_uppercase() {
        MONOISOTOPIC_F32[(aa - b'A') as usize]
    } else {
        0.0
    }
}

/// Parsed FASTA database.
#[derive(Clone)]
pub struct Fasta {
    /// `(accession, sequence)` for every protein kept.
    pub targets: Vec<(Arc<str>, String)>,
    decoy_tag: String,
    generate_decoys: bool,
}

impl Fasta {
    /// Parse FASTA text. With `generate_decoys`, proteins whose accession contains
    /// `decoy_tag` are dropped (decoys are generated later instead).
    pub fn parse(contents: &str, decoy_tag: &str, generate_decoys: bool) -> anyhow::Result<Fasta> {
        let mut targets = Vec::new();
        let mut last_id = "";
        let mut s = String::new();
        let push =
            |id: &str, seq: String, targets: &mut Vec<(Arc<str>, String)>| -> anyhow::Result<()> {
                let acc = id
                    .split_ascii_whitespace()
                    .next()
                    .ok_or_else(|| anyhow::anyhow!("FASTA entry with an empty header"))?;
                let acc: Arc<str> = Arc::from(acc);
                if !acc.contains(decoy_tag) || !generate_decoys {
                    targets.push((acc, seq));
                }
                Ok(())
            };
        for line in contents.lines() {
            if line.is_empty() {
                continue;
            }
            let line = line.trim();
            if let Some(id) = line.strip_prefix('>') {
                if !s.is_empty() {
                    push(last_id, std::mem::take(&mut s), &mut targets)?;
                }
                last_id = id;
            } else {
                s.push_str(line);
            }
        }
        if !s.is_empty() {
            push(last_id, s, &mut targets)?;
        }
        Ok(Fasta {
            targets,
            decoy_tag: decoy_tag.to_string(),
            generate_decoys,
        })
    }

    fn digest(&self, enzyme: &EnzymeParameters) -> Vec<Digest> {
        self.targets
            .par_iter()
            .flat_map_iter(|(protein, sequence)| {
                let is_tagged = protein.contains(self.decoy_tag.as_str());
                let generate = self.generate_decoys;
                enzyme
                    .digest(sequence, protein.clone())
                    .into_iter()
                    .filter_map(move |mut d| {
                        if is_tagged {
                            if generate {
                                None
                            } else {
                                d.decoy = true;
                                Some(d)
                            }
                        } else {
                            Some(d)
                        }
                    })
            })
            .collect()
    }
}

/// Where a peptide sits in its protein.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Hash, Default)]
pub enum Position {
    Nterm,
    Cterm,
    Full,
    #[default]
    Internal,
}

#[derive(Clone, Debug)]
struct Digest {
    decoy: bool,
    sequence: String,
    protein: Arc<str>,
    missed_cleavages: u8,
    position: Position,
}

struct DigestGroup {
    reference: Digest,
    proteins: Vec<Arc<str>>,
}

fn group_digests(mut digests: Vec<Digest>) -> Vec<DigestGroup> {
    let mut groups = Vec::new();
    if digests.is_empty() {
        return groups;
    }
    digests.sort_unstable_by(|a, b| {
        a.position
            .cmp(&b.position)
            .then(a.decoy.cmp(&b.decoy))
            .then(a.sequence.cmp(&b.sequence))
    });
    let mut curr = DigestGroup {
        reference: digests[0].clone(),
        proteins: Vec::new(),
    };
    for d in digests {
        if d.decoy == curr.reference.decoy
            && d.position == curr.reference.position
            && d.sequence == curr.reference.sequence
        {
            curr.proteins.push(d.protein);
        } else {
            curr.proteins.sort_unstable();
            groups.push(curr);
            curr = DigestGroup {
                proteins: vec![d.protein.clone()],
                reference: d,
            };
        }
    }
    groups.push(curr);
    groups
}

struct Enzyme {
    /// Residues to cleave at; `None` cleaves only at the protein C-terminus ("$").
    cleave: Option<[bool; 256]>,
    skip_suffix: [bool; 256],
    c_terminal: bool,
}

struct EnzymeParameters {
    missed_cleavages: u8,
    min_len: usize,
    max_len: usize,
    /// `None` = non-specific digestion (empty `cleave_at`).
    enzyme: Option<Enzyme>,
}

#[derive(Clone)]
struct Site {
    range: std::ops::Range<usize>,
    missed_cleavages: u8,
}

impl Enzyme {
    fn new(cleave: &str, skip_suffix: &str, c_terminal: bool) -> anyhow::Result<Option<Self>> {
        anyhow::ensure!(
            cleave == "$" || cleave.bytes().all(|b| VALID_AA.contains(&b) || b == b'?'),
            "enzyme cleavage sequence contains non-amino-acid characters: {cleave}"
        );
        anyhow::ensure!(
            skip_suffix.bytes().all(|b| VALID_AA.contains(&b)),
            "enzyme cleavage restriction contains non-amino-acid characters: {skip_suffix}"
        );
        let mut skip = [false; 256];
        for b in skip_suffix.bytes() {
            skip[b as usize] = true;
        }
        Ok(match cleave.replace('?', "").as_str() {
            "" if cleave.is_empty() => None,
            "$" => Some(Enzyme {
                cleave: None,
                skip_suffix: [false; 256],
                c_terminal: true,
            }),
            set => {
                let mut arr = [false; 256];
                for b in set.bytes() {
                    arr[b as usize] = true;
                }
                Some(Enzyme {
                    cleave: Some(arr),
                    skip_suffix: skip,
                    c_terminal,
                })
            }
        })
    }

    fn cleavage_sites(&self, sequence: &str) -> Vec<Site> {
        let bytes = sequence.as_bytes();
        let mut sites = Vec::new();
        let mut left = 0;
        let matches: Vec<(usize, usize)> = match &self.cleave {
            None => vec![(bytes.len(), bytes.len())],
            Some(set) => bytes
                .iter()
                .enumerate()
                .filter(|(_, b)| set[**b as usize])
                .map(|(i, _)| (i, i + 1))
                .collect(),
        };
        for (start, end) in matches {
            let right = if self.c_terminal { end } else { start };
            if bytes
                .get(right)
                .is_some_and(|b| self.skip_suffix[*b as usize])
            {
                continue;
            }
            sites.push(Site {
                range: left..right,
                missed_cleavages: 0,
            });
            left = right;
        }
        sites.push(Site {
            range: left..bytes.len(),
            missed_cleavages: 0,
        });
        sites
    }
}

impl EnzymeParameters {
    fn cleavage_sites(&self, sequence: &str) -> Vec<Site> {
        match &self.enzyme {
            Some(e) => e.cleavage_sites(sequence),
            None => {
                let mut v = Vec::new();
                for len in self.min_len..=self.max_len {
                    for i in 0..=sequence.len().saturating_sub(len) {
                        v.push(Site {
                            range: i..i + len,
                            missed_cleavages: 0,
                        });
                    }
                }
                v
            }
        }
    }

    fn digest(&self, sequence: &str, protein: Arc<str>) -> Vec<Digest> {
        let n = sequence.len();
        let mut sites = self.cleavage_sites(sequence);
        let missed = if self.enzyme.is_none() {
            0
        } else {
            self.missed_cleavages
        };
        if missed > 0 {
            let mut extra = Vec::new();
            for cleavage in 1..=(1 + missed) {
                for win in sites.windows(cleavage as usize) {
                    extra.push(Site {
                        range: win[0].range.start..win[cleavage as usize - 1].range.end,
                        missed_cleavages: cleavage - 1,
                    });
                }
            }
            sites.append(&mut extra);
        }
        let mut seen: HashSet<&str> = HashSet::new();
        let mut digests = Vec::new();
        for site in &sites {
            let (start, end) = (site.range.start, site.range.end);
            let Some(pep) = sequence.get(start..end) else {
                continue;
            };
            let len = pep.len();
            let position = match (start == 0, end == n) {
                (true, true) => Position::Full,
                (true, false) => Position::Nterm,
                (false, true) => Position::Cterm,
                (false, false) => Position::Internal,
            };
            if len >= self.min_len && len <= self.max_len && len > 0 && seen.insert(pep) {
                digests.push(Digest {
                    decoy: false,
                    sequence: pep.to_string(),
                    protein: protein.clone(),
                    missed_cleavages: site.missed_cleavages,
                    position,
                });
            }
        }
        digests
    }
}

/// Where a modification may sit.
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ModSpecificity {
    PeptideN(Option<u8>),
    PeptideC(Option<u8>),
    ProteinN(Option<u8>),
    ProteinC(Option<u8>),
    Residue(u8),
}

impl ModSpecificity {
    /// Parse Sage's notation: `M`, `^` / `^Q` (peptide N-term), `$` (peptide
    /// C-term), `[` (protein N-term), `]` (protein C-term). `None` if invalid.
    pub fn parse(s: &str) -> Option<Self> {
        if s.len() > 2 {
            return None;
        }
        let rest = |r: &str| r.chars().next().map(|c| c as u8);
        if let Some(r) = s.strip_prefix('^') {
            return Some(Self::PeptideN(rest(r)));
        }
        if let Some(r) = s.strip_prefix('$') {
            return Some(Self::PeptideC(rest(r)));
        }
        if let Some(r) = s.strip_prefix('[') {
            return Some(Self::ProteinN(rest(r)));
        }
        if let Some(r) = s.strip_prefix(']') {
            return Some(Self::ProteinC(rest(r)));
        }
        let c = s.chars().next()?;
        VALID_AA
            .contains(&(c as u8))
            .then_some(Self::Residue(c as u8))
    }
}

fn parse_mods<V>(input: &HashMap<String, V>) -> Vec<(ModSpecificity, V)>
where
    V: Clone,
{
    let mut out: Vec<(ModSpecificity, V)> = Vec::new();
    for (s, v) in input {
        match ModSpecificity::parse(s) {
            Some(m) => out.push((m, v.clone())),
            None => log::error!("skipping invalid modification string: {s:?}"),
        }
    }
    // Deterministic order; on a duplicate specificity the later key wins, as in a map.
    out.sort_by_key(|a| a.0);
    out.dedup_by(|later, earlier| {
        if later.0 == earlier.0 {
            std::mem::swap(later, earlier);
            true
        } else {
            false
        }
    });
    out
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum ModSite {
    Nterm,
    Cterm,
    Sequence(u32),
}

/// A digested, modified peptide (target or decoy).
#[derive(Clone, Debug, PartialEq)]
pub struct Peptide {
    pub decoy: bool,
    pub sequence: Arc<[u8]>,
    /// Per-residue modification mass (0.0 = unmodified).
    pub modifications: Vec<f32>,
    pub nterm: Option<f32>,
    pub cterm: Option<f32>,
    /// f32 monoisotopic mass including modifications (used for the mass filter).
    pub monoisotopic: f32,
    pub missed_cleavages: u8,
    pub position: Position,
    /// Sorted source protein accessions.
    pub proteins: Vec<Arc<str>>,
}

impl Peptide {
    fn from_group(g: DigestGroup) -> Option<Peptide> {
        let seq = g.reference.sequence;
        if !seq.is_ascii() {
            return None;
        }
        let mut mass = H2O_F32;
        for &c in seq.as_bytes() {
            let m = monoisotopic_f32(c);
            if m == 0.0 {
                return None;
            }
            mass += m;
        }
        Some(Peptide {
            decoy: g.reference.decoy,
            modifications: vec![0.0; seq.len()],
            sequence: Arc::from(seq.into_bytes().into_boxed_slice()),
            nterm: None,
            cterm: None,
            monoisotopic: mass,
            missed_cleavages: g.reference.missed_cleavages,
            position: g.reference.position,
            proteins: g.proteins,
        })
    }

    fn modification_mass(&self) -> f32 {
        self.modifications.iter().sum::<f32>()
            + self.nterm.unwrap_or(0.0)
            + self.cterm.unwrap_or(0.0)
    }

    fn apply_site(&mut self, site: ModSite, mass: f32) {
        match site {
            ModSite::Nterm => {
                if self.nterm.is_none() {
                    self.nterm = Some(0.0 + mass);
                }
            }
            ModSite::Cterm => {
                if self.cterm.is_none() {
                    self.cterm = Some(0.0 + mass);
                }
            }
            ModSite::Sequence(i) => {
                if self.modifications[i as usize] == 0.0 {
                    self.modifications[i as usize] += mass;
                }
            }
        }
    }

    fn first(&self) -> u8 {
        *self.sequence.first().unwrap_or(&0)
    }
    fn last_idx(&self) -> u32 {
        self.sequence.len().saturating_sub(1) as u32
    }
    fn last(&self) -> u8 {
        *self.sequence.last().unwrap_or(&0)
    }

    /// Sites a terminal specificity targets on this peptide (shared by static and
    /// variable mods). Residue specificities are handled by the callers.
    fn terminal_site(&self, target: ModSpecificity) -> Option<ModSite> {
        use ModSpecificity::*;
        let prot_n = matches!(self.position, Position::Nterm | Position::Full);
        let prot_c = matches!(self.position, Position::Cterm | Position::Full);
        match target {
            PeptideN(None) => Some(ModSite::Nterm),
            PeptideN(Some(r)) if r == self.first() => Some(ModSite::Sequence(0)),
            PeptideC(None) => Some(ModSite::Cterm),
            PeptideC(Some(r)) if r == self.last() => Some(ModSite::Sequence(self.last_idx())),
            ProteinN(None) if prot_n => Some(ModSite::Nterm),
            ProteinN(Some(r)) if prot_n && r == self.first() => Some(ModSite::Sequence(0)),
            ProteinC(None) if prot_c => Some(ModSite::Cterm),
            ProteinC(Some(r)) if prot_c && r == self.last() => {
                Some(ModSite::Sequence(self.last_idx()))
            }
            _ => None,
        }
    }

    fn push_sites(&self, acc: &mut Vec<(ModSite, f32)>, target: ModSpecificity, mass: f32) {
        if let ModSpecificity::Residue(r) = target {
            acc.extend(
                self.sequence
                    .iter()
                    .enumerate()
                    .filter(|(_, &aa)| aa == r)
                    .map(|(i, _)| (ModSite::Sequence(i as u32), mass)),
            );
        } else if let Some(site) = self.terminal_site(target) {
            acc.push((site, mass));
        }
    }

    fn apply_static(&mut self, target: ModSpecificity, mass: f32) {
        if let ModSpecificity::Residue(r) = target {
            for i in 0..self.sequence.len() {
                if self.sequence[i] == r && self.modifications[i] == 0.0 {
                    self.modifications[i] = mass;
                }
            }
        } else if let Some(site) = self.terminal_site(target) {
            self.apply_site(site, mass);
        }
    }

    /// Expand variable modifications (up to `max_mods` sites), then apply static
    /// modifications to every form and add the modification mass.
    fn apply(
        self,
        variable: &[(ModSpecificity, f32)],
        fixed: &[(ModSpecificity, f32)],
        max_mods: usize,
    ) -> Vec<Peptide> {
        let mut modified = vec![self.clone()];
        if !variable.is_empty() {
            let mut sites = Vec::new();
            for &(target, mass) in variable {
                self.push_sites(&mut sites, target, mass);
            }
            for n in 1..=max_mods {
                for_each_combination(sites.len(), n, |idx| {
                    let mut used = HashSet::new();
                    if !idx.iter().all(|&i| used.insert(sites[i].0)) {
                        return;
                    }
                    let mut p = self.clone();
                    for &i in idx {
                        p.apply_site(sites[i].0, sites[i].1);
                    }
                    modified.push(p);
                });
            }
        }
        for p in modified.iter_mut() {
            for &(target, mass) in fixed {
                p.apply_static(target, mass);
            }
            p.monoisotopic += p.modification_mass();
        }
        modified
    }

    /// Decoy with the interior residues (and their mods) reversed.
    fn reverse(&self) -> Peptide {
        let mut p = self.clone();
        p.decoy = !self.decoy;
        let n = p.sequence.len().saturating_sub(1);
        if n > 1 {
            let mut s = p.sequence.to_vec();
            s[1..n].reverse();
            p.sequence = Arc::from(s.into_boxed_slice());
            p.modifications[1..n].reverse();
        }
        p
    }

    fn sort_key_cmp(&self, other: &Self) -> Ordering {
        self.monoisotopic
            .total_cmp(&other.monoisotopic)
            .then_with(|| self.sequence.cmp(&other.sequence))
            .then_with(|| {
                self.modifications
                    .partial_cmp(&other.modifications)
                    .unwrap_or(Ordering::Equal)
            })
            .then_with(|| {
                self.nterm
                    .partial_cmp(&other.nterm)
                    .unwrap_or(Ordering::Equal)
            })
            .then_with(|| {
                self.cterm
                    .partial_cmp(&other.cterm)
                    .unwrap_or(Ordering::Equal)
            })
    }
}

/// Call `f` with every `k`-subset of `0..n`, in lexicographic order.
fn for_each_combination(n: usize, k: usize, mut f: impl FnMut(&[usize])) {
    if k == 0 || k > n {
        return;
    }
    let mut idx: Vec<usize> = (0..k).collect();
    loop {
        f(&idx);
        let mut i = k;
        while i > 0 && idx[i - 1] == n - k + i - 1 {
            i -= 1;
        }
        if i == 0 {
            return;
        }
        idx[i - 1] += 1;
        for j in i..k {
            idx[j] = idx[j - 1] + 1;
        }
    }
}

/// Digestion parameters (Sage database block semantics).
#[derive(Debug, Clone)]
pub struct DigestParams {
    pub missed_cleavages: u8,
    pub min_len: usize,
    pub max_len: usize,
    /// Residues to cleave at; "" = non-specific, "$" = no cleavage.
    pub cleave_at: String,
    /// Residues that block cleavage when they follow the site.
    pub restrict: String,
    pub c_terminal: bool,
    pub peptide_min_mass: f32,
    pub peptide_max_mass: f32,
    pub static_mods: HashMap<String, f32>,
    pub variable_mods: HashMap<String, Vec<f32>>,
    pub max_variable_mods: usize,
    pub generate_decoys: bool,
}

/// Digest `fasta` into the mass-sorted, deduplicated target (+ reversed decoy)
/// peptide list.
pub fn digest(fasta: &Fasta, params: &DigestParams) -> anyhow::Result<Vec<Peptide>> {
    let enzyme = EnzymeParameters {
        missed_cleavages: params.missed_cleavages,
        min_len: params.min_len,
        max_len: params.max_len,
        enzyme: Enzyme::new(&params.cleave_at, &params.restrict, params.c_terminal)?,
    };
    let fixed = parse_mods(&params.static_mods);
    let variable: Vec<(ModSpecificity, f32)> = parse_mods(&params.variable_mods)
        .into_iter()
        .flat_map(|(t, ms)| ms.into_iter().map(move |m| (t, m)))
        .collect();
    let max_mods = params.max_variable_mods.max(1);

    let groups = group_digests(fasta.digest(&enzyme));
    let targets: HashSet<Vec<u8>> = groups
        .iter()
        .filter(|g| !g.reference.decoy)
        .map(|g| g.reference.sequence.as_bytes().to_vec())
        .collect();

    let mut peptides: Vec<Peptide> = groups
        .into_par_iter()
        .filter_map(Peptide::from_group)
        .flat_map_iter(|p| {
            let mut out = Vec::new();
            for p in p.apply(&variable, &fixed, max_mods) {
                if p.monoisotopic < params.peptide_min_mass
                    || p.monoisotopic > params.peptide_max_mass
                {
                    continue;
                }
                if params.generate_decoys {
                    let d = p.reverse();
                    if !d.decoy || !targets.contains(&d.sequence[..]) {
                        out.push(d);
                    }
                }
                if !p.decoy || !targets.contains(&p.sequence[..]) {
                    out.push(p);
                }
            }
            out
        })
        .collect();

    peptides.par_sort_unstable_by(|a, b| a.sort_key_cmp(b));
    peptides.dedup_by(|remove, keep| {
        if remove.monoisotopic == keep.monoisotopic
            && remove.sequence == keep.sequence
            && remove.modifications == keep.modifications
            && remove.nterm == keep.nterm
            && remove.cterm == keep.cterm
        {
            keep.proteins.extend(remove.proteins.iter().cloned());
            keep.decoy &= remove.decoy;
            true
        } else {
            false
        }
    });
    peptides
        .par_iter_mut()
        .for_each(|p| p.proteins.sort_unstable());
    Ok(peptides)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params() -> DigestParams {
        DigestParams {
            missed_cleavages: 0,
            min_len: 2,
            max_len: 50,
            cleave_at: "KR".into(),
            restrict: "P".into(),
            c_terminal: true,
            peptide_min_mass: 0.0,
            peptide_max_mass: 1e6,
            static_mods: HashMap::new(),
            variable_mods: HashMap::new(),
            max_variable_mods: 2,
            generate_decoys: false,
        }
    }

    #[test]
    fn trypsin_sites() {
        let fasta = Fasta::parse(
            ">P1 desc\nMADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGN\n",
            "rev_",
            false,
        )
        .unwrap();
        let peps = digest(&fasta, &params()).unwrap();
        let mut seqs: Vec<String> = peps
            .iter()
            .map(|p| String::from_utf8(p.sequence.to_vec()).unwrap())
            .collect();
        seqs.sort();
        assert_eq!(
            seqs,
            ["LPPGWEK", "MADEEK", "MSR", "SSGR", "VYYFNHITNASQWERPSGN"]
        );
        assert!(peps.iter().all(|p| &*p.proteins[0] == "P1"));
    }

    #[test]
    fn variable_mods_and_decoys() {
        let fasta = Fasta::parse(">P1\nAMCMGK\n", "rev_", true).unwrap();
        let mut p = params();
        p.generate_decoys = true;
        p.static_mods.insert("C".into(), 57.02146);
        p.variable_mods.insert("M".into(), vec![15.9949]);
        let peps = digest(&fasta, &p).unwrap();
        // 1 + 2 + 1 modforms, each with a reversed decoy.
        assert_eq!(peps.iter().filter(|p| !p.decoy).count(), 4);
        assert_eq!(peps.iter().filter(|p| p.decoy).count(), 4);
        assert!(peps
            .iter()
            .all(|p| p.modifications[2] == 57.02146 || p.decoy));
    }

    #[test]
    fn combinations_are_lexicographic() {
        let mut v = Vec::new();
        for_each_combination(4, 2, |c| v.push(c.to_vec()));
        assert_eq!(
            v,
            vec![
                vec![0, 1],
                vec![0, 2],
                vec![0, 3],
                vec![1, 2],
                vec![1, 3],
                vec![2, 3]
            ]
        );
    }
}
