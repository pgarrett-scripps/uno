//! Protein grouping for group-aware ID and quant.
//!
//! Highly conserved proteins (isoforms, orthologs, paralogs) share most or all of
//! their peptides. uno's greedy inference otherwise reports a single arbitrary
//! "leading" protein and silently drops its indistinguishable twins — and, worse
//! for quant, would attribute a shared peptide's intensity (the SUM over every
//! protein that produces it) to that one leading protein.
//!
//! This groups proteins that are **indistinguishable** by the observed data — i.e.
//! that map to the identical set of observed peptides — and reports them as one
//! group. For quant it also marks each observed peptide **group-unique** (maps only
//! within a single group) vs **cross-group shared** (ambiguous, its intensity is a
//! blend across groups). Only group-unique peptides are trustworthy for quant;
//! shared ones are excluded. Principle: quantify at the level the data can resolve.

use std::collections::{BTreeSet, HashMap, HashSet};

use crate::database::PeptideDb;
use crate::score::PeptideHit;

/// A set of indistinguishable proteins (identical observed-peptide sets).
#[derive(Debug, Clone)]
pub struct ProteinGroup {
    /// Representative accession (a confident member if any, else the first).
    pub leading: String,
    /// All indistinguishable member accessions (includes `leading`), sorted.
    pub members: Vec<String>,
    pub decoy: bool,
}

/// Build protein groups from the observed (matched) peptides.
///
/// Returns `(groups, pep_group)` where `pep_group[peptide_idx] = (decoy, leading)`
/// only for **group-unique** peptides — those that map within a single group and
/// are therefore usable for that group's quant. Cross-group-shared peptides are
/// absent from `pep_group` (excluded from quant).
pub fn build_groups(
    db: &PeptideDb,
    observed: &[PeptideHit],
    confident: &HashSet<String>,
) -> (Vec<ProteinGroup>, HashMap<u32, (bool, String)>) {
    // Per protein: the set of observed peptides mapping to it (target/decoy kept
    // separate). Only proteins that carry observed evidence appear here.
    let mut prot_peps: HashMap<(bool, &str), BTreeSet<u32>> = HashMap::new();
    for h in observed {
        for prot in &db.peptides[h.peptide_idx as usize].proteins {
            prot_peps
                .entry((h.decoy, prot.as_ref()))
                .or_default()
                .insert(h.peptide_idx);
        }
    }

    // Group proteins by identical observed-peptide signature.
    let mut sig_group: HashMap<(bool, Vec<u32>), usize> = HashMap::new();
    let mut groups: Vec<ProteinGroup> = Vec::new();
    let mut prot_group: HashMap<(bool, &str), usize> = HashMap::new();
    for ((decoy, acc), peps) in &prot_peps {
        let sig = (*decoy, peps.iter().copied().collect::<Vec<_>>());
        let gid = *sig_group.entry(sig).or_insert_with(|| {
            groups.push(ProteinGroup {
                leading: String::new(),
                members: Vec::new(),
                decoy: *decoy,
            });
            groups.len() - 1
        });
        groups[gid].members.push((*acc).to_string());
        prot_group.insert((*decoy, *acc), gid);
    }
    for g in &mut groups {
        g.members.sort();
        g.leading = g
            .members
            .iter()
            .find(|m| confident.contains(*m))
            .cloned()
            .unwrap_or_else(|| g.members[0].clone());
    }

    // Groups each observed peptide touches (deduped).
    let mut pep_gids: HashMap<u32, Vec<usize>> = HashMap::new();
    for h in observed {
        let mut set: HashSet<usize> = HashSet::new();
        for prot in &db.peptides[h.peptide_idx as usize].proteins {
            if let Some(&gid) = prot_group.get(&(h.decoy, prot.as_ref())) {
                set.insert(gid);
            }
        }
        // Sorted, so `gids[0]` is a defined element rather than whichever one the
        // HashSet happened to yield first. The `max_by` below already tie-breaks
        // deterministically; this covers the `unwrap_or(&gids[0])` fallback, which
        // is inert today only because its result is discarded at the `confident`
        // check a few lines down.
        pep_gids.entry(h.peptide_idx).or_insert_with(|| {
            let mut v: Vec<usize> = set.into_iter().collect();
            v.sort_unstable();
            v
        });
    }
    // Unique-peptide count per group (peptides touching exactly one group).
    let mut uniq_count: HashMap<usize, usize> = HashMap::new();
    for gids in pep_gids.values() {
        if gids.len() == 1 {
            *uniq_count.entry(gids[0]).or_default() += 1;
        }
    }
    // Assign each peptide to ONE group: unique peptides to their group; shared
    // (razor) peptides to the group with the most unique evidence (deterministic
    // tie-break by leading accession). One group per peptide -> no double-counting.
    // Tested universe = confident target groups + all decoy groups.
    let mut pep_group: HashMap<u32, (bool, String)> = HashMap::new();
    for h in observed {
        let gids = &pep_gids[&h.peptide_idx];
        if gids.is_empty() {
            continue;
        }
        let gid = *gids
            .iter()
            .filter(|&&g| h.decoy || confident.contains(&groups[g].leading))
            .max_by(|&&a, &&b| {
                uniq_count
                    .get(&a)
                    .copied()
                    .unwrap_or(0)
                    .cmp(&uniq_count.get(&b).copied().unwrap_or(0))
                    .then(groups[b].leading.cmp(&groups[a].leading))
            })
            .unwrap_or(&gids[0]);
        let leading = &groups[gid].leading;
        if h.decoy || confident.contains(leading) {
            pep_group.insert(h.peptide_idx, (h.decoy, leading.clone()));
        }
    }

    (groups, pep_group)
}

/// Map each group's `(decoy, leading accession)` to its full member list (for output).
///
/// Keyed on the decoy flag as well as the accession, and it must be: uno's decoy
/// peptides carry the SAME protein accessions as the targets they were shuffled
/// from (see `database.rs`), so a target group and a decoy group routinely share a
/// leading accession. Keyed on the accession alone, `collect()` silently kept
/// whichever of the two came later in `groups` -- and `groups` is built by
/// iterating a `HashMap`, so that was decided by the process hash seed.
///
/// This was the last source of run-to-run difference in `--out-proteins`: two runs
/// agreed on every protein, score, peptide count and q-value, and disagreed on the
/// group members of ~1% of rows, because those rows were being handed a decoy
/// group's membership at random.
pub fn members_by_leading(groups: &[ProteinGroup]) -> HashMap<(bool, String), Vec<String>> {
    groups
        .iter()
        .map(|g| ((g.decoy, g.leading.clone()), g.members.clone()))
        .collect()
}
