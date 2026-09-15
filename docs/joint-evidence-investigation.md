# Joint evidence: investigation and next experiment

Pooling evidence has a measurable signal worth pursuing. An offline check of
saved protein scores improved accession coverage in all four disjoint five-run
IonStar subsets, without running new searches or changing production code.
This supports combining evidence across files, but does not validate the proposed
peptide-level model.

## What the saved results show

The table counts human/E. coli accessions after expanding same-species protein
groups. Each ranking uses the largest complete score-threshold set with empirical
yeast-entrapment FDP at most 1%; these are development evaluation thresholds,
not calibrated cohort q-values. The separate-search baseline ranks accessions by
their minimum per-run q-value and tightens the union threshold to the same
empirical error criterion.

| Frozen five-run subset | Separate-search union | Sum of two strongest protein scores | Sum of all five protein scores |
| --- | ---: | ---: | ---: |
| 1 | 1,820 | 2,093 | 2,079 |
| 2 | 1,790 | 1,773 | 1,891 |
| 3 | 1,614 | 1,699 | 1,773 |
| 4 | 1,387 | 1,529 | 1,470 |

Summing all scores gains 6–14% across these subsets. In subset 1 it includes
200 accessions absent from every individual run's q ≤ 1% list, while losing
220 from that unadjusted union. Pooling therefore rescues weak evidence but does
not preserve every original call. Top-two pooling gains 15% over the matched-error
baseline in subset 1, but slightly loses in subset 2; its apparent pilot advantage
does not establish it as the better general rule.

These scores have already undergone per-run inference and competition, so sums
are ranking diagnostics, not probability calculations. The rules were fixed
before examining subsets 2–4, but all subsets belong to the same previously used
development dataset. Each threshold contains only 2–5 yeast accessions, making
the empirical error estimates coarse. No claim of externally validated improvement
or protein-group FDR control follows from these accession counts.

## Why recurrence alone is insufficient

Among 4,643,765 saved, unthresholded peptide-candidate rows, 56.9% of distinct
present-species candidates and 56.5% of yeast candidates recur in all five files.
Even requiring the same KOTH consensus group across all five leaves 5,740 yeast
and 28,605 shuffled-decoy candidate sequences. These are candidate assignments,
not accepted identifications. A persistent real peak can repeatedly match the
same wrong peptide, so repeated measurements cannot simply be treated as
independent confirmations of sequence identity.

The coordinate join found no unmatched candidates and excluded two candidate
rows affected by one rounded-coordinate collision. Sequence/decoy identity is
valid for this export: the audit rejects within-run duplicates, and found none.
A general exporter must retain modification-aware peptide and observation IDs.

## Where the implementation loses evidence

`long_input.rs::load_consensus_with_representative` creates one feature per
consensus; `score.rs::run_ml` then chooses one best feature per peptide.
Its additional group statistics summarize all mass-compatible matches, rather
than a coherent set of observations. The final `PeptideHit` exposes only peptide
index, feature index and decoy status, losing the ML score and calibrated errors.

The existing `UNO_DUMP_ML` diagnostic is unsuitable for collecting batch evidence:
each evaluated calibration and run overwrites the same path, so the file can
describe the last evaluated calibration rather than the selected winner.
Evidence must travel with `FinalizedCandidate` and be exported after selection.

## Recommended bounded prototype

1. Export the selected calibration's per-run peptide candidates, including scores,
   mass/RT residuals, original observation IDs and target/decoy identity. Preserve
   the existing per-run outputs exactly and reuse the working batch cache.
2. First pool these ranked candidates before protein inference. Count each unique
   peptide once per protein; complementary peptides from different files can add
   evidence without pretending repeated observations are new peptides. Compare
   this baseline against a separate, capped recurrence contribution using coherent
   original observations; keep singleton candidates and exclude MBR-only cells
   from direct identification evidence.
3. Recompute protein inference and target/decoy competition for the cohort, with
   one common grouping policy. Do not multiply per-run q-values or interpret ML
   ranks as probabilities. Validate duplicate-file invariance, absent-run handling,
   correlated wrong candidates, and target/decoy symmetry. Any new learned model
   must hold related peptide observations together across runs during validation.
4. Compare both prototypes with the union at matched empirical error, reporting
   rescues, losses and runtime. Freeze the method before testing another dataset;
   measure quantification separately after cohort identification is established.

There is a relevant baseline in
[ms1combine](https://github.com/markmipt/ms1searchpy/blob/master/ms1searchpy/combine.py):
it combines per-run peptide-match ML tables, orders run-specific score bins by
their decoy fractions, and reruns protein inference. This is a useful comparator
for pooled inference, rather than evidence that multiplying replicate match
probabilities is valid.

Reproduce with `python3 benchmark/joint_evidence_audit.py` followed by
`python3 benchmark/joint_evidence_validation.py`. Results and input hashes are in
`benchmark/results/lfq_matrix_fixes/joint_evidence_audit/`; the independent integer
implementation reproduces all four pilot rankings at all three audited thresholds.
Only these audit scripts and this note were added; no paper or search behavior
changed.
