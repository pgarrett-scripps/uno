# KOTH consensus matrix algorithm audit — 2026-09-12

The current matrix route loses information before identification and has
reproducible grouping problems. The evidence does not support treating the
identification deficit as just an overly strict recurrence threshold.

Scope: read the local KOTH and UNO implementations, inspect the existing fresh
IonStar benchmark, and reproduce grouping failure mechanisms using the actual
KOTH consensus functions. No production algorithm or paper results changed;
no identification searches were rerun for this audit.

KOTH `koth_ff/src/lfq/consensus.rs` is byte-identical between the local repository
and the benchmark source snapshot. SHA-256:
`15414204b5f9a45e30962b33c5237bbc7d2f5a03417a639fad92d540c39c0747`.

## Findings, in implementation priority order

1. **Final group filters run before fragmented groups are reunited.**
   `consensus.rs:151` applies minimum group size and minimum seed score;
   `build_consensus` calls it at lines 280 and 289, before `merge_consensus`.
   The first sweep closes a group whenever an intervening feature fails its
   RT/IM window, even if a later mass-sorted feature belongs to the earlier group.
   A discarded fragment cannot contribute to the final merge.

   Reproduction with the benchmark's minimum group size of one: an allowed
   member with score 0.70 and a nearby member from another run with score 0.90
   should form a two-run group whose seed clears 0.75. An unrelated intervening
   feature separates them; the 0.70 fragment disappears before merging, leaving
   one contributing run. UNO's minimum of two then rejects the surviving row.
   With KOTH's default group size of two, even two high-quality singleton
   fragments can both disappear. The benchmark explicitly used KOTH group size
   one, so that second variant does not explain its results.

   Fix: assemble the final groups before applying group-level size/seed filters.
   Keep member eligibility distinct from final group eligibility.

2. **The merge can collapse features beyond its own pairwise tolerances.**
   `merge_consensus` uses transitive union-find connections. A–B and B–C links
   merge all three even if A and C cannot match. There is no final cluster-span
   check. Its pairwise grouping windows are already twice the LFQ extraction
   windows: 20 ppm and 2.84 minutes for the 20-run IonStar reference span.
   The synthetic test merges endpoints 30 ppm and three minutes apart under
   pairwise limits of 20 ppm and two minutes.

   Each resulting cluster retains only one seed mass/RT and one feature per run.
   Consequently, overmerging can remove an alternative identification coordinate
   and misrepresent how many runs support the same precursor. The code establishes
   this possibility; its frequency and contribution to lost benchmark IDs have
   not been measured. The observed reduction from 534,773 pre-merge groups to
   156,267 final rows is not itself proof of erroneous merging.

   Fix: constrain final clusters relative to a representative and/or their
   full coordinate span, and retain conflicting members as separate hypotheses.
   Separate grouping tolerances from extraction tolerances.

3. **The matrix removes feature quality needed by UNO's early calibration.**
   KOTH emits one seed's coordinates and combined score, but does not preserve
   its isotope count, scan count, or separate shape scores in the matrix schema.
   UNO `input.rs:220` sets those descriptors to `None`. Its mass-anchor stage
   (`score.rs:454`) therefore skips the three-isotope eligibility filter used
   with single-file KOTH inputs. Appending extraction-quality summaries later
   in the ML feature vector does not restore that early filter.

   Existing 20-run IonStar logs show bootstrap random-match parameter `p=0.537`
   and 36 confident target proteins for the matrix; the 20 single-file searches
   have median `p=0.3195` and median 159.5 confident proteins. Those different
   inputs prevent attributing the whole gap to isotope information alone, but
   the matrix starts its calibration with a materially weaker anchor set.

   Fix: preserve seed and constituent feature descriptors, including genuine
   original isotope counts. Use them to choose calibration anchors before ML.
   Do not substitute the LFQ grid's isotope count for the original feature's
   isotope count: those are different measurements.

4. **The per-run merge winner is chosen using another feature's score.**
   `merge_consensus` starts with the best-seed group's run slots and only fills
   empty slots from other groups. Group ordering uses each group's seed score,
   not the score of the competing feature in that run. A group with a 0.99 seed
   in run A and a 0.51 feature in run B retains that 0.51 run-B feature even when
   another merged group contains a 0.90 run-B feature. The executable diagnostic
   reproduces this exactly.

   Fix: retain the actual per-run member quality and compare competing members
   directly. In this benchmark `detected_use_grid=true`, so the chosen feature's
   original intensity is not directly exported; do not attribute the observed
   quantitative ratios to this defect without a separate test. It directly
   affects feature selection and the alternative direct-intensity mode.

5. **Current consistency evidence mixes recurrence, extraction background,
   and candidate multiplicity.** UNO's base classifier counts positive LFQ
   intensities and uses their sum. With unmasked matrix extraction nearly every
   cell is positive, which is not independent evidence of detection. Peptide
   consistency statistics at `score.rs:1429` also summarize all its candidate
   mass matches, including RT-inconsistent alternatives, rather than a coherent
   set of observed precursors. Collapsing many original observations into one
   consensus feature changes what those counts mean. These are modeling concerns,
   not a demonstrated classifier rule requiring presence in every run.

   The experimental quality fractions use all cohort runs as denominator;
   condition-specific or load-specific features can therefore have low fractions
   despite strong evidence in their relevant samples. A suitable comparison
   would preserve best-observation quality and absolute confident direct support
   alongside recurrence, with held-out evaluation of any learned cohort effects.

## What the completed ablations establish

At the existing empirical 1% entrapment threshold, the standard 20-run matrix
has 942 present-species protein groups. Allowing singleton rows reduces that to
751; suppressing recurrence and intensity variation while retaining exactly the
baseline rows yields 1,074. The latter is a compound ablation: it also changes
positive-cell count, intensity ranking, and intensity-based tie-breaking.
It does not isolate the recurrence feature. At internal 1% q, present-species
accessions are 938, 739, and 931, respectively.

The default UNO floor removes 40,581 of 156,267 matrix rows (26.0%). Before
grouping, KOTH admits 3,827,321 of 6,403,171 original features at its 0.5 member
score threshold; roughly 40.2% are excluded. These are feature counts, not
counts of true precursors or explanations of how many proteins were lost.

The separate-search union's 2,550 present-species accessions have estimated
entrapment FDP 5.78%, versus 938 and 0.87% for the matrix at internal 1% q.
Nevertheless, even the stricter union at per-file q <= 0.002 retains 1,649
accessions with estimated FDP 0.25%. These are discrete, noisy entrapment
estimates, but indicate that differing error rates alone do not explain the gap.

The 5-, 10-, and 20-run cohorts also select different alignment references.
Their changing depth therefore cannot be assigned solely to increasing run
count or a consistency penalty.

## Recommended bounded follow-up

First correct group filtering order and cluster bounds, preserve constituent
quality, and evaluate one fixed cohort/reference/seed with stage-by-stage
precursor retention counts. Compare protein coverage at matched empirical error
rates. Test recurrence-only ablation separately from intensity suppression.

For a union-like identification route, retain the original per-run mass/RT
observations behind each aligned row. A shared database and cached peptide
predictions can support one batch invocation with run-specific calibration;
combine evidence and estimate cohort error after scoring. A single consensus
coordinate is not equivalent to searching all constituent observations, and
the current results do not guarantee how close a corrected matrix will get.

## Reproduction

`benchmark/matrix_consensus_audit.py` compiles the production consensus functions
with identity-alignment fixtures and no external Rust dependencies. All four
diagnostics passed, confirming the current problematic behaviors. They assert
those behaviors for diagnosis; corrected code should change their assertions.

```sh
python3 benchmark/matrix_consensus_audit.py \
  /home/patrick-garrett/Repos/koth_rust/koth_ff/src/lfq/consensus.rs \
  --output benchmark/results/matrix_20260912/algorithm-audit/consensus-reproductions.json
```

Supporting existing-run measurements are saved in
`benchmark/results/matrix_20260912/algorithm-audit/study-diagnostics.json`;
source data are the cohort `align_report.json`, `consensus_features.tsv`,
search logs, and `summary.json` in the same study directory.
