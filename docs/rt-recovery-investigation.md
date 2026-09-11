# Broader RT calibration: recovery and mechanism investigation

**Investigation complete: the persistent broad RT error floor contributed to two selected failures. Removing it recovered identifications, but did not establish a reliable replacement.** Separately, the recovery-only prototype reproduced three known rescues and preserved every tested successful V3 result; its new-draw screen had no recovery opportunities. V3 remains the default.

## Recovery-only implementation

The new candidate first runs complete V3, including V3's existing evidence-retention retry. If V3 reports any target at q≤1%, it returns the complete original result. Only after a complete failure does it run the unchanged, previously tested broad RT algorithm. It accepts a broad result only when it reports target discoveries; otherwise the complete V3 result is retained. Protein lists and peptide–feature evidence stay together. RT predictions are reused, while failed searches can incur an additional full scoring search.

This changes when broad RT runs. It does not alter the broad fit, its error floor, the classifier, decoys or final reporting threshold. It also cannot stabilize counts among runs that already succeed. Source and binary are isolated in the [candidate archive](../benchmark/results/rt_recovery_investigation/candidate/manifest.json), with the [exact protocol](../benchmark/results/rt_recovery_investigation/candidate/PROTOCOL.md) and [frozen study plan](../benchmark/results/rt_recovery_investigation/PLAN.md).

## Executed results

| Development replay | V3 targets | Recovery targets | Complete output matches |
|---|---:|---:|---|
| 1ng rep4, seed1 | 0 | 272 | Previous broad RT |
| 200pg rep1, seed8 | 0 | 208 | Previous broad RT |
| 200pg rep3, seed5 | 0 | 410 | Previous broad RT |
| 200pg rep3, seed0 | 220 | 220 | V3 |

All four fresh candidate executions passed exact-output and activation checks. These cases were selected from known prior results and verify implementation correspondence, not independent efficacy. [Replay verification](../benchmark/results/rt_recovery_investigation/development-verification.json).

The prospective screen ran **both engines fresh on six low-input files at seeds 9 and 10: 24 searches, 12 pairs**. V3 already reported targets in all 12 cases. Recovery therefore never activated, and every complete candidate output was byte-identical to V3. Both arms reported a summed 4,614 human and 11 yeast targets, pooled entrapment FDP 0.811%, and median relative target-count range 54.000% across the same six files. There were no losses or changed IDs. Two seeds per known file are correlated algorithmic draws, not biological replication.

The formal screen failed only its required-improvement criterion: empty lists stayed 0→0 and count variability was unchanged. This is **absence of a recovery test opportunity**, not evidence that an activated rescue failed. All preservation/error/depth guards passed; the both-empty preservation check had zero applicable cases. Under the frozen stopping rule, the conditional primary/raw/additional panels were not launched. [Formal report](../benchmark/results/rt_recovery_investigation/screen-report.md) and [full results](../benchmark/results/rt_recovery_investigation/screen-report.json).

![Every new-draw pair](../benchmark/results/rt_recovery_investigation/new-draw-counts.png)

[PDF figure](../benchmark/results/rt_recovery_investigation/new-draw-counts.pdf). Lines connect matching seeds. Panel scales differ; all data points are shown. The figure displays exact preservation and the remaining differences between seeds.

Error remains a separate concern: three of 12 individual lists exceeded 1% empirical entrapment FDP in both arms; 200pg rep1 seed10 had 120 human and 5 yeast IDs, estimated FDP 14.18%. The pooled result does not establish per-run nominal FDR. This was an unchanged V3 output, not an effect of broad recovery.

For context, applying the new selection policy retrospectively to all 30 old V3/RT pairs gives 5→2 empty lists, three rescues and no newly empty lists, preserving all 25 V3-positive outputs. Pooled FDP changes 0.610%→0.622%; the rescues alone contain 888 human and 2 yeast IDs, FDP 0.767%, with one rescue individually at 1.669%. Relative count range changes 125.066%→84.483% across the same six files. This is explicitly **post hoc selection**, not 30 fresh hybrid executions. [Independent retrospective analysis](../benchmark/results/rt_recovery_investigation/selection_analysis/findings.md).

## Why the global RT replacement could regress

The [source and trace audit](../benchmark/results/rt_recovery_investigation/mechanism_audit/findings.md) separated three coupled mechanisms in the old broad RT implementation:

1. Its uncertainty floor affects normalized RT penalties and training evidence throughout the search. On the problematic 200pg rep3 seeds 0 and 8, it exceeded 1 minute. The nominal floor 0/floor 1 retries consequently used identical effective floors, reducing their intended diversity; the 1.5× candidate widened the floor again.
2. The initial affine RT map stays frozen through later stages, preventing later confident matches from correcting it. An inherited minimum-anchor gate still determines whether later updates run.
3. All 38 archived selected fits could geometrically intersect only six of eight RT strata within the fitting rectangle. Thus six supported strata did not establish coverage across the full predicted-RT range. The normalized slope cap and the uncertainty estimate are also coupled.

These verified implementation facts identify mechanisms to isolate; they do not by themselves prove what caused the zero-ID results. The separate floor diagnostic changes only the applied broad error floor to zero while retaining the original broad map, freezing policy and V3 retry floors. It tests exactly the two previously problematic cases and is diagnostic-only, outside the recovery screen and its adoption gates.

## Isolated floor experiment

The floor-only candidate completed both prespecified cases, 200pg rep3 seeds 0 and 8, against verified original broad-RT outputs. Restoring its single changed assignment reproduces the entire original broad-RT score source byte-for-byte. Initial map and computed-floor diagnostic messages were identical; only the applied broad floor changed. No parameters were tuned after the runs.

| Seed | Original broad RT targets | Floor removed: human / yeast | Entrapment depth at 1% / 5% | Estimated FDP |
|---|---:|---:|---:|---:|
| 0 | 0 | 140 / 1 | 139 / 206 | 2.43% |
| 8 | 0 | 222 / 0 | 515 / 970 | 0% observed |

The corresponding actual stage RT sigmas decreased from about 1.76 to 0.56 minutes and from 1.95 to 0.51 minutes. Both candidates succeeded before the seeded floor retries or evidence-retention pass; the 1.5× sigma alternatives still returned zero targets. Thus this directly demonstrates a consequential floor effect in the two selected cases, but does not directly measure restored retry diversity.

This is not yet an adoption result. Seed 0 still has fewer targets than V3's 220 and estimated FDP above 1%. No observed yeast in seed 8 does not prove zero error. These deliberately selected known failures cannot establish general stability or preserve the benefits seen on other files. The floor correction and the recovery-only policy were tested **separately**; their combination has not been validated. [Diagnostic protocol](../benchmark/results/rt_recovery_investigation/floor_diagnostic/PROTOCOL.md), [paired findings](../benchmark/results/rt_recovery_investigation/floor_diagnostic/diagnostic-report.md) and [complete provenance](../benchmark/results/rt_recovery_investigation/floor_diagnostic/finalized-manifest.json).

## Validation and status

Recovery candidate: 36 required-data core tests and 5 CLI tests passed, along with release/help checks. Five focused report-preservation tests passed. The [independent screen audit](../benchmark/results/rt_recovery_investigation/independent-screen-audit-notes.md) verified all 24 fresh executions, 438 unique artifact hashes, 25 frozen-plan entries, 56 candidate source files, 13 artifacts and 153 external dependency files. Every successful V3 result was preserved in the four development checks and 12 new pairs; the three known empty cases were rescued.

The floor diagnostic passed 31 required-data core tests and release/help checks. An additional [root verification](../benchmark/results/rt_recovery_investigation/floor-and-preservation-verification.json) independently recomputed its metrics, verified the one-assignment source difference and 309 hashes, and confirmed the main sources, external dependencies and default/recovery binaries were preserved. Overall execution: 16 fresh recovery-candidate searches, 12 fresh V3 searches, and two fresh floor-only searches, plus six verified archived reference outputs. Neither experiment was promoted.

Recovery candidate SHA256: `fc31d06f6c62e2d2bdf0c994a71f8d6816498dd80b0d1445cc3d5edef00795fb`; floor-only candidate: `ad4223e92ea9aaf297195154ffe0ed76e3a8034886dcff8180922ae879e828c7`. Default V3 remains `50b4227bca68a58283368d7189dd65dbbd54854f340d05d2ca355824c4ad1738`. No main engine or paper figures were replaced. The concrete prototypes and the newly isolated mechanism are retained; the evidence does not yet justify a default change.
