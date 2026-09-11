# Confidence-weighted startup calibration experiment

**Complete: reject this prototype and retain V3.** Thirty fresh candidate searches were compared with verified V3 results on six low-input files, five seeds per file. One additional identical-seed execution checked reproducibility. All searches completed successfully as processes; “empty” below means no reported target identifications at q≤1%.

| Measure across 30 paired runs | V3 | Soft startup |
|---|---:|---:|
| Empty reported target lists | 5 | 21 |
| Total reported targets, summed across runs | 10,626 | 4,328 |
| Reported human / yeast targets | 10,607 / 19 | 4,316 / 12 |
| Pooled entrapment FDP | 0.610% | 0.946% |

There were **16 newly empty runs and no rescues**. All five seeds of 1ng rep2, 1ng rep4, 200pg rep2 and 200pg rep3 produced empty candidate target lists. Their median reported human counts and entrapment-based identification depths fell to zero. These sums count repeated identifications across runs, not unique proteins across the study.

Four of the five prespecified gates failed: new empty runs, median identification-depth loss, median reported-human loss, and improvement in failures or reported-count variability. The pooled FDP increase of 0.337 percentage points passed the fixed 0.5-point guard, but does not offset the identification collapse. Four of nine nonempty candidate lists had individual entrapment FDP above 1%; FDP for the 21 empty lists is undefined.

The common-file median relative target-count range was 108.996% in both arms **only for the two files with positive medians in both arms**, 1ng rep1 and 200pg rep1. The other four files have undefined candidate relative ranges because all counts are zero. This is not evidence of unchanged or improved overall stability. The [formal report](../benchmark/results/soft_bootstrap/soft-stage-a-report.md) and [full metrics](../benchmark/results/soft_bootstrap/soft-stage-a-report.json) retain exclusions and all paired results.

![All five seeds per file, including empty outputs](../benchmark/results/soft_bootstrap/stage-a-reported-counts.png)

The [PDF figure](../benchmark/results/soft_bootstrap/stage-a-reported-counts.pdf) shows all seeds, with bars for medians and lines spanning each range. Paper figures were not changed by this experiment.

## What was tested

The frozen prototype estimated a bounded, monotone target/decoy score-to-error curve from the original decoy population. It transformed confidence into protein budgets, distributed those budgets over peptide–feature matches, and used weighted robust mass and affine retention-time (RT) fits for startup. Support checks could reject a fit and execute the original V3 startup instead. Final reported q≤1% was unchanged.

These were heuristic calibration weights, including confidence squaring and crowding adjustments; they were **not validated posterior error probabilities**. Mass and RT startup changed together. The experiment therefore tests this specific combined implementation, not every possible soft-confidence method. The [exact protocol](../benchmark/results/soft_bootstrap/candidate/PROTOCOL.md), [prospective plan](../benchmark/results/soft_bootstrap/PLAN.md), [frozen provenance](../benchmark/results/soft_bootstrap/plan-manifest.json) and [separate worktree record](../benchmark/results/soft_bootstrap/worktree.json) preserve the design. No algorithm or threshold was retuned after scientific results were observed.

## What we learned

Soft calibration installed in **20 runs, all of which returned zero target identifications**. The remaining ten runs always fell back to V3 and matched their V3 output files byte-for-byte. Installed RT fits only exceeded constant-fit support by factors of 1.058–1.084 despite initial feature effective sample sizes of 1,707–2,985. Their median fitted mass sigma was 1.060 ppm versus 0.466 ppm in the corresponding hard bootstrap. Large nominal support did not ensure useful alignment.

Source inspection identifies a plausible failure mechanism: protein-level confidence does not establish that each peptide was matched to the correct feature. Startup selects the closest mass match per peptide, which can be an accidental feature. The RT fitter then receives the original broad-window weights without refiltering or reweighting matches around the fitted mass residual. Those ambiguities can persist despite confident protein evidence. This is consistent with contaminated calibration anchors, but the combined experiment does not isolate causality. [Mechanism findings and parameter ranges](../benchmark/results/soft_bootstrap/mechanism-stage-a/mechanism-findings.md) retain the supporting traces.

The useful design lesson is to assess confidence in the peptide–feature correspondence and joint mass/RT consistency before treating a smooth protein-confidence curve as calibration evidence. This study does not establish that a looser threshold or another fitted curve would solve the instability.

## Verification and decision

The frozen implementation passed 35 required-data core tests, five CLI tests and release/help checks. Five report checks and five diagnostic-parser checks also passed. The [independent audit](../benchmark/results/soft_bootstrap/independent-final-audit-notes.md) recomputed all five gates, checked all 30 baseline metric rows and output bytes against archived V3, and verified 398 unique artifact hashes and 13 frozen-plan entries.

The [identical-seed repeat](../benchmark/results/soft_bootstrap/repeat-verification.json) produced exactly the same complete output as its Stage-A run and V3: 655 reported targets. Both candidate executions fell back, so this verifies that executed fallback path only; it does not establish repeatability of an installed soft fit or add an independent efficacy replicate.

The [preservation audit](../benchmark/results/soft_bootstrap/preservation-audit.json) verified unchanged main sources, external dependencies, default V3 binary and all 223 frozen candidate artifacts/source/dependency hashes. Candidate SHA256 is `d42fc56ca415482581fa13b42bfab48fbd9bc94274cd1caaa5f884e149c0cb02`; retained V3 is `50b4227bca68a58283368d7189dd65dbbd54854f340d05d2ca355824c4ad1738`.

Under the frozen stopping rule, Stage A failed and Stages B/C were not launched. The prototype remains archived in its separate worktree; it was not promoted. These correlated seeds on known files provide an exploratory engineering comparison, not independent biological confirmation.
