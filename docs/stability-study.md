# UNO seed stability study

Status: **V3 passes the frozen stability and calibration criteria.** All 174 results (87 matched pairs) were verified. V2 remains rejected; historical manuscript results remain archived separately.

The subsequent [RT-calibration and peptide-identity-fold experiments](stability-followup.md) each completed 30 searches and were rejected for regressions. V3 remains the default.

The default `target/release/uno` was rebuilt and is byte-identical to the frozen validated V3 binary. [Release verification](../benchmark/results/stability_study/validation/v3-final/release-verification.json) records its SHA-256, tests, source match, and figure checks.

## Final finding

| Panel | Paired runs | Empty reported 1% lists, baseline → V3 | Pooled reported-list FDP, baseline → V3 |
|---|---:|---:|---:|
| Low-input biosaur2 | 54 | 15 → 5 | 0.464% → 0.526% |
| Primary biosaur2 | 18 | 0 → 0 | 0.945% → 0.945% |
| Raw 200 pg | 15 | 8 → 3 | 0.556% → 0.572% |

All **64 previously successful outputs are byte-identical**. No paired depth decreases occurred at either entrapment threshold. Low-input median relative 1%-depth range fell from 123.821% to 99.147% (about 20% lower), comparing the same four defined files; 200 pg replicates 1 and 3 have zero baseline medians and are excluded from that range statistic. Median file depth remained 637 / 1134.5 at the 1% / 5% entrapment thresholds.

On unseen draws, low-input failures fell from 6/12 to 3/12 and raw failures from 4/6 to 2/6. These are new random draws on known files, not independent biological samples. Among rescued runs only, pooled FDP was 0.937% for low-input biosaur2 and 0.588% for raw. Five of ten rescued low-input lists and one of five rescued raw lists individually exceeded 1% empirical FDP. Thus the result supports fewer failures with acceptable pooled calibration under the frozen engineering guard; it does not establish per-run nominal FDR or eliminate seed variation.

[Combined decision and full results](../benchmark/results/stability_study/overall-v3.md) · [All reported-identification draws](../benchmark/results/stability_study/figures-v3/lowinput_internal_q1_human.png) · [1% entrapment-depth comparison](../benchmark/results/stability_study/figures-v3/lowinput_depth_1pct.png) · [5% entrapment-depth comparison](../benchmark/results/stability_study/figures-v3/lowinput_depth_5pct.png)

Each figure also has a PDF beside its PNG. All nine seeds per low-input file are shown, including failures; filled circles are development seeds 0–6 and hollow diamonds are unseen seeds 7–8.

## Implemented changes

- **Recovery only after complete original failure.** First run the original empty-cut and final-candidate selection. If it reports any target at q≤1%, preserve that ranking. Only a zero-discovery result triggers a second pass with evidence retention; accept that retry only if it produces actual reported discoveries, otherwise keep the full original result. Clear and drain both evidence stores for each pass, including early returns. This reuses the same random seeds and RT predictions, but can substantially increase calibration/scoring time on failed searches.
- **Retain training evidence after a zero-discovery cut sweep.** Previously, if all deviation cuts had zero proteins at the conservative 1% target-decoy threshold, UNO returned an empty score list. This also removed the existing 25%-FDR training positives and prevented ML recovery. The selector now retains the cut with the largest existing training set when all 1% counts are zero; nonzero count rankings and first-cut ties are unchanged. The FDR thresholds and correction are unchanged.
- **Protect existing discoveries.** A training-only ranking must not displace a candidate with actual reported 1% identifications. Final selection first prefers a positive discovery count, then uses the existing score criterion and first-candidate tie rule. If all candidates remain unconfident, the ranked evidence is retained without claiming confident discoveries.
- **Keep candidate evidence together.** Protein lists and peptide-to-feature maps now belong to the same finalized candidate. Previously, a later losing ML candidate could overwrite the map used for grouping and quantification. A selected candidate that never ran ML has an empty map, preserving the existing no-ML limitation rather than borrowing another candidate’s evidence.
- **Independent random-stream controls.** Optional `scoring.rt_seed`, `cv_seed`, and `classifier_seed` override the shared scoring seed; absent overrides preserve old behavior. `--seed` still sets the parent scoring seed and decoy seed.

## What the ablations showed

Astral 200 pg replicate 1, all other seed settings held at zero. Depths are human proteins at the entrapment 1% / 5% thresholds. These are separate from UNO’s reported 1% lists.

| Seed intervention | Seed 1 | Seed 2 |
|---|---:|---:|
| Decoy generation | 0 / 0 | 0 / 0 |
| RT/RANSAC sampling | 0 / 0 | 97 / 513 |
| CV fold assignment | 46 / 162 | 161 / 457 |
| Classifier sampling | 389 / 656 | 432 / 695 |

The seed-zero control was 389 / 769 and byte-identical to its archive. Decoy generation and RT sampling can each trigger the early failure; CV alone also changes quality substantially. The failing RT draw had a *smaller* inlier residual, so residual size alone does not identify a useful calibration. Decoy changes also alter downstream peptide indices and CV assignments: these interventions are not independent estimates of variance components. Full commands, hashes and actual reported-list FDP are in [the ablation archive](../benchmark/results/stability_study/ablations/findings.md).

## Validation and safeguards

Phase 1 compared 30 low-input searches and 18 primary searches per binary. All primary pairs were byte-identical. Low-input empty outputs fell from 6 to 0; median relative depth range on the same five defined files fell from 1.081 to 0.629. Actual reported-list pooled FDP changed from 0.545% to 0.622%. However, two previously nonempty reported 1% lists became empty because a training-only list won by soft depth. That functional regression was rejected despite passing the initial aggregate criteria.

Candidate v2 adds the discovery-preservation safeguard. Its frozen [validation plan](../benchmark/results/stability_study/v2-plan.md) specifies six low-input files at seeds 0–6, six primary files at seeds 0–2, and three raw 200 pg files at seeds 0–2. Seeds 0–4 are now development evidence; 5–6 are unseen random draws on the same files. This is not a holdout of biological samples. Verified baseline outputs may be reused, but every v2 candidate search is executed afresh.

Acceptance requires halving low-input zero-depth searches, reducing relative variability, preserving median depth, no newly empty reported 1% lists, and no pooled reported-list FDP increase above 0.5 percentage points in each panel. Every individual depth loss above 10% is retained and flagged. Pooled FDP and ranges are descriptive because repeated seeds/files are correlated; this finite-panel tolerance is not a proof of nominal FDR control.

V2 passed its complete biosaur2 evaluation: zero entrapment-depth searches fell from 9/42 to 0/42, but two results still had no actual reported 1% identifications. All 18 primary pairs were byte-identical. **V2 failed the separate raw guard:** pooled reported-list FDP increased from 0.426% to 1.591%. In raw 200 pg replicate 3, seed 2, newly trained recovery candidates displaced an unchanged successful original candidate, reducing entrapment-1% depth from 732 to 356. V2 therefore remains rejected.

The frozen [v3 plan](../benchmark/results/stability_study/v3-plan.md) adds strict fallback behavior and requires unchanged scoring rows for every baseline-positive search. It includes all six low-input files at seeds 0–8 (new draws 7–8), primary seeds 0–2, and raw seeds 0–4 (new raw draws 3–4). All candidate searches are fresh; only verified matching baseline outputs may be copied. Calibration guards remain unchanged.

The implementation passes 24 core tests plus 3 CLI seed-precedence tests. Reproducibility harnesses were tested for empty outputs, conservative tied-score handling, checked resume, tampered hashes, undefined ranges, group-only changes, and baseline reuse without execution.

## Provenance and remaining limits

The frozen pre-study, rejected v2, and accepted v3 binaries, source snapshots, model/input hashes, dependency fingerprints, commands, environments, logs, output hashes, plans and tests are under `benchmark/results/stability_study/`. No old paper output is overwritten or relabeled as a new engine run. Interrupted searches preserve their old logs; successful results are reused only after validation.

Some older nonzero-seed archives differ from fresh frozen-baseline reruns. Three controls on one discrepancy—4 threads/cache off, 20 threads/cache on, and repeated 4 threads/cache on—were byte-identical to the fresh baseline. Thus those controls do not explain the historical mismatch. Older studies incompletely recorded the environment and did not recheck mutable binary paths at each queued launch; the cause remains unresolved. The new matched comparison does not rely on those archived values. Sequoia’s `nthread(1)` field is not enforced by its current backend, but no causal thread effect was established here, so no unsupported threading change was promoted.

Pinning a favorable seed can reproduce one output but does not cure seed sensitivity. These changes target the failure path and preserve discovery confidence; ordinary differences in rankings and depth across seeds can remain. Before using the updated engine for manuscript numbers, rerun and version-pin that manuscript’s full benchmark inventory.

## Reproduce the final report

The completed v3 biosaur2 panels can be rechecked with:

```sh
python3 benchmark/stability_report.py \
  --studies benchmark/results/stability_study/lowinput_v3 benchmark/results/stability_study/primary_v3 \
  --protocol benchmark/results/stability_study/v3-plan.md \
  --lowinput-seed-count 9 --protect-discoveries \
  --output benchmark/results/stability_study/report-v3
```

The report evaluates biosaur2 panels only: overall acceptance additionally requires the separate raw guard and baseline-positive score-preservation audits. The exact search commands and baseline-copy origins are recorded in each study’s `study.json` and run stamps. Standalone figures are generated by `benchmark/plot_stability.py` in the paper analysis environment; manuscript figures are kept separate.
