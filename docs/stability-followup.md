# UNO stability follow-up: RT calibration and peptide-identity folds

**Completed 2026-09-10: both candidates were rejected; accepted V3 remains the default.** Each isolated candidate completed 30 fresh searches against 30 verified V3 references: six low-input files at seeds 0, 1, 2, 5, and 8. All 60 experimental searches completed successfully. Scientific regressions, rather than execution errors, caused rejection. Neither candidate advanced to primary/raw controls or reserved-seed confirmation.

| Metric across the same 30 cases | Accepted V3 | RT calibration candidate | Peptide-identity folds candidate |
| --- | ---: | ---: | ---: |
| Empty actual reported q≤1% target lists | 5 | 4 | 6 |
| Rescued / newly empty runs versus V3 | — | 2 / 1 | 0 / 1 |
| Median within-file relative entrapment-1% depth range | 116.44% | 111.45% | 120.25% |
| Median within-file relative reported-q1 target-count range | 125.07% | 138.17% | 147.09% |
| Pooled reported-list entrapment FDP | 0.610% | 0.609% | 0.672% |
| Runs with >10% loss at either entrapment threshold | — | 16 | 12 |
| Decision | Keep | Reject | Reject |

Relative range means (maximum−minimum)/median across five seeds, followed by the median across the same six positive-median files. It can exceed 100%. The reported-count range is descriptive; the frozen improvement criterion uses entrapment-1% depth range or empty reported lists. These distinct measures moved in opposite directions for RT. Both candidates passed the pooled FDP guard, but both failed the no-new-empty and per-file median-depth guards. Pooled FDP is descriptive across correlated reruns; the largest individual FDP increased from 2.55% to 5.80% for RT and 4.48% for folds.

RT calibration has a useful recovery signal, but this implementation is unsuitable as a global replacement. It rescued 1 ng rep4 seed1 (0→220 reported targets) and 200 pg rep3 seed5 (0→118), with 338 human and no yeast entries across those two lists. It also made 200 pg rep1 seed2 newly empty (144→0). Median entrapment-1% depth decreased by 5.9%, 23.5%, 29.4%, and 77.5% in four files. The rescued 1 ng case selected the original RT seed's stream, so extra seed attempts alone do not explain the rescue. The new failure passed average validation support despite only 3/49 early-RT validation anchors matching within 0.5 minutes. Validation coverage and average support were insufficient safeguards. Weak early-bin support also occurred in the 200 pg rescue, so it is not a uniquely causal discriminator. The contributions of grouping, training split, positive-slope proposals, and stream selection were not separately ablated.

Peptide-identity folds removed database-index dependence in fold assignment, as verified by unit tests, but did not improve the complete screen. Initial calibration and first-pass training-count logs matched V3 in all 30 pairs; rounded logs do not prove identity of every internal value. One V3 recovery was lost: 1 ng rep4 seed8 changed from 211 accepted targets to zero. Median entrapment-1% depth fell 21.1% in 1 ng rep2 and 75.7% in rep4. Core seeds 0–2 alone looked better for entrapment-depth variability, but stress seeds 5 and 8 reversed that finding. Stable fold identities do not guarantee stable classifier results; training order, negative subsampling, rank ties, and upstream selection remain possible couplings.

The RT candidate grouped identical sequences and shared features, reserved anchor groups for validation, and selected among three fixed positive-slope RANSAC streams using equal-weight support across four RT bins against a training-only constant control. It fell back to V3 when eligibility or selection failed. All 36 calibration passes in this screen selected an experimental fit; none used that fallback. The fold candidate changed only fold assignment to a fixed hash of stripped peptide sequence and effective CV seed. Modified forms and identical target/decoy sequences shared folds. The candidates were developed independently in isolated source copies, frozen before their searches, and never combined or retuned after results.

The prespecified adoption rules protected individual files and successful runs, rather than accepting an overall improvement that concealed new failures. The data suggest investigating RT recovery further, but do not establish that a recovery-only variant would be accurate or stable; no such third variant was tested in this round.

Validation: RT passed 30 core tests, folds passed 28, and the acceptance-report implementation passed 11 tests covering coverage, provenance, null metrics, per-file guards, overlap semantics, and rejected screening. Both release builds and CLI checks passed. Baseline reuse verified binary, scientific inputs, environment, commands, and output/log hashes; copied timings are unavailable. The default binary and all 35 checked Rust/Cargo source files still match accepted V3. No manuscript results were replaced. The reserved draws remain unused in this experiment.

- [RT full results and gates](../benchmark/results/stability_followup/rt-stage-a-report.md) · [mechanism](../benchmark/results/stability_followup/rt_lowinput/rt-screen-notes.md) · [machine-readable details](../benchmark/results/stability_followup/rt-stage-a-report.json)
- [Fold full results and gates](../benchmark/results/stability_followup/folds-stage-a-report.md) · [mechanism](../benchmark/results/stability_followup/folds_lowinput/mechanism.md)
- [Frozen plan](../benchmark/results/stability_study/followup_experiments/screening-plan-v1.md) · [pre-launch amendment](../benchmark/results/stability_study/followup_experiments/screening-plan-amendment-2026-09-10.md) · [RT specification](../benchmark/results/stability_study/followup_experiments/rt-arm-specification-v1.md)
- [Default-preservation audit](../benchmark/results/stability_followup/default-preservation.json) · [prior accepted V3 study](stability-study.md)

Plans, source snapshots, candidate patches/binaries, commands, test/build logs, run outputs, and diagnostic tables are archived under `benchmark/results/stability_followup/`. The RT validation is conditional on mass-selected anchors; correlated reruns on these acquisition files are not independent biological validation. Rejecting these candidates does not establish that all RT or identity-based-fold approaches will fail.
