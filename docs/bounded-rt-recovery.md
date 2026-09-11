# Bounded RT recovery evaluation

Adopted into the main source and rebuilt default executable on 2026-09-10 after
all prospectively specified accuracy, preservation, speed and repeatability
requirements passed. The fixed protocol is in
[`PLAN.md`](../benchmark/results/bounded_rt_recovery/PLAN.md).

## What changes

UNO first completes V3, including its evidence-retention retry. Any successful
V3 result is returned unchanged. Only a complete failure triggers one additional
broad RT-calibration pass, with the persistent broad RT-error floor removed.
The pass keeps the existing broad map and freezing policy and does not add a
second broad retry or extra seeds. A failed recovery returns the original V3
result. The extra pass still includes its existing internal scoring candidates;
it is not a single model fit.

This change targets empty results. It does not alter successful V3 rankings or
eliminate differences among successful runs with different random seeds.

## Completed development-case evaluation

All 30 known low-input file/seed pairs passed the frozen accuracy requirements.
Five V3-empty results became one candidate-empty result. All 25 successful V3
outputs and the one unrecovered output remained byte-identical.

| Recovered case | Human / yeast reported IDs | Entrapment FDP estimate |
|---|---:|---:|
| 1 ng replicate 4, seed 1 | 367 / 1 | 0.927% |
| 200 pg replicate 1, seed 8 | 146 / 1 | 2.331% |
| 200 pg replicate 3, seed 5 | 369 / 1 | 0.922% |
| 200 pg replicate 3, seed 8 | 222 / 0 | 0.000% |

The recovered lists pooled to 1,104 human and three yeast IDs, giving an
entrapment FDP estimate of 0.925%. The pooled estimate across all runs changed
from 0.610% to 0.639%. One recovered list exceeded 1%; all were below the fixed
5% individual cap. Finite entrapment estimates do not establish nominal per-run
FDR, and these are development cases rather than independent confirmation.

## Controlled timing and exact repeats

All 16 prespecified serialized executions passed. Two rounds per case used
V3/candidate then candidate/V3 order, with unchanged warm RT cache, no extra
warm-up searches and no other agent-controlled scientific work. The complete
command was timed, including parsing and scoring.

| Case | V3 median seconds | Candidate median seconds | Candidate / V3 |
|---|---:|---:|---:|
| Successful low-input case | 56.74 | 57.32 | 1.010 |
| Successful primary case | 53.19 | 53.25 | 1.001 |
| Failed 1 ng case | 91.69 | 122.81 | 1.339 |
| Failed 200 pg case | 74.71 | 96.62 | 1.293 |

Both within-version repeats were byte-identical in every case. Successful
outputs matched across versions; recovered outputs matched the development
panel. All observed within-version timing spreads were below the fixed 10%
limit. Two timing repetitions are engineering checks, not confidence intervals
or a universal performance guarantee.

## Fresh-seed confirmation

The same six acquisitions at previously untested seeds 11–15 passed all frozen
accuracy requirements, with both engines executed fresh for all 30 pairs.
Empty results fell from five to two. All 25 successful V3 outputs and both
unrecovered outputs remained byte-identical.

All three rescues were on 200 pg replicate 1: seed 11 returned 234 human / zero
yeast IDs, seed 12 returned 218 / zero, and seed 14 returned 127 / one. The
recovered lists pooled to 579 human / one yeast, FDP 0.588%. Seed 14's individual
estimate was 2.680%, above 1% but below the frozen 5% cap. Overall FDP changed
from 0.633% to 0.631%. The remaining empty cases were 200 pg replicate 2 seed 11
and replicate 3 seed 15.

The median across files of relative target-count ranges remained 0.888 in both
versions. Recovery improved one acquisition's failures without changing the
other files or the median. This supports fewer empty runs, not a claim that
general variability among successful runs is solved. These are fresh random
draws on familiar acquisitions, not biological holdouts.

## Controls and acceptance

All 24 control pairs passed independently: 12 primary-feature pairs, six raw
input pairs and six higher-input feature pairs. Every output was byte-identical
between V3 and the candidate. The raw panel included one empty result that
neither version recovered. No panel had a new empty result or identification
loss. All five efficacy panels and the controlled timing panel passed their
prospectively fixed requirements; the candidate is eligible for promotion.

Across the 84 efficacy pairs, seven of 11 V3-empty runs recovered, leaving four
empty results. All 73 V3-positive outputs and all four unrecovered outputs were
preserved. The known and fresh low-input evidence must still be interpreted
separately: four of five development failures and three of five fresh-seed
failures recovered. The controls contributed no additional recovery benefit.

## Verified adoption

The reviewed four-file patch was applied to the main checkout, preserving the
other audited sources and dependencies. All 37 required-data core tests and five
CLI tests passed. The default release executable was rebuilt in the main target
directory; help and runtime-library checks passed.

Four main-build replays exactly matched the frozen candidate's complete output:
an ordinary success (1 ng replicate 1, seed 0), the 368-ID recovery (1 ng
replicate 4, seed 1), the 222-ID recovery (200 pg replicate 3, seed 8), and a failed
recovery preserving the empty V3 output (200 pg replicate 2, seed 8). Activation
logs and input/source/binary hashes also passed verification. These are build
equivalence checks, not additional independent efficacy evidence.

- Frozen evaluated binary SHA256:
  `5ba14ed0e2ac4d7ab90501ce2548d6ad3f0519116be2816b40d3fa49b3422cd0`.
- Rebuilt default `target/release/uno` SHA256:
  `c376954dbd7c725c8cd280b9c392ba3854fd9c831eb423601a3eb7ebc5539eda`.

These are separate builds with different build locations and runtime-library
paths. The four changed source files match the frozen candidate exactly, and
scientific output equivalence is verified above. An independent promotion audit
also confirmed all 52 other audited main source files were unchanged. The prior
V3 executable is retained in the promotion evidence. Paper results and figures
were not replaced.

## Evidence

- [Known-case report](../benchmark/results/bounded_rt_recovery/known-report.json)
- [Independent known-case audit](../benchmark/results/bounded_rt_recovery/independent_review/known-efficacy-audit.json)
- [Recovery-path diagnostics](../benchmark/results/bounded_rt_recovery/known-diagnostics.json)
- [Controlled timing report](../benchmark/results/bounded_rt_recovery/timing/report.json)
- [Independent timing audit](../benchmark/results/bounded_rt_recovery/independent_timing_audit/)
- [Fresh-seed report](../benchmark/results/bounded_rt_recovery/fresh-report.json)
- [Complete acceptance report](../benchmark/results/bounded_rt_recovery/complete-report.json)
- [Independent complete audit](../benchmark/results/bounded_rt_recovery/independent_review/complete-efficacy-audit.json)
- [Verified main-build adoption](../benchmark/results/bounded_rt_recovery/promotion-validation/promotion-result.json)
- [Independent promotion audit](../benchmark/results/bounded_rt_recovery/independent_promotion_audit/)
- [Frozen source and binary registry](../benchmark/results/bounded_rt_recovery/candidate-registry.json)
