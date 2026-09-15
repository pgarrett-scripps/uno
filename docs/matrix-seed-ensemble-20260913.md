# Matrix seed ensemble diagnostic — 13 September 2026

Combining three saved matrix searches does not close the IonStar gap. Gains vary by cohort, scoring arm, and aggregation rule; the experiment supports testing search stability but does not establish randomness as the cause of the separate-file advantage. All matrix ensembles use only matrix search outputs.

## Matched minimum-q comparison

Each entry is the median of the three one-seed results → the ensemble containing all three seeds. Both routes rank accessions by their minimum saved internal q-value, then independently select the deepest whole-tied-score cutoff at ≤1% empirical entrapment FDP. The matrix columns therefore have different one-seed baselines from the manuscript’s binomial-score ranking.

| Cohort | Separate-file union | Matrix | Matrix + cell covariates |
|---|---:|---:|---:|
| Astral 200pg | 410 → 557 | 706 → 740 | 0 → 0 |
| Astral 1ng | 1,209 → 1,221 | 819 → 722 | 849 → 880 |
| Astral 10ng | 1,458 → 1,279 | 1,306 → 1,357 | 1,199 → 1,254 |
| Astral 100ng | 1,352 → 1,558 | 1,687 → 1,774 | 1,721 → 1,779 |
| Astral all | 1,275 → 1,455 | 1,252 → 1,479 | 1,194 → 1,042 |
| Ionstar 5 | 1,892 → 2,063 | 1,375 → 1,271 | 1,384 → 1,497 |
| Ionstar 10 | 1,846 → 1,935 | 1,418 → 1,275 | 1,349 → 1,337 |
| Ionstar 20 | 1,911 → 1,870 | 656 → 877 | 890 → 1,030 |

The zero median for Astral 200 pg with cell covariates reflects minimum-q ties: two of the three searches have no complete q-value block satisfying the 1% entrapment cutoff. It does not mean those searches identify no proteins when ranked by binomial score.

## Sensitivity to aggregation

Using maximum binomial score for both routes, the IonStar 20-file matrix changes from 782 to 969 accessions; adding cell covariates changes 897 to 844. The separate-file union changes from 1,642 to 1,520 under this same ranking. Thus the IonStar deficit remains with either aggregation rule. For Astral 100 ng, the matrix changes 1,718 to 1,711 and the cell-covariate arm changes 1,723 to 1,918; the separate-file union changes 1,279 to 1,451.

Pooling individually thresholded lists without recalculating error overstates the benefit. For example, the three IonStar 20-file matrix lists contain 1,402 unique present-species accessions at their individual binomial-score cutoffs, but their union has 1.46% empirical FDP. Re-ranking the combined results by maximum score and restoring ≤1% leaves 969.

## Scope and verification

All seven nonempty subsets of seeds 0, 1, and 2 were evaluated in all eight cohorts, separately for both matrix arms and the independent-file route. Minimum q was the primary diagnostic and maximum score a sensitivity analysis; no best-performing seed combination was selected. All 144 saved protein tables match the original hashes, all original baseline counts and accession sets were reproduced, and an independent evaluator agrees at all 1,008 evaluated cutoffs (1%, 3%, 5%). Whole-tie rejection and exact-boundary tests passed.

This is a post hoc diagnostic, not a new UNO ensemble implementation or a calibrated cohort q-value method. Seeds vary RT calibration, CV folds, and classifier training together. Three matrix searches are compared with 3N separate searches when both routes use all seeds; compute is not matched, and three seeds cannot isolate the cause of the remaining gap. No new searches were needed.

[Complete results](../benchmark/results/matrix_seed_ensemble_20260913/report.json) · [Validation](../benchmark/results/matrix_seed_ensemble_20260913/validation.json) · [Reproduction script](../benchmark/matrix_seed_ensemble.py)
