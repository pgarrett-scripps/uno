# RT-first IM CNN experiment — 14 September 2026

**Soft CNN IM improved six of nine paired searches, worsened two, and tied one.** Median depth increased in all three acquisitions, but the result did not meet the predefined requirement of at least seven improved pairs. This is a promising, variable benefit; it does not justify claiming either a consistent gain or that IM is unhelpful in general. Following this experiment, the bundled model was enabled by default at the user’s request; `--no-im` restores the off arm. This software choice does not change the experimental consistency result.

Median present-species accessions across seeds 0–2, at **1% empirical entrapment FDP**:

| Acquisition | IM off | Tree IM, soft feature | CNN IM, soft feature | CNN median change |
|---|---:|---:|---:|---:|
| HYE A2 | 263 | 210 | 280 | +6.5% |
| HYE B2 | 69 | 194 | 76 | +10.1% |
| HYE C2 | 381 | 378 | 384 | +0.8% |

The CNN paired deltas, in run order A2/B2/C2 and seed order 0/1/2, were `+16, +176, +14; +7, −3, 0; −109, +3, +67`. The large C2 seed-0 loss (193→84) is retained. The tree soft-feature control improved four pairs, worsened four, and tied one. No search outcome was used to change the model, scoring weight, or acceptance threshold.

These are retrospective accession counts, not UNO’s internally accepted q-value counts. The report separately preserves internal q≤1% target/entrapment counts and all 1%, 3%, and 5% empirical thresholds. It expands same-species protein groups, excludes mixed groups, and includes whole tied-score blocks. All HYE species are present; absent Arabidopsis provides entrapment.

## What changed

The new IM model uses the same trainable sequence CNN, optimizer, normalization, validation, and fine-tuning machinery as RT. Charge/4 is appended to the 11 global sequence descriptors; IM uses separate weights and a UIMC schema. `uno im train`, `inspect`, `predict`, and `evaluate` provide a reusable workflow, including sequence-disjoint validation and parent-lineage checks. See [usage and training](im-models.md).

IM is now an additional classifier feature after RT and mass have fixed the ML shortlist. A label-free per-charge median ratio calibrates observed/predicted mobility symmetrically across target and decoy candidates. The classifier receives absolute relative IM error and availability; no candidate is discarded for IM error. IM does not alter RT fitting or the initial shortlist. Missing IM skips prediction. The archived experiment used explicit opt-in; the current search enables the bundled model by default. Joint-member matrix scoring does not yet use it.

## Prediction validation

The model was trained from scratch on 123,565 public observations from 100,000 sequences, with 12,284 observations from 10,000 other sequences for validation. A fixed 20-epoch run selected epoch 19. All charge states of each peptide stay together. The independent public test contains 45,929 observations from 37,101 unseen sequences: MAE 0.01987 in proportional inverse-mobility units, RMSE 0.02865, and R² 0.96381.

The compatible tree trained on a larger public split is slightly better on that exact public test (MAE 0.01911; R² 0.96865). The CNN therefore provides a consistent sequence-training architecture, not an established overall accuracy improvement. Reversing internal residue order in 1,000 held-out peptides while preserving composition and termini increased CNN MAE from 0.01897 to 0.02429, confirming use of sequence order.

On the separate A1 calibration diagnostic, median relative errors after same-set per-charge scaling were CNN/tree: charge 2, 1.127%/1.160%; charge 3, 2.427%/2.798%; charge 4, 3.334%/2.963%. These same-set calibration errors are optimistic diagnostics, not independent accuracy estimates. A1 labels did not train or select the model and A1 was excluded from searches. Public training peptides can occur in independent search acquisitions; peptide-disjoint generalization is assessed by the separate public test. Training covers charges 2–4 and lengths 7–30 with fixed carbamidomethyl C.

## Cost and verification

The first IM stage requested 318,925–416,597 unique sequence/charge predictions per search. Its candidates were 21.2–25.5% of the mass-matched reduced PFM pool. Prediction is deduplicated within a stage; different RT-calibration candidates can still repeat predictions. Observed median runtime was 73.6 s with IM off and 149.6 s with CNN IM; median paired overhead was 70.0 s. These end-to-end times include concurrency and unequal cache warmth, so they are not a controlled inference-speed measurement.

All 77 Rust tests and eight Python workflow tests passed, including actual CLI training, model adaptation, invalid data and artifact rejection, leakage guards, missing IM, and Rust/PyTorch parity. Maximum Rust/PyTorch discrepancy on 1,000 real public test observations was 2.51e−7. An independent Fraction-based calculation reproduced all 81 threshold counts, verified disjoint training/validation/test sequences, and checked that the initial RT shortlist counts match across paired arms.

## Reproduction and provenance

The fixed study driver is `benchmark/im_cnn_study.py`; `benchmark/im_cnn_prepare.py`, `im_cnn_diagnostics.py`, and `im_cnn_check.py` cover data preparation and verification. The [machine-readable report](../benchmark/im_cnn_20260914.json) includes all 27 searches, model metrics, hashes, and the frozen plan. Working artifacts, full logs, source archive, and the exact search binary are retained in `benchmark/results/im_cnn_20260914`; the reusable checkpoint is `models/im_cnn_public_20260914.uimc`.

The prior IM experiment was not modified. A subsequent SI section documents this study and the published-model comparison. This experiment compares arms using current-source search behavior and preserves the earlier RT model, KOTH features, quality thresholds, FASTA, and seeds within each triple. The exact benchmark binary is archived separately so later builds cannot erase its provenance. Post-freeze changes add tests, strict Python UIMC-version checking, shared RT/IM evaluation metrics, and archive-aware reporting; search scoring and model predictions are unchanged.

## Published IM2Deep comparison

A separate matched calibration experiment used the official `im2deep==2.0.2`
`IM2DeepUni.ckpt` single-output TIMS checkpoint, with explicit fixed
carbamidomethyl cysteine. Both models received the same 20% sequence calibration
partition and per-charge median observed/predicted mapping; the remaining 80%
provided errors. Every uno training/validation sequence was excluded, with all
charges kept together. Predictions were never fitted to test labels.

| Cohort | Test observations | uno median error | IM2Deep median error | uno P95 | IM2Deep P95 |
|---|---:|---:|---:|---:|---:|
| Public reserved test | 36,805 | 1.29% | 1.15% | 5.47% | 4.23% |
| A1 MS/MS assignments | 18,143 | 1.31% | 1.17% | 4.98% | 4.36% |

Errors are absolute relative to each observed mobility value. IM2Deep was more
accurate on both cohorts. Its exact training membership cannot be excluded,
including potential public Meier overlap, so this is not a comparison certified
as unseen by both models. The benefit of the uno model is its shared trainable
RT/IM workflow; these results do not establish superior prediction accuracy.
IM2Deep was not an arm in the search-depth experiment above.

Source: [matched comparison report](../benchmark/im_published_20260914.json),
[reproduction script](../benchmark/im_published_comparison.py), and per-observation
predictions under `benchmark/results/im_published_20260914/`.
Published models: [TIMS2Rescore](https://doi.org/10.1021/acs.jproteome.4c00609),
[IM2Deep multiconformer extension](https://doi.org/10.1021/acs.analchem.5c01142).
The manuscript SI contains the training protocol and generated comparison tables.
