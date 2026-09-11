# Retention-time predictor benchmark

This benchmark compares uno's shipped sequence CNN with DeepLC 3.1.13, the RT
predictor enabled by `ms1searchpy -deeplc 1`. It is deliberately scoped to
unmodified peptides (fixed carbamidomethyl cysteine is allowed) because the
current uno model has no variable-modification encoding.

The three checked-in sets cover the peptide-disjoint Chronologer holdout,
unmodified ProteomeTools/Prosit holdout peptides, and DeepLC's public
`unmod.csv` example dataset. Every peptide is 7--30 canonical residues and is
absent from the Chronologer train and validation partitions used for the uno
model. Repeated experimental observations are reduced to their median label.

`prepare_sets.py` rebuilds `sets/*.tsv` from the public files documented in
`data/rt/README.md`. It ranks sequences with SHA-256 and retains at most 10,000
per source, so the checked-in sets are deterministic and reviewable. This proves
sequence disjointness from uno's training data. The exact peptide membership of
the packaged DeepLC model's training corpus is not available here, so the same
exclusion cannot be established for DeepLC.

## Run

```bash
cargo build --release -p uno-rt-train --bin rt_eval
.venv/bin/python benchmark/rt_predictor/run.py
cargo build --release -p uno-rt-train --bin rt_throughput
.venv/bin/python benchmark/rt_predictor/prepare_library.py   # synthetic tryptic library
.venv/bin/python benchmark/rt_predictor/throughput.py        # results/throughput.tsv
cd paper/analysis
uv run --quiet scripts/figS2_rt_predictor.py
uv run --quiet scripts/figS3_rt_throughput.py
```

The sequence hash assigns 20% of each set to calibration and 80% to testing.
Both models receive the same calibration: a least-squares affine map from the
model's native output scale to the observed scale, fitted on the 20%. DeepLC is
run with `calibrate=False` and mapped by that same fit. Its own piecewise
`calibrate_preds` (what ms1searchpy uses in-search) is deliberately not used:
with ~200 calibration peptides it extrapolates poorly at the elution-range ends,
which would measure the calibrator rather than the model.
Only the held-out 80% contributes accuracy. Normalized MAE is MAE divided by
the test-label 5th--95th percentile range, which makes datasets with different
RT units comparable. Throughput times prediction after calibration; it excludes
calibration time, which is reported separately in `results/summary.tsv`.

The output includes exact sample counts, MAE, normalized MAE, delta-t 95, R²,
Spearman correlation, calibration time, and prediction throughput.
`results/predictions_<set>.tsv` holds the per-peptide test predictions of both
models, which the paper's Figure S2 uses for residual and error-distribution
panels. Hardware,
package versions, model hash, and git revision are recorded in
`results/environment.json`.

## Throughput

Throughput depends on the peptides, not the experimental set, so it is measured
separately on `sets/synthetic_tryptic.tsv`: an in-silico tryptic digest of the
human entries in `../data/entrap_human_yeast.fasta` (7--30 canonical residues,
unique, SHA-256 ordered, 100,000 peptides). `throughput.py` times both models
on the first N peptides for N = 100 ... 100,000, five repeats after a warm-up
at the largest size, inference only. uno is timed in-process by the
`rt_throughput` binary; DeepLC by `make_preds(calibrate=False)` on a prebuilt
DataFrame (its feature encoding is included, since every DeepLC prediction pays
it). Results go to `results/throughput.tsv` and
`results/throughput_environment.json`.
