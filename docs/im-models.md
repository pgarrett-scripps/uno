# Sequence-aware IM models

`uno im train` uses the same sequence CNN and training engine as `uno rt train`:
AA embeddings, masked residual dilated convolutions, mean/max pooling, and an
MLP regression head. IM adds charge/4 after the 11 RT global features. All
network weights are trainable, so the model can learn sequence order rather
than only amino-acid composition. The IM model has separate weights and an
explicit `UIMC` format/schema; an RT artifact cannot silently be used as IM.

Search uses the bundled model by default:

```sh
uno --features run.features.tsv --fasta proteins.fasta \
  --seed 0
```

Use `--no-im` to disable IM or `--im-model adapted.uimc` to override the bundled
weights. These flags are mutually exclusive.

RT remains responsible for calibration, bootstrap, and candidate separation.
Before calibration, search precomputes raw IM predictions for every unique
target and decoy sequence in the digested library at each charge present in
valid IM features. These predictions are cached on disk and in the loaded model
for reuse across runs, seeds, and alternative calibration attempts. Only the
existing mass/RT shortlist consumes them during ML rescoring.
A per-charge median observed/predicted ratio calibrates the
scale using targets and decoys equally, without identity labels. At least 20
valid candidates are required for that charge. The absolute relative residual
and an availability flag enter the existing classifier; there is **no IM
cutoff**, no change to the RT shortlist, and no IM-only recovery of failed RT
bootstraps. Predictions can still affect the final classifier ranking/FDR.

Missing IM requires no predictions and has no effect on search. IM is currently
available for single-run/representative-feature rescoring, not joint aligned
member evidence. Disabling ML rescoring also disables this soft IM path. The
legacy `protein.im_window` field is not used by this scorer. The bundled UIMC
model is loaded by default. Legacy tree loading checks the exact descriptor count;
the historical 70-feature `im_im2deep.gbdt` is incompatible with the current
100-feature tree encoder and is rejected.

## Prediction cache

IM caching is enabled by default and reuses **Parquet + Zstd** files across runs.
The default is the OS temporary directory: on Unix, a private `uno-im-<uid>`
directory under `std::env::temp_dir()` (normally `/tmp`, or `TMPDIR` when set).
On Windows, it uses `uno-im` under the user's OS temp directory. Cache files
survive process exit but can be removed by OS temp cleanup; missing entries are
recomputed as needed. `UNO_IM_CACHE_DIR` overrides this root when persistent
storage is wanted; `XDG_CACHE_HOME` no longer selects the IM cache location.
`UNO_IM_CACHE=0` disables disk and memory caching and library precomputation,
returning to uncached prediction of the shortlist. For example:

```sh
UNO_IM_CACHE_DIR=benchmark/results/current/im-cache uno --features run.features.tsv --fasta proteins.fasta
```

Under the selected root, `parquet-v2/<model-fingerprint>/` contains immutable
Parquet chunks with `sequence`, `charge`, and `prediction` columns. Each new
batch writes a same-directory temporary file, closes and syncs it, then
atomically publishes a new `.parquet` file. This avoids rewriting the full
cache for each batch and lets concurrent processes publish independently.
Readers ignore incomplete `.tmp` files and skip corrupt chunks. Normal write
errors remove temporary files; an abruptly killed process can leave a `.tmp`
file for OS cleanup. Legacy TSV caches are ignored and left untouched.

Cache identity includes a SHA-256 of the exact model bytes and prediction-schema
version; each row is keyed by sequence and charge. Float64 predictions preserve
their exact bits, including signed zero. Calibration and residuals are recomputed
from each shortlist. A new model uses a separate cache; new shuffled-decoy
sequences and new charges fill only missing entries. Unavailable disk storage
falls back to memory reuse. The first complete library
prediction costs more than a single shortlist prediction, while later searches
reuse the covered pairs. Query batches are bounded, but the loaded cache occupies
memory proportional to its stored sequence/charge pairs.

Logs report library pair counts, cache misses, and elapsed time, followed by
per-shortlist prediction lookup timing. RT uses its independent `UNO_RT_CACHE`
and `UNO_RT_CACHE_DIR` controls. Neither cache needs to be discarded just because
KOTH features change when the model and prediction semantics are unchanged.

### Cache verification

Cache tests verify the Parquet schema and Zstd codec, exact floating-point bits,
charge and model separation, memory and fresh-model disk hits, disabled caching,
partial/corrupt files, concurrent publication, unavailable or evicted storage,
private temporary directories, and independent recalibration of different
shortlists.

The measurements below describe the earlier **TSV cache implementation**. They
are retained as historical verification, not as Parquet performance measurements.

A serialized C2/seed1 check with four threads compared uncached, cold-cache,
and warm-cache searches using the same binary, inputs, and warm RT cache.
Protein, precursor LFQ, and protein abundance files matched byte for byte in
all three cases; the protein output also matched the archived experiment.
The full library contained 2,004,278 unique sequences at five observed charges
(10,021,390 pairs). Initial IM cache population took 385.1 seconds; the warm
library check/load took 8.7 seconds with zero misses. Shortlist lookups took
0.24 seconds versus 12.9 seconds of uncached prediction. Total search wall
times were 43.2 seconds uncached, 413.0 seconds cold, and 32.4 seconds warm.
These are single-run observations, not a repeated performance benchmark.

The verification driver, binary, input/output hashes, timings, logs, and
`report.json` are archived in `benchmark/results/im_cache_verification/`.

## Train or adapt

Install the same optional Python package used for RT:

```sh
python3 -m venv .venv-rt
.venv-rt/bin/python -m pip install ./crates/uno-rt-train/python
```

Use independently assigned sequences and observed mobility, with TSV headers
`sequence`, `charge`, `im`, and optionally `run_id`. Canonical sequences of
1–50 residues are supported; cysteine is implicitly carbamidomethylated and
variable modifications are rejected. Charge must be an integer from 1 to 8;
IM must be finite and positive. Duplicate sequence/charge/run observations
are reduced to their median. Use observations from a compatible mobility
method/scale; a descriptive method name does not normalize incompatible data.

```sh
uno im train --data labels.tsv --units inverse-mobility --method "TIMS method A" \
  --init models/im_cnn_public_20260914.uimc --out adapted.uimc \
  --python .venv-rt/bin/python
```

Use `--from-scratch` instead of `--init` for fresh weights. There is no implicit
bundled IM initialization. Default training settings match RT: 40 epochs,
Adam, smooth-L1 loss, deterministic peptide splits, early stopping, and
full-network fine-tuning at a lower learning rate than scratch training.
Fine-tuning affine-calibrates the parent on training observations only and
retains the calibrated epoch-zero parent if adaptation does not improve
validation MAE. This global affine initialization is separate from search's
label-free per-charge calibration.

`--validation`, `--holdout-run`, `--seed`, `--epochs`, and `--config` work like
RT training. All charges and replicates of a sequence stay in the same split;
validation overlap with the parent's recorded training/selection lineage is
rejected. Architecture overrides require scratch training; IM requires
`n_global=12`. The Python API is `uno_rt.im.train`; its unit strings use
underscores (`inverse_mobility`, `proportional_inverse_mobility`).

## Inspect, predict, and evaluate

These commands use the Rust search predictor and need no Python:

```sh
uno im inspect adapted.uimc
uno im predict --model adapted.uimc --data sequence_charge.tsv --out predictions.tsv
uno im evaluate --model adapted.uimc --data independent.tsv --units inverse-mobility
```

Prediction input needs `sequence` and `charge` headers. Output contains each
unique sequence/charge in sorted order. Evaluation input additionally needs
`im`. Evaluation rejects sequences used in training or checkpoint selection.
Use `--calibration-data calibration.tsv` to fit per-charge multiplicative
scales on separate observations; calibration/evaluation sequences must be
disjoint, including across charges. At least two calibration observations per
charge are required. Reports include MAE, RMSE, R², and 95th-percentile absolute
error; constant labels produce null R². `--out-predictions` retains observed
and predicted values. Existing model and prediction files are never overwritten.

The experimental public model predicts **proportional inverse mobility**, not
absolute 1/K0 or CCS. Its labels are the prior experiment's measured public
CCS labels converted to `CCS * sqrt(reduced_mass) / charge / 1000`; a scale
calibration is required for measured 1/K0. It was trained on charges 2–4 and
lengths 7–30; predictions outside that training domain are extrapolations.

## Reproduce the experiment

`benchmark/im_cnn_prepare.py` selects a deterministic subset of the existing
public training data: 100,000 training sequences, 10,000 validation sequences,
and 37,101 independent test sequences. All their charge observations are kept.
The model used 20 epochs with seed 0 and selected epoch 19 by validation MAE.
Source identities, hashes, the training log, model, source archive and every
search receipt are in `benchmark/results/im_cnn_20260914`.

`benchmark/im_cnn_study.py` implements the frozen 27-search comparison: IM off,
the previous compatible tree as a soft feature, and the new CNN as a soft
feature, across A2/B2/C2 and seeds 0–2. The existing RT model, KOTH features,
quality filters, and HYE-plus-absent-Arabidopsis FASTA are identical within
pairs. A1 is used only for prediction diagnostics. Public training peptides
can occur in the independent search acquisitions; the peptide-disjoint public
test is the separate check of unseen-sequence prediction performance.

The primary outcome is retrospective accession depth at 1% empirical
entrapment FDP, with whole tied-score blocks and same-species group expansion.
All nine paired CNN outcomes are retained. Improved prediction error alone
does not establish improved protein identification.

The [completed 27-search comparison](https://github.com/pgarrett-scripps/uno-paper/blob/main/docs/im-cnn-results-20260914.md) found six CNN gains, two losses, and one tie, with higher medians in all three files; the predefined consistency rule was not met.
