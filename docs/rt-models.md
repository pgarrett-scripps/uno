# Custom retention-time CNNs

Uno embeds `crates/uno-core/models/rt_cnn_chrono_global_v3_ch32.urtc` (symlinked at `models/`) in the search binary. A
normal search needs no model path, Python, network access, or runtime download.
Use `--rt-model custom.urtc` to select a custom CNN. Search still calibrates its
predictions to each run; it does not train or overwrite reusable model files.

Training is an explicit, optional PyTorch workflow. The supported architecture
is the sequence CNN; MLP and GBDT RT training are deprecated research tools.
Legacy model loading remains available during migration, with a warning. [IM training](im-models.md) uses the same CNN/training engine with a charge
input, separate weights, and explicit mobility units.

## Train for a different column or separation

1. Install the optional package in a Python environment. The search binary has
   no Python dependency. A GPU is optional; the trainer automatically selects
   CUDA when available, otherwise CPU.

   ```sh
   python3 -m venv .venv-rt
   .venv-rt/bin/python -m pip install ./crates/uno-rt-train/python
   ```

2. Prepare a headered TSV of independently assigned peptide sequences and
   observed RTs. Unidentified raw features are not training labels.

   ```text
   sequence	rt
   PEPTIDEK	12.4
   ACDEFGHIK	18.2
   ```

   Use actual tab characters. Required columns are `sequence` and `rt`;
   `run_id` is optional but must be populated when present. Extra columns are
   ignored. Specify `--units minutes`, `seconds`, `irt`, or `hi` explicitly.
   Duplicate observations of the same peptide in the same run are reduced to
   their median RT. Sequences must contain 1–50 uppercase canonical amino acids;
   cysteine is implicitly carbamidomethylated. Variable modifications are not
   supported by this encoding and are rejected rather than silently discarded.

   Train one model per method. A method name records the column/separation; it
   does not make heterogeneous gradients interchangeable. Multiple raw-time
   runs in one input are rejected: normalize them to a common iRT/HI scale
   first. Without `run_id`, the file is treated as one run; the user is
   responsible for not silently pooling incompatible gradients.

3. Adapt the bundled CNN. The default is full-network fine-tuning with a lower
   learning rate than training from scratch.

   ```sh
   uno rt train --data peptides.tsv --units minutes --method "C18 method A" \
     --out custom.urtc --python .venv-rt/bin/python
   ```

   `--init existing.urtc` adapts another CNN; its dimensions, feature count and
   normalization are read from the artifact. `--from-scratch` initializes a
   fresh CNN. These flags are mutually exclusive. No `--arch` choice is needed.
   Output files are never overwritten; choose a new path for each experiment.

   By default, 20% of distinct peptide sequences are held out deterministically
   (`--seed 0`). Every replicate of a sequence stays in the same split. Supply
   `--validation validation.tsv` to choose a separate validation file, or
   `--holdout-run runB` for normalized multi-run data. Run holdout also removes
   those peptides from the other training runs. Explicit train/validation
   overlap is rejected. At least two distinct peptides are required per split;
   this minimum is a technical requirement, not evidence that the data are
   sufficient for a useful model.

   Warm-start predictions are affine-calibrated using training observations
   only, then re-expressed on the new target normalization. Validation compares
   the adapted CNN against that calibrated parent. Early stopping keeps the
   best checkpoint, including epoch zero if training does not improve it.

4. Inspect and evaluate independently, using the same Rust predictor as search.

   ```sh
   uno rt inspect custom.urtc
   uno rt evaluate --model custom.urtc --data heldout.tsv --units minutes
   ```

   The evaluation file must be unused in training **and** model selection.
   The artifact records hashes of peptides used in both stages and known parent
   lineage. Evaluation rejects overlaps. Legacy/pretrained models do not record
   their complete training corpora, so independence cannot be certified for
   those unknown ancestors. The trainer likewise rejects validation peptides
   present in the starting model's recorded lineage; use fresh validation data
   when adapting an existing custom model.

   For another gradient, different label scales, or a legacy model with unknown
   units, provide separate calibration observations in the evaluation units:

   ```sh
   uno rt evaluate --model custom.urtc --data heldout.tsv --units minutes \
     --calibration-data calibration.tsv \
     --baseline models/rt_cnn_chrono_global_v3_ch32.urtc \
     --out-predictions predictions.tsv
   ```

   Calibration and evaluation peptides must be disjoint. The baseline and custom
   model each fit their affine mapping on calibration observations, never on
   evaluation labels. Minutes/seconds conversions need no fitted calibration.
   Reports include MAE, RMSE, R², and the 95th percentile absolute error. A
   constant-label set has undefined R², reported as null. Prediction output is
   headered `sequence`, `observed_rt`, `predicted_rt` and never overwrites a file.

5. Use the custom model in a search:

   ```sh
   uno --features run.features.tsv --fasta proteins.fasta --rt-model custom.urtc
   ```

   Evaluate identification depth and FDR on independent searches before adopting
   a model. A lower validation RT error alone does not establish better search
   performance. The existing dataset-specific
   `crates/uno-rt-train/scripts/validate_rt_in_search.sh` is a research harness,
   not a portable quality certification. Training never replaces Uno's shipped
   default model.

## API and advanced configuration

The CLI calls the installable `uno_rt` package; Python users call the same entry
point directly:

```python
from uno_rt import train

report = train(
    "peptides.tsv", "custom.urtc",
    units="minutes", method="C18 method A",
    epochs=40, seed=0,
)
```

The Rust API retains `RtModel::load(path)` and `predict_many(&sequences)`;
`RtModel::bundled()` constructs the portable default and `cnn_info()` exposes
its architecture/metadata when available. Per-run calibration remains the
search scorer's responsibility; direct predictions are on the model's scale.

`--config advanced.json` (or a Python `config` dict) accepts only:

```json
{
  "lr": 0.0001,
  "batch_size": 512,
  "patience": 6,
  "val_fraction": 0.2,
  "device": "cpu",
  "threads": 4
}
```

For `--from-scratch` only, an `architecture` object can override `emb`, `ch`,
`k`, `blocks`, `head`, or `n_global`. Ordinary users need none of these knobs.
The default is emb32/ch32, kernel 5, dilations [1,2,4], head64 and 11 global
features. Supported global feature prefixes have 0, 5 or 11 entries. Seeded
CPU runs are repeatable within the same runtime; different devices/versions
are not promised bit-identical results. Training stores only the current
minibatch on the accelerator.

## Model format and external predictors

URTC v3 contains architecture and output normalization in its binary header,
then JSON metadata and float32 weights. Metadata records feature schema,
sequence encoding, RT units, method name, timestamps, software versions,
training settings, dataset and parent hashes, peptide lineage and the complete
validation report. No Python pickle or executable code is loaded by search.
The v3 layout extends the v2 header after `n_global` with a little-endian u32
JSON byte count followed by UTF-8 JSON, before the unchanged weight layout.
Rust and Python still read v1/v2 CNNs, including the historical 5-feature model.

The Python wheel includes the same pretrained weights as the Rust binary.
`bundled_cnn_is_portable_and_matches_training_package` verifies they match;
update both copies together when deliberately changing the shipped default.
The CLI passes its embedded weights to the trainer explicitly, so its default
does not depend on which Python package version happens to be installed.

Keep using `--rt-library predictions.tsv` for an external predictor. It is
mutually exclusive with `--rt-model` and replaces the bundled predictions.
Use `--dump-candidate-peptides` to obtain all sequences, including decoys.
Direct loading of arbitrary PyTorch/ONNX architectures is not supported.

## Developer checks

```sh
cargo test -p uno -- rt_command_tests
cargo test -p uno-core -- cnn_workflow_tests
cargo build -p uno
UNO_TEST_BINARY=target/debug/uno .venv-rt/bin/python -m unittest discover \
  -s crates/uno-rt-train/python/tests -v
```

The integration test trains a small CNN, checks Rust/PyTorch parity for v1/v2/v3,
checks artifact inspection and evaluation, and rejects leaked evaluation data.
