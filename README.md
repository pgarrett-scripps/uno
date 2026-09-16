# uno — MS1-only search engine (Rust)

A Rust MS1-only protein identification engine, in the spirit of DirectMS1 /
ms1searchpy. `uno` accepts a single-run biosaur2 or KOTH feature TSV (`--features`),
finds features with embedded KOTH (`--raw`), or searches a KOTH aligned intensity
matrix across several runs (`--matrix`). It assigns peptide/protein identities
and estimates FDR with target/decoy competition.

See [matrix input](docs/matrix-input.md) for multi-run searching and the
experimental confidence-aware matrix schema. Matrix input retains cross-run
intensities but omits some quality descriptors from individual feature tables.

Digestion is reused from [sage-core](../sage); retention-time prediction is a
pure-Rust learned model, with an on-data calibration step per run.

## Experimental quantitative options (LFQ branch)

This branch adds [uno-quant](crates/uno-quant/README.md), a Rust
DirectMS1Quant-style contrast analyzer, and opt-in `--out-quant-evidence` exports
from independent uno searches. The existing KOTH-matrix route remains the second
quantitative option. Native-reference parity is tested; biological accuracy and
error calibration are pending. Development is isolated on `codex/uno-lfq`.

## Workspace layout

| crate | purpose |
| --- | --- |
| `uno-core` | search, scoring, and FDR library |
| `uno` | command-line binary |
| `uno-quant` | experimental quantitative contrast analysis |
| `uno-rt-train/python` | shared optional PyTorch sequence-CNN trainer for `uno rt train` and `uno im train` |
| `uno-rt-train` | legacy RT/IM tree-model research tools |
| `uno-lightgbm-sys` | FFI bindings to a vendored LightGBM |

## Prerequisites

- A recent stable Rust toolchain (`rustup`, edition 2021).
- **A sibling `sage` checkout.** The workspace has a path dependency
  `sage-core = { path = "../sage/crates/sage" }`, so `sage` must be cloned next
  to this repo before it will build.
- **KOTH 0.8.0 candidate source** at
  `../koth_rust/.worktrees/release-0.8.0` (commit
  `dce83d5f2f01281d9464cbdc147d8759d78e0834`) and the sibling `Sequoia`
  checkout. This LFQ branch enforces KOTH version `=0.8.0`; the source remains a
  local path dependency. See [the LFQ migration](docs/matrix-input.md#koth-070-compatibility)
  for new settings and schema-3 support.

```
Repos/
├── sage/    # https://github.com/lazear/sage
├── koth_rust/ # https://github.com/tacular-omics/koth_rust, tag v0.7.0
├── Sequoia/ # sequoia-boost dependency
└── uno/     # this repo
```

## Build

```sh
cargo build --release
```

## Run

```sh
cargo run --release -p uno -- \
  --features path/to/run.features.tsv \
  --fasta path/to/proteins.fasta \
  --config example_config.toml \
  --out-proteins proteins.tsv
```

See `example_config.toml` for all tunables (digest, mass tolerance, `top_k`) and
their documented defaults. Search uses the embedded pretrained CNN by default;
`--rt-model custom.urtc` selects a custom model. Training uses PyTorch; search,
model inspection, and evaluation run entirely in Rust.

See [Custom RT models](docs/rt-models.md) for `uno rt train`, `uno rt inspect`,
`uno rt evaluate`, the data format, and the Python training API.
`uno --help` lists every flag, including the external predictor
RT-library workflow (`--dump-candidate-peptides` → `--rt-library`).

## Experimental sequence-aware IM

[IM models](docs/im-models.md) use the same sequence CNN and training workflow
as RT, with charge added as an input. `uno im train` supports training from
scratch and full-network adaptation; inspection, prediction, and independent
evaluation use the Rust search predictor.

The bundled IM CNN contributes a soft classifier feature by default after the
mass/RT shortlist. RT remains the main separator; there is no IM cutoff. Use
`--no-im` to disable it or `--im-model adapted.uimc` to use custom weights.
Missing IM leaves the search unchanged. See the model guide for supported chemistry, units, and limitations.

## Search reproducibility

`--seed N` controls decoy generation and the default RT, cross-validation, and
classifier random streams. Reuse the same seed for repeat searches;
test multiple seeds to measure whether the result is stable. Keep the inputs,
model, config, software version, and experimental environment settings fixed.

Optional `scoring.rt_seed`, `scoring.cv_seed`, and `scoring.classifier_seed`
settings override their individual streams, including when `--seed` is supplied.
Leave them absent for shared-seed behavior; see `example_config.toml` for examples.
Changing decoys also changes downstream peptide indexing, so an individual seed
intervention is not a pure variance-component estimate.

The [stability study](https://github.com/pgarrett-scripps/uno-paper/blob/main/docs/stability-study.md) records the low-input failure
mechanism, paired evaluation, and calibration checks. Historical manuscript
outputs remain archived separately from new engine experiments.

The current search includes [bounded RT recovery](https://github.com/pgarrett-scripps/uno-paper/blob/main/docs/bounded-rt-recovery.md):
when the usual search and its evidence-retention retry find no target proteins
at 1% FDR, UNO attempts one additional broad RT calibration. Successful results
are returned unchanged. Validation showed fewer empty low-input runs with
0–1% normal-path runtime overhead; variation between successful seeds remains.

## Repository contents

- `crates/` — the Rust workspace.
- `models/` — sequence CNNs and historical RT/IM tree models.
- [uno-paper](https://github.com/pgarrett-scripps/uno-paper) — the separate manuscript repository; `paper/` is a compatibility link to the sibling checkout.
- [uno-lfq-paper](https://github.com/pgarrett-scripps/uno-lfq-paper) — LFQ manuscript and workflow studies.
- Study analysis, figure/table generators, and evidence live with their respective papers.
- `tests/fixtures/` — software regression-test configurations.

Build artifacts (`target/`), the Python virtualenv (`.venv/`), RT training data
(`data/`), and scratch outputs (`experiments/`) are intentionally not tracked.

## License

MIT © Patrick Garrett
