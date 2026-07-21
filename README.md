# uno — MS1-only search engine (Rust)

A Rust MS1-only protein identification engine, in the spirit of DirectMS1 /
ms1searchpy. `uno` does **not** read raw data or find features — it consumes a
**single-run biosaur2 feature file** (`*.features.tsv`), the same input
ms1searchpy uses, assigns peptide/protein identities, and controls FDR with
target/decoy competition. Single-run input is required: it carries the
per-feature quality columns that separate target from decoy, which cross-run
alignment strips out.

Digestion is reused from [sage-core](../sage); retention-time prediction is a
pure-Rust learned model, with an on-data calibration step per run.

## Workspace layout

| crate | purpose |
| --- | --- |
| `uno-core` | search, scoring, and FDR library |
| `uno` | command-line binary |
| `uno-rt-train` | trains/distills the RT (and IM/CCS) models |
| `uno-lightgbm-sys` | FFI bindings to a vendored LightGBM |

## Prerequisites

- A recent stable Rust toolchain (`rustup`, edition 2021).
- **A sibling `sage` checkout.** The workspace has a path dependency
  `sage-core = { path = "../sage/crates/sage" }`, so `sage` must be cloned next
  to this repo before it will build.

```
Repos/
├── sage/    # https://github.com/lazear/sage
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
  --rt-model models/rt_cnn_chrono_global_v3_ch32.urtc \
  --config example_config.toml \
  --out-proteins proteins.tsv
```

See `example_config.toml` for all tunables (digest, mass tolerance, `top_k`) and
their documented defaults. `uno --help` lists every flag, including the DeepLC
RT-library workflow (`--dump-candidate-peptides` → `--rt-library`).

## Repository contents

- `crates/` — the Rust workspace.
- `models/` — distilled RT/IM models used at search time.
- `paper/` — the manuscript (Typst source + PDF).
- `figures/` — figures for the paper.

Build artifacts (`target/`), the Python virtualenv (`.venv/`), RT training data
(`data/`), and scratch outputs (`experiments/`) are intentionally not tracked.

## License

MIT © Patrick Garrett
