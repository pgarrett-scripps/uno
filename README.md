# uno — MS1-only search engine (Rust)

A Rust MS1-only protein identification + quantitation engine, in the spirit of
DirectMS1 / ms1searchpy. `uno` does **not** read raw data or find features — it
consumes an aligned LFQ matrix (from [koth](../koth_rust)), assigns
peptide/protein identities, controls FDR with target/decoy competition, and
rolls the intensities up into peptide/protein quant matrices.

Digestion is reused from [sage-core](../sage); retention-time prediction is
distilled from DeepLC into a pure-Rust model, with an on-data calibration/
fine-tuning step per run.

## Workspace layout

| crate | purpose |
| --- | --- |
| `uno-core` | search, scoring, FDR, and quant-rollup library |
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
  --intensity-matrix path/to/intensity_matrix.tsv \
  --fasta path/to/proteins.fasta \
  --rt-model models/rt_cnn_chrono_ch32.urtc \
  --config example_config.toml \
  --out-proteins proteins.tsv \
  --out-peptides peptides.tsv
```

See `example_config.toml` for all tunables (digest, mass tolerance, `top_k`,
on-data RT fine-tuning) and their documented defaults. `uno --help` lists every
flag, including quant outputs (`--out-protein-quant`, `--out-peptide-quant`),
the peptide-centric hill second pass (`--hills`), and DeepLC RT-library
workflow (`--dump-candidate-peptides` → `--rt-library`).

## Repository contents

- `crates/` — the Rust workspace.
- `models/` — distilled RT/IM models used at search time.
- `paper/` — the manuscript (Typst source + PDF).
- `figures/` — figures for the paper.

Build artifacts (`target/`), the Python virtualenv (`.venv/`), RT training data
(`data/`), and scratch outputs (`experiments/`) are intentionally not tracked.

## License

MIT © Patrick Garrett
