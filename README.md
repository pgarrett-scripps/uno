# uno — MS1-only search engine (Rust)

A Rust MS1-only protein identification engine, in the spirit of DirectMS1 /
ms1searchpy. `uno` accepts a single-run biosaur2 or KOTH feature TSV (`--features`),
finds features with embedded KOTH (`--raw`), or searches a KOTH aligned intensity
matrix across several runs (`--matrix`). It assigns peptide/protein identities
and estimates FDR with target/decoy competition.

See [matrix input](docs/matrix-input.md) for multi-run searching and the
experimental confidence-aware matrix schema. Matrix input retains cross-run
intensities but omits some quality descriptors from individual feature tables.

Digestion follows [Sage](https://github.com/lazear/sage)'s database builder (ported into `uno-core`); retention-time prediction is a
pure-Rust learned model, with an on-data calibration step per run.

## Experimental quantitative options

`uno` includes `uno quant` ([details](crates/uno-core/src/quant/directms1/README.md)), a Rust
DirectMS1Quant-style contrast analyzer in `uno-core`, and opt-in `--out-quant-evidence` exports
from independent uno searches. The existing KOTH-matrix route remains the second
quantitative option. Native-reference parity is tested; biological accuracy and
error calibration are pending.

## Workspace layout

| directory | crates.io package | purpose |
| --- | --- | --- |
| `crates/uno-core` | `uno-ms-core` | search, scoring, FDR, and quantitative contrast library |
| `crates/uno` | `uno-ms` | command-line binary `uno` |
| `crates/uno-lightgbm-sys` | `uno-ms-lightgbm-sys` | FFI bindings to a vendored LightGBM |
| `crates/uno-rt-train` | not published | legacy RT/IM tree-model research tools (`uno-rt-train` binary) |
| `crates/uno-rt-train/python` | not published | optional PyTorch sequence-CNN trainer for `uno rt train` and `uno im train` |

## Install

uno is published on crates.io as `uno-ms` (the `uno` name there belongs to an
unrelated crate); the installed command is still `uno`:

```sh
cargo install uno-ms
```

This needs Rust 1.97 or newer. Bruker `.d` and Thermo `.raw` input (through
[koth-ms](https://github.com/pgarrett-scripps/koth), pure Rust, no .NET runtime)
are built in; `cargo install uno-ms --no-default-features` drops `.raw`.
**Platform:** `cargo install uno-ms` currently supports Linux x86_64 only.
LightGBM is linked from a vendored prebuilt `lib_lightgbm.so`.

## Build from source

```sh
git clone https://github.com/pgarrett-scripps/uno
cd uno
cargo build --release   # binary at target/release/uno
```

All dependencies are fetched by Cargo; no sibling checkouts are needed.

## Citation

If you use uno, please cite it using the metadata in [CITATION.cff](CITATION.cff).

## Run

```sh
cargo run --release -p uno-ms -- \
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

The current search includes bounded RT recovery:
when the usual search and its evidence-retention retry find no target proteins
at 1% FDR, UNO attempts one additional broad RT calibration. Successful results
are returned unchanged. Validation showed fewer empty low-input runs with
0–1% normal-path runtime overhead; variation between successful seeds remains.

## Repository contents

- `crates/` — the Rust workspace.
- `models/` — notes on the bundled RT and IM CNNs (the weights ship in `crates/uno-core/models/`).
- `docs/` — model, matrix-input and quantification notes.
- `tools/quant-reference/` — the DirectMS1Quant reference used for parity checks.
- `tests/fixtures/` — software regression-test configurations.

Build artifacts (`target/`), the Python virtualenv (`.venv/`), RT training data
(`data/`), and scratch outputs (`experiments/`) are intentionally not tracked.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) and the [Code of Conduct](CODE_OF_CONDUCT.md).
Report security issues privately as described in [SECURITY.md](SECURITY.md).

## License

MIT © Patrick Garrett
