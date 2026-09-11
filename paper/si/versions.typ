// Software and hardware used for every measurement in this manuscript.
// TRANSCRIBED BY HAND from the environment itself, 2026-09-11:
//   uno            `git -C . rev-parse HEAD`, `sha256sum target/release/uno`
//   path deps      `git -C ../{Sequoia,sage,koth_rust} log -1`, plus `git diff`
//                  where the dependency carries local modifications
//   rust           `rustc --version`, `cargo --version`
//   comparator     `importlib.metadata.version(...)` inside the ms1searchpy
//                  virtualenv; patch state from `git -C ../ms1searchpy diff`
//   analysis env   `uv run python -c "import matplotlib, numpy, pandas"` in
//                  paper/analysis (matplotlib is pinned in its pyproject.toml)
//   host           `nproc` / /proc/cpuinfo / `free -g` / `uname -sr`
// Every uno row in this manuscript was produced by the single commit below;
// benchmark/results/publication_rerun/manifest.json records the binary hash
// beside every search.
#table(
  columns: 3,
  align: (left, left, left),
  table.header([Component], [Version], [Role]),
  [uno],
  [commit `0710225`, binary SHA-256 `c376954d`],
  [every uno search reported here],
  [rustc / cargo], [1.97.1], [uno build (`cargo build --release`)],
  [sage-core],
  [commit `0639176` plus a local patch, diff SHA-256 `4a45766d`],
  [digestion and peptide mass, path dependency],
  [sequoia-boost],
  [commit `5c2772e`],
  [pure-Rust gradient boosting, path dependency],
  [koth_ff],
  [0.3.0, commit `c80a133`],
  [built-in feature detection (`uno --raw`), path dependency],
  [ms1searchpy],
  [commit `9df4a94`, reported version 2.8.12, plus the comparator patch
    `benchmark/ms1searchpy_comparator.patch` (SHA-256 `f8d0c455`)],
  [reference MS1-only engine],
  [DeepLC], [3.1.13], [ms1searchpy retention-time model],
  [LightGBM], [4.6.0], [ms1searchpy rescoring],
  [pyteomics], [4.7.5], [ms1searchpy FDR and dependency],
  [biosaur2], [0.2.28], [feature detection, both engines],
  [ThermoRawFileParser], [1.4.5], [`.raw` to mzML conversion],
  [Python (comparator)], [3.12.3], [ms1searchpy environment],
  [NumPy / SciPy], [1.26.4 / 1.16.2], [ms1searchpy dependencies],
  [Python (analysis)], [3.14.0], [entrapment scoring, figures, tables],
  [Matplotlib / NumPy / pandas],
  [3.11.0 / 2.5.1 / 3.0.5],
  [figure and table generation],
  [Typst], [0.14.2], [manuscript build],
  table.hline(),
  [CPU], [Intel Core i7-12700H, 14 cores / 20 threads], [all timings],
  [Memory], [62 GB], [all peak-RSS figures],
  [Operating system], [Ubuntu 24.04.4 LTS, Linux 7.0.0-30-generic], [all measurements],
)
