// Software and hardware used for every measurement in this manuscript.
// Sources: benchmark/README.md, benchmark/run.sh, `rustc --version`,
// `importlib.metadata.version(...)` inside the ms1searchpy virtualenv,
// `nproc` / /proc/cpuinfo / `free -g` on the benchmark host.
#table(
  columns: 3,
  align: (left, left, left),
  table.header([Component], [Version], [Role]),
  [uno], [0.1.0, commit `7c80199`], [MS1-only search and quantitation],
  [rustc], [1.91.0], [uno build (`cargo build --release`)],
  [ms1searchpy], [2.8.12], [reference MS1-only engine],
  [DeepLC], [3.1.13], [ms1searchpy retention-time model],
  [LightGBM], [4.6.0], [ms1searchpy rescoring],
  [pyteomics], [4.7.5], [ms1searchpy dependency],
  [biosaur2], [0.2.28], [feature detection, both engines],
  // koth version from koth_rust/Cargo.toml, commit from `git log -1`.
  [koth], [0.1.0, commit `eac09b7`], [feature detection and cross-run alignment,
  quantitation only],
  [ThermoRawFileParser], [1.4.5], [`.raw` to mzML conversion],
  [Python], [3.12.3], [ms1searchpy environment, scoring scripts],
  [NumPy / SciPy], [1.26.4 / 1.16.2], [scoring scripts],
  table.hline(),
  [CPU], [Intel Core i7-12700H, 14 cores / 20 threads], [all timings],
  [Memory], [62 GB], [all peak-RSS figures],
)
