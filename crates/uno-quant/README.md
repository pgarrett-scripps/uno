# uno-quant (experimental)

Rust implementation of the two-condition DirectMS1Quant workflow, pinned to
ms1searchpy `9df4a941a901aac657b1c9ab5229cb49671e3c7d`. It is separate from the
uno identification CLI; no existing search behavior is changed by default.

## Native reference inputs

```sh
cargo build --release -p uno-quant
./target/release/uno-quant \
  --s1 runA1_proteins_full.tsv runA2_proteins_full.tsv \
  --s2 runB1_proteins_full.tsv runB2_proteins_full.tsv \
  --out new-quant-directory
```

As in upstream, each input stem needs `_proteins.tsv`, `_PFMs.tsv`, and
`_PFMs_ML.tsv` companions. `--all-proteins`, `--all-pfms`, and explicit
`--allowed-proteins`/`--allowed-peptides` lists select the corresponding native
policies. Default normalization is the sum of the 1,000 largest intensities;
`--normalization 1` selects median and `0` disables normalization. Fold changes
are S2/S1. Inputs need at least two replicates per condition, adequate variation
for calibration, and decoys in every used intensity quintile.

The output directory must not exist. `peptides.tsv` contains intermediate
statistics, resampling counts, and normalized/minimum-filled sample values;
these sample columns are NOT raw observed intensities. `proteins.tsv` includes
all target protein results; `significant.tsv` contains reference-rule calls.
`manifest.json` records input/output/executable hashes, options, calibration,
normalization factors, sample order, and the experimental status.

`score`, `FDR_pass`, and `BH_pass` reproduce the reference definitions. They
are not newly calibrated q-values. Reference protein group labels are opaque;
compare group membership rather than numeric labels.

## Uno per-run evidence (experimental adaptation)

Export each single-run uno search without changing its identification settings:

```sh
./target/release/uno --features run.tsv --fasta database.fasta \
  --out-proteins run.proteins.tsv --out-quant-evidence run.evidence.tsv
./target/release/uno-quant --samples samples.tsv --out new-uno-quant-directory
```

`samples.tsv` has tab-separated columns `sample`, `condition`, `evidence`, and
`proteins`. Conditions must be `S1` and `S2`; paths resolve relative to the
manifest. Include at least two independent files per condition. Evidence files
are schema 1 and preserve feature/peptide indices, charge, native mobility,
observed mass/RT, intensity, explicit decoy labels, and the identity evidence type.
The `run` field retains the input loader's quantitative column label (for example
`intensityApex`); the sample manifest assigns the biological run name explicitly.
The exporter supports single-run input only; use uno's existing matrix path for
the separate KOTH-LFQ matrix option.

The adapter admits proteins identified in any run at the specified protein
ID threshold (`--max-protein-qvalue`, default 0.01) and their prefix-matched decoy
counterparts. Uno stores target and decoy proteins under the same accession;
the authoritative `decoy` flag creates a separate quantitative namespace using
`--decoy-prefix`. No peptide confidence is inferred from that label. It uses uno's
RT/mass-gated assignments instead of inventing an upstream PFM-ML qpred score.
This is an explicit methodological adaptation requiring independent biological
calibration; the union is not claimed to have cohort-wide 1% identification FDR.
Zero intensity means missing signal. Multiple observations of a precursor in
one run reduce to its greatest positive intensity; retained assignments may
still require signal-ownership validation before scientific use.

## Reference parity

```sh
uv sync --project tools/quant-reference --locked
uv run --project tools/quant-reference python tools/quant-reference/check_parity.py \
  --out target/parity-new
uv run --project tools/quant-reference python tools/quant-reference/check_parity.py \
  --out target/parity-missing-new --missing --structural
cargo test -p uno-quant
```

The harness verifies the unmodified reference source hash, creates deterministic
native inputs, and compares peptide values, exact resampling counts, protein
effects/scores/calls, and group partitions. Numeric tolerances are absolute 1e-8
plus relative 1e-6. It covers all three normalizations, absolute/relative effect
thresholds, partial/whole-condition missingness, duplicate observations, multiple
charges, and shared peptides through explicit test options. These are numerical
fixtures, not evidence for accuracy or FDR in biological samples.

NumPy's CPU-dependent quicksort tie ordering affects reference resampling when
fold changes are tied. The reference runner sets
`NPY_DISABLE_CPU_FEATURES=AVX512F,AVX2`; Rust reproduces NumPy 1.26's scalar
ordering and CPython's integer-seeded `random.choice` stream. This qualification
is part of the parity contract, not a change to the upstream Python source.
The first unrestricted missing-data run is retained as a failed diagnostic.

Nonpositive native intensities, empty samples, missing decoy pools, repeated
sample files, and degenerate calibration panels fail explicitly instead of
creating confidence values. FASTA gene annotation, experimental protein shifts,
and multi-condition DirectMS1QuantMulti are not implemented. Cross-platform
floating-point parity and biological calibration remain validation tasks.
