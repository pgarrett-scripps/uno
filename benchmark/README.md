# uno vs ms1searchpy — benchmark

Apples-to-apples comparison of **uno** and **ms1searchpy** on **identical input**,
scored by **yeast entrapment** on **depth and time**.

## Protocol

**Held identical**
- **Input:** one biosaur2 single-run feature file (`data/run10.features.tsv`) fed to
  *both* engines. This is the whole point — same features, same run, with the
  per-feature quality columns (`cosine_score`, `isotope_score`, `ppm_error`,
  `nScans`, …) intact. No aligned/consensus matrices.
- **Database:** one human+yeast FASTA (`data/entrap_human_yeast.fasta`,
  20 647 human + 6 066 yeast). Each engine generates its **own** decoys internally.
- **Digest:** peptides 7–30, 0 missed cleavages, Carbamidomethyl@C
  (`uno.toml` aligns uno to ms1searchpy's defaults, so only the *engines* differ).
- **RT:** both use a learned sequence RT predictor — uno its CNN, ms1searchpy
  **DeepLC** (`-deeplc 1`, always).

**Measured**
- **Depth:** yeast-entrapment FDP. The sample is human (HeLa); yeast is absent, so
  any accepted `_YEAST` protein is a false positive. Scale-corrected:
  `FDP(k) = yeast·(20647/6066)/human` at rank `k` (targets ranked by score desc),
  `q` = running-min. Reported at **1%** and **5%**. Same scorer for both engines
  (`entrapment.py`, schema auto-detected).
- **Time:** `user-CPU` seconds (load-immune) + wall-clock + peak RSS, via
  `/usr/bin/time -v`.

**Fairness controls**
- ms1searchpy is **non-deterministic** (random decoy shuffle + random HP search),
  so it is run **N=3 replicates** (median + range); uno is byte-deterministic → 1 run.
- Runs are **serialized** — no concurrent load corrupting the timings.

## Datasets

**Canonical (paper-faithful): the original DirectMS1 data.** Five high-loading 5-min HeLa
runs from **PXD015669** (Ivanov et al., DirectMS1, Anal Chem 2020), processed through the
papers' own **biosaur2** feature-finder (no koth post-processing — biosaur2-only is how
ms1searchpy is normally run). All HeLa searched against the human+yeast entrapment DB.

| name | features | source |
|---|---|---|
| `hela5min_09746/09766/09734/09778/09950` | 60k–87k | PXD015669, Q Exactive HF, 5-min, no FAIMS |
| `faims_LUM2_05433_VG` | 101k | **PXD023977** (ms1searchpy/DirectMS1Quant "2000 proteins"), Lumos + FAIMS |

The FAIMS run is one file with 3 CV fractions (−50/−65/−80 V); uno reads it as one run
(ignores CV) and still reaches ~2000 @1%.

**Legacy (koth-processed, not paper datasets):** `run10` (211k, medium HeLa), `deephela`
(305k, 240-min — off-design long gradient), `09744` (8.8k, low-loading PXD015669). These
carry koth-extra quality columns (proven inert). Kept for continuity; the canonical set
above is the proper benchmark.

### Rebuilding a dataset from raw

```sh
# raw -> mzML (docker; pwiz/thermorawfileparser images cached)
docker run --rm -v $PWD:/data quay.io/biocontainers/thermorawfileparser:1.4.5--ha8f3691_0 \
  ThermoRawFileParser -i /data/RUN.raw -o /data -f 2
# mzML -> features (biosaur2)
biosaur2 RUN.mzML -o data/NAME.features.tsv
```
Raw files: `https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/02/PXD015669/`

## Run

```sh
./run.sh                       # all datasets, N_MS1=3 reps each
DATASETS="run10 09744" ./run.sh   # subset
N_MS1=5 ./run.sh               # more ms1 replicates
```

Requires: `cargo build --release` (uno), and `ms1searchpy` + `deeplc 3.1.13` (see
`MS1=` in `run.sh`). Writes `results/summary.tsv` (with a `dataset` column) and per-run
logs under `results/`.

## Files

| path | what |
|---|---|
| `data/` | shared inputs (gitignored): the biosaur2 features + the entrapment FASTA |
| `uno.toml` | uno digest aligned to ms1searchpy |
| `entrapment.py` | unified entrapment scorer (reads uno *or* ms1 protein output) |
| `run.sh` | runs both engines on identical input; captures depth + time + RSS |
| `results/` | outputs + `summary.tsv` (gitignored) |

## Adding a dataset

Drop `data/<name>.features.tsv` (biosaur2, HeLa, with quality columns) into `data/`.
`run.sh` auto-discovers every `data/*.features.tsv`; the FASTA and everything else stay
the same.
