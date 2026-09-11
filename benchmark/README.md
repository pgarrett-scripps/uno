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
- **Both engines are randomized, and both are measured as spreads.** uno's
  byte-identical output is a property of its default seed, not of the algorithm:
  `uno --seed N` perturbs the decoy shuffle, RANSAC sampling, CV folds and
  booster subsampling together (seed 0 = the shipped constants, so published
  numbers still reproduce). `seed_sweep.sh` runs uno over N seeds,
  `ms1_extra_reps.sh` runs ms1searchpy over N replicates, and
  `spread_table.py` scores both into `results/spread.tsv`. Comparing one locked
  uno run against three unlocked ms1searchpy runs measures the locking, not the
  engines. `ms1_seeded.py` shows the converse: ms1searchpy also becomes
  reproducible once its generator, its hyperparameter search, LightGBM, and
  `PYTHONHASHSEED` are all pinned.
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

### Variability in the manuscript

The six primary files use the complete draw inventory in `results/spread.tsv`
(11 uno seeds per file, 3–4 comparator runs). The inventory drives depth,
protein-list overlap, and per-run calibration. Resource plots retain only the
original serialized searches in `results/summary.tsv`.

`variability_extension.py` measures seeds 0, 1, 2 for both feature sources on
all 18 primary/Astral files. It reuses archived seed-zero outputs and the
primary biosaur2 sweep, then runs the 60 missing searches two at a time. Run
from `paper/analysis` with `uv run ../../benchmark/variability_extension.py`.
Successful outputs carry command and SHA-256 provenance; `--summarize` validates
and re-scores them. New outputs stay in `results/variability_extension/`.
Concurrent timings are excluded from resource claims. Within-file seeds and
between-file acquisition variation are kept separate in the Astral figure.

### Engine stability study

`stability_study.py` compares two frozen binaries on matched feature files and
seeds, records commands and SHA-256 hashes, and resumes only matching outputs.
It writes per-run depth, reported-list entrapment FDP, paired differences, and
within-file protein overlap. Empty lists and undefined FDP stay distinct.

```sh
python3 benchmark/stability_study.py \
  --baseline /path/to/baseline/uno --candidate /path/to/candidate/uno \
  --output benchmark/results/new-study/lowinput \
  --datasets lowinput --seeds 0-4 --jobs 2 --threads 4 \
  --env UNO_RT_CACHE_DIR=/tmp/uno-rt-cache
```

Run `--datasets primary --seeds 0-2` in a separate output directory for the
primary guard panel. Optional `--baseline-study DIR` copies only verified matching
baseline results; missing seeds run afresh and every candidate search stays fresh.
For the accepted V3 validation, specify `--seeds 0-8 --held-out-seeds 7-8` and pass
`--lowinput-seed-count 9 --protect-discoveries` to the report.
`stability_report.py --studies DIR/lowinput DIR/primary
--protocol PATH/TO/FROZEN-PLAN.md --output DIR/report` verifies completed
outputs and evaluates the documented acceptance checks. `plot_stability.py`
renders every low-input seed in standalone figures using the paper analysis
environment. See [the study protocol and findings](../docs/stability-study.md).

### The publication re-run

`publication_rerun.py` re-executes the complete uno side of the manuscript --
156 searches: six primary files x biosaur2 x seeds 0-10, six primary files x
`--raw` x seeds 0-2, and twelve PXD066701 Astral files x both feature sources x
seeds 0-2 -- against one pinned binary, so every reported uno number comes from
a single commit. Seed-0 searches run serialized, because their `/usr/bin/time -v`
output is what the resource panels report; the remaining seeds run two at a
time and their timings are recorded but not reported.

Outputs go to the canonical paths the figure scripts already read. The
superseded tables are copied first into `results/archive_<old build>/` with
their hashes and depths, and `results/publication_rerun/cases.tsv` then reports,
per search, the new depths beside the old ones and whether the protein table is
byte-identical to the superseded build's.

```sh
python3 benchmark/publication_rerun.py --plan   # inventory, no searching
python3 benchmark/publication_rerun.py          # run it (restartable)
python3 benchmark/refresh_summaries.py          # rebuild every derived summary
```

`refresh_summaries.py` re-scores what is on disk -- `results/spread.tsv`, the
seed-sweep and variability-extension summaries, the uno rows of the serialized
resource cohort, `raw_mode/summary.tsv` and the Astral input-series summary. It
never recomputes a comparator number from a new search: ms1searchpy did not
change, so its rows are carried through unchanged.

The comparator itself is pinned by `ms1searchpy_comparator.patch`, the diff
between the installed ms1searchpy checkout and its last commit; it is what gives
ms1searchpy the per-feature quality columns uno reads.
