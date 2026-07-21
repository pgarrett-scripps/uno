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

## Run

```sh
./run.sh                # N_MS1=3 by default
N_MS1=5 ./run.sh        # more ms1 replicates
```

Requires: `cargo build --release` (uno), and `../.venv` with `ms1searchpy` +
`deeplc 3.1.13`. Writes `results/summary.tsv` and per-run logs under `results/`.

## Files

| path | what |
|---|---|
| `data/` | shared inputs (gitignored): the biosaur2 features + the entrapment FASTA |
| `uno.toml` | uno digest aligned to ms1searchpy |
| `entrapment.py` | unified entrapment scorer (reads uno *or* ms1 protein output) |
| `run.sh` | runs both engines on identical input; captures depth + time + RSS |
| `results/` | outputs + `summary.tsv` (gitignored) |

## Adding a dataset

Drop `data/<name>.features.tsv` (biosaur2) + point `FEATURES` in `run.sh` at it (or
parameterize). The FASTA and everything else stay the same.
