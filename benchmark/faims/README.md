# FAIMS benchmark: uno vs ms1searchpy (DirectMS1 PXD023977)

DirectMS1/ms1searchpy's own FAIMS "2000-protein" benchmark: 3 replicate 200 ng
HeLa Lumos+FAIMS runs. Apples-to-apples — SAME biosaur2 features, SAME target-only
human+E.coli DB (both engines make their own decoys), honest E. coli-entrapment FDP.

## Data (download once)
ProteomeXchange PXD023977 (PRIDE):
  https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/03/PXD023977/
Grab: LUM2_054{33,34,35}_VG.features.tsv, sprot_human_and_ecoli_shuffled.fasta.
Kept at /mnt/labdata/pgarrett/PXD023977_FAIMS (override with DATA=...).

## Run
    ./run_faims.sh                 # DATA=... MS1=... OUT=... to override

## Result (re-measured 2026-07-28, honest E. coli-entrapment, mean of 3 reps)
| metric | uno  | ms1searchpy |
|--------|------|-------------|
| @1% FDR| 1591 | 1474  (+8%   uno) |
| @5% FDR| 2060 | 1915  (+7.6% uno) |

Per replicate (LUM2_05433 / 05434 / 05435):

| metric | uno              | ms1searchpy      |
|--------|------------------|------------------|
| @1%    | 1642 / 1460 / 1672 | 1662 / 1205 / 1554 |
| @5%    | 2112 / 1863 / 2204 | 2007 / 1878 / 1861 |

uno wins both means and 4 of the 6 replicate x threshold cells. ms1searchpy is
deeper on 05433 @1% (1662 vs 1642) and on 05434 @5% (1878 vs 1863). ms1searchpy
still swings rep-to-rep (1205-1662 @1%) while uno is byte-deterministic: two
independent invocations gave identical uno numbers. uno ~2060 @5% = DirectMS1's
headline "2000". Complements the Orbitrap (no-FAIMS) result in ../README.md.

### 2026-07-28 re-measurement note (the 2026-07-27 numbers did not reproduce)

Rebuilt uno (`cargo build --release`, commit 7c80199) and re-ran the whole
benchmark. Commands:

    ./benchmark/faims/run_faims.sh
    python paper/scripts/figures/figS2_faims.py   # scores both engines, prints the table

Saved outputs: `benchmark/faims/faims_out/` (per-engine protein tables, engine
logs, `run_faims_stdout.txt`, and `figS2_stdout.txt` which is the table above).

What changed and why:

- **uno reproduced exactly.** 1591 @1% and 2060 @5%, byte-identical across two
  independent invocations. No change.
- **ms1searchpy did NOT reproduce.** Measured 1474 @1% (was recorded 1432) and
  1915 @5% (was 1883); per-replicate 1662/1205/1554 @1% (was 1659/1136/1503).
  Root cause: `run_faims.sh` used relative `MS1`/`ROOT`/`OUT` paths while the
  ms1searchpy step runs inside `cd "$OUT"`, so **the committed script never
  actually ran ms1searchpy** and the recorded numbers came from an ad-hoc
  invocation that was not captured. The path bug is fixed in this commit.
- **The margin is therefore +8% / +7.6%, not +11% / +9%**, and the claim "uno
  beats ms1searchpy on every replicate" was never true: ms1searchpy wins
  05433 @1% and 05434 @5%. Report the per-replicate table, not just the mean.
