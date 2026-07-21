# Benchmark results

**run10** — biosaur2 single-run features (`run10.features.tsv`, 211 064 features), human+yeast
entrapment FASTA (20 647 human / 6 066 yeast), digest aligned (7–30, 0 MC, Carbamidomethyl@C).
Same entrapment scorer for both engines. ms1searchpy run with **DeepLC** (`-deeplc 1 -ad 1`).
Serialized; ms1 3 replicates (non-deterministic), uno 1 (deterministic).

_Generated 2026-07-21. `entrapment.py` on each engine's protein output; `/usr/bin/time -v` for
user-CPU (load-immune) + wall + peak RSS._

| tool | depth @1% | depth @5% | wall (s) | user-CPU (s) | peak RSS (GB) |
|---|---|---|---|---|---|
| **uno** (deterministic) | **1309** | 1887 | **46** | **436** | **3.9** |
| ms1 rep1 | 1220 | 1941 | 352 | 1998 | 9.5 |
| ms1 rep2 | 1155 | 1651 | 360 | 2060 | 9.5 |
| ms1 rep3 | 1079 | 1534 | 397 | 2012 | 9.6 |
| **ms1 median** | **1155** | **1651** | **360** | **2012** | **9.5** |

## Verdict — uno wins on every axis

- **Depth @1%:** 1309 vs ms1 median 1155 (**+13%**), above ms1's best replicate (1220).
- **Depth @5%:** 1887 vs median 1651 (**+14%**). ms1's best rep (1941) edges it — that's variance.
- **Speed:** **~7.8× faster** wall (46 s vs 360 s), **~4.6× less** user-CPU.
- **Memory:** **~2.4× less** (3.9 vs 9.5 GB).
- **Determinism:** uno is exact; ms1 swings **1079–1220 @1%** and **1534–1941 @5%** across
  identical reruns (random decoy shuffle + random HP search).

## Notes

- **Honest metric.** ms1's *own* internal-FDR count was 1319 @1% (rep1), but by the shared
  **yeast-entrapment** metric it is 1220 — its target/decoy FDR is looser than the entrapment
  truth. uno's 1309 already is the honest-entrapment number.
- **Feature usage (A/B tested).** uno's classifier *reads* the biosaur quality columns
  (`cosine_score, isotope_score, ppm_error, neutron_offset, combined_score`), the same ones
  ms1 feeds LightGBM — but they are **not** what drives uno's depth. Stripping all of them
  changes honest depth by ~0 @1% (1309 → 1306) and only +6.5% @5% (1887 → 1772); the ML
  ensemble still jumps 737 → ~1240 confident without them. uno's depth comes from its
  **ensemble band-averaging + mass/RT/composition** features, not the quality columns.
  Single-run input still matters — for a fair same-input comparison and the small @5% gain —
  just not as a big depth lever.
- **RT.** uno used its CNN; ms1 used DeepLC (second-stage calibrated RT σ ≈ 4.2 min).

Reproduce: `./run.sh` (writes `results/summary.tsv`).
