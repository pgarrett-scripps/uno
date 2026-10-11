# URTC retention-time CNN

For the supported CLI/API workflow, see [Custom RT models](../docs/rt-models.md).
`uno rt train` adapts this CNN and exports a self-describing URTC v3 artifact;
search embeds the shipped model and needs no Python. Historical benchmarks below
use earlier scoring configurations and should not be treated as current results.

A sequence-aware RT model for uno, trained on **experimental iRT** (not distilled
from DeepLC): learned AA embedding → residual dilated Conv1d stack → masked
mean+max pool → MLP head. Pure-Rust inference in [`crates/uno-core/src/rt_cnn.rs`](../crates/uno-core/src/rt_cnn.rs)
(`URTC` binary, auto-detected by file magic in `RtModel::load`).

Replaces the positional one-hot MLP (`rt_posmlp.bin`): a convolution shares
weights across positions and reads a learned embedding, so it generalizes far
better than absolute-position one-hot. AA-only — C is implicitly carbamidomethyl
(uno's fixed mod); no variable mods.

## Files
- `rt_cnn_chrono_global_v3_ch32.urtc` — **SHIPPED DEFAULT** (2026-07-09): ch32 on the
  Chronologer corpus **plus 11 global physicochemical features** (URTC v2 format,
  N_GLOBAL=11). Best in-domain val MAE (0.447). ⚠️ On the run10 honest benchmark it
  scores 248 @1% — LOWER than the 5-feature variant (401) and the no-global baseline
  (259); shipped by explicit product decision (val-MAE / broader-generalization
  rationale), not because it wins run10. The build computes 11 global features and also loads the historical
  5-feature prefix models. See "Global features" below.
- `rt_cnn_chrono_global_ch32.urtc` — 5-feature variant (URTC v2, N_GLOBAL=5). Best
  run10 depth of any model tested: **401 @1%** (+55% over baseline). Superseded as
  default by v3 above; loads with the current build through the supported 5-feature prefix.
- `rt_cnn_chrono_ch32.urtc` — original plain ch32, no global feats (URTC v1, still
  loadable): Chronologer corpus, 259 @1% on run10. ~40k params, 60k pred/s.
- `rt_cnn_chrono_global_ch64.urtc` — ch64 + 5 global feats. Best val MAE among the
  v2 set (0.423) but REGRESSES in-search (229 @1%): ch64 overfits, loses under domain
  shift. Do NOT use; uses the supported N_GLOBAL=5 feature prefix.
- `rt_cnn_chrono_ch64.urtc` — plain ch64: +12% depth (760 vs 670) under the OLD config below, 4x slower.
- `rt_cnn_chrono_ch16.urtc` — same data, ch=16: fastest (205k pred/s) but loses depth.
- `rt_cnn_ch32.urtc` / `rt_cnn_ch64.urtc` — OLD, Prosit-iRT-only. Do NOT use: they
  overfit synthetic iRT and HALVE in-search depth (546 vs 760 for ch64). Kept for history.

## In-search benchmark (run10, honest entrapment — THE metric that matters)
ProteomeTools holdout MAE is NOT predictive of search depth (the old Prosit CNN won
holdout but lost 2x in search). Always judge with `scripts/validate_rt_in_search.sh`.
Honest = yeast-covered entrapment library; emp FDR from yeast IDs.

| RT model                  | honest proteins | emp FDR | run10 PWL-MAE | pred/s | s/1M |
|---------------------------|-----------------|---------|---------------|--------|------|
| DeepLC generic (Python)   |       827       |  3.0%   |    0.76 min   | ~1.7k  | ~600 |
| chrono ch64               |       760       |  2.7%   |    2.04 min   |  15k   |  66  |
| **chrono ch32 (default)** |     **670**     |  1.2%   |    2.04 min   |  60k   |  21  |
| chrono ch16               |       476       |  1.7%   |    2.13 min   | 205k   |   5  |
| OLD Prosit ch64           |       546       |  1.5%   |    2.65 min   |  15k   |  66  |

chrono ch32 = 81% of DeepLC's depth, pure-Rust, ~40x faster, no Python. Retrain/extend:
see [RT models](../docs/rt-models.md). Only the default model ships in this repository.

## Global features (URTC v2) — the current lever
DeepLC's edge over uno's CNN on run10 (residual σ 1.44 vs 5.07 min) traces partly
to its explicit global physicochemical descriptors. Appending 5 of them to uno's
`[mean;max]` pool before the head closes much of the gap. `N_GLOBAL = 5`:
`[len/50, (D+E)/len, (K+R)/len, His/len, (FWYLIVM)/len]` — must stay in lockstep
between `global_feats` (rt_cnn.rs) and `feat_global` (train_rt.py). v1 models carry
n_global=0 and load unchanged; v2 stores the count in the header.

A/B on run10 (entrap_db_C; max_len=30, offsets[-1,0,1], top_k=50, ondata_rt;
honest yeast-entrapment FDR — NOT comparable to the 670-row table above, which used
a different config):

| RT model            | val MAE (iRT) | honest @1% | honest @5% |
|---------------------|---------------|------------|------------|
| ch32 (no global)    |     ~0.50     |    259     |    690     |
| **ch32 + global**   |     0.456     |  **401**   |    707     |
| ch64 + global       |     0.423     |    229     |    568     |

The gain is NOT a tighter global RT sigma (two-stage rt_sig pins at the 3.50 clamp
in both) — it's that more accurate per-peptide RT predictions win the target/decoy
deviation competition, so the grid picks a tighter cut (dev_cut 6→3). ch64 again
proves holdout/val MAE is anti-predictive of search depth.

**On run10 depth, 5 features beat 11.** The v3 set (N_GLOBAL=11: + KD-hydropathy sum,
GRAVY, aromatic FWY, Pro, N-term KD, C-term KD) scored 248 @1% — below the 5-feature
401 and even the 259 baseline — despite the best val MAE of all (0.447). Likely cause:
the hydropathy sum ≈ the HI label, so the net leans on it and predicts a smooth
composition-only RT that matches many candidates equally, eroding the per-peptide
discrimination the target/decoy competition needs.

Full A/B (run10, same config; honest yeast-entrapment FDR):

| RT model                | val MAE (iRT) | honest @1% | honest @5% |
|-------------------------|---------------|------------|------------|
| ch32 (0 global)         |     ~0.50     |    259     |    690     |
| ch32 + 5 global         |     0.456     |  **401**   |    707     |
| ch32 + 11 global (v3)   |   **0.447**   |    248     |    648     |
| ch64 + 5 global         |     0.423     |    229     |    568     |

**v3 (11 global) is the SHIPPED default** by product decision despite the lower run10
depth — it has the best in-domain val MAE and its features are domain-invariant
physics (rationale: broader-dataset generalization). If optimizing run10 depth
specifically, use `rt_cnn_chrono_global_ch32.urtc` (401).
To push RT further, don't pile on global scalars — ablate + add order-aware/atom
encoding one change at a time.

## Benchmark (ProteomeTools/Prosit holdout, peptide-disjoint split) — HISTORICAL, misleading
Aligned = after affine calibration to true iRT (how uno calibrates per-run).
Replicate noise floor on this data is ~0.28 iRT MAE.

| model            | aligned MAE | dt95  | pred/s (14 core) | 1M peptides |
|------------------|-------------|-------|------------------|-------------|
| cnn_ch64         | 1.98        | 4.93  | 15,300           | ~65 s       |
| **cnn_ch32**     | **2.21**    | 5.67  | **67,000**       | **~15 s**   |
| posMLP (real iRT)| 3.52        | 9.64  | 233,000          | ~4 s        |
| DeepLC (no calib)| 4.58        | 12.0  | ~1,700           | ~10 min     |

The CNN beats DeepLC ~2x on MAE/dt95 and runs ~40x faster (DeepLC is TF/CPU).
Rust inference is parity-exact with PyTorch (max |diff| ~1e-4 iRT).

NOTE: the DeepLC comparison is in-domain (both vs ProteomeTools) and DeepLC is run
without per-dataset calibration; uno's downstream per-run calibration is exactly
the affine alignment used here, so this is the relevant comparison for uno.

## Retrain
Data prep + splits: see [`../data/rt/README.md`](../data/rt/README.md).
```
# in the GPU venv (torch + cuda)
python crates/uno-rt-train/python/train_rt.py --arch cnn --ch 32 --blocks 1,2,4 \
  --train data/rt/prosit/train.tsv --val data/rt/prosit/val.tsv \
  --out models/rt_cnn_ch32.urtc --epochs 60 --bs 4096 --patience 10
```
Parity check (torch oracle vs Rust): `urtc_torch_predict.py` + `rt_eval` dump,
diff with `compare_rt.py`. Head-to-head vs DeepLC: `deeplc_predict.py`.

## Architecture notes
- Trained with masking after every conv → PAD positions stay zero → Rust runs
  over the real peptide length (no fixed-50 padding work). Keep this if retraining.
- Inference is memory-bandwidth-bound on the conv weights (ch²·k per position),
  NOT compute-bound: `target-cpu=native` gave ~0 speedup, so the portable build
  is fine. Smaller `ch` is the speed lever (ch32 is 4x faster than ch64).
- Still improving slowly at 60 epochs; more epochs / LR tuning could gain a little.
