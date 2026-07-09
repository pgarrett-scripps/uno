# URTC retention-time CNN

A sequence-aware RT model for uno, trained on **experimental iRT** (not distilled
from DeepLC): learned AA embedding → residual dilated Conv1d stack → masked
mean+max pool → MLP head. Pure-Rust inference in [`crates/uno-core/src/rt_cnn.rs`](../crates/uno-core/src/rt_cnn.rs)
(`URTC` binary, auto-detected by file magic in `RtModel::load`).

Replaces the positional one-hot MLP (`rt_posmlp.bin`): a convolution shares
weights across positions and reads a learned embedding, so it generalizes far
better than absolute-position one-hot. AA-only — C is implicitly carbamidomethyl
(uno's fixed mod); no variable mods.

## Files
- `rt_cnn_chrono_ch32.urtc` — **DEPLOYED DEFAULT** (2026-06-27): trained on the
  Chronologer harmonized corpus (1.5M real-LC peptides, HI space). The
  speed/accuracy sweet spot — see in-search table below. ~40k params, 60k pred/s.
- `rt_cnn_chrono_ch64.urtc` — same data, ch=64: +12% depth (760 vs 670), 4x slower.
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
see [`../crates/uno-rt-train/PRETRAIN_PLAN.md`](../crates/uno-rt-train/PRETRAIN_PLAN.md).

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
