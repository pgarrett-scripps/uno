#!/usr/bin/env python3
"""figS2_faims.png -- FAIMS benchmark (PXD023977), uno vs ms1searchpy.

DirectMS1/ms1searchpy's own FAIMS "2000-protein" benchmark: three replicate
200 ng HeLa Lumos+FAIMS runs (LUM2_05433/05434/05435). Apples-to-apples --
the SAME biosaur2 feature files, the SAME target-only human + E. coli FASTA
(each engine makes its own decoys), and an honest E. coli-entrapment protein
FDP for both.

PROVENANCE -- freshly measured, not transcribed. Run 2026-07-28 from the repo
root with the release binary (`cargo build --release`):

    ./benchmark/faims/run_faims.sh

which, per replicate R in {05433, 05434, 05435}, executes

    target/release/uno --features $DATA/LUM2_${R}_VG.features.tsv \
        --fasta faims_out/hye_target.fasta \
        --rt-model models/rt_cnn_chrono_global_v3_ch32.urtc \
        --config benchmark/uno.toml --out-proteins faims_out/uno_${R}.tsv

    ms1searchpy LUM2_${R}_VG.features.tsv -d hye_target.fasta -deeplc 1 -ad 1
        (from /home/patrick-garrett/Repos/FeatureFinder/.venv/bin/ms1searchpy;
         -deeplc 1 -ad 1 is mandatory in this lab -- achrom RT is never used)

$DATA = /mnt/labdata/pgarrett/PXD023977_FAIMS. Depth is scored by
``benchmark/faims/score_entrap.py``, imported here so the figure and the
benchmark's own console output cannot drift: targets ranked by score, and

    FDP(k) = entrap_in_top_k * (n_human_db / n_ecoli_db) / target_in_top_k

with n_human_db = 20193 and n_ecoli_db = 4431 counted from the target-only
FASTA the script builds. Reported depth is the largest accepted human count
whose FDP is still <= the threshold.

Panels
    a  accepted human proteins at 1% entrapment FDP, per replicate.
    b  the same at 5%.
"""
from __future__ import annotations

import os
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import style  # noqa: E402

FAIMS = style.BENCH / "faims"
OUT = FAIMS / "faims_out"
REPS = ["05433", "05434", "05435"]

# DB sizes for the entrapment scale factor, counted from the target-only FASTA
# that run_faims.sh writes (see the console line "DB: 20193 human + 4431 ecoli").
N_TARGET_DB, N_ENTRAP_DB = 20193, 4431
os.environ.setdefault("N_TARGET_DB", str(N_TARGET_DB))
os.environ.setdefault("N_ENTRAP_DB", str(N_ENTRAP_DB))

sys.path.insert(0, str(FAIMS))
import score_entrap as se  # noqa: E402  (reads the env vars above at import)


def measure():
    """{engine: {rep: {0.01: n, 0.05: n}}} from the run outputs on disk."""
    missing, out = [], {"uno": {}, "ms1searchpy": {}}
    for r in REPS:
        paths = {
            "uno": (OUT / f"uno_{r}.tsv", se.rows_uno),
            "ms1searchpy": (OUT / f"LUM2_{r}_VG.features_proteins_full.tsv", se.rows_ms1),
        }
        for engine, (p, reader) in paths.items():
            if not p.exists():
                missing.append(str(p))
                continue
            out[engine][r] = se.curve(list(reader(str(p))))
    if missing:
        raise SystemExit(
            "FAIMS run outputs missing -- run benchmark/faims/run_faims.sh first:\n  "
            + "\n  ".join(missing))
    return out


def panel(ax, data, thr, title):
    x = np.arange(len(REPS))
    w = 0.36
    for k, engine in enumerate(("uno", "ms1searchpy")):
        col, _ = style.ENGINE[engine]
        v = [data[engine][r][thr] for r in REPS]
        off = (k - 0.5) * w
        mean = float(np.mean(v))
        ax.bar(x + off, v, width=w, color=col, edgecolor="white", zorder=3,
               label=f"{engine} — mean {mean:.0f}")
        for xi, vi in zip(x + off, v):
            ax.annotate(f"{vi}", xy=(xi, vi), xytext=(0, 2),
                        textcoords="offset points", ha="center", va="bottom",
                        fontsize=6.4, color=style.C["neutral"])
        # dashed guide at the 3-replicate mean, labelled in the legend
        ax.axhline(mean, color=col, lw=0.9, ls="--", zorder=2)
    ax.legend(loc="upper center", fontsize=6.8, ncol=1)
    ax.set_xticks(x)
    ax.set_xticklabels([f"LUM2_{r}" for r in REPS], fontsize=6.8)
    ax.set_ylabel("human proteins accepted")
    ax.set_title(title, pad=4)
    ax.set_ylim(0, max(data[e][r][thr] for e in data for r in REPS) * 1.24)
    return {e: [data[e][r][thr] for r in REPS] for e in data}


def main():
    data = measure()
    fig, axes = plt.subplots(1, 2, figsize=(style.W2, 2.8))
    a = panel(axes[0], data, 0.01, "1% entrapment FDP")
    b = panel(axes[1], data, 0.05, "5% entrapment FDP")
    for ax, letter in zip(axes, "ab"):
        style.panel_label(ax, letter, x=-0.14)
    fig.tight_layout(w_pad=2.2)
    style.save(fig, "figS2_faims", tight=False)

    for name, d in (("@1%", a), ("@5%", b)):
        for e, v in d.items():
            print(f"{name} {e:12s} " + " ".join(f"{r}={n}" for r, n in zip(REPS, v))
                  + f"  mean={np.mean(v):.0f}")


if __name__ == "__main__":
    main()
