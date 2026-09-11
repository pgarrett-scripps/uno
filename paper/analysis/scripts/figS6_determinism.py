#!/usr/bin/env python3
"""Fig S6 -- randomization spread, measured the same way for both engines.

WHAT CHANGED AND WHY (2026-09-09). This figure used to plot ms1searchpy's three
replicates against a single uno point labelled "deterministic". That comparison
was not fair: uno's byte-identical output came from hardcoded constants (booster
seed 42, a fixed RANSAC LCG constant, a fixed decoy-shuffle salt, a fixed CV
fold hash), not from an algorithm with less randomness in it. Locking those
constants is something either engine can do, so a locked run against an unlocked
one measures the locking, not the engine.

uno now exposes `--seed`, which perturbs all four stages together (seed 0 = the
shipped constants). Both engines are therefore drawn the same way here: one open
marker per independent random draw, a connector spanning min..max, and a tick at
the median. uno's draws are seeds 0..N; ms1searchpy's are reruns.

Panels:
  a  depth @1% entrapment FDR, per dataset: uno's seed draws vs ms1searchpy's
     reruns, both as spreads
  b  same at 5% entrapment FDR
  c  where the previously published single uno run (seed 0) falls inside uno's
     own seed distribution, as a percentile per dataset -- the diagnostic for
     whether a reported single number was a lucky draw

DATA
  benchmark/results/spread.tsv, written by benchmark/spread_table.py, which
  scores every available uno seed run (benchmark/seed_sweep.sh) and every
  ms1searchpy replicate directory with benchmark/entrapment.py. Regenerate it
  before rebuilding this figure; this script does not re-score protein tables.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import style  # noqa: E402
import seed_data
from _assets import record

SPREAD = style.BENCH / "results" / "spread.tsv"

LABEL = {
    "hela5min_09746": "09746",
    "hela5min_09766": "09766",
    "hela5min_09734": "09734",
    "hela5min_09778": "09778",
    "hela5min_09950": "09950",
    "faims_LUM2_05433_VG": "05433\n(FAIMS)",
    "run10": "run10",
}

# A spread is only a spread with more than one draw; a dataset carrying a single
# run of an engine is plotted as a lone marker and excluded from panel c.
MIN_DRAWS = 3


def load() -> pd.DataFrame:
    if not SPREAD.exists():
        raise SystemExit(
            f"missing {SPREAD}\nrun: python3 benchmark/spread_table.py")
    df = pd.DataFrame(seed_data.rows())
    for col in ("draw", "depth_1pct", "depth_5pct"):
        df[col] = df[col].astype(int)
    missing = set(df["dataset"]) - set(LABEL)
    if missing:
        raise SystemExit(f"unlabelled datasets in {SPREAD}: {sorted(missing)}")
    return df


def draws(df: pd.DataFrame, dataset: str, engine: str, col: str) -> np.ndarray:
    sel = df[(df.dataset == dataset) & (df.engine == engine)][col]
    return np.asarray(sel, dtype=float)


def main() -> int:
    df = load()
    datasets = [d for d in LABEL if d in set(df["dataset"])]
    x = np.arange(len(datasets), dtype=float)

    fig, axes = plt.subplots(1, 3, figsize=(style.W2, 2.9))

    # -- a, b: both engines as spreads over independent draws -----------------
    for pi, (col, ylab) in enumerate([
        ("depth_1pct", "proteins @ 1% entrapment FDR"),
        ("depth_5pct", "proteins @ 5% entrapment FDR"),
    ]):
        ax = axes[pi]
        for i, d in enumerate(datasets):
            for engine, off, first_label, first_med in [
                ("uno", -0.16, "uno seeds", "uno median"),
                ("ms1searchpy", 0.16, "ms1searchpy reruns", "ms1searchpy median"),
            ]:
                v = draws(df, d, engine, col)
                if v.size == 0:
                    continue
                lab = first_label if (pi == 0 and i == 0) else None
                mlab = first_med if (pi == 0 and i == 0) else None
                if v.size == 1:
                    ax.plot(x[i] + off, v[0], ls="none", marker="o", mfc="white",
                            mec=style.C[engine], mew=0.9, ms=3.4, zorder=4,
                            label=lab)
                else:
                    style.spread_marker(ax, x[i] + off, v, style.C[engine],
                                        half_width=0.12, jitter=0.022, label=lab,
                                        median_label=mlab)
        ax.set_xticks(x)
        ax.set_xticklabels([LABEL[d].replace("\n", " ") for d in datasets],
                           fontsize=6, rotation=45, ha="right")
        ax.set_ylabel(ylab)
        ax.set_xlim(-0.6, len(datasets) - 0.4)
        ax.margins(y=0.18)
        style.panel_label(ax, "abc"[pi], x=-0.19)
    axes[0].legend(loc="upper left", fontsize=5.6)

    # -- c: is the published seed-0 run a typical draw? -----------------------
    # Percentile of the seed-0 depth within uno's own seed distribution. 50 is a
    # median draw; near 100 means the published single number sat at the top of
    # the distribution and overstated what a fresh seed would give.
    ax = axes[2]
    labels, pcts = [], []
    for d in datasets:
        v = draws(df, d, "uno", "depth_1pct")
        seed0 = df[(df.dataset == d) & (df.engine == "uno") & (df.draw == 0)]
        if v.size < MIN_DRAWS or seed0.empty:
            continue
        s0 = float(seed0["depth_1pct"].iloc[0])
        labels.append(LABEL[d].replace("\n", " "))
        pcts.append(100.0 * float((v < s0).sum()) / (v.size - 1))
    if labels:
        y = np.arange(len(labels), dtype=float)
        ax.barh(y, pcts, 0.6, color=style.C["uno"], edgecolor="white",
                linewidth=0.4, zorder=2)
        ax.axvline(50, color=style.C["neutral"], lw=0.9, ls="--", zorder=3)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=6.5)
        ax.set_xlim(0, 100)
        ax.set_xlabel("other sampled seeds below\nthe default run (%)",
                      fontsize=6.5)
        ax.invert_yaxis()
        for yi, p in zip(y, pcts):
            ax.text(p - 2 if p >= 80 else p + 2, yi, f"{p:.0f}", va="center",
                    ha="right" if p >= 80 else "left", fontsize=5.6,
                    color="white" if p >= 80 else "black")
    else:
        ax.text(0.5, 0.5, "no dataset has enough uno seeds yet",
                ha="center", va="center", transform=ax.transAxes, fontsize=7)
        ax.set_xticks([])
        ax.set_yticks([])
    style.panel_label(ax, "c", x=-0.19)

    style.save(fig, "figS6_determinism")
    for suffix, ext in (("", "png"), ("-pdf", "pdf")):
        record(f"fig.figs6-determinism{suffix}", f"figures/figS6_determinism.{ext}",
               kind="figure", inputs=["../benchmark/results/spread.tsv"],
               desc="Primary seed/rerun depth distributions and default-seed ranks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
