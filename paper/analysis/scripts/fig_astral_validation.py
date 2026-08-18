#!/usr/bin/env python3
"""External Orbitrap Astral validation on three independent 10 ng HeLa runs.

uno is deterministic and was run once per raw file. ms1searchpy is stochastic;
all ten serialized searches per file are shown. Identification depth is scored
with the same human/yeast entrapment procedure used throughout the manuscript.

DATA: ~/Data/PXD066701_Astral/results/astral_summary.tsv, generated from primary
protein tables by analysis/astral/summarize.py.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import style  # noqa: E402

SUMMARY = Path.home() / "Data/PXD066701_Astral/results/astral_summary.tsv"
STEMS = [
    "CW0710_OldFAIMS_10ng__RA4_rep1_ch1_DIA10ng",
    "CW0710_OldFAIMS_10ng__RA4_rep4_ch2_DIA10ng",
    "CW0710_OldFAIMS_10ng__RA4_rep6_ch2_DIA10ng",
]
LABELS = ["Astral 1", "Astral 2", "Astral 3"]


def wall_seconds(value: str) -> float:
    fields = [float(x) for x in str(value).split(":")]
    return fields[0] * 60 + fields[1] if len(fields) == 2 else fields[0] * 3600 + fields[1] * 60 + fields[2]


def depth_panel(ax, df, column: str, threshold: str, letter: str) -> None:
    x = np.arange(3, dtype=float)
    width = 0.34
    uno = [float(df[(df.dataset == s) & (df.engine == "uno")][column].iloc[0]) for s in STEMS]
    ax.bar(x - width / 2, uno, width, color=style.C["uno"], edgecolor="none",
           label="uno" if letter == "a" else None)
    for i, stem in enumerate(STEMS):
        reps = df[(df.dataset == stem) & (df.engine == "ms1searchpy")][column].astype(float)
        if len(reps) != 10:
            raise SystemExit(f"{stem}: expected 10 ms1searchpy runs, found {len(reps)}")
        med = style.spread_marker(
            ax, x[i] + width / 2, reps, style.C["ms1searchpy"], jitter=0.018,
            half_width=width / 2, ms=2.8,
            label="ms1searchpy runs" if (letter == "a" and i == 0) else None,
            median_label="ms1searchpy median" if (letter == "a" and i == 0) else None,
        )
        ax.text(x[i] - width / 2, uno[i] + 35, f"{uno[i]:.0f}", ha="center", va="bottom",
                fontsize=6.5, color=style.C["uno"])
        ax.text(x[i] + width / 2, max(reps.max(), med) + 35, f"{med:.0f}", ha="center",
                va="bottom", fontsize=6.5, color=style.C["ms1searchpy"])
    ax.set_xticks(x, LABELS)
    ax.set_ylabel(f"proteins @ {threshold} entrapment FDR")
    ax.set_ylim(0, max(uno) * 1.15)
    ax.set_xlim(-0.55, 2.55)
    style.panel_label(ax, letter, x=-0.15)


def resource_panel(ax, df) -> None:
    uno = df[df.engine == "uno"].copy()
    ms1 = df[df.engine == "ms1searchpy"].copy()
    for block, engine, label in [(uno, "uno", "uno"), (ms1, "ms1searchpy", "ms1searchpy")]:
        color, marker = style.ENGINE[engine]
        ax.scatter([wall_seconds(x) for x in block.wall], block.rss_kib.astype(float) / 1048576,
                   s=18 if engine == "uno" else 12, marker=marker,
                   facecolors=color if engine == "uno" else "white", edgecolors=color,
                   linewidths=0.8, label=label, zorder=3)
    ax.set_xlabel("wall-clock (s)")
    ax.set_ylabel("peak RSS (GiB)")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.legend(loc="lower right")
    style.panel_label(ax, "c", x=-0.15)


def main() -> int:
    df = pd.read_csv(SUMMARY, sep="\t")
    expected = {(s, "uno") for s in STEMS} | {(s, "ms1searchpy") for s in STEMS}
    if not expected <= set(zip(df.dataset, df.engine)):
        raise SystemExit("Astral summary is incomplete")

    fig, axes = plt.subplots(1, 3, figsize=(style.W2, 2.65), gridspec_kw={"width_ratios": [1, 1, 0.86]})
    depth_panel(axes[0], df, "depth_1pct_entrapment_q", "1%", "a")
    depth_panel(axes[1], df, "depth_5pct_entrapment_q", "5%", "b")
    resource_panel(axes[2], df)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.34, 1.015),
               ncol=3, fontsize=6.7)
    style.save(fig, "fig_astral_validation")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
