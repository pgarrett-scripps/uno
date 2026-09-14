#!/usr/bin/env python3
"""Fig 2 -- uno vs ms1searchpy across the six DirectMS1 benchmark datasets.

Panels (per FIGURE_MANIFEST.md):
  a  protein depth at 1% entrapment FDR
  b  wall-clock (s) versus peak resident set size (GiB), one point per search,
     the same resource view as panel b of fig_astral_validation.py

Both engines are shown as individual random draws, full ranges, and medians.
Depths use benchmark/results/spread.tsv; resources retain the original
serialized searches in benchmark/results/summary.tsv. Concurrent seed-sweep
timings are deliberately excluded from the resource comparison.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import style  # noqa: E402
import seed_data  # noqa: E402
from _assets import record  # noqa: E402

SUMMARY = style.BENCH / "results" / "summary.tsv"

# dataset key -> two-line axis label. Order = order of appearance in summary.tsv
# (five 5-min HeLa runs from PXD015669, then the FAIMS run from PXD023977).
LABEL = {
    "hela5min_09746": "09746",
    "hela5min_09766": "09766",
    "hela5min_09734": "09734",
    "hela5min_09778": "09778",
    "hela5min_09950": "09950",
    "faims_LUM2_05433_VG": "05433\n(FAIMS)",
}


def depth_panel(ax, datasets):
    rows = seed_data.rows()
    x = np.arange(len(datasets), dtype=float)
    for i, ds in enumerate(datasets):
        for engine, offset in (("uno", -0.19), ("ms1searchpy", 0.19)):
            values = [int(r["depth_1pct"]) for r in rows
                      if r["dataset"] == ds and r["engine"] == engine]
            style.spread_marker(ax, x[i] + offset, values, style.C[engine],
                                jitter=0.022, half_width=0.13, ms=3,
                                label=engine if i == 0 else None)
            ax.text(x[i] + offset, max(values) + 85, f"n={len(values)}",
                    fontsize=5.7, ha="center", color=style.C[engine])
    ax.set_xticks(x, [LABEL[d] for d in datasets], fontsize=6.5)
    ax.set_ylabel("proteins at 1% entrapment FDR")
    ax.set_xlim(-0.65, len(datasets) - 0.35)
    ax.set_ylim(bottom=0)
    ax.margins(y=0.17)
    style.panel_label(ax, "a", x=-0.16)


def resource_panel(ax, df, datasets, letter):
    """Panel b: wall-clock vs peak RSS, one marker per search (all replicates)."""
    for tool, engine in [("uno", "uno"), ("ms1", "ms1searchpy")]:
        block = df[(df.tool == tool) & (df.dataset.isin(datasets))]
        color, marker = style.ENGINE[engine]
        ax.scatter(block.wall_s.astype(float),
                   block.peak_rss_mb.astype(float) / 1024,
                   s=18 if tool == "uno" else 12, marker=marker,
                   facecolors=color if tool == "uno" else "white",
                   edgecolors=color, linewidths=0.8, zorder=3,
                   label="uno" if tool == "uno" else "ms1searchpy (3 runs each)")
    ax.set_xlabel("wall time (s)")
    ax.set_ylabel("peak resident memory (GiB)")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.legend(loc="lower right")
    style.panel_label(ax, letter, x=-0.15)


def main() -> int:
    df = pd.read_csv(SUMMARY, sep="\t")
    datasets = [d for d in LABEL if d in set(df["dataset"])]
    missing = sorted(set(df["dataset"]) - set(LABEL))
    if missing:
        raise SystemExit(f"summary.tsv has datasets not in LABEL: {missing}")

    fig, (ax_a, ax_b) = plt.subplots(
        1, 2, figsize=(style.W2, 2.7), gridspec_kw={"width_ratios": [1.45, 1]})

    depth_panel(ax_a, datasets)
    ax_a.set_xlabel("dataset")
    ax_a.legend(loc="upper left", ncol=1, fontsize=7,
                bbox_to_anchor=(0.0, 1.02))
    resource_panel(ax_b, df, datasets, "b")

    style.save(fig, "fig2_headtohead")
    record("fig.fig2-headtohead", "figures/fig2_headtohead.png", kind="figure",
           inputs=["../benchmark/results/summary.tsv", "../benchmark/results/spread.tsv"],
           desc="uno versus ms1searchpy on six DirectMS1 datasets: depth at 1% "
                "entrapment FDR and wall-clock versus peak memory")
    record("fig.fig2-headtohead-pdf", "figures/fig2_headtohead.pdf", kind="figure",
           inputs=["../benchmark/results/summary.tsv", "../benchmark/results/spread.tsv"],
           desc="Vector twin of the head-to-head figure")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
