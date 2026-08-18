#!/usr/bin/env python3
"""Fig 2 -- uno vs ms1searchpy across the six DirectMS1 benchmark datasets.

Panels (per FIGURE_MANIFEST.md):
  a  protein depth at 1% entrapment FDR
  b  protein depth at 5% entrapment FDR
  c  ms1searchpy's replicate range as a percentage of its own median,
     (max-min)/median; uno's is 0 by construction (identical reruns are
     byte-identical), drawn as a zero baseline
  d  wall-clock (s)
  e  peak resident set size (GiB)

uno is deterministic -> one bar per dataset. ms1searchpy is not -> all three
replicates are plotted as open circles with a solid tick at their median, so the
run-to-run spread is visible rather than hidden behind a summary statistic.
Panel c is the scale-free form of that spread; the per-dataset absolute view and
the run10 close-up live in the SI (figS6_determinism.py).

DATA: benchmark/results/summary.tsv, written by benchmark/run.sh
      (uno --features ... --fasta entrap_human_yeast.fasta --rt-model
       rt_cnn_chrono_global_v3_ch32.urtc --config uno.toml   vs
       ms1searchpy <features> -d entrap_human_yeast.fasta -deeplc 1 -ad 1),
      depths scored by benchmark/entrapment.py (human+yeast entrapment FDP).
No value in this figure is computed anywhere but that TSV.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import style  # noqa: E402

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


def bars_and_reps(ax, df, datasets, col, ylab, scale, letter, first,
                  rotate=False):
    """One bar-vs-replicates panel (a, b, d, e)."""
    x = np.arange(len(datasets), dtype=float)
    bw = 0.34
    uno = np.array([
        df[(df.dataset == d) & (df.tool == "uno")][col].astype(float).iloc[0] * scale
        for d in datasets
    ])
    ax.bar(x - bw / 2, uno, bw, color=style.C["uno"], edgecolor="none",
           label="uno (deterministic)" if first else None, zorder=2)
    for i, d in enumerate(datasets):
        reps = df[(df.dataset == d) & (df.tool == "ms1")][col].astype(float) * scale
        if len(reps) != 3:
            raise SystemExit(f"{d}: expected 3 ms1 replicates, got {len(reps)}")
        style.spread_marker(
            ax, x[i] + bw / 2, reps.to_numpy(), style.C["ms1searchpy"],
            half_width=bw / 2,
            label="ms1searchpy replicates" if (first and i == 0) else None,
            median_label="ms1searchpy median" if (first and i == 0) else None,
        )
    ax.set_xticks(x)
    if rotate:
        ax.set_xticklabels([LABEL[d].replace("\n", " ") for d in datasets],
                           fontsize=6, rotation=45, ha="right")
    else:
        ax.set_xticklabels([LABEL[d] for d in datasets], fontsize=6.5)
    ax.set_ylabel(ylab)
    ax.set_xlim(-0.6, len(datasets) - 0.4)
    ax.margins(y=0.14)
    style.panel_label(ax, letter, x=-0.13)


def spread_panel(ax, df, datasets, letter):
    """Panel c: replicate range as % of median, 1% and 5% FDR side by side."""
    x = np.arange(len(datasets), dtype=float)
    # FDR threshold is encoded by hatch; the colour stays the engine's, so
    # nothing here can be mistaken for a second tool.
    for col, lab, off, hatch in [
        ("depth_1pct", "@ 1% FDR", -0.19, None),
        ("depth_5pct", "@ 5% FDR", 0.19, "///"),
    ]:
        vals = []
        for d in datasets:
            r = df[(df.dataset == d) & (df.tool == "ms1")][col].astype(float).to_numpy()
            vals.append(100.0 * (r.max() - r.min()) / np.median(r))
        ax.bar(x + off, vals, 0.34, color=style.C["ms1searchpy"], edgecolor="white",
               hatch=hatch, linewidth=0.4, label=f"ms1searchpy {lab}", zorder=2)
    ax.axhline(0, color=style.C["uno"], lw=2.0, zorder=3,
               label="uno: 0% (byte-identical)")
    ax.set_xticks(x)
    ax.set_xticklabels([LABEL[d].replace("\n", " ") for d in datasets],
                       fontsize=6, rotation=45, ha="right")
    ax.set_ylabel("replicate range,\n(max − min) / median (%)")
    ax.set_xlim(-0.6, len(datasets) - 0.4)
    ax.set_ylim(0, None)
    ax.margins(y=0.20)
    ax.legend(loc="upper left", fontsize=6.3)
    style.panel_label(ax, letter, x=-0.13)


def main() -> int:
    df = pd.read_csv(SUMMARY, sep="\t")
    datasets = [d for d in LABEL if d in set(df["dataset"])]
    missing = sorted(set(df["dataset"]) - set(LABEL))
    if missing:
        raise SystemExit(f"summary.tsv has datasets not in LABEL: {missing}")

    fig = plt.figure(figsize=(style.W2, 5.4))
    gs = fig.add_gridspec(2, 6)
    ax_a = fig.add_subplot(gs[0, :3])
    ax_b = fig.add_subplot(gs[0, 3:])
    ax_c = fig.add_subplot(gs[1, :2])
    ax_d = fig.add_subplot(gs[1, 2:4])
    ax_e = fig.add_subplot(gs[1, 4:])

    bars_and_reps(ax_a, df, datasets, "depth_1pct",
                  "proteins @ 1% entrapment FDR", 1.0, "a", first=True)
    bars_and_reps(ax_b, df, datasets, "depth_5pct",
                  "proteins @ 5% entrapment FDR", 1.0, "b", first=False)
    spread_panel(ax_c, df, datasets, "c")
    bars_and_reps(ax_d, df, datasets, "wall_s", "wall-clock (s)", 1.0, "d",
                  first=False, rotate=True)
    bars_and_reps(ax_e, df, datasets, "peak_rss_mb", "peak RSS (GiB)", 1 / 1024,
                  "e", first=False, rotate=True)

    ax_a.legend(loc="upper left", ncol=1, fontsize=7,
                bbox_to_anchor=(0.0, 1.02))
    for ax in (ax_d, ax_e):
        ax.set_ylim(bottom=0)
        ax.set_xlabel("dataset")
    ax_c.set_xlabel("dataset")

    style.save(fig, "fig2_headtohead")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
