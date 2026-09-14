#!/usr/bin/env python3
"""Independent Orbitrap Astral validation across sample input (PXD066701).

Twelve 11-min dual-column, fixed-CV HeLa files at 200 pg, 1 ng, 10 ng and
100 ng (three per level, the first three by file name). uno is measured over
three seeds per file. Self-calibrated ms1searchpy searches were run three
times per file (ten at 10 ng) and, because its RT self-calibration failed on
every file, three more times per file with DeepLC calibrated on an external
MS/MS peptide list (`-pl`, analysis/astral/RT_CALIBRATION.md). Depth is the
common human/yeast entrapment scorer.

Panels: (a) depth at 1% entrapment FDR versus input; (b) wall-clock versus
peak RSS for every serialized search (uno and self-calibrated
ms1searchpy at all inputs; the 10 ng `-pl` runs were concurrent and are
excluded, the other `-pl` runs were serialized and are shown).

DATA: ~/Data/PXD066701_Astral/results/input_series_summary.tsv from
analysis/astral/summarize_input_series.py.
"""
from __future__ import annotations
import os, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import style  # noqa: E402
from _assets import record  # noqa: E402

SUMMARY = Path.home() / "Data/PXD066701_Astral/results/input_series_summary.tsv"
EXTENSION = style.BENCH / "results/variability_extension/summary.tsv"
LEVELS = [0.2, 1.0, 10.0, 100.0]
TICKS = ["200 pg", "1 ng", "10 ng", "100 ng"]
PL_COLOR = style.C["accent2"]
SERIES = [("uno", style.C["uno"], "uno"),
          ("ms1searchpy", style.C["ms1searchpy"], "ms1searchpy, self-calibrated"),
          ("ms1searchpy-pl", PL_COLOR, "ms1searchpy, external RT calibration")]
OFFS = {"uno": 0.62, "ms1searchpy": 1.0, "ms1searchpy-pl": 1.6}  # multiplicative x offsets


def wall_seconds(v: str) -> float:
    f = [float(x) for x in str(v).split(":")]
    return f[0] * 60 + f[1] if len(f) == 2 else f[0] * 3600 + f[1] * 60 + f[2]


def depth_panel(ax, df, col, threshold, letter):
    for eng, color, label in SERIES:
        meds = []
        for i, lv in enumerate(LEVELS):
            block = df[(df.engine == eng) & (df.input_ng == lv)]
            files = sorted(block.dataset.unique())
            if len(files) != 3:
                raise ValueError(f"need three acquisition files at {lv}: {eng}")
            file_medians = []
            for j, dataset in enumerate(files):
                v = block[block.dataset == dataset][col].astype(float).values
                if len(v) < 3:
                    raise ValueError(f"need at least three searches: {dataset}, {eng}")
                x = lv * OFFS[eng] * [0.88, 1.0, 1.12][j]
                ax.plot([x, x], [v.min(), v.max()], color=color, lw=0.8, zorder=2)
                offsets = np.linspace(-0.025, 0.025, len(v))
                ax.plot(x * (1 + offsets), np.sort(v), ls="none", marker="o", ms=2.4,
                        mfc="white", mec=color, mew=0.7, zorder=4,
                        label=label if i == 0 and j == 0 else None)
                med = float(np.median(v))
                ax.plot([x / 1.05, x * 1.05], [med, med], color=color, lw=1.5, zorder=5)
                file_medians.append(med)
            meds.append(float(np.median(file_medians)))
        ax.plot([lv * OFFS[eng] for lv in LEVELS], meds, color=color, lw=1.2, zorder=3,
                label=("median" if (eng == "uno" and letter == "a") else None))
    ax.set_xscale("log")
    ax.set_xticks(LEVELS, TICKS)
    ax.set_xlim(0.08, 220)
    ax.set_ylabel(f"proteins at {threshold} entrapment FDR")
    ax.set_xlabel("HeLa input")
    ax.set_ylim(bottom=-45)
    ax.set_yticks([v for v in ax.get_yticks() if v >= 0])
    style.panel_label(ax, letter, x=-0.18)


def resource_panel(ax, df):
    for eng, color, label in SERIES:
        block = df[df.engine == eng]
        if eng == "ms1searchpy-pl":
            block = block[block.input_ng != 10.0]  # those ran concurrently
        marker = {"uno": "o", "ms1searchpy": "s", "ms1searchpy-pl": "D"}[eng]
        ax.scatter([wall_seconds(w) for w in block.wall], block.rss_kib.astype(float) / 1048576,
                   s=16 if eng == "uno" else 11, marker=marker,
                   facecolors=color if eng == "uno" else "white", edgecolors=color,
                   linewidths=0.8, zorder=3)
    ax.set_xlabel("wall time (s)")
    ax.set_ylabel("peak resident memory (GiB)")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    style.panel_label(ax, "b", x=-0.18)


def main() -> int:
    df = pd.read_csv(SUMMARY, sep="\t")
    df = df[df.engine.isin([s[0] for s in SERIES])]
    fig, axes = plt.subplots(1, 2, figsize=(style.W2, 2.5), gridspec_kw={"width_ratios": [1.3, 1]})
    extra = pd.read_csv(EXTENSION, sep="\t")
    extra = extra[(extra.arm == "biosaur2") & extra.astral_stem.notna()].copy()
    extra["dataset"] = extra.astral_stem
    extra["engine"] = "uno"
    extra["input_ng"] = extra.dataset.str.extract(r"OldFAIMS_(\d+(?:ng|pg))__")[0].map(
        {"200pg": 0.2, "1ng": 1.0, "10ng": 10.0, "100ng": 100.0})
    depth_df = pd.concat([df[df.engine != "uno"], extra], ignore_index=True)
    depth_panel(axes[0], depth_df, "depth_1pct", "1%", "a")
    resource_panel(axes[1], df)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 1.06), ncol=4, fontsize=6.4)
    fig.subplots_adjust(wspace=0.35)
    style.save(fig, "fig_astral_validation")
    inputs = [os.path.relpath(SUMMARY, style.PAPER), os.path.relpath(EXTENSION, style.PAPER)]
    record("fig.fig-astral-validation", "figures/fig_astral_validation.png", kind="figure", inputs=inputs,
           desc="Independent Orbitrap Astral validation across input: depth at 1% entrapment FDR "
                "(uno, self-calibrated ms1searchpy, ms1searchpy with external RT calibration) and wall-clock versus peak memory")
    record("fig.fig-astral-validation-pdf", "figures/fig_astral_validation.pdf", kind="figure", inputs=inputs,
           desc="Vector twin of the Astral validation figure")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
