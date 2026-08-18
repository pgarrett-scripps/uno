#!/usr/bin/env python3
"""Fig S4 -- the biosaur2 per-feature quality columns are a non-lever.

Panels:
  a  ablation: protein depth on run10 with vs. without the five quality columns
     (cosine_score, isotope_score, combined_score, ppm_error, neutron_offset)
     stripped from the SAME feature file, for both engines, at 1% and 5%
     entrapment FDR
  b  how much target/decoy signal each column actually carries, as |AUC - 0.5|
     (0 = none), next to the two things that do separate: retention time alone
     and the full ML discriminant

DATA: paper/figdata/quality_columns.tsv -- see that file's header for
per-row provenance. Its values are transcribed from benchmark/RESULTS.md
(section "Feature analysis"); the depth rows that are still recomputable from
protein tables on disk were recomputed with benchmark/entrapment.py and agree.
The per-column AUC diagnostics were throwaway scripts removed after use, so
RESULTS.md is their only surviving record.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import style  # noqa: E402

TSV = style.FIGDATA / "quality_columns.tsv"
ENGINE_LABEL = {"uno": "uno", "ms1": "ms1searchpy"}
ENGINE_COLOR = {"uno": style.C["uno"], "ms1": style.C["ms1searchpy"]}

# panel b, bottom-to-top: the five quality columns, then the real discriminants
BROWS = [
    ("auc", "neutron_offset"),
    ("auc", "ppm_error"),
    ("auc", "cosine_score"),
    ("auc", "isotope_score"),
    ("auc", "combined_score"),
    ("separator", "RT (single feature)"),
    ("separator", "full ML classifier"),
]


def main() -> int:
    df = style.read_tsv(TSV)

    def val(block, item, engine, metric):
        s = df[(df.block == block) & (df["item"] == item)
               & (df.engine == engine) & (df.metric == metric)]["value"]
        if len(s) != 1:
            raise SystemExit(f"expected 1 row for {block}/{item}/{engine}/{metric}, got {len(s)}")
        return float(s.iloc[0])

    fig, axes = plt.subplots(1, 2, figsize=(style.W2, 2.9))

    # -- a: ablation ---------------------------------------------------------
    ax = axes[0]
    metrics = [("depth_1pct", "1% entrapment FDR"), ("depth_5pct", "5% entrapment FDR")]
    series = [("uno", "with", None), ("uno", "without", "///"),
              ("ms1", "with", None), ("ms1", "without", "///")]
    x = np.arange(len(metrics), dtype=float)
    w = 0.2
    for k, (eng, item, hatch) in enumerate(series):
        vals = [val("ablation", item, eng, m) for m, _ in metrics]
        b = ax.bar(x + (k - 1.5) * w, vals, w * 0.92, color=ENGINE_COLOR[eng],
                   edgecolor="white", hatch=hatch, linewidth=0.5, zorder=2)
        ax.bar_label(b, fmt="%d", fontsize=6.5, padding=2, rotation=90)
    # legend: colour = engine, hatch = ablation arm
    for eng in ("uno", "ms1"):
        ax.bar(np.nan, np.nan, color=ENGINE_COLOR[eng], label=ENGINE_LABEL[eng])
    for item, hatch in (("with quality columns", None), ("without quality columns", "///")):
        ax.bar(np.nan, np.nan, color=style.C["light"], edgecolor="white", hatch=hatch,
               label=item)
    ax.set_xticks(x)
    ax.set_xticklabels([lab for _, lab in metrics], fontsize=7)
    ax.set_ylabel("proteins (run10, entrapment FDR)")
    ax.set_ylim(0, 2750)
    ax.set_xlim(-0.55, 1.55)
    ax.legend(loc="upper left", fontsize=6.5, ncol=2, columnspacing=1.0)
    style.panel_label(ax, "a", x=-0.19)

    # -- b: per-column target/decoy signal -----------------------------------
    ax = axes[1]
    y = np.arange(len(BROWS), dtype=float)
    h = 0.36
    for k, eng in enumerate(("uno", "ms1")):
        vals = []
        for block, item in BROWS:
            if block == "auc":
                vals.append(abs(val(block, item, eng, "auc") - 0.5))
            else:
                vals.append(val(block, item, eng, "auc_dist"))
        ax.barh(y + (0.5 - k) * h, vals, h, color=ENGINE_COLOR[eng], edgecolor="none",
                label=ENGINE_LABEL[eng], zorder=2)
    ax.axhline(4.5, color=style.C["neutral"], lw=0.7, ls=":", zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels([item for _, item in BROWS], fontsize=7)
    ax.set_xlabel("target/decoy separation, |AUC − 0.5|")
    ax.set_xlim(0, 0.20)
    ax.set_ylim(-0.6, len(BROWS) - 0.4)
    ax.text(0.198, 4.62, "actually discriminating", fontsize=6.5, style="italic",
            ha="right", va="bottom", color=style.C["neutral"])
    ax.text(0.198, 4.38, "biosaur2 quality columns", fontsize=6.5, style="italic",
            ha="right", va="top", color=style.C["neutral"])
    ax.legend(loc="lower right", fontsize=7)
    style.panel_label(ax, "b", x=-0.36)

    style.save(fig, "figS4_quality_columns")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
