#!/usr/bin/env python3
"""TOC / abstract graphic for the JPR submission (SUBMISSION_TODO B7).

ACS "For Table of Contents Only" graphic: the artwork must fit 3.25 in wide by
1.75 in tall (the canvas is exactly that size; savefig cropping is disabled so
the submitted image keeps those dimensions) at 300 dpi.

Layout, colours and type sizes reproduce `paper/toc_graphic.pptx`, the slide the
graphic was designed on. That slide is drawn at 2x the ACS box (6.5 in x 3.5 in)
with type at 2x printed size, so this script works in *slide* units throughout
and halves every length and point size on the way out; SLIDE_SCALE is the only
conversion. Keeping slide coordinates makes the two artefacts diffable.

Left column is schematic and carries no data: an MS1 feature map for one run
(the feature blobs are the slide's own hand-placed positions, not a random
scatter), feeding a binomial protein score with an illustrative 1% FDR cutoff.
Three feature colours mark the three scored target proteins.

The right-hand bar chart is the only part with numbers in it, and they are the
manuscript's: per-dataset median depth at 1% entrapment FDR over every search
draw in benchmark/results/spread.tsv -- the same medians Figure 2a reports. The
slide itself carried uno's single default-seed depth against ms1searchpy's
3-replicate median, which overstated uno on 09746 and 09766; this script reads
the file so the TOC graphic cannot drift from the figure it advertises.
"""

from __future__ import annotations

import statistics
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

import style  # noqa: E402
from _assets import record  # noqa: E402

# a TOC graphic is a fixed-size canvas; cropping it breaks the ACS dimensions
plt.rcParams["savefig.bbox"] = "standard"

# ---------------------------------------------------------------------------
# canvas: work in the slide's inches, emit at half size
# ---------------------------------------------------------------------------
SLIDE_W, SLIDE_H = 6.5, 3.5
SLIDE_SCALE = 0.5
W_IN, H_IN = SLIDE_W * SLIDE_SCALE, SLIDE_H * SLIDE_SCALE  # 3.25 x 1.75

# palette of the slide (deliberately not style.C: the TOC graphic is a cover
# image, and its uno/comparator pair is warm/neutral rather than Wong blue)
INK = "#1F2A37"
MUTED = "#5B6672"
CARD = "#F2F4F6"
UNO = "#D94F1E"
MS1 = "#9AA5B1"
AXIS = "#C7CFD8"
DOT = "#D3D9E0"
PROTEIN = ("#D94F1E", "#1B7F79", "#3B5BA5")
DECOY_EDGE = "#8F9AA6"
GRID = "#E6E9ED"
TICK = "#888888"


def pt(size: float) -> float:
    """Slide point size -> point size on the half-scale canvas."""
    return size * SLIDE_SCALE


def blob(ax, x, y, w, h, facecolor, edgecolor=None, lw=0.0, zorder=3):
    """A stadium-cornered rectangle, the slide's roundRect at max radius."""
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0,rounding_size={h / 2}",
        facecolor=facecolor, edgecolor=edgecolor or "none",
        linewidth=lw, zorder=zorder))


def card(ax, x, y, w, h):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0,rounding_size=0.06",
        facecolor=CARD, edgecolor="none", zorder=1))


# ---------------------------------------------------------------------------
# left panel: MS1 feature map (schematic; slide coordinates, x, y, w, h, colour)
# ---------------------------------------------------------------------------
FEATURES = [
    (1.728, 1.069, 0.113, 0.022, DOT), (1.945, 1.286, 0.080, 0.022, DOT),
    (2.348, 1.168, 0.056, 0.022, DOT), (1.948, 1.206, 0.090, 0.022, DOT),
    (1.554, 0.964, 0.072, 0.022, DOT), (2.306, 0.976, 0.098, 0.022, DOT),
    (2.352, 1.038, 0.091, 0.022, DOT), (1.115, 1.297, 0.089, 0.022, DOT),
    (1.350, 0.922, 0.077, 0.022, DOT), (1.247, 1.350, 0.066, 0.022, DOT),
    (1.698, 1.024, 0.058, 0.022, DOT), (1.394, 0.897, 0.087, 0.022, DOT),
    (1.220, 1.338, 0.085, 0.022, DOT), (1.429, 1.132, 0.102, 0.022, DOT),
    (1.623, 1.283, 0.113, 0.022, DOT), (0.907, 1.341, 0.125, 0.022, DOT),
    (1.579, 1.247, 0.111, 0.022, DOT), (0.951, 1.080, 0.071, 0.022, DOT),
    (1.326, 0.948, 0.108, 0.022, DOT), (0.929, 1.352, 0.082, 0.022, DOT),
    (1.622, 0.959, 0.058, 0.022, DOT), (1.180, 1.179, 0.072, 0.022, DOT),
    (1.763, 1.271, 0.111, 0.022, DOT), (1.937, 0.875, 0.114, 0.022, DOT),
    (0.691, 0.852, 0.106, 0.022, DOT), (1.409, 1.004, 0.067, 0.022, DOT),
    (1.883, 0.880, 0.053, 0.022, DOT), (1.161, 1.097, 0.054, 0.022, DOT),
    (1.952, 0.957, 0.089, 0.022, DOT), (2.225, 1.043, 0.078, 0.022, DOT),
    (1.632, 1.071, 0.072, 0.022, DOT), (0.742, 1.264, 0.129, 0.022, DOT),
    (1.168, 0.895, 0.105, 0.022, DOT), (1.737, 0.941, 0.129, 0.022, DOT),
    (1.027, 1.079, 0.109, 0.022, DOT), (2.242, 0.909, 0.087, 0.022, DOT),
    (1.080, 1.051, 0.116, 0.022, DOT), (1.650, 0.962, 0.078, 0.022, DOT),
    (2.396, 1.121, 0.114, 0.022, DOT), (1.898, 0.857, 0.096, 0.022, DOT),
    (0.716, 0.961, 0.119, 0.022, DOT), (2.078, 1.084, 0.063, 0.022, DOT),
    (2.185, 1.091, 0.101, 0.022, DOT), (2.085, 1.044, 0.070, 0.022, DOT),
    (2.441, 1.345, 0.122, 0.022, DOT), (1.024, 0.854, 0.086, 0.022, DOT),
    (1.726, 1.069, 0.115, 0.022, DOT), (2.194, 1.088, 0.118, 0.022, DOT),
    (2.338, 0.853, 0.165, 0.032, PROTEIN[0]),
    (0.955, 0.893, 0.091, 0.032, PROTEIN[0]),
    (1.377, 1.339, 0.153, 0.032, PROTEIN[0]),
    (0.804, 1.159, 0.107, 0.032, PROTEIN[0]),
    (1.467, 1.331, 0.111, 0.032, PROTEIN[0]),
    (1.047, 0.899, 0.124, 0.032, PROTEIN[0]),
    (2.117, 1.226, 0.096, 0.032, PROTEIN[1]),
    (1.035, 0.900, 0.141, 0.032, PROTEIN[1]),
    (0.712, 0.874, 0.167, 0.032, PROTEIN[1]),
    (1.564, 1.083, 0.128, 0.032, PROTEIN[1]),
    (2.005, 1.030, 0.127, 0.032, PROTEIN[1]),
    (2.280, 1.067, 0.127, 0.032, PROTEIN[2]),
    (1.070, 1.190, 0.158, 0.032, PROTEIN[2]),
    (0.822, 1.019, 0.112, 0.032, PROTEIN[2]),
]

# illustrative protein scores: three targets above the cutoff, two decoys below
SCORE_ROWS = [
    (2.30, 1.4308, "12.4", PROTEIN[0], None, True),
    (2.45, 1.0500, "9.1", PROTEIN[1], None, True),
    (2.60, 0.7846, "6.8", PROTEIN[2], None, True),
    (2.75, 0.2423, "2.1", "#FFFFFF", "decoy", False),
    (2.90, 0.1385, "1.2", "#FFFFFF", "decoy", False),
]
SCORE_X = 0.44
CUTOFF_X = 0.9708

# ---------------------------------------------------------------------------
# right panel: the manuscript's depth medians
# ---------------------------------------------------------------------------
SPREAD = style.PAPER.parent / "benchmark" / "results" / "spread.tsv"
DATASETS = [
    ("hela5min_09746", "09746"), ("hela5min_09766", "09766"),
    ("hela5min_09734", "09734"), ("hela5min_09778", "09778"),
    ("hela5min_09950", "09950"), ("faims_LUM2_05433_VG", "05433"),
]
ENGINES = [("uno", "uno", UNO), ("ms1searchpy", "ms1searchpy", MS1)]


def median_depths() -> dict[tuple[str, str], float]:
    """Median depth at 1% entrapment FDR per (engine, dataset), from spread.tsv."""
    draws: dict[tuple[str, str], list[int]] = defaultdict(list)
    for row in style.read_tsv(SPREAD).itertuples(index=False):
        draws[(row.engine, row.dataset)].append(int(row.depth_1pct))
    return {k: statistics.median(v) for k, v in draws.items()}


HEADLINE = [("uno", UNO), (": fast MS1-only protein identification in Rust", INK)]
SERIF = ["Liberation Serif", "Nimbus Roman", "DejaVu Serif"]


def headline(fig, ax, x=0.28, y=0.30, right=6.23, size=pt(17)) -> None:
    """Two-colour headline, laid out run by run and shrunk to fit the box.

    Matplotlib has no rich text, so the runs are separate artists whose
    positions come from measuring the one before it; the slide's 17 pt is an
    upper bound, since the substituted serif is not metric-compatible with the
    Cambria the slide was set in.
    """
    for attempt in range(12):
        artists, cursor = [], x
        for text, colour in HEADLINE:
            t = ax.text(cursor, y, text, ha="left", va="center",
                        fontsize=size, fontweight="bold", color=colour,
                        family=SERIF)
            artists.append(t)
            box = t.get_window_extent(fig.canvas.get_renderer())
            cursor = ax.transData.inverted().transform(
                (box.x1, box.y1))[0]
        if cursor <= right or attempt == 11:
            return
        for t in artists:
            t.remove()
        size *= (right - x) / (cursor - x)


def feature_map(ax) -> None:
    card(ax, 0.28, 0.56, 2.45, 1.05)
    ax.text(1.505, 0.69, "MS1 feature map, one run", ha="center", va="center",
            fontsize=pt(11), fontweight="bold", color=INK, zorder=3)
    ax.plot([0.66, 0.66], [0.81, 1.45], color=AXIS, lw=pt(1.0),
            solid_capstyle="butt", zorder=2)
    ax.plot([0.66, 2.61], [1.45, 1.45], color=AXIS, lw=pt(1.0),
            solid_capstyle="butt", zorder=2)
    ax.text(0.31, 0.885, "m/z", ha="left", va="center", fontsize=pt(9),
            color=MUTED)
    ax.text(2.61, 1.535, "retention time", ha="right", va="center",
            fontsize=pt(9), color=MUTED)
    for x, y, w, h, colour in FEATURES:
        blob(ax, x, y, w, h, colour)


def score_panel(ax) -> None:
    card(ax, 0.28, 1.85, 2.45, 1.37)
    ax.text(1.505, 1.975, "binomial protein score", ha="center", va="center",
            fontsize=pt(11), fontweight="bold", color=INK, zorder=3)
    ax.text(1.505, 2.135, "−log₁₀ P(X ≥ k | n, p)",
            ha="center", va="center", fontsize=pt(9.5), color=MUTED, zorder=3)
    for top, width, label, colour, note, above in SCORE_ROWS:
        blob(ax, SCORE_X, top, width, 0.09, colour,
             edgecolor=None if above else DECOY_EDGE,
             lw=0.0 if above else pt(1.0))
        ax.text(2.63, top + 0.045, label, ha="right", va="center",
                fontsize=pt(8.5), fontweight="bold" if above else "normal",
                color=INK if above else MUTED, zorder=4)
        if note:
            ax.text(1.0608, top + 0.045, note, ha="left", va="center",
                    fontsize=pt(8.5), color=MUTED, zorder=4)
    ax.plot([CUTOFF_X, CUTOFF_X], [2.25, 3.03], color=INK, lw=pt(1.25),
            linestyle=(0, (3.2, 2.2)), zorder=5)
    ax.text(1.505, 3.125, "dashed line: 1% FDR cutoff", ha="center",
            va="center", fontsize=pt(8.5), color=MUTED)


def depth_chart(fig, depths) -> None:
    """Clustered columns in the slide's chart rectangle (3.05, 0.56, 3.17 x 2.18)."""
    # leave room inside that rectangle for the title, tick labels and legend
    left, top, width, height = 3.05, 0.56, 3.17, 2.18
    plot = (left + 0.36, top + 0.26, width - 0.42, height - 0.26 - 0.52)
    ax = fig.add_axes([
        plot[0] / SLIDE_W,
        1.0 - (plot[1] + plot[3]) / SLIDE_H,
        plot[2] / SLIDE_W,
        plot[3] / SLIDE_H,
    ])

    x = np.arange(len(DATASETS))
    bar_w = 1.0 / (len(ENGINES) + 0.55)  # OOXML gapWidth 55%, overlap 0
    for i, (key, label, colour) in enumerate(ENGINES):
        vals = [depths[(key, ds)] for ds, _ in DATASETS]
        ax.bar(x + (i - 0.5) * bar_w, vals, bar_w, label=label, color=colour,
               linewidth=0, zorder=3)

    ax.set_ylim(0, 2500)
    ax.set_yticks(range(0, 2501, 500))
    ax.set_xticks(x)
    ax.set_xticklabels([short for _, short in DATASETS])
    ax.set_xlim(-0.5, len(DATASETS) - 0.5)
    ax.set_title("proteins at 1% entrapment FDR", fontsize=pt(12), color=INK,
                 pad=pt(7))
    ax.yaxis.grid(True, color=GRID, lw=pt(0.75), zorder=0)
    ax.set_axisbelow(True)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(TICK)
        ax.spines[side].set_linewidth(pt(1.0))
    ax.tick_params(colors=TICK, labelcolor=MUTED, labelsize=pt(10),
                   width=pt(1.0), length=pt(2.5), pad=pt(2))
    ax.legend(
        handles=[Line2D([], [], marker="s", linestyle="none", markersize=pt(6),
                        markerfacecolor=c, markeredgecolor="none", label=lab)
                 for _, lab, c in ENGINES],
        loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2,
        fontsize=pt(11), frameon=False, handletextpad=0.4, columnspacing=1.6,
        labelcolor=INK, borderpad=0, borderaxespad=0)


def main() -> int:
    depths = median_depths()
    missing = [k for k in ((e, d) for e, _, _ in ENGINES for d, _ in DATASETS)
               if k not in depths]
    if missing:
        raise SystemExit(f"spread.tsv is missing draws for {missing}")

    fig = plt.figure(figsize=(W_IN, H_IN))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, SLIDE_W)
    ax.set_ylim(SLIDE_H, 0)  # slide coordinates run downwards
    ax.set_axis_off()

    headline(fig, ax)

    feature_map(ax)
    # feature map -> protein score
    ax.add_patch(FancyArrowPatch(
        (1.505, 1.65), (1.505, 1.81), arrowstyle="-|>", color=UNO,
        linewidth=pt(2.25), mutation_scale=pt(11), shrinkA=0, shrinkB=0,
        zorder=4))
    score_panel(ax)

    depth_chart(fig, depths)
    ax.text(4.635, 2.90, "6 DirectMS1 datasets, identical inputs",
            ha="center", va="center", fontsize=pt(11), color=MUTED)
    ax.text(4.635, 3.08, "3.9–6.6× faster, half the memory",
            ha="center", va="center", fontsize=pt(11), fontweight="bold",
            color=INK)

    style.save(fig, "fig_toc_graphic", tight=False)
    desc = ("ACS TOC graphic: what an MS1-only search is, and the per-dataset "
            "median depth at 1% entrapment FDR behind the claim")
    record("fig.fig-toc-graphic", "figures/fig_toc_graphic.png", kind="figure",
           inputs=["../benchmark/results/spread.tsv"], desc=desc)
    record("fig.fig-toc-graphic-pdf", "figures/fig_toc_graphic.pdf",
           kind="figure", inputs=["../benchmark/results/spread.tsv"],
           desc="Vector twin of the TOC graphic")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
