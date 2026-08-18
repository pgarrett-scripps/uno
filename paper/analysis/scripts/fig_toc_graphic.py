#!/usr/bin/env python3
"""TOC / abstract graphic for the JPR submission (SUBMISSION_TODO B7).

ACS "For Table of Contents Only" graphic: the artwork must fit 3.25 in wide by
1.75 in tall (the canvas is exactly that size; savefig cropping is disabled so
the submitted image keeps those dimensions) at 300 dpi.

A conceptual graphic, not a results plot: what MS1-only search is. One run's
MS1 feature map (no MS/MS acquired), each feature carrying only mass, RT and
intensity, is matched by an engine against an in-silico peptide database on
precursor mass and predicted RT, and identity is read out at the protein level
under FDR control. Schematic only; no data (the feature cloud is a seeded
random scatter, drawn to look like an LC-MS feature map).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Circle, FancyBboxPatch  # noqa: E402

import style  # noqa: E402

# a TOC graphic is a fixed-size canvas; cropping it breaks the ACS dimensions
plt.rcParams["savefig.bbox"] = "standard"

W_IN, H_IN = 3.25, 1.75
NOTE = "#5A7A8C"


def arrow(ax, x0, y0, x1, y1):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="-|>", color=style.C["neutral"],
                                linewidth=0.9, shrinkA=0, shrinkB=0,
                                mutation_scale=7), zorder=4)


def main() -> int:
    fig = plt.figure(figsize=(W_IN, H_IN))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W_IN)
    ax.set_ylim(0, H_IN)
    ax.set_axis_off()

    # ---------------- headline ----------------------------------------------
    ax.text(W_IN / 2, 1.64, "protein identification from MS1 features alone",
            ha="center", va="center", fontsize=7.6, fontweight="bold",
            color=style.C["neutral"])

    # ---------------- left: one run's MS1 feature map ------------------------
    px0, px1, py0, py1 = 0.14, 1.06, 0.42, 1.36
    ax.plot([px0, px0, px1], [py1, py0, py0], color=style.C["decoy"], lw=0.8,
            solid_capstyle="butt")
    rng = np.random.default_rng(7)
    n = 90
    fx = px0 + 0.04 + rng.random(n) * (px1 - px0 - 0.08)
    fy = py0 + 0.04 + rng.random(n) * (py1 - py0 - 0.10)
    ax.scatter(fx, fy, s=rng.random(n) * 2.2 + 0.5, color=style.C["decoy"],
               alpha=0.55, linewidth=0, zorder=2)
    # the one feature being searched: highlighted, with its isotope envelope
    hx, hy = px0 + 0.36, py0 + 0.33
    for i, ih in enumerate((0.115, 0.075, 0.042)):
        ax.plot([hx + i * 0.035, hx + i * 0.035], [hy, hy + ih],
                color=style.C["uno"], lw=1.1, solid_capstyle="butt", zorder=3)
    ax.add_patch(Circle((hx + 0.035, hy + 0.05), 0.085, facecolor="none",
                        edgecolor=style.C["uno"], linewidth=0.8, zorder=3))
    ax.text(px0 - 0.045, (py0 + py1) / 2, "m/z", ha="center", va="center",
            fontsize=5.4, color=NOTE, rotation=90)
    ax.text((px0 + px1) / 2, py0 - 0.055, "RT", ha="center", va="center",
            fontsize=5.4, color=NOTE)
    ax.text((px0 + px1) / 2, 0.20, "one run's MS1 feature map",
            ha="center", va="center", fontsize=5.8, color=style.C["neutral"])
    ax.text((px0 + px1) / 2, 0.09, "mass, RT and intensity only",
            ha="center", va="center", fontsize=5.4, style="italic", color=NOTE)

    # "no MS/MS" mark, top-right corner of the map
    cx, cy, r = px1 - 0.13, py1 - 0.14, 0.105
    for i, sh in enumerate((0.055, 0.09, 0.045, 0.07)):
        sx = cx - 0.055 + i * 0.037
        ax.plot([sx, sx], [cy - 0.05, cy - 0.05 + sh], color=style.C["decoy"],
                lw=0.8, solid_capstyle="butt", zorder=3)
    ax.add_patch(Circle((cx, cy), r, facecolor="none",
                        edgecolor=style.C["neutral"], linewidth=0.9, zorder=4))
    d = r / np.sqrt(2)
    ax.plot([cx - d, cx + d], [cy + d, cy - d], color=style.C["neutral"],
            lw=0.9, zorder=5, solid_capstyle="butt")
    ax.text(cx, cy - r - 0.05, "no MS/MS", ha="center", va="center",
            fontsize=5.0, style="italic", color=NOTE, zorder=4)

    # ---------------- middle: the search engine ------------------------------
    bx, bw = 1.63, 0.92
    # in-silico peptide database feeding the engine
    for i in range(3):
        ax.add_patch(FancyBboxPatch(
            (bx - 0.17 + i * 0.02, 1.24 - i * 0.045), 0.30, 0.032,
            boxstyle="round,pad=0,rounding_size=0.012",
            facecolor="white", edgecolor=style.C["decoy"], linewidth=0.7,
            zorder=3 - i * 0.1))
    ax.text(bx, 1.34, "in-silico peptide database", ha="center", va="center",
            fontsize=5.4, color=style.C["neutral"])
    ax.text(bx, 1.44, "protein FASTA → digest", ha="center", va="center",
            fontsize=5.0, style="italic", color=NOTE)
    arrow(ax, bx, 1.10, bx, 0.97)

    ax.add_patch(FancyBboxPatch(
        (bx - bw / 2, 0.52), bw, 0.45,
        boxstyle="round,pad=0,rounding_size=0.035",
        facecolor="#E4EFF7", edgecolor=style.C["uno"], linewidth=1.0,
        zorder=2))
    ax.text(bx, 0.86, "uno", ha="center", va="center", fontsize=8.5,
            fontweight="bold", color=style.C["uno"], zorder=3)
    ax.text(bx, 0.66, "match by mass and\npredicted RT",
            ha="center", va="center", fontsize=5.4, color=style.C["neutral"],
            zorder=3, linespacing=1.35)
    ax.text(bx, 0.40, "MS1 search engine", ha="center", va="center",
            fontsize=5.4, style="italic", color=NOTE)

    arrow(ax, px1 + 0.05, 0.745, bx - bw / 2, 0.745)

    # ---------------- right: FDR-controlled protein list ---------------------
    gx, gw = 2.63, 0.92
    arrow(ax, bx + bw / 2, 0.745, gx - gw / 2 - 0.05, 0.745)
    for i, frac in enumerate((0.95, 0.72, 0.50)):
        yb = 1.00 - i * 0.21
        ax.add_patch(FancyBboxPatch(
            (gx - gw / 2, yb), gw * frac, 0.13,
            boxstyle="round,pad=0,rounding_size=0.025",
            facecolor="#E6F4EF", edgecolor=style.C["target"], linewidth=0.9,
            zorder=2))
        ax.text(gx - gw / 2 + 0.05, yb + 0.065, "protein", ha="left",
                va="center", fontsize=5.2, color=style.C["neutral"], zorder=3)
    ax.text(gx, 0.475, "⋮", ha="center", va="center", fontsize=6,
            color=style.C["neutral"])
    ax.text(gx, 0.20, "protein-level identity", ha="center", va="center",
            fontsize=5.8, color=style.C["neutral"])
    ax.text(gx, 0.09, "target/decoy FDR control", ha="center", va="center",
            fontsize=5.4, style="italic", color=NOTE)

    style.save(fig, "fig_toc_graphic", tight=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
