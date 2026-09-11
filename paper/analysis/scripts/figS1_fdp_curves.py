#!/usr/bin/env python3
"""Per-dataset calibration across all primary random draws.

Thin step curves show individual uno seeds and ms1searchpy reruns. Bold curves
are pointwise medians, evaluated only where all runs have a nonempty accepted
list. Open circles retain every draw at reported 1% and 5%. Inputs are protein
outputs inventoried by benchmark/results/spread.tsv; seed_data validates their
depths against that inventory. Writes figdata/entrapment_curves.tsv for stats.
"""

from __future__ import annotations

import os
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

import style  # noqa: E402
import entrapment_data as ed  # noqa: E402
import seed_data
from _assets import record  # noqa: E402

XLIM = (1.5e-3, 0.32)
# FDP is a count ratio, so it can be exactly 0 (no yeast accepted at all), which
# a log axis cannot show. Zero is drawn on the floor tick, labelled "0".
FLOOR = 1.65e-3
YLIM = (1.35e-3, 0.42)
TICKS = [0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2]
TICKLBL = ["0.2", "0.5", "1", "2", "5", "10", "20"]
# drop the 0.2 y-tick: it would collide with the "0" floor tick
YTICKS = [FLOOR] + TICKS[1:]
YTICKLBL = ["0"] + TICKLBL[1:]


def draw(ax, ds, runs, show_x, show_y):
    ax.set_xscale("log")
    ax.set_yscale("log")
    diag = np.array(XLIM)
    ax.fill_between(diag, diag, YLIM[1] * 1.4, color=style.C["light"],
                    alpha=0.55, lw=0, zorder=0)
    ax.plot(diag, diag, ls=(0, (4, 2)), lw=0.9, color=style.C["neutral"], zorder=2)
    for op in (0.01, 0.05):
        ax.axvline(op, color=style.C["neutral"], lw=0.5, alpha=0.35, zorder=1)

    def clamp(v):
        return np.maximum(np.asarray(v, float), FLOOR)

    grid = np.unique(np.r_[np.geomspace(*XLIM, 250), 0.01, 0.05])
    for engine, label in (("uno", "uno"), ("ms1", "ms1searchpy")):
        matrix = []
        for run in runs[engine]:
            xs, ys = run.calibration()
            xs, ys = np.asarray(xs), np.asarray(ys)
            ax.step(xs, clamp(ys), where="post", color=style.C[label], lw=0.6,
                    alpha=0.25, zorder=3)
            indices = np.searchsorted(xs, grid, side="right") - 1
            sampled = np.where(indices >= 0, ys[np.maximum(indices, 0)], np.nan)
            matrix.append(sampled)
            for op in (0.01, 0.05):
                human, _, fdp = run.at_reported(op)
                if human:
                    ax.plot(op, max(fdp, FLOOR), "o", ms=2.5, mfc="white",
                            mec=style.C[label], mew=0.6, zorder=6)
        matrix = np.asarray(matrix)
        # Draw medians only where all searches have a nonempty accepted list.
        valid = np.all(np.isfinite(matrix), axis=0)
        ax.plot(grid[valid], clamp(np.median(matrix[:, valid], axis=0)),
                color=style.C[label], lw=1.6, zorder=5, label=label)
    ax.text(0.97, 0.05, f"n={len(runs['uno'])} / {len(runs['ms1'])}",
            transform=ax.transAxes, ha="right", fontsize=6)

    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    ax.set_xticks(TICKS)
    ax.set_xticklabels(TICKLBL)
    ax.set_yticks(YTICKS)
    ax.set_yticklabels(YTICKLBL)
    ax.tick_params(labelbottom=show_x, labelleft=show_y)
    ax.minorticks_off()
    ax.set_title(ed.NICE[ds], pad=3)
    if show_x:
        ax.set_xlabel("engine-reported q (%)")
    if show_y:
        ax.set_ylabel("entrapment FDP (%)")


def main():
    data = seed_data.load()
    fig, axes = plt.subplots(2, 3, figsize=(style.W2, 4.7))
    for i, ds in enumerate(ed.DATASETS):
        r, c = divmod(i, 3)
        draw(axes[r][c], ds, data[ds], show_x=(r == 1), show_y=(c == 0))
        style.panel_label(axes[r][c], "abcdef"[i], x=-0.14 if c == 0 else -0.07,
                          y=1.08, fontsize=9)

    h, lab = axes[0][0].get_legend_handles_labels()
    order = [lab.index("uno"), lab.index("ms1searchpy")]
    fig.legend([h[i] for i in order], [lab[i] for i in order], ncol=2,
               loc="lower center", bbox_to_anchor=(0.5, 0.038), fontsize=8)
    fig.text(0.5, 0.010, "shaded = engine's reported FDR is optimistic "
             "(empirical FDP above the dashed identity line)",
             ha="center", fontsize=7, color=style.C["neutral"], style="italic")
    fig.tight_layout(rect=(0, 0.085, 1, 1))
    fig.subplots_adjust(hspace=0.34, wspace=0.14)
    style.save(fig, "figS1_fdp_curves", tight=False)

    # The SI table and every cal.* stat read this; it is written here because
    # this figure is now the only consumer of the same loaded runs.
    ed.write_table(style.FIGDATA / "entrapment_curves.tsv", data)

    inputs = seed_data.inputs()
    record("fig.figs1-fdp-curves", "figures/figS1_fdp_curves.png", kind="figure",
           inputs=inputs,
           desc="Per-dataset reported q against entrapment FDP, all six primary datasets")
    record("fig.figs1-fdp-curves-pdf", "figures/figS1_fdp_curves.pdf", kind="figure",
           inputs=inputs, desc="Vector twin of the per-dataset FDP curve figure")



if __name__ == "__main__":
    main()
