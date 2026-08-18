#!/usr/bin/env python3
"""fig3_entrapment.png -- protein FDR calibration at the operating points.

(a) Entrapment FDP of each engine's own reported 1% protein list.
(b) Entrapment FDP of each engine's own reported 5% protein list.

The full calibration and depth/FDP curves are intentionally left to Figure S1.
The main-text figure answers the simpler question the calibration claim rests on:
does an engine's reported list meet its nominal operating point?

The two axes in (a) are DIFFERENT quantities and are labelled as such: the
depth numbers quoted in the manuscript (benchmark/RESULTS_canonical.md) are the
entrapment numbers, not the engines' self-reported q.

Data: benchmark/results/uno_*.tsv and
      benchmark/results/ms1_<dataset>_rep{1,2,3}/*_proteins_full.tsv
Scoring: benchmark/entrapment.py (copied verbatim into entrapment_data.py)

Run:  .venv/bin/python paper/scripts/figures/fig3_entrapment.py
"""

from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

import style  # noqa: E402
import entrapment_data as ed  # noqa: E402

def step_lookup(xs, ys, grid):
    """Value of the step function (xs -> ys) at each grid point; NaN before xs[0]."""
    xs = np.asarray(xs, float)
    ys = np.asarray(ys, float)
    idx = np.searchsorted(xs, grid, side="right") - 1
    out = np.where(idx >= 0, ys[np.clip(idx, 0, len(ys) - 1)], np.nan)
    return out


def fdp_at_reported_q(run, q):
    xs, ys = run.calibration()
    return float(step_lookup(xs, ys, np.array([q]))[0])


def operating_point_panel(ax, data, q):
    x = np.arange(len(ed.DATASETS), dtype=float)
    uno = np.array([fdp_at_reported_q(data[ds]["uno"], q)
                    for ds in ed.DATASETS]) * 100
    ms1 = np.array([[fdp_at_reported_q(run, q) for run in data[ds]["ms1"]]
                    for ds in ed.DATASETS]) * 100

    # A faint connector makes the engine comparison readable without implying
    # that the three stochastic comparator runs are biological replicates.
    ms1_med = np.median(ms1, axis=1)
    for i in range(len(x)):
        ax.plot([x[i] - 0.12, x[i] + 0.12], [uno[i], ms1_med[i]],
                color=style.C["neutral"], alpha=0.35, lw=0.7, zorder=1)
    ax.scatter(x - 0.12, uno, s=30, marker="D", color=style.C["uno"],
               zorder=4, label="uno")
    for rep in range(ms1.shape[1]):
        ax.scatter(x + 0.12, ms1[:, rep], s=25, marker="o", facecolors="white",
                   edgecolors=style.C["ms1searchpy"], linewidths=1.0, zorder=3,
                   label="ms1searchpy runs" if rep == 0 else None)
    ax.scatter(x + 0.12, ms1_med, s=34, marker="_",
               color=style.C["ms1searchpy"], linewidths=2.0, zorder=5,
               label="ms1searchpy median")

    nominal = q * 100
    ax.axhline(nominal, color=style.C["neutral"], lw=0.9,
               ls=(0, (4, 2)), zorder=0)
    ax.text(0.99, nominal, f" nominal {nominal:g}%", transform=ax.get_yaxis_transform(),
            va="bottom", ha="right", fontsize=6.5, color=style.C["neutral"])
    ax.fill_between([-0.5, len(x) - 0.5], nominal, max(np.max(ms1), np.max(uno)) * 1.18,
                    color=style.C["light"], alpha=0.45, lw=0, zorder=-1)

    ax.set_xlim(-0.5, len(x) - 0.5)
    ax.set_ylim(0, max(nominal * 1.45, np.max(ms1) * 1.15, np.max(uno) * 1.15))
    ax.set_xticks(x)
    ax.set_xticklabels([ed.SHORT[ds].replace("LUM2_", "") for ds in ed.DATASETS],
                       rotation=35, ha="right")
    ax.set_title(f"engine-reported {nominal:g}% protein list")
    ax.set_ylabel("entrapment FDP (%)")
    return uno, ms1


def main():
    data = ed.load_all()
    fig, axes = plt.subplots(1, 2, figsize=(style.W2, 2.85))
    uno1, ms1_1 = operating_point_panel(axes[0], data, 0.01)
    uno5, ms1_5 = operating_point_panel(axes[1], data, 0.05)
    style.panel_label(axes[0], "a", x=-0.16)
    style.panel_label(axes[1], "b", x=-0.16)
    axes[1].legend(loc="lower left", frameon=False)
    fig.subplots_adjust(wspace=0.30, bottom=0.24)
    style.save(fig, "fig3_entrapment")

    # the SI table cites exactly these numbers
    ed.write_table(style.FIGDATA / "entrapment_curves.tsv", data)

    print(f"  median empirical FDP at reported q=1%%: uno {np.median(uno1):.2f}%, "
          f"ms1 {np.median(np.median(ms1_1, axis=1)):.2f}%")
    print(f"  median empirical FDP at reported q=5%%: uno {np.median(uno5):.2f}%, "
          f"ms1 {np.median(np.median(ms1_5, axis=1)):.2f}%")


if __name__ == "__main__":
    main()
