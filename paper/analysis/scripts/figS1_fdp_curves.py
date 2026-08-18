#!/usr/bin/env python3
"""figS1_fdp_curves.png -- per-dataset protein FDR calibration, all six datasets.

fig4a pooled the six DirectMS1 benchmark datasets; this breaks the same curve
out one small panel per dataset so a reader can see which dataset drives the
pooled behaviour. Axes are identical to fig4a:

    x = engine-REPORTED protein q-value (uno's own `q_value` column;
        ms1searchpy's own pyteomics target/decoy FDR, reproduced exactly)
    y = entrapment-corrected empirical FDP (yeast entrapment, DB-size
        corrected, running-min q) -- scoring copied verbatim from
        benchmark/entrapment.py

A curve above the dashed y = x line means the engine's reported FDR is
OPTIMISTIC on that dataset.

Data: benchmark/results/uno_*.tsv,
      benchmark/results/ms1_<dataset>_rep{1,2,3}/*_proteins_full.tsv

Run:  .venv/bin/python paper/scripts/figures/figS1_fdp_curves.py
"""

from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

import style  # noqa: E402
import entrapment_data as ed  # noqa: E402

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
    ax.axvline(0.01, color=style.C["neutral"], lw=0.5, alpha=0.35, zorder=1)

    def clamp(v):
        return np.maximum(np.asarray(v, float), FLOOR)

    for i, m in enumerate(runs["ms1"]):
        xs, ys = m.calibration()
        ax.plot(xs, clamp(ys), color=style.C["ms1searchpy"], lw=0.85, alpha=0.9,
                zorder=3, label="ms1searchpy (3 reps)" if i == 0 else None)
    u = runs["uno"]
    xs, ys = u.calibration()
    ax.plot(xs, clamp(ys), color=style.C["uno"], lw=1.6, zorder=5, label="uno")

    # the reported-1% readout for this dataset
    _, _, e_uno = u.at_reported(0.01)
    ax.plot([0.01], [max(e_uno, FLOOR)], marker="o", ms=4.6, mfc="white",
            mec=style.C["uno"], mew=1.3, zorder=7)

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
    data = ed.load_all()
    fig, axes = plt.subplots(2, 3, figsize=(style.W2, 4.7))
    for i, ds in enumerate(ed.DATASETS):
        r, c = divmod(i, 3)
        draw(axes[r][c], ds, data[ds], show_x=(r == 1), show_y=(c == 0))
        style.panel_label(axes[r][c], "abcdef"[i], x=-0.14 if c == 0 else -0.07,
                          y=1.08, fontsize=9)

    h, lab = axes[0][0].get_legend_handles_labels()
    order = [lab.index("uno"), lab.index("ms1searchpy (3 reps)")]
    fig.legend([h[i] for i in order], [lab[i] for i in order], ncol=2,
               loc="lower center", bbox_to_anchor=(0.5, 0.038), fontsize=8)
    fig.text(0.5, 0.010, "shaded = engine's reported FDR is optimistic "
             "(empirical FDP above the dashed identity line)",
             ha="center", fontsize=7, color=style.C["neutral"], style="italic")
    fig.tight_layout(rect=(0, 0.085, 1, 1))
    fig.subplots_adjust(hspace=0.34, wspace=0.14)
    style.save(fig, "figS1_fdp_curves", tight=False)

    for ds in ed.DATASETS:
        u = data[ds]["uno"]
        h, y, e = u.at_reported(0.01)
        ms = [m.at_reported(0.01)[2] for m in data[ds]["ms1"]]
        print(f"{ed.SHORT[ds]:6s} uno rep-q 1% -> FDP {e * 100:5.2f}%   "
              f"ms1 reps -> " + ", ".join(f"{v * 100:.2f}%" for v in ms))


if __name__ == "__main__":
    main()
