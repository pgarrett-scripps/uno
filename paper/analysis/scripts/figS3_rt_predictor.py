#!/usr/bin/env python3
"""Fig S2 -- sequence-disjoint RT-predictor benchmark: uno CNN vs DeepLC.

Three unmodified-peptide test sets (Chronologer holdout, ProteomeTools/Prosit
holdout, DeepLC's public unmod.csv), each split 20% calibration / 80% test by a
sequence hash. Both models are calibrated on the same 20%; only the 80% is
scored. Errors are expressed as a percentage of the test set's 5th--95th
percentile RT range so the three RT scales are comparable.

  a-c  residual (predicted - observed) against observed RT, one column per
       test set: median line and 10th--90th percentile band over 12
       equal-count bins of observed RT. Shows WHERE each model errs.
  d-f  cumulative distribution of absolute error; the 95th percentile
       (delta-t95) is marked on each curve. Shows the whole error
       distribution, not one summary number.
Throughput is a separate figure (figS3_rt_throughput.py): it depends on the
peptides, not on which experimental set they came from.

DATA: benchmark/rt_predictor/results/summary.tsv (metrics) and
      benchmark/rt_predictor/results/predictions_<set>.tsv (per-peptide test
      predictions), both written by benchmark/rt_predictor/run.py.
"""

from __future__ import annotations

import os
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

import style  # noqa: E402
from _assets import record  # noqa: E402

BENCHMARK_DIR = style.BENCH / "rt_predictor"
RESULTS = BENCHMARK_DIR / "results"
SUMMARY = RESULTS / "summary.tsv"

# (display name, x-axis label, scale factor applied to the set's RT units)
SETS = {
    "chronologer_holdout": ("Chronologer holdout", "observed RT (hydrophobic index)", 1.0),
    "prosit_unmodified_holdout": ("ProteomeTools holdout", "observed RT (iRT)", 1.0),
    "deeplc_unmodified": ("DeepLC unmodified", "observed RT (min)", 1 / 60),
}
MODELS = (("uno CNN", "uno_cnn", style.C["uno"]),
          ("DeepLC 3.1.13", "deeplc", style.C["ms1searchpy"]))
N_BINS = 12


def load_predictions(name: str):
    d = style.read_tsv(RESULTS / f"predictions_{name}.tsv")
    rng = float(np.percentile(d["observed"], 95) - np.percentile(d["observed"], 5))
    return d, rng


def binned(x, y, n_bins):
    """Median and 10th/90th percentiles of ``y`` in ``n_bins`` equal-count bins of ``x``."""
    order = np.argsort(x, kind="stable")
    edges = np.array_split(order, n_bins)
    cx = np.array([np.median(x[i]) for i in edges])
    med = np.array([np.median(y[i]) for i in edges])
    lo = np.array([np.percentile(y[i], 10) for i in edges])
    hi = np.array([np.percentile(y[i], 90) for i in edges])
    return cx, med, lo, hi


def main() -> int:
    summary = style.read_tsv(SUMMARY).set_index(["dataset", "model"])
    missing = set(SETS) - set(summary.index.get_level_values(0))
    if missing:
        raise RuntimeError(f"missing predictor benchmark datasets: {sorted(missing)}")

    fig = plt.figure(figsize=(style.W2, 4.9))
    gs = fig.add_gridspec(2, 3, hspace=0.62, wspace=0.42,
                          left=0.085, right=0.985, top=0.9, bottom=0.09)
    ax_res = [fig.add_subplot(gs[0, k]) for k in range(3)]
    ax_cdf = [fig.add_subplot(gs[1, k]) for k in range(3)]

    # ---- rows 1-2: residual structure and error distribution, per set ----
    band_lim = 0.0
    preds = {}
    for k, (name, (title, xlabel, scale)) in enumerate(SETS.items()):
        d, rng = load_predictions(name)
        preds[name] = (d, rng)
        obs = d["observed"].to_numpy(float)
        ax = ax_res[k]
        for label, col, color in MODELS:
            resid = 100 * (d[col].to_numpy(float) - obs) / rng
            cx, med, lo, hi = binned(obs * scale, resid, N_BINS)
            ax.fill_between(cx, lo, hi, color=color, alpha=0.18, lw=0, zorder=1)
            ax.plot(cx, med, color=color, lw=1.4, zorder=3)
            band_lim = max(band_lim, float(np.max(np.abs([lo, hi]))))
        ax.axhline(0, color=style.C["neutral"], lw=0.7, ls=(0, (4, 2)), zorder=2)
        ax.set_title(f"{title} (n = {len(d):,})", fontsize=7.5, pad=3)
        ax.set_xlabel(xlabel)
        if k == 0:
            ax.set_ylabel("predicted − observed\n(% of RT range)")
        style.panel_label(ax, "abc"[k], x=-0.30 if k == 0 else -0.18, y=1.04)

        ax = ax_cdf[k]
        for label, col, color in MODELS:
            ae = np.sort(100 * np.abs(d[col].to_numpy(float) - obs) / rng)
            frac = np.arange(1, len(ae) + 1) / len(ae)
            ax.plot(ae, frac, color=color, lw=1.4, zorder=3)
            p95 = float(np.percentile(ae, 95))
            ax.plot([p95], [0.95], marker="o", ms=3.6, mfc="white", mec=color, mew=1.0,
                    zorder=4)
            ax.text(p95, 0.95 - 0.06, f"{p95:.0f}", ha="center", va="top", fontsize=6,
                    color=color)
        ax.axhline(0.95, color=style.C["neutral"], lw=0.6, ls=(0, (2, 2)), zorder=2)
        ax.text(0.5, 0.955, "95%", ha="left", va="bottom", fontsize=6,
                color=style.C["neutral"])
        ax.set_xlim(0, 30)
        ax.set_ylim(0, 1.0)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_xlabel("|error| (% of RT range)")
        if k == 0:
            ax.set_ylabel("fraction of test peptides")
        style.panel_label(ax, "def"[k], x=-0.30 if k == 0 else -0.18, y=1.04)

    lim = 5 * np.ceil(band_lim / 5)
    for ax in ax_res:
        ax.set_ylim(-lim, lim)

    handles = [Line2D([], [], color=c, lw=1.6, label=l) for l, _, c in MODELS]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.99), ncol=2,
               fontsize=7, handlelength=1.6, columnspacing=2.0, borderaxespad=0.0)

    style.save(fig, "figS3_rt_predictor", tight=False)
    provenance = [
        SUMMARY,
        *(RESULTS / f"predictions_{name}.tsv" for name in SETS),
        RESULTS / "environment.json",
        BENCHMARK_DIR / "run.py",
        *(sorted((BENCHMARK_DIR / "sets").glob("*.tsv"))),
        style.REPO / "models" / "rt_cnn_chrono_global_v3_ch32.urtc",
    ]
    inputs = [pathlib.Path(os.path.relpath(path, style.PAPER)).as_posix()
              for path in provenance]
    record("fig.figs3-rt-predictor", "figures/figS3_rt_predictor.png", kind="figure",
           inputs=inputs,
           desc="Uno-training-disjoint unmodified-peptide RT accuracy and throughput")
    record("fig.figs3-rt-predictor-pdf", "figures/figS3_rt_predictor.pdf", kind="figure",
           inputs=inputs,
           desc="Vector twin of the RT predictor benchmark figure")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
