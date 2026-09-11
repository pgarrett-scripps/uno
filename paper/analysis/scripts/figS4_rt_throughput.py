#!/usr/bin/env python3
"""Fig S3 -- RT-prediction throughput of the uno CNN vs DeepLC 3.1.13.

Throughput depends on the peptides, not on which experimental set they came
from, so both models are timed on one synthetic tryptic library
(benchmark/rt_predictor/sets/synthetic_tryptic.tsv) at batch sizes from 100 to
100,000 peptides, five repeats each after a warm-up call. Inference only: no
calibration.

  a  wall time per batch against batch size (log-log). The flat part at small
     batches is each engine's fixed per-call overhead; the slope-1 part is
     steady-state cost per peptide.
  b  peptides per second against batch size. Lines are medians of the five
     repeats; the shaded band is the min-max range.

DATA: benchmark/rt_predictor/results/throughput.tsv, written by
      benchmark/rt_predictor/throughput.py.
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
DATA = RESULTS / "throughput.tsv"
MODELS = (("uno CNN", style.C["uno"]), ("DeepLC 3.1.13", style.C["ms1searchpy"]))


def main() -> int:
    d = style.read_tsv(DATA)
    fig, (ax_t, ax_r) = plt.subplots(1, 2, figsize=(style.W2, 2.7))
    rate_at_max = {}
    for label, color in MODELS:
        g = d[d["model"] == label].groupby("n")
        n = np.array(sorted(g.groups))
        sec = np.array([g.get_group(k)["seconds"].to_numpy(float) for k in n])
        rate = n[:, None] / sec
        for ax, y in ((ax_t, sec), (ax_r, rate)):
            ax.fill_between(n, y.min(axis=1), y.max(axis=1), color=color, alpha=0.2, lw=0)
            ax.plot(n, np.median(y, axis=1), color=color, lw=1.4, marker="o", ms=3.2,
                    mfc="white", mec=color, mew=1.0, zorder=3)
        rate_at_max[label] = float(np.median(rate[-1]))
        up = label != "uno CNN"
        ax_r.text(n[-1], rate_at_max[label] * (1.25 if up else 0.8),
                  f"{rate_at_max[label]:,.0f}/s", ha="right",
                  va="bottom" if up else "top", fontsize=6.3, color=color)

    for ax in (ax_t, ax_r):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("peptides per batch")
        ax.yaxis.set_minor_formatter(plt.NullFormatter())
    ax_t.set_ylabel("wall time per batch (s)")
    ax_r.set_ylabel("throughput (peptides / s)")
    speedup = rate_at_max["uno CNN"] / rate_at_max["DeepLC 3.1.13"]
    ax_t.set_title("inference only, same workstation", fontsize=7.5, pad=3)
    ax_r.set_title(f"{speedup:.0f}× at 100,000 peptides", fontsize=7.5, pad=3)
    style.panel_label(ax_t, "a", x=-0.22, y=1.03)
    style.panel_label(ax_r, "b", x=-0.22, y=1.03)
    handles = [Line2D([], [], color=c, lw=1.6, marker="o", ms=3.2, mfc="white", mec=c,
                      label=l) for l, c in MODELS]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.0), ncol=2,
               fontsize=7, handlelength=1.8, columnspacing=2.0, borderaxespad=0.0)
    fig.tight_layout(w_pad=1.6, rect=(0, 0, 1, 0.9))
    style.save(fig, "figS4_rt_throughput", tight=False)

    provenance = [DATA, RESULTS / "throughput_environment.json",
                  BENCHMARK_DIR / "throughput.py", BENCHMARK_DIR / "prepare_library.py",
                  BENCHMARK_DIR / "sets" / "synthetic_tryptic.tsv",
                  style.REPO / "models" / "rt_cnn_chrono_global_v3_ch32.urtc"]
    inputs = [pathlib.Path(os.path.relpath(p, style.PAPER)).as_posix() for p in provenance]
    record("fig.figs4-rt-throughput", "figures/figS4_rt_throughput.png", kind="figure",
           inputs=inputs, desc="RT-prediction throughput of the uno CNN and DeepLC vs batch size")
    record("fig.figs4-rt-throughput-pdf", "figures/figS4_rt_throughput.pdf", kind="figure",
           inputs=inputs, desc="Vector twin of the RT throughput figure")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
