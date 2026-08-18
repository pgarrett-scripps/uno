#!/usr/bin/env python3
"""Fig S6 -- run-to-run reproducibility: ms1searchpy's replicate spread vs uno's
single deterministic value. SI detail behind main-text Section 3.2; the
scale-free spread panel ((max-min)/median) is panel c of main-text Fig 2.

Panels:
  a  depth @1% entrapment FDR, per dataset: ms1searchpy's three replicates
     (open circles + median tick) against uno's one value (filled diamond)
  b  same at 5% entrapment FDR
  c  run10 close-up: the three individual ms1searchpy replicates and uno,
     @1% and @5%

DATA
  six benchmark datasets : benchmark/results/summary.tsv (written by
      benchmark/run.sh; depths scored with benchmark/entrapment.py)
  run10                  : scored HERE, at run time, straight from the protein
      tables with benchmark/entrapment.py --
        uno   benchmark/results/uno_proteins.tsv
        ms1   benchmark/results/ms1_rep{1,2,3}/run10.features_proteins_full.tsv
      (run10 predates summary.tsv and is not in it. Values are recomputed rather
      than transcribed from benchmark/RESULTS.md; they agree with that table:
      uno 1309/1887, ms1 1220/1941, 1155/1651, 1079/1534.)

uno's determinism claim is a property of the implementation (seeded shuffles,
fixed fold assignment, order-invariant tiebreaks), not a measured spread: only
one uno run exists per dataset.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import style  # noqa: E402

SUMMARY = style.BENCH / "results" / "summary.tsv"
RUN10_UNO = style.BENCH / "results" / "uno_proteins.tsv"
RUN10_MS1 = [style.BENCH / "results" / f"ms1_rep{i}" / "run10.features_proteins_full.tsv"
             for i in (1, 2, 3)]

LABEL = {
    "hela5min_09746": "09746",
    "hela5min_09766": "09766",
    "hela5min_09734": "09734",
    "hela5min_09778": "09778",
    "hela5min_09950": "09950",
    "faims_LUM2_05433_VG": "05433\n(FAIMS)",
    "run10": "run10",
}


def _entrapment():
    """Import benchmark/entrapment.py as a module (it is not a package)."""
    path = style.BENCH / "entrapment.py"
    spec = importlib.util.spec_from_file_location("entrapment", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def score(path: Path) -> dict[float, int]:
    """Entrapment depth at 1% and 5% for one protein table."""
    e = _entrapment()
    out, _, _ = e.curve(e.load(str(path)))
    return out


def collect() -> pd.DataFrame:
    """Long table: dataset, tool, rep, depth_1pct, depth_5pct."""
    df = pd.read_csv(SUMMARY, sep="\t")[
        ["dataset", "tool", "rep", "depth_1pct", "depth_5pct"]
    ]
    rows = []
    u = score(RUN10_UNO)
    rows.append(dict(dataset="run10", tool="uno", rep=1,
                     depth_1pct=u[0.01], depth_5pct=u[0.05]))
    for i, p in enumerate(RUN10_MS1, start=1):
        m = score(p)
        rows.append(dict(dataset="run10", tool="ms1", rep=i,
                         depth_1pct=m[0.01], depth_5pct=m[0.05]))
    return pd.concat([df, pd.DataFrame(rows)], ignore_index=True)


def main() -> int:
    df = collect()
    datasets = [d for d in LABEL if d in set(df["dataset"])]
    x = np.arange(len(datasets), dtype=float)

    fig, axes = plt.subplots(1, 3, figsize=(style.W2, 2.9))

    # -- a, b: absolute depth, ms1 spread vs uno's single value ---------------
    for pi, (col, ylab) in enumerate([
        ("depth_1pct", "proteins @ 1% entrapment FDR"),
        ("depth_5pct", "proteins @ 5% entrapment FDR"),
    ]):
        ax = axes[pi]
        for i, d in enumerate(datasets):
            reps = df[(df.dataset == d) & (df.tool == "ms1")][col].astype(float)
            style.spread_marker(
                ax, x[i] + 0.16, reps.to_numpy(), style.C["ms1searchpy"],
                half_width=0.15,
                label="ms1searchpy replicates" if (pi == 0 and i == 0) else None,
                median_label="ms1searchpy median" if (pi == 0 and i == 0) else None,
            )
            uno = float(df[(df.dataset == d) & (df.tool == "uno")][col].iloc[0])
            ax.plot(x[i] - 0.16, uno, marker="D", ms=4, color=style.C["uno"],
                    ls="none", zorder=6,
                    label="uno (1 run, deterministic)" if (pi == 0 and i == 0) else None)
        ax.set_xticks(x)
        ax.set_xticklabels([LABEL[d].replace("\n", " ") for d in datasets],
                           fontsize=6, rotation=45, ha="right")
        ax.set_ylabel(ylab)
        ax.set_xlim(-0.6, len(datasets) - 0.4)
        ax.margins(y=0.18)
        style.panel_label(ax, "abc"[pi], x=-0.19)
    axes[0].legend(loc="upper left", fontsize=6)

    # -- c: run10 close-up ----------------------------------------------------
    ax = axes[2]
    r10 = df[df.dataset == "run10"]
    cats = ["uno", "rep 1", "rep 2", "rep 3"]
    xr = np.arange(4, dtype=float)
    colors = [style.C["uno"]] + [style.C["ms1searchpy"]] * 3
    for col, off, hatch, lab in [
        ("depth_1pct", -0.19, None, "@ 1% FDR"),
        ("depth_5pct", 0.19, "///", "@ 5% FDR"),
    ]:
        vals = [float(r10[(r10.tool == "uno")][col].iloc[0])] + [
            float(r10[(r10.tool == "ms1") & (r10.rep == i)][col].iloc[0]) for i in (1, 2, 3)
        ]
        b = ax.bar(xr + off, vals, 0.34, color=colors, edgecolor="white",
                   hatch=hatch, linewidth=0.4, zorder=2)
        ax.bar_label(b, fmt="%d", fontsize=5.6, padding=1, rotation=90)
        ax.bar(np.nan, np.nan, color=style.C["light"], hatch=hatch,
               edgecolor="white", label=lab)
    ax.set_xticks(xr)
    ax.set_xticklabels(cats, fontsize=6.5)
    ax.set_ylabel("proteins (run10)")
    ax.set_xlabel("run10: uno vs 3 reruns")
    ax.set_ylim(0, 2500)
    ax.legend(loc="upper right", fontsize=6)
    style.panel_label(ax, "c", x=-0.19)

    style.save(fig, "figS6_determinism")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
