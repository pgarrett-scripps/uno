#!/usr/bin/env python3
"""Fig S3 -- resource use vs input size: wall-clock, user-CPU and peak RSS
against the number of biosaur2 features in the input file.

Panels:
  a  wall-clock (s)
  b  user CPU (s)          -- load-immune, the fair total-work measure
  c  peak resident set (GiB)

DATA
  six benchmark datasets : benchmark/results/summary.tsv
  input size             : data rows (lines - 1) counted HERE in
                           benchmark/data/<dataset>.features.tsv
  run10 (extra point)    : parsed HERE from the /usr/bin/time -v logs
                           benchmark/results/uno.time and
                           benchmark/results/ms1_rep{1,2,3}/ms1.time,
                           with benchmark/data/run10.features.tsv row count.
                           run10 predates summary.tsv; it is included because it
                           is the largest single-run input measured and extends
                           the x-range from 60k to 211k features.

uno is deterministic -> one point per dataset. ms1searchpy is run three times ->
median point with a min..max bar. Both engines were serialized on the same
machine by benchmark/run.sh so the timings are directly comparable.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import style  # noqa: E402

SUMMARY = style.BENCH / "results" / "summary.tsv"
DATA = style.BENCH / "data"

DATASETS = ["hela5min_09746", "hela5min_09766", "hela5min_09734",
            "hela5min_09778", "hela5min_09950", "faims_LUM2_05433_VG"]
RUN10 = "run10"


def n_features(dataset: str) -> int:
    """Data rows in the biosaur2 feature file (lines minus the header)."""
    path = DATA / f"{dataset}.features.tsv"
    with path.open("rb") as f:
        n = sum(1 for _ in f)
    return n - 1


def parse_time(path: Path) -> dict[str, float]:
    """Pull user-CPU (s), elapsed wall (s) and peak RSS (MB) out of a
    `/usr/bin/time -v` log -- the same three fields benchmark/run.sh parses."""
    txt = path.read_text(errors="replace")
    user = float(re.search(r"User time \(seconds\): ([\d.]+)", txt).group(1))
    rss = float(re.search(r"Maximum resident set size \(kbytes\): (\d+)", txt).group(1)) / 1024
    el = re.search(r"Elapsed \(wall clock\) time.*?:\s*([\d:.]+)\s*$", txt, re.M).group(1)
    parts = [float(p) for p in el.split(":")]
    wall = parts[0] * 3600 + parts[1] * 60 + parts[2] if len(parts) == 3 else parts[0] * 60 + parts[1]
    return {"user_cpu_s": user, "wall_s": wall, "peak_rss_mb": rss}


def collect() -> pd.DataFrame:
    df = pd.read_csv(SUMMARY, sep="\t")
    rows = []
    for d in DATASETS:
        sub = df[df.dataset == d]
        rows.append(dict(dataset=d, tool="uno", **{
            c: float(sub[sub.tool == "uno"][c].iloc[0])
            for c in ("user_cpu_s", "wall_s", "peak_rss_mb")}))
        for _, r in sub[sub.tool == "ms1"].iterrows():
            rows.append(dict(dataset=d, tool="ms1",
                             **{c: float(r[c]) for c in ("user_cpu_s", "wall_s", "peak_rss_mb")}))
    rows.append(dict(dataset=RUN10, tool="uno",
                     **parse_time(style.BENCH / "results" / "uno.time")))
    for i in (1, 2, 3):
        rows.append(dict(dataset=RUN10, tool="ms1",
                         **parse_time(style.BENCH / "results" / f"ms1_rep{i}" / "ms1.time")))
    out = pd.DataFrame(rows)
    out["n_features"] = out["dataset"].map(n_features)
    return out


PANELS = [
    ("wall_s", "wall-clock (s)", 1.0),
    ("user_cpu_s", "user CPU (s)", 1.0),
    ("peak_rss_mb", "peak RSS (GiB)", 1 / 1024),
]


def main() -> int:
    df = collect()
    order = sorted(df["dataset"].unique(), key=lambda d: df[df.dataset == d]["n_features"].iloc[0])
    xs = np.array([df[df.dataset == d]["n_features"].iloc[0] for d in order], dtype=float) / 1e3

    fig, axes = plt.subplots(1, 3, figsize=(style.W2, 2.5))
    for pi, (col, ylab, scale) in enumerate(PANELS):
        ax = axes[pi]
        uno = np.array([df[(df.dataset == d) & (df.tool == "uno")][col].iloc[0]
                        for d in order]) * scale
        med, lo, hi = [], [], []
        for d in order:
            v = df[(df.dataset == d) & (df.tool == "ms1")][col].to_numpy(float) * scale
            m, a, b = style.median_range(v)
            med.append(m); lo.append(m - a); hi.append(b - m)
        ax.errorbar(xs, med, yerr=[lo, hi], color=style.C["ms1searchpy"], marker="s",
                    ms=4, lw=1.1, capsize=2.2, elinewidth=0.8,
                    label="ms1searchpy (median of 3)" if pi == 0 else None, zorder=3)
        ax.plot(xs, uno, color=style.C["uno"], marker="o", ms=4, lw=1.1,
                label="uno (deterministic)" if pi == 0 else None, zorder=4)
        ax.set_xlabel("input features (×10³)")
        ax.set_ylabel(ylab)
        ax.set_ylim(bottom=0)
        ax.set_xlim(0, xs.max() * 1.08)
        style.panel_label(ax, "abc"[pi], x=-0.28)
    axes[0].legend(loc="upper left", fontsize=7)

    style.save(fig, "figS3_resources")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
