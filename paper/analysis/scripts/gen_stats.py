#!/usr/bin/env python3
"""Write stats.json at the manuscript root: every number the prose states.

STATUS: skeleton. The manuscript's numbers are still typed into the prose from
the pre-scaffold era, and every identification number is provisional until the
C1 re-run lands (see ../../SUBMISSION_TODO.md). As sentences are migrated onto
`#s("id")`, the values they state get declared here, derived from the benchmark
outputs -- benchmark/results/summary.tsv and friends -- or from the committed
extracts under ../../figdata/, never typed as constants (a constant typed here
is a number nobody can reproduce).

One script per project: it writes a single file, so splitting it across scripts
would mean each clobbering the last. Import from other modules if it grows.
"""
from __future__ import annotations

import csv
import statistics
from pathlib import Path

from _stats import Stats

HERE = Path(__file__).resolve().parent
# Lives outside the manuscript, at the uno repository root; declared as an
# input so `just check-stats` notices when a re-run moves it.
SUMMARY = HERE.parent.parent.parent / "benchmark" / "results" / "summary.tsv"


def main() -> int:
    with SUMMARY.open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    uno_wall = [float(r["wall_s"]) for r in rows if r["tool"] == "uno"]
    ms1_wall = [float(r["wall_s"]) for r in rows if r["tool"] == "ms1"]
    # Per-dataset factor: the ms1searchpy median replicate over uno's single
    # (deterministic) run, matching how Section 3.1 states the speed gap.
    factors = []
    for d in {r["dataset"] for r in rows}:
        uno = [float(r["wall_s"]) for r in rows
               if r["dataset"] == d and r["tool"] == "uno"]
        ms1 = [float(r["wall_s"]) for r in rows
               if r["dataset"] == d and r["tool"] == "ms1"]
        if uno and ms1:
            factors.append(statistics.median(ms1) / uno[0])

    st = Stats()
    st.add("bench.uno_wall_min", min(uno_wall), fmt=".0f", unit="s",
           desc="Fastest uno single-run search wall time across the six datasets",
           sign="+", between=(1, 600))
    st.add("bench.uno_wall_max", max(uno_wall), fmt=".0f", unit="s",
           desc="Slowest uno single-run search wall time across the six datasets",
           sign="+", between=(1, 600))
    st.add("bench.ms1_wall_min", min(ms1_wall), fmt=".0f", unit="s",
           desc="Fastest ms1searchpy replicate wall time across the six datasets",
           sign="+", between=(1, 3600))
    st.add("bench.ms1_wall_max", max(ms1_wall), fmt=".0f", unit="s",
           desc="Slowest ms1searchpy replicate wall time across the six datasets",
           sign="+", between=(1, 3600))
    # Guarded above 1: the prose calls uno "several times cheaper", which stops
    # being true the day this dips toward parity.
    st.add("bench.speed_factor_min", min(factors), fmt=".1f", unit="x",
           desc="Smallest per-dataset ms1searchpy-median/uno wall-time factor",
           sign="+", between=(1, 50))
    st.add("bench.speed_factor_max", max(factors), fmt=".1f", unit="x",
           desc="Largest per-dataset ms1searchpy-median/uno wall-time factor",
           sign="+", between=(1, 50))

    return st.write(inputs=["../benchmark/results/summary.tsv"])


if __name__ == "__main__":
    raise SystemExit(main())
