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

from _stats import Stats


def main() -> int:
    st = Stats()
    return st.write()


if __name__ == "__main__":
    raise SystemExit(main())
