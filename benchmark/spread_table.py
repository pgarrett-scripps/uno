#!/usr/bin/env python3
"""Combined randomization-spread table for both engines.

The comparison the manuscript needs is like for like: uno's depth over N
independent seeds (`uno --seed`, benchmark/seed_sweep.sh) against ms1searchpy's
depth over N independent runs (benchmark/run.sh replicates + ms1_extra_reps.sh).
uno's seed 0 and ms1searchpy's rep 1 are ordinary draws here, not privileged
values.

Writes results/spread.tsv: engine, dataset, draw, depth_1pct, depth_5pct.
"""
from __future__ import annotations

import glob
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import entrapment as e  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results", "spread.tsv")


def depths(path: str) -> tuple[int, int]:
    d, _, _ = e.curve(e.load(path))
    return d[0.01], d[0.05]


def main() -> int:
    rows: list[tuple[str, str, int, int, int]] = []

    # uno: seed 0 is the shipped default run; seeds 1..N from the sweep.
    # run10 predates summary.tsv and is stored under its own filenames
    # (results/uno_proteins.tsv, results/ms1_rep{1,2,3}/) -- normalize both to
    # the dataset name "run10" so it joins the rest of the table.
    for path in sorted(glob.glob(os.path.join(HERE, "results", "uno_*.tsv"))):
        base = os.path.basename(path)
        if base.startswith("uno_raw") or not os.path.getsize(path):
            continue
        m = re.fullmatch(r"uno_(.+)\.tsv", base)
        if not m:
            continue
        ds = "run10" if m.group(1) == "proteins" else m.group(1)
        rows.append(("uno", ds, 0, *depths(path)))
    for d in sorted(glob.glob(os.path.join(HERE, "results", "ms1_rep*"))):
        m = re.fullmatch(r"ms1_rep(\d+)", os.path.basename(d))
        pf = os.path.join(d, "run10.features_proteins_full.tsv")
        if m and os.path.exists(pf) and os.path.getsize(pf):
            rows.append(("ms1searchpy", "run10", int(m.group(1)), *depths(pf)))
    for path in sorted(glob.glob(os.path.join(HERE, "results", "seed_sweep", "uno_*_seed*.tsv"))):
        m = re.fullmatch(r"uno_(.+)_seed(\d+)\.tsv", os.path.basename(path))
        if m and os.path.getsize(path):
            rows.append(("uno", m.group(1), int(m.group(2)), *depths(path)))

    # ms1searchpy: one directory per replicate, from run.sh and ms1_extra_reps.sh.
    for d in sorted(glob.glob(os.path.join(HERE, "results", "ms1_*_rep*"))):
        m = re.fullmatch(r"ms1_(.+)_rep(\d+)", os.path.basename(d))
        if not m:
            continue
        ds, rep = m.group(1), int(m.group(2))
        pf = os.path.join(d, f"{ds}.features_proteins_full.tsv")
        if os.path.exists(pf) and os.path.getsize(pf):
            rows.append(("ms1searchpy", ds, rep, *depths(pf)))

    rows.sort()
    with open(OUT, "w") as fh:
        fh.write("engine\tdataset\tdraw\tdepth_1pct\tdepth_5pct\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print(f"wrote {OUT} ({len(rows)} rows)")

    # Per-engine, per-dataset summary to stdout.
    import statistics as st
    keyed: dict[tuple[str, str], list[int]] = {}
    for eng, ds, _, d1, _d5 in rows:
        keyed.setdefault((eng, ds), []).append(d1)
    print(f"\n{'dataset':<24}{'engine':<13}{'n':>3}{'min':>7}{'med':>7}{'max':>7}{'range':>7}{'%med':>7}")
    for (eng, ds), v in sorted(keyed.items(), key=lambda kv: (kv[0][1], kv[0][0])):
        med = st.median(v)
        rng = max(v) - min(v)
        print(f"{ds:<24}{eng:<13}{len(v):>3}{min(v):>7}{med:>7.0f}{max(v):>7}{rng:>7}"
              f"{100 * rng / med:>7.1f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
