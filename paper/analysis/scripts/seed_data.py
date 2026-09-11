"""The frozen inventory of primary random draws, shared by figures and stats.

Only the six validation files are included; run10 remains development data.
Every row must identify an existing protein output and an unambiguous draw.
"""
import csv
import os

import entrapment_data as ed
import style

SPREAD = style.BENCH / "results/spread.tsv"


def rows():
    with SPREAD.open() as handle:
        data = [r for r in csv.DictReader(handle, delimiter="\t") if r["dataset"] in ed.DATASETS]
    seen = set()
    for r in data:
        key = r["dataset"], r["engine"], int(r["draw"])
        if key in seen:
            raise ValueError(f"duplicate draw: {key}")
        seen.add(key)
    for ds in ed.DATASETS:
        for engine in ("uno", "ms1searchpy"):
            if sum(r["dataset"] == ds and r["engine"] == engine for r in data) < 3:
                raise ValueError(f"need at least three draws: {ds} {engine}")
    return data


def protein_path(r):
    ds, draw = r["dataset"], int(r["draw"])
    if r["engine"] == "uno":
        return (ed.RESULTS / f"uno_{ds}.tsv" if draw == 0 else
                ed.RESULTS / f"seed_sweep/uno_{ds}_seed{draw}.tsv")
    if r["engine"] != "ms1searchpy":
        raise ValueError(f"unknown engine: {r['engine']}")
    return ed.ms1_path(ds, draw)


def inputs():
    return [os.path.relpath(p, style.PAPER) for p in [SPREAD, *(protein_path(r) for r in rows())]]


def load():
    data = {ds: {"uno": [], "ms1": []} for ds in ed.DATASETS}
    for r in rows():
        engine = "uno" if r["engine"] == "uno" else "ms1"
        run = ed.Run(r["dataset"], engine, int(r["draw"]), str(protein_path(r)))
        for cutoff, col in ((0.01, "depth_1pct"), (0.05, "depth_5pct")):
            if run.depth_at_entrapment(cutoff) != int(r[col]):
                raise ValueError(f"spread.tsv disagrees with protein output: {r}")
        data[r["dataset"]][engine].append(run)
    return data
