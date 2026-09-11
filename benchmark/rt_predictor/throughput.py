#!/usr/bin/env python3
"""Steady-state RT-prediction throughput of the uno CNN and DeepLC 3.1.13.

Both predictors are timed on the same synthetic tryptic library
(``sets/synthetic_tryptic.tsv``, see prepare_library.py) at a series of batch
sizes, after a warm-up call at the largest size, with REPEATS repeats per size.
Only inference is timed: no calibration, and for DeepLC the timed call is
``make_preds(calibrate=False)`` on an already-built DataFrame, which includes
its own feature encoding because that is part of every DeepLC prediction. Both
engines use every core of the same workstation with their default threading.

Run from the repository root after ``cargo build --release -p uno-rt-train
--bin rt_throughput``::

    .venv/bin/python benchmark/rt_predictor/throughput.py
"""

from __future__ import annotations

import csv
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from deeplc import DeepLC

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
LIBRARY = HERE / "sets" / "synthetic_tryptic.tsv"
RESULTS = HERE / "results"
MODEL = ROOT / "models" / "rt_cnn_chrono_global_v3_ch32.urtc"
RT_THROUGHPUT = ROOT / "target" / "release" / "rt_throughput"
SIZES = (100, 300, 1_000, 3_000, 10_000, 30_000, 100_000)
REPEATS = 5


def read_library() -> list[str]:
    with LIBRARY.open() as handle:
        return [line.strip() for line in handle
                if line.strip() and not line.startswith("#") and line.strip() != "sequence"]


def modifications(seq: str) -> str:
    return "|".join(f"{i}|Carbamidomethyl" for i, aa in enumerate(seq, 1) if aa == "C")


def uno_rows() -> list[dict[str, object]]:
    if not RT_THROUGHPUT.is_file():
        raise SystemExit(f"missing {RT_THROUGHPUT}; build rt_throughput in release mode first")
    proc = subprocess.run(
        [str(RT_THROUGHPUT), str(MODEL), str(LIBRARY), str(REPEATS),
         ",".join(str(n) for n in SIZES)],
        cwd=ROOT, check=True, text=True, capture_output=True,
    )
    rows = []
    for line in proc.stdout.splitlines()[1:]:
        n, rep, seconds = line.split("\t")
        rows.append({"model": "uno CNN", "n": int(n), "repeat": int(rep),
                     "seconds": float(seconds)})
    return rows


def deeplc_rows(seqs: list[str]) -> list[dict[str, object]]:
    frame = pd.DataFrame({"seq": seqs, "modifications": [modifications(s) for s in seqs]})
    predictor = DeepLC(verbose=False, batch_num=100_000)
    predictor.make_preds(seq_df=frame.iloc[:max(SIZES)].reset_index(drop=True),
                         calibrate=False)  # warm-up
    rows = []
    for n in SIZES:
        batch = frame.iloc[:n].reset_index(drop=True)
        for rep in range(REPEATS):
            start = time.perf_counter()
            predictor.make_preds(seq_df=batch, calibrate=False)
            rows.append({"model": "DeepLC 3.1.13", "n": n, "repeat": rep,
                         "seconds": time.perf_counter() - start})
    return rows


def main() -> None:
    seqs = read_library()
    if len(seqs) < max(SIZES):
        raise SystemExit(f"library has {len(seqs)} peptides; run prepare_library.py")
    rows = uno_rows() + deeplc_rows(seqs)
    for row in rows:
        row["peptides_per_s"] = row["n"] / row["seconds"]
    out = RESULTS / "throughput.tsv"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: f"{v:.8g}" if isinstance(v, float) else v for k, v in row.items()})

    import deeplc
    env = {
        "python": sys.version, "platform": platform.platform(),
        "cpu_count": os.cpu_count(), "deeplc": getattr(deeplc, "__version__", "unknown"),
        "numpy": np.__version__, "pandas": pd.__version__,
        "model": str(MODEL.relative_to(ROOT)), "sizes": list(SIZES), "repeats": REPEATS,
        "thread_env": {k: os.environ.get(k) for k in (
            "RAYON_NUM_THREADS", "OMP_NUM_THREADS", "TF_NUM_INTRAOP_THREADS",
            "TF_NUM_INTEROP_THREADS")},
    }
    (RESULTS / "throughput_environment.json").write_text(json.dumps(env, indent=2) + "\n")
    summary = pd.DataFrame(rows).groupby(["model", "n"])["peptides_per_s"].median()
    print(summary.to_string())
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
