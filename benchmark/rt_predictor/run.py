#!/usr/bin/env python3
"""Benchmark the shipped uno CNN against DeepLC 3.1.13.

Both predictors receive the same deterministic 20% calibration / 80% test
split and the SAME calibration: a least-squares affine map from the model's
native output scale to the observed scale, fitted on the 20%.  This is a model
comparison, so neither predictor gets a calibration the other does not.
(DeepLC's own piecewise ``calibrate_preds`` is what ms1searchpy uses in-search;
on the smallest set here it extrapolates poorly with ~200 calibration peptides,
which would measure the calibrator rather than the model.)  Accuracy is
measured only on the held-out 80%.  Prediction throughput excludes calibration.

Run from the repository root after ``cargo build --release -p uno-rt-train
--bin rt_eval`` and with DeepLC 3.1.13 plus pandas/numpy available::

    .venv/bin/python benchmark/rt_predictor/run.py
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import platform
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
from deeplc import DeepLC

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
SETS = HERE / "sets"
RESULTS = HERE / "results"
MODEL = ROOT / "models" / "rt_cnn_chrono_global_v3_ch32.urtc"
RT_EVAL = ROOT / "target" / "release" / "rt_eval"
SPEED_RE = re.compile(r"speed\s+:\s+(\d+) pred/s parallel")


def read_set(path: Path) -> tuple[list[str], np.ndarray]:
    seqs: list[str] = []
    labels: list[float] = []
    with path.open() as handle:
        for row in csv.reader((line for line in handle if not line.startswith("#")), delimiter="\t"):
            if not row or row[0] == "sequence":
                continue
            seqs.append(row[0])
            labels.append(float(row[1]))
    return seqs, np.asarray(labels, dtype=float)


def labeled(path: Path) -> bool:
    """True when the set carries observed RTs, i.e. a second column."""
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields and fields[0] == "sequence":
                return len(fields) > 1
            return len(fields) > 1
    return False


def calibration_mask(seqs: list[str]) -> np.ndarray:
    return np.asarray([
        int.from_bytes(hashlib.sha256(f"uno-rt-cal-v1\0{seq}".encode()).digest()[:8], "big") % 5 == 0
        for seq in seqs
    ])


def modifications(seq: str) -> str:
    # Same fixed-C mapping used by ms1searchpy/DeepLC for alkylated peptides.
    return "|".join(f"{i}|Carbamidomethyl" for i, aa in enumerate(seq, 1) if aa == "C")


def affine_fit(y: np.ndarray, pred: np.ndarray) -> tuple[float, float]:
    design = np.column_stack((pred, np.ones_like(pred)))
    slope, intercept = np.linalg.lstsq(design, y, rcond=None)[0]
    return float(slope), float(intercept)


def ranks(values: np.ndarray) -> np.ndarray:
    return pd.Series(values).rank(method="average").to_numpy(dtype=float)


def metrics(y: np.ndarray, pred: np.ndarray) -> dict[str, float]:
    err = np.abs(y - pred)
    spread = float(np.percentile(y, 95) - np.percentile(y, 5))
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    return {
        "mae": float(err.mean()),
        "normalized_mae_pct": 100.0 * float(err.mean()) / spread,
        "dt95": float(np.percentile(err, 95)),
        "r2": 1.0 - ss_res / ss_tot,
        "spearman": float(np.corrcoef(ranks(y), ranks(pred))[0, 1]),
    }


def cnn_predictions(set_path: Path, seqs: list[str]) -> tuple[np.ndarray, float]:
    if not RT_EVAL.is_file():
        raise SystemExit(f"missing {RT_EVAL}; build rt_eval in release mode first")
    with tempfile.TemporaryDirectory(prefix="uno-rt-benchmark-") as tmp:
        pred_path = Path(tmp) / "predictions.tsv"
        proc = subprocess.run(
            [str(RT_EVAL), str(MODEL), str(set_path), str(pred_path)],
            cwd=ROOT,
            check=True,
            text=True,
            capture_output=True,
        )
        match = SPEED_RE.search(proc.stdout)
        if not match:
            raise RuntimeError(f"could not parse rt_eval throughput:\n{proc.stdout}\n{proc.stderr}")
        by_seq: dict[str, float] = {}
        with pred_path.open() as handle:
            for line in handle:
                seq, value = line.rstrip().split("\t")
                by_seq[seq] = float(value)
    return np.asarray([by_seq[seq] for seq in seqs]), float(match.group(1))


def deeplc_predictions(seqs: list[str], labels: np.ndarray, cal: np.ndarray) -> tuple[np.ndarray, float, float]:
    """Uncalibrated DeepLC output mapped to the observed scale by the same affine
    fit uno receives. Throughput is timed on the uncalibrated test-set pass."""
    frame = pd.DataFrame({
        "seq": seqs,
        "modifications": [modifications(seq) for seq in seqs],
    })
    predictor = DeepLC(verbose=False, batch_num=100_000)
    cal_frame = frame.loc[cal].reset_index(drop=True)
    raw_cal = np.asarray(predictor.make_preds(seq_df=cal_frame, calibrate=False), dtype=float)
    start = time.perf_counter()
    slope, intercept = affine_fit(labels[cal], raw_cal)
    calibration_s = time.perf_counter() - start

    test_frame = frame.loc[~cal].reset_index(drop=True)
    start = time.perf_counter()
    raw_test = np.asarray(predictor.make_preds(seq_df=test_frame, calibrate=False), dtype=float)
    prediction_s = time.perf_counter() - start
    return slope * raw_test + intercept, calibration_s, len(raw_test) / prediction_s


def result_row(dataset: str, model: str, n: int, n_cal: int, n_test: int,
               values: dict[str, float], calibration_s: float,
               predictions_per_s: float) -> dict[str, object]:
    return {
        "dataset": dataset,
        "model": model,
        "n_total": n,
        "n_calibration": n_cal,
        "n_test": n_test,
        **values,
        "calibration_s": calibration_s,
        "predictions_per_s": predictions_per_s,
    }


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    # Only the labeled accuracy sets. sets/ also holds the unlabeled synthetic
    # library that throughput.py times, and a bare glob picks it up and dies on
    # its missing RT column.
    set_paths = [p for p in sorted(SETS.glob("*.tsv")) if labeled(p)]
    if not set_paths:
        raise SystemExit(f"no labeled benchmark sets under {SETS}; run prepare_sets.py")

    rows: list[dict[str, object]] = []
    for set_path in set_paths:
        dataset = set_path.stem
        seqs, labels = read_set(set_path)
        cal = calibration_mask(seqs)
        if cal.sum() < 100 or (~cal).sum() < 100:
            raise RuntimeError(f"{dataset}: calibration/test split is unexpectedly small")
        print(f"{dataset}: {len(seqs)} total, {cal.sum()} calibration, {(~cal).sum()} test", flush=True)

        cnn_raw, cnn_speed = cnn_predictions(set_path, seqs)
        start = time.perf_counter()
        slope, intercept = affine_fit(labels[cal], cnn_raw[cal])
        cnn_cal_s = time.perf_counter() - start
        cnn_test = slope * cnn_raw[~cal] + intercept
        rows.append(result_row(
            dataset, "uno CNN", len(seqs), int(cal.sum()), int((~cal).sum()),
            metrics(labels[~cal], cnn_test), cnn_cal_s, cnn_speed,
        ))

        dlc_test, dlc_cal_s, dlc_speed = deeplc_predictions(seqs, labels, cal)
        rows.append(result_row(
            dataset, "DeepLC 3.1.13", len(seqs), int(cal.sum()), int((~cal).sum()),
            metrics(labels[~cal], dlc_test), dlc_cal_s, dlc_speed,
        ))

        # Per-peptide test predictions, so the figure can show error
        # distributions rather than only the summary statistics.
        test_seqs = [seq for seq, is_cal in zip(seqs, cal) if not is_cal]
        with (RESULTS / f"predictions_{dataset}.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sequence", "observed", "uno_cnn", "deeplc"])
            for seq, obs, a, b in zip(test_seqs, labels[~cal], cnn_test, dlc_test):
                writer.writerow([seq, f"{obs:.8g}", f"{a:.8g}", f"{b:.8g}"])

    fieldnames = list(rows[0])
    summary = RESULTS / "summary.tsv"
    with summary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: f"{value:.8g}" if isinstance(value, float) else value
                             for key, value in row.items()})

    import deeplc  # imported here only to record the installed package metadata

    environment = {
        "python": sys.version,
        "platform": platform.platform(),
        "processor": platform.processor(),
        "deeplc": getattr(deeplc, "__version__", "unknown"),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "model": str(MODEL.relative_to(ROOT)),
        "model_sha256": hashlib.sha256(MODEL.read_bytes()).hexdigest(),
        "rt_eval": str(RT_EVAL.relative_to(ROOT)),
        "git_head": subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True,
            text=True, capture_output=True,
        ).stdout.strip(),
        "thread_env": {key: os.environ.get(key) for key in (
            "RAYON_NUM_THREADS", "OMP_NUM_THREADS", "TF_NUM_INTRAOP_THREADS",
            "TF_NUM_INTEROP_THREADS",
        )},
    }
    (RESULTS / "environment.json").write_text(json.dumps(environment, indent=2) + "\n")
    print(f"wrote {summary}")


if __name__ == "__main__":
    main()
