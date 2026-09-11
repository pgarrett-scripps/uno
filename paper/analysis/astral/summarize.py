#!/usr/bin/env python3
"""Summarize the external Astral benchmark from primary search outputs."""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

PAPER = Path(__file__).resolve().parents[2]
ROOT = PAPER.parent
DATA = Path("/home/patrick-garrett/Data/PXD066701_Astral")
sys.path.insert(0, str(ROOT / "benchmark"))
sys.path.insert(0, str(PAPER / "analysis" / "scripts"))

import entrapment  # noqa: E402
from entrapment_data import ms1_reported_q  # noqa: E402

STEMS = (
    "CW0710_OldFAIMS_10ng__RA4_rep1_ch1_DIA10ng",
    "CW0710_OldFAIMS_10ng__RA4_rep4_ch2_DIA10ng",
    "CW0710_OldFAIMS_10ng__RA4_rep6_ch2_DIA10ng",
)


def resources(path: Path):
    text = path.read_text()
    pats = {
        "user_cpu_s": r"User time \(seconds\): ([0-9.]+)",
        "wall": r"Elapsed \(wall clock\) time.*: ([0-9:.]+)",
        "rss_kib": r"Maximum resident set size \(kbytes\): ([0-9]+)",
    }
    out = {}
    for key, pat in pats.items():
        hit = re.search(pat, text)
        out[key] = hit.group(1) if hit else ""
    return out


def target_rows(path: Path):
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if str(row.get("decoy", "0")).lower() not in ("0", "false"):
                continue
            accession = row.get("accession", row.get("dbname", ""))
            if "_HUMAN" not in accession and "_YEAST" not in accession:
                continue
            yield row, accession


def reported_fdp(path: Path, engine: str, threshold: float):
    step = ms1_reported_q(path) if engine.startswith("ms1searchpy") else None
    human = yeast = 0
    for row, accession in target_rows(path):
        if engine == "uno":
            q = float(row["q_value"])
        else:
            q = step(float(row["score"]))
        if q <= threshold:
            if "_YEAST" in accession:
                yeast += 1
            else:
                human += 1
    fdp = yeast * entrapment.SCALE / human if human else 0.0
    return human, yeast, fdp


def one_row(stem: str, engine: str, rep: int, output: Path, timing: Path):
    depth, _, _ = entrapment.curve(entrapment.load(output))
    h1, y1, fdp1 = reported_fdp(output, engine, 0.01)
    h5, y5, fdp5 = reported_fdp(output, engine, 0.05)
    res = resources(timing)
    features = DATA / "features" / f"{stem}.features.tsv"
    with features.open() as handle:
        feature_count = sum(1 for _ in handle) - 1
    return {
        "dataset": stem,
        "engine": engine,
        "rep": rep,
        "features": feature_count,
        "depth_1pct_entrapment_q": depth[0.01],
        "depth_5pct_entrapment_q": depth[0.05],
        "reported_1pct_human": h1,
        "reported_1pct_yeast": y1,
        "reported_1pct_raw_fdp": f"{fdp1:.6f}",
        "reported_5pct_human": h5,
        "reported_5pct_yeast": y5,
        "reported_5pct_raw_fdp": f"{fdp5:.6f}",
        **res,
    }


def main():
    rows = []
    for stem in STEMS:
        uno = DATA / "results" / "uno" / f"{stem}.tsv"
        uno_time = DATA / "logs" / f"{stem}.uno.time"
        if uno.exists() and uno_time.exists():
            rows.append(one_row(stem, "uno", 1, uno, uno_time))
        for rep_dir in sorted((DATA / "results" / "ms1" / stem).glob("rep*")):
            rep = int(rep_dir.name.removeprefix("rep"))
            output = rep_dir / f"{stem}.features_proteins_full.tsv"
            timing = DATA / "logs" / f"{stem}.ms1.rep{rep:02d}.time"
            if output.exists() and timing.exists():
                rows.append(one_row(stem, "ms1searchpy", rep, output, timing))
        # External-RT-calibration arm: ms1searchpy with `-pl` on the DIA-NN
        # MS/MS list (analysis/astral/RT_CALIBRATION.md). These runs were
        # launched three-way concurrent, so wall/RSS are recorded but are not
        # comparable to the serialized production timings.
        pl_root = DATA / "results" / "rt_calibration" / "pl" / stem
        for rep_dir in sorted(pl_root.glob("rep*")) if pl_root.exists() else []:
            rep = int(rep_dir.name.removeprefix("rep"))
            output = rep_dir / f"{stem}.features_proteins_full.tsv"
            timing = rep_dir / "stderr.log"
            if output.exists() and timing.exists():
                rows.append(one_row(stem, "ms1searchpy-pl", rep, output, timing))
    if not rows:
        raise SystemExit("no completed production results")
    out = DATA / "results" / "astral_summary.tsv"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {out} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
