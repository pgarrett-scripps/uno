#!/usr/bin/env python3
"""Rebuild every derived benchmark summary from the re-run protein tables.

Run after benchmark/publication_rerun.py. Nothing here searches anything: it
re-scores the tables on disk and rewrites the summaries the manuscript's figure
and table generators read.

  results/spread.tsv                     primary seed/rerun inventory (both engines)
  results/seed_sweep/summary.tsv         the uno seed sweep's own summary
  results/variability_extension/summary.tsv   18 files x 2 arms x 3 seeds
  results/summary.tsv                    serialized resource cohort (uno rows only;
                                         the ms1searchpy rows are left untouched,
                                         the comparator did not change)
  raw_mode/summary.tsv                   built-in feature detection vs biosaur2
  ~/Data/PXD066701_Astral/results/input_series_summary.tsv   Astral input series

The comparator's numbers are never recomputed from a new search here; only uno
rows are refreshed, and only from the tables the re-run produced.
"""
from __future__ import annotations

import csv
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

BENCH = Path(__file__).resolve().parent
ROOT = BENCH.parent
sys.path.insert(0, str(BENCH))
import entrapment  # noqa: E402

RESULTS = BENCH / "results"
RUN_DIR = RESULTS / "publication_rerun"
EXTENSION = RESULTS / "variability_extension"
SWEEP = RESULTS / "seed_sweep"
PRIMARY = ["hela5min_09746", "hela5min_09766", "hela5min_09734",
           "hela5min_09778", "hela5min_09950", "faims_LUM2_05433_VG"]


def run(*cmd: str, cwd: Path = ROOT) -> None:
    print(f"\n$ {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, cwd=cwd, check=True)


def sidecar_stamps() -> int:
    """Mirror the re-run's provenance next to the outputs that live elsewhere.

    benchmark/variability_extension.py validates its own directory by a stamp
    beside each output; the re-run writes one stamp per case in
    results/publication_rerun. Copy the ones that belong to the extension
    directory across so `--summarize` verifies against the current binary
    rather than against the superseded one.
    """
    n = 0
    for stamp in sorted(RUN_DIR.glob("*.json")):
        if stamp.name == "manifest.json":
            continue
        meta = json.loads(stamp.read_text())
        out = Path(meta["command"][-1])
        if out.parent != EXTENSION:
            continue
        shutil.copy2(stamp, out.with_suffix(".json"))
        n += 1
    return n


def seed_sweep_summary() -> None:
    rows = ["dataset\tseed\tdepth_1pct\tdepth_5pct"]
    for ds in PRIMARY:
        for seed in range(11):
            path = RESULTS / f"uno_{ds}.tsv" if seed == 0 else SWEEP / f"uno_{ds}_seed{seed}.tsv"
            if not path.exists():
                continue
            d, _, _ = entrapment.curve(entrapment.load(path))
            rows.append(f"{ds}\t{seed}\t{d[0.01]}\t{d[0.05]}")
    (SWEEP / "summary.tsv").write_text("\n".join(rows) + "\n")
    print(f"wrote {SWEEP / 'summary.tsv'} ({len(rows) - 1} rows)")


def resource_cohort() -> None:
    """Refresh the uno rows of results/summary.tsv from the serialized re-run.

    The ms1searchpy rows are copied through unchanged: the comparator was not
    rebuilt, so re-timing it would only add scheduler noise to a published
    number. uno's depth, CPU, wall and peak RSS all come from the seed-0
    searches, which the re-run executed one at a time.
    """
    path = RESULTS / "summary.tsv"
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    fields = list(rows[0])
    for row in rows:
        if row["tool"] != "uno":
            continue
        ds = row["dataset"]
        out, timing = RESULTS / f"uno_{ds}.tsv", RESULTS / f"uno_{ds}.time"
        if not (out.exists() and timing.exists()):
            raise SystemExit(f"missing re-run output for {ds}")
        d, _, _ = entrapment.curve(entrapment.load(out))
        text = timing.read_text()
        user = re.search(r"User time \(seconds\): ([\d.]+)", text).group(1)
        wall = [float(x) for x in re.search(
            r"Elapsed \(wall clock\) time.*: ([\d:.]+)", text).group(1).split(":")]
        rss = int(re.search(r"Maximum resident set size \(kbytes\): (\d+)", text).group(1))
        row.update(depth_1pct=d[0.01], depth_5pct=d[0.05], user_cpu_s=f"{float(user):.1f}",
                   wall_s=f"{wall[-1] + 60 * wall[-2] + (3600 * wall[-3] if len(wall) == 3 else 0):.1f}",
                   peak_rss_mb=f"{rss / 1024:.0f}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"refreshed the uno rows of {path}")


def main() -> int:
    print(f"copied {sidecar_stamps()} provenance stamps into {EXTENSION}")
    run(sys.executable, "benchmark/spread_table.py")
    seed_sweep_summary()
    run(sys.executable, "benchmark/variability_extension.py", "--summarize")
    resource_cohort()
    run(sys.executable, "benchmark/raw_mode/summarize.py")
    run(sys.executable, "paper/analysis/astral/summarize_input_series.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
