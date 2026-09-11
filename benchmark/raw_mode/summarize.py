#!/usr/bin/env python3
"""Built-in feature detection (`uno --raw`) versus the biosaur2 feature tables.

One row per (dataset, arm). The biosaur2 arm is the manuscript's primary uno
search (benchmark/run.sh, ~/Data/PXD066701_Astral for the Astral files); the
raw arm is `uno --raw <mzML>` with the bundled koth profile, same FASTA, RT
model and config (run_primary.sh, paper/analysis/astral/run_raw_mode.sh).

Depth is scored by benchmark/entrapment.py, the scorer used everywhere else.
`shared_1pct` is the number of human proteins accepted at 1% entrapment FDR by
BOTH arms; `jaccard_1pct` is shared / union. Wall and peak RSS come from
`/usr/bin/time -v`; the raw arm's wall time includes feature detection.

Writes benchmark/raw_mode/summary.tsv (tracked); protein tables stay on disk.
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
BENCH = HERE.parent
sys.path.insert(0, str(BENCH))
import entrapment  # noqa: E402

ASTRAL = Path.home() / "Data/PXD066701_Astral"
PRIMARY = [
    "hela5min_09746", "hela5min_09766", "hela5min_09734",
    "hela5min_09778", "hela5min_09950", "faims_LUM2_05433_VG",
]
# The 12-file PXD066701 input series (paper/analysis/astral/README.md):
# 200 pg x3, 1 ng x3, 10 ng x3, 100 ng x3, first three by file name per level. Both arms use the fixed uno binary
# (low-input collapse retry, 2026-09-04): results/uno_fixed and results/uno_raw_fixed.
ASTRAL_STEMS = [
    "CW0710_OldFAIMS_200pg__RA9_rep1_ch1_DIA200pg",
    "CW0710_OldFAIMS_200pg__RA9_rep2_ch2_DIA200pg",
    "CW0710_OldFAIMS_200pg__RA9_rep3_ch1_DIA200pg",
    "CW0710_OldFAIMS_1ng__RA7_rep1_ch1_DIA200pg",
    "CW0710_OldFAIMS_1ng__RA7_rep2_ch2_DIA200pg",
    "CW0710_OldFAIMS_1ng__RA7_rep4_ch2_DIA200pg",
    "CW0710_OldFAIMS_10ng__RA4_rep1_ch1_DIA10ng",
    "CW0710_OldFAIMS_10ng__RA4_rep4_ch2_DIA10ng",
    "CW0710_OldFAIMS_10ng__RA4_rep6_ch2_DIA10ng",
    "CW0710_OldFAIMS_100ng__RA1_rep2_ch2_DIA10ng",
    "CW0710_OldFAIMS_100ng__RA1_rep3_ch1_DIA10ng",
    "CW0710_OldFAIMS_100ng__RA1_rep4_ch2_DIA10ng",
]
def _astral_name(stem: str) -> str:
    level = stem.split("OldFAIMS_")[1].split("__")[0]
    rep = stem.split("_rep")[1].split("_")[0]
    return f"astral_{level}_rep{rep}"


ASTRAL_NAME = {s: _astral_name(s) for s in ASTRAL_STEMS}


def accepted(path: Path, thr: float) -> set[str]:
    """Human accessions accepted at entrapment q <= thr (entrapment.py rule)."""
    rows = []
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if str(row.get("decoy", "0")) not in ("0", "False", "false"):
                continue
            name = row["accession"]
            if "_YEAST" not in name and "_HUMAN" not in name:
                continue
            rows.append((float(row["binomial_score"]), "_YEAST" in name, name))
    rows.sort(key=lambda r: (-r[0], not r[1]))
    human = yeast = 0
    fdp = []
    for _, is_y, _ in rows:
        yeast += is_y
        human += not is_y
        fdp.append(yeast * entrapment.SCALE / human if human else 0.0)
    q = [0.0] * len(fdp)
    run = float("inf")
    for i in range(len(fdp) - 1, -1, -1):
        run = min(run, fdp[i])
        q[i] = run
    last = max((i for i, v in enumerate(q) if v <= thr), default=-1)
    return {n for _, is_y, n in rows[: last + 1] if not is_y}


def timing(path: Path) -> dict[str, str]:
    text = path.read_text()
    out = {}
    m = re.search(r"Elapsed \(wall clock\) time.*: ([0-9:.]+)", text)
    if m:
        parts = [float(x) for x in m.group(1).split(":")]
        out["wall_s"] = f"{parts[-1] + 60 * parts[-2] + (3600 * parts[-3] if len(parts) == 3 else 0):.1f}"
    m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
    if m:
        out["peak_rss_mb"] = f"{int(m.group(1)) / 1024:.0f}"
    m = re.search(r"loaded (\d+) features", text)
    out["features"] = m.group(1) if m else ""
    return out


def feature_count_tsv(path: Path) -> str:
    with path.open() as fh:
        return str(sum(1 for _ in fh) - 1)


def pair(dataset, bio_tsv, bio_time, bio_features, raw_tsv, raw_time):
    for p in (bio_tsv, raw_tsv, raw_time):
        if not p.exists():
            return []
    d_bio, _, _ = entrapment.curve(entrapment.load(bio_tsv))
    d_raw, _, _ = entrapment.curve(entrapment.load(raw_tsv))
    a1, b1 = accepted(bio_tsv, 0.01), accepted(raw_tsv, 0.01)
    shared, union = len(a1 & b1), len(a1 | b1)
    t_bio = timing(bio_time) if bio_time.exists() else {}
    t_raw = timing(raw_time)
    rows = []
    for arm, depth, t, feats in (
        ("biosaur2", d_bio, t_bio, feature_count_tsv(bio_features)),
        ("raw", d_raw, t_raw, t_raw.get("features", "")),
    ):
        rows.append({
            "dataset": dataset, "arm": arm, "features": feats,
            "depth_1pct": depth[0.01], "depth_5pct": depth[0.05],
            "shared_1pct": shared, "union_1pct": union,
            "jaccard_1pct": f"{shared / union:.4f}" if union else "",
            "wall_s": t.get("wall_s", ""), "peak_rss_mb": t.get("peak_rss_mb", ""),
        })
    return rows


def main() -> int:
    res = BENCH / "results"
    rows = []
    for ds in PRIMARY:
        rows += pair(ds, res / f"uno_{ds}.tsv", res / f"uno_{ds}.time",
                     BENCH / "data" / f"{ds}.features.tsv",
                     res / f"uno_raw_{ds}.tsv", res / f"{ds}.uno_raw.time")
    for stem in ASTRAL_STEMS:
        rows += pair(ASTRAL_NAME[stem], ASTRAL / "results/uno_fixed" / f"{stem}.tsv",
                     ASTRAL / "logs" / f"{stem}.uno_fixed.time",
                     ASTRAL / "features" / f"{stem}.features.tsv",
                     ASTRAL / "results/uno_raw_fixed" / f"{stem}.tsv",
                     ASTRAL / "logs" / f"{stem}.uno_raw_fixed.time")
    if not rows:
        raise SystemExit("no completed raw-mode pairs found")
    out = HERE / "summary.tsv"
    with out.open("w", newline="") as fh:
        fh.write("# benchmark/raw_mode/summary.tsv -- generated by benchmark/raw_mode/summarize.py\n")
        fh.write("# biosaur2 arm: benchmark/results/uno_<ds>.tsv (run.sh) and ~/Data/PXD066701_Astral/results/uno_fixed/\n")
        fh.write("# raw arm: benchmark/results/uno_raw_<ds>.tsv (run_primary.sh) and ~/Data/PXD066701_Astral/results/uno_raw_fixed/\n")
        fh.write("# depth = human proteins at entrapment q (benchmark/entrapment.py); shared/union over the 1% lists\n")
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out} ({len(rows)} rows, {len(rows) // 2} datasets)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
