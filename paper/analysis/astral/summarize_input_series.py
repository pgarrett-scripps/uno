#!/usr/bin/env python3
"""Summarize the PXD066701 input-series subset (200 pg, 1, 10, 100 ng).

One row per search: uno (x1), stock ms1searchpy (x3; x10 at 10 ng), and
ms1searchpy with the external DIA-NN calibration list (`-pl`, x3). Depth is the
common entrapment scorer at 1% and 5%. Writes
$DATA/results/input_series_summary.tsv. Timing columns are copied from
`/usr/bin/time -v` but the 10 ng `-pl` runs were concurrent (see
RT_CALIBRATION.md); every other run here was serialized.
"""
from __future__ import annotations
import csv, re, sys
from pathlib import Path
PAPER = Path(__file__).resolve().parents[2]; ROOT = PAPER.parent
DATA = Path("/home/patrick-garrett/Data/PXD066701_Astral")
sys.path.insert(0, str(ROOT / "benchmark"))
import entrapment  # noqa: E402

LEVELS = {"200pg": 0.2, "1ng": 1.0, "10ng": 10.0, "100ng": 100.0}
# 2026-09-04: the paper reports three files per level, the first three by file
# name (README.md amendment). The two extra 1 ng files that had already been
# searched are written to input_series_dropped.tsv, not to the paper summary.
DROPPED = {"CW0710_OldFAIMS_1ng__RA7_rep5_ch1_DIA200pg", "CW0710_OldFAIMS_1ng__RA7_rep6_ch2_DIA200pg"}
STEMS = sorted(p.stem for p in (DATA / "features").glob("CW0710_OldFAIMS_*.features.tsv"))
STEMS = [s.replace(".features", "") for s in STEMS]

def level(stem):
    return re.search(r"OldFAIMS_(\d+(?:ng|pg))__", stem).group(1)

def resources(path):
    t = path.read_text()
    g = lambda pat: (re.search(pat, t) or [None, ""])[1]
    return {"wall": g(r"Elapsed \(wall clock\) time.*: ([0-9:.]+)"),
            "rss_kib": g(r"Maximum resident set size \(kbytes\): ([0-9]+)"),
            "rt_sigma_min": g(r"Second-stage calibrated RT sigma: ([0-9.]+)")}

def row(stem, engine, rep, output, timing):
    depth, _, _ = entrapment.curve(entrapment.load(output))
    return {"dataset": stem, "level": level(stem), "input_ng": LEVELS[level(stem)], "engine": engine, "rep": rep,
            "depth_1pct": depth[0.01], "depth_5pct": depth[0.05], **resources(timing)}

rows = []
for stem in STEMS:
    if level(stem) not in LEVELS: continue
    # "uno" = the fixed binary (results/uno_fixed, low-input collapse retry, 2026-09-04);
    # "uno-frozen" = the binary the series was preregistered with (results/uno).
    for eng, sub, log in (("uno", "uno_fixed", "uno_fixed"), ("uno-frozen", "uno", "uno")):
        uno = DATA / "results" / sub / f"{stem}.tsv"; t = DATA / "logs" / f"{stem}.{log}.time"
        if uno.exists() and t.exists(): rows.append(row(stem, eng, 1, uno, t))
    for d in sorted((DATA / "results" / "ms1" / stem).glob("rep*")):
        rep = int(d.name[3:]); out = d / f"{stem}.features_proteins_full.tsv"
        t = DATA / "logs" / f"{stem}.ms1.rep{rep:02d}.time"
        if not t.exists(): t = DATA / "logs" / f"{stem}.ms1.rep{rep}.time"
        if out.exists() and t.exists(): rows.append(row(stem, "ms1searchpy", rep, out, t))
    pl_dirs = list(sorted((DATA / "results" / "ms1_pl" / stem).glob("rep*"))) or \
              list(sorted((DATA / "results" / "rt_calibration" / "pl" / stem).glob("rep*")))
    for d in pl_dirs:
        rep = int(d.name[3:]); out = d / f"{stem}.features_proteins_full.tsv"
        t = DATA / "logs" / f"{stem}.ms1pl.rep{rep}.time"
        if not t.exists(): t = d / "stderr.log"
        if out.exists() and t.exists(): rows.append(row(stem, "ms1searchpy-pl", rep, out, t))
for name, keep in (("input_series_summary.tsv", lambda r: r["dataset"] not in DROPPED),
                   ("input_series_dropped.tsv", lambda r: r["dataset"] in DROPPED)):
    out = DATA / "results" / name; sel = [r for r in rows if keep(r)]
    with out.open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(sel)
    print(f"wrote {out} ({len(sel)} rows)")
