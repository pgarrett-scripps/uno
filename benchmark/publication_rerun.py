#!/usr/bin/env python3
"""Re-execute every uno search the manuscript reports, on one pinned binary.

Why this exists
---------------
Every uno protein table behind the Technical Note was produced by build
``574df1fd`` (2026-09-09). Three engine changes have landed since: `--seed`,
V3 candidate selection, and bounded RT recovery (docs/stability-study.md,
docs/bounded-rt-recovery.md). Both engine changes are documented to leave
successful runs byte-identical and to move only runs that previously returned
nothing -- but that was verified on the study panels, not on the manuscript's
own file/seed inventory. So the inventory is re-executed here rather than
argued about, and this script records which outputs actually changed.

The inventory (156 searches, the complete uno side of the paper)
---------------------------------------------------------------
  * 6 primary DirectMS1 files x biosaur2 features x seeds 0-10   (66)
  * 6 primary files x ``--raw`` (in-process koth) x seeds 0-2    (18)
  * 12 PXD066701 Astral files x {biosaur2, raw} x seeds 0-2      (72)

Outputs are written to the canonical paths the figure scripts already read, so
nothing downstream needs repointing; the superseded tables are copied into
``results/archive_<old build>/`` first, with their hashes.

Seed-0 searches are the ones whose ``/usr/bin/time -v`` output feeds the
resource panels, so they run **serialized** (phase A). The remaining seeds run
two at a time (phase B); their timings are recorded but flagged concurrent, as
in benchmark/variability_extension.py.

Usage (from the repository root)::

    python3 benchmark/publication_rerun.py --plan        # inventory only
    python3 benchmark/publication_rerun.py               # run everything
    python3 benchmark/publication_rerun.py --summarize   # verify + re-score

Completed cases are skipped when their provenance stamp matches this binary and
their output hash still agrees, so an interrupted run can simply be restarted.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path

BENCH = Path(__file__).resolve().parent
ROOT = BENCH.parent
sys.path.insert(0, str(BENCH))
sys.path.insert(0, str(BENCH / "raw_mode"))

import entrapment  # noqa: E402
from summarize import ASTRAL, ASTRAL_NAME, ASTRAL_STEMS, PRIMARY  # noqa: E402

UNO = ROOT / "target/release/uno"
FASTA = BENCH / "data/entrap_human_yeast.fasta"
MODEL = ROOT / "models/rt_cnn_chrono_global_v3_ch32.urtc"
CONFIG = BENCH / "uno.toml"
RESULTS = BENCH / "results"
EXTENSION = RESULTS / "variability_extension"
SWEEP = RESULTS / "seed_sweep"
RUN_DIR = RESULTS / "publication_rerun"
PRIMARY_MZML = Path.home() / "Data/primary_raw/mzML"

PRIMARY_SEEDS = tuple(range(11))      # biosaur2 arm: the 11-draw spread
EXTENSION_SEEDS = (0, 1, 2)           # raw arm and every Astral file


def digest(path: Path) -> str:
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def raw_mzml(dataset: str) -> Path:
    """Primary-dataset mzML, mapped through raw_mode/urls.txt as run_primary.sh does."""
    for line in (BENCH / "raw_mode/urls.txt").read_text().splitlines():
        parts = line.split()
        if parts and parts[0] == dataset:
            return PRIMARY_MZML / f"{Path(parts[1]).stem}.mzML"
    raise KeyError(dataset)


@dataclass
class Case:
    dataset: str          # paper-facing name, e.g. hela5min_09746 / astral_1ng_rep1
    stem: str             # Astral file stem, "" for the primary files
    arm: str              # biosaur2 | raw
    seed: int
    source: Path          # features .tsv or mzML
    out: Path             # protein table, at the path the figure scripts read
    timing: Path          # /usr/bin/time -v capture
    log: Path
    serialized: bool
    old: dict = field(default_factory=dict)

    @property
    def key(self) -> str:
        return f"{self.dataset}.{self.arm}.seed{self.seed}"

    @property
    def stamp(self) -> Path:
        return RUN_DIR / f"{self.key}.json"

    def command(self) -> list[str]:
        return [str(UNO),
                "--features" if self.arm == "biosaur2" else "--raw", str(self.source),
                "--fasta", str(FASTA), "--rt-model", str(MODEL), "--config", str(CONFIG),
                *(["--seed", str(self.seed)] if self.seed else []),
                "--out-proteins", str(self.out)]


def inventory() -> list[Case]:
    cases: list[Case] = []
    for ds in PRIMARY:
        for seed in PRIMARY_SEEDS:
            cases.append(Case(
                ds, "", "biosaur2", seed, BENCH / f"data/{ds}.features.tsv",
                RESULTS / f"uno_{ds}.tsv" if seed == 0 else SWEEP / f"uno_{ds}_seed{seed}.tsv",
                RESULTS / f"uno_{ds}.time" if seed == 0 else SWEEP / f"uno_{ds}_seed{seed}.time",
                RESULTS / f"uno_{ds}.log" if seed == 0 else SWEEP / f"uno_{ds}_seed{seed}.log",
                seed == 0))
        for seed in EXTENSION_SEEDS:
            cases.append(Case(
                ds, "", "raw", seed, raw_mzml(ds),
                RESULTS / f"uno_raw_{ds}.tsv" if seed == 0 else EXTENSION / f"{ds}.raw.seed{seed}.tsv",
                RESULTS / f"{ds}.uno_raw.time" if seed == 0 else EXTENSION / f"{ds}.raw.seed{seed}.time",
                RESULTS / f"{ds}.uno_raw.log" if seed == 0 else EXTENSION / f"{ds}.raw.seed{seed}.log",
                seed == 0))
    for stem in ASTRAL_STEMS:
        name = ASTRAL_NAME[stem]
        for arm, sub, tag in (("biosaur2", "uno_fixed", "uno_fixed"),
                              ("raw", "uno_raw_fixed", "uno_raw_fixed")):
            source = (ASTRAL / f"features/{stem}.features.tsv" if arm == "biosaur2"
                      else ASTRAL / f"mzML/{stem}.mzML")
            for seed in EXTENSION_SEEDS:
                cases.append(Case(
                    name, stem, arm, seed, source,
                    ASTRAL / f"results/{sub}/{stem}.tsv" if seed == 0
                    else EXTENSION / f"{name}.{arm}.seed{seed}.tsv",
                    ASTRAL / f"logs/{stem}.{tag}.time" if seed == 0
                    else EXTENSION / f"{name}.{arm}.seed{seed}.time",
                    ASTRAL / f"logs/{stem}.{tag}.log" if seed == 0
                    else EXTENSION / f"{name}.{arm}.seed{seed}.log",
                    seed == 0))
    return cases


def commit() -> str:
    return subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                          capture_output=True, text=True, check=True).stdout.strip()


def fingerprints(cases: list[Case]) -> dict[str, str]:
    paths = {UNO, FASTA, MODEL, CONFIG} | {c.source for c in cases}
    return {str(p): digest(p) for p in sorted(paths)}


def depths(path: Path) -> tuple[int, int]:
    out, _, _ = entrapment.curve(entrapment.load(path))
    return out[0.01], out[0.05]


def archive(cases: list[Case], build: str) -> Path:
    """Copy every superseded output out of the way, recording its hash.

    Keyed by the *old* build so a second re-run cannot silently overwrite the
    first archive. Files that live outside the repository (the Astral trees)
    are archived here too, under their absolute path with separators flattened.
    """
    dest = RESULTS / f"archive_{build[:8]}"
    dest.mkdir(parents=True, exist_ok=True)
    manifest_path = dest / "manifest.json"
    if manifest_path.exists():
        return dest
    manifest = {}
    for case in cases:
        entry = {}
        for role in ("out", "timing", "log"):
            src = getattr(case, role)
            if not src.exists():
                continue
            rel = (src.relative_to(ROOT) if src.is_relative_to(ROOT)
                   else Path("external") / str(src).lstrip("/").replace("/", "__"))
            target = dest / rel
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, target)
            entry[role] = {"path": str(src), "archived": str(target.relative_to(dest))}
            if role == "out":
                entry[role]["sha256"] = digest(src)
                d1, d5 = depths(src)
                entry[role]["depth_1pct"], entry[role]["depth_5pct"] = d1, d5
        if entry:
            manifest[case.key] = entry
    manifest_path.write_text(json.dumps(
        {"superseded_build_sha256": build, "archived_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
         "cases": manifest}, indent=2) + "\n")
    print(f"archived {len(manifest)} previous cases into {dest}")
    return dest


def load_old(cases: list[Case], dest: Path) -> None:
    manifest = json.loads((dest / "manifest.json").read_text())["cases"]
    for case in cases:
        case.old = manifest.get(case.key, {}).get("out", {})


def execute(case: Case, inputs: dict[str, str], head: str, summarize: bool) -> str:
    if case.stamp.exists() and case.out.exists():
        meta = json.loads(case.stamp.read_text())
        if meta["inputs"] == inputs and meta["output_sha256"] == digest(case.out):
            return f"reused  {case.key}"
        if summarize:
            raise SystemExit(f"stale output for {case.key}: rerun without --summarize")
    if summarize:
        raise SystemExit(f"missing run: {case.key}")
    case.out.parent.mkdir(parents=True, exist_ok=True)
    case.timing.parent.mkdir(parents=True, exist_ok=True)
    partial = case.out.with_suffix(".partial")
    with case.log.open("w") as log, case.timing.open("w") as timing:
        subprocess.run(["/usr/bin/time", "-v", *[a.replace(str(case.out), str(partial))
                                                 for a in case.command()]],
                       stdout=log, stderr=timing, cwd=ROOT, check=True)
    partial.replace(case.out)
    d1, d5 = depths(case.out)
    case.stamp.write_text(json.dumps(
        {"command": case.command(), "commit": head, "inputs": inputs,
         "output_sha256": digest(case.out), "serialized": case.serialized,
         "depth_1pct": d1, "depth_5pct": d5}, indent=2) + "\n")
    return f"ran     {case.key}  @1%={d1} @5%={d5}"


def report(cases: list[Case], head: str, build: str) -> None:
    rows = []
    for case in cases:
        meta = json.loads(case.stamp.read_text())
        old = case.old
        rows.append({
            "dataset": case.dataset, "astral_stem": case.stem, "arm": case.arm,
            "seed": case.seed, "serialized": int(case.serialized),
            "depth_1pct": meta["depth_1pct"], "depth_5pct": meta["depth_5pct"],
            "prev_depth_1pct": old.get("depth_1pct", ""), "prev_depth_5pct": old.get("depth_5pct", ""),
            "identical_to_prev": int(bool(old) and old.get("sha256") == meta["output_sha256"]),
            "had_prev": int(bool(old)),
            "sha256": meta["output_sha256"],
            "path": str(case.out.relative_to(ROOT)) if case.out.is_relative_to(ROOT) else str(case.out),
        })
    out = RUN_DIR / "cases.tsv"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    changed = [r for r in rows if r["had_prev"] and not r["identical_to_prev"]]
    (RUN_DIR / "manifest.json").write_text(json.dumps(
        {"commit": head, "uno_sha256": build, "superseded_build": "574df1fd188b3f2bf9e3baa610f25461b9faf92fa564a2a25cd5e661ce9c05e7",
         "cases": len(rows), "unchanged": sum(r["identical_to_prev"] for r in rows),
         "changed": [r["dataset"] + "." + r["arm"] + ".seed" + str(r["seed"]) for r in changed],
         "finished_at": time.strftime("%Y-%m-%dT%H:%M:%S%z")}, indent=2) + "\n")
    print(f"\nwrote {out}: {len(rows)} cases, "
          f"{sum(r['identical_to_prev'] for r in rows)} byte-identical to the superseded build")
    if changed:
        print(f"{len(changed)} changed:")
        for r in changed:
            print(f"  {r['dataset']:<20} {r['arm']:<9} seed{r['seed']:<3} "
                  f"1%: {r['prev_depth_1pct']:>5} -> {r['depth_1pct']:<5} "
                  f"5%: {r['prev_depth_5pct']:>5} -> {r['depth_5pct']}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", action="store_true", help="print the inventory and exit")
    parser.add_argument("--summarize", action="store_true", help="verify completed runs only")
    parser.add_argument("--jobs", type=int, default=2, help="concurrency for the seeded phase")
    parser.add_argument("--only", default="", help="substring filter on <dataset>.<arm>.seed<N>")
    args = parser.parse_args()

    cases = inventory()
    if args.only:
        cases = [c for c in cases if args.only in c.key]
    missing = sorted({str(c.source) for c in cases if not c.source.is_file()})
    if missing:
        raise SystemExit("missing inputs:\n  " + "\n  ".join(missing))
    if args.plan:
        for c in cases:
            print(f"{c.key:<40} {'serial' if c.serialized else 'pooled'}  {c.source}")
        print(f"{len(cases)} searches "
              f"({sum(c.serialized for c in cases)} serialized, {sum(not c.serialized for c in cases)} pooled)")
        return 0

    RUN_DIR.mkdir(parents=True, exist_ok=True)
    build, head = digest(UNO), commit()
    print(f"uno {build[:12]} at {head[:12]}; {len(cases)} searches")
    dest = archive(cases, "574df1fd188b3f2bf9e3baa610f25461b9faf92fa564a2a25cd5e661ce9c05e7")
    load_old(cases, dest)
    prints = fingerprints(cases)
    jobs = {c.key: {str(p): prints[str(p)] for p in (UNO, FASTA, MODEL, CONFIG, c.source)} for c in cases}

    serial = [c for c in cases if c.serialized]
    pooled = [c for c in cases if not c.serialized]
    for i, case in enumerate(serial, 1):
        print(f"[A {i}/{len(serial)}] {execute(case, jobs[case.key], head, args.summarize)}", flush=True)
    if pooled:
        with ThreadPoolExecutor(max_workers=1 if args.summarize else args.jobs) as pool:
            futures = {pool.submit(execute, c, jobs[c.key], head, args.summarize): c for c in pooled}
            for i, future in enumerate(as_completed(futures), 1):
                print(f"[B {i}/{len(pooled)}] {future.result()}", flush=True)
    report(cases, head, build)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
