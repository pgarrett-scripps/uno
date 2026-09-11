#!/usr/bin/env python3
"""Measure search-seed variability for Astral and raw-mode comparisons.

Prespecified extension: seeds 0, 1, 2 for each of the 18 raw-mode files and
both feature sources. Seed 0 is the archived default; primary biosaur2 seeds
1 and 2 are reused from the existing 11-seed sweep. Other searches are run
with the current seeded binary, same model/config/database, two at a time.
Concurrent timings are recorded for provenance but excluded from resource
comparisons. Outputs and the complete summary stay inside benchmark/results.

Run with the paper analysis environment: uv run ../../benchmark/variability_extension.py
from paper/analysis. --summarize only validates and re-scores completed outputs.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import entrapment
from raw_mode.summarize import ASTRAL, ASTRAL_NAME, ASTRAL_STEMS, PRIMARY, accepted

BENCH = Path(__file__).resolve().parent
ROOT = BENCH.parent
OUT = BENCH / "results/variability_extension"
UNO = ROOT / "target/release/uno"
FASTA = BENCH / "data/entrap_human_yeast.fasta"
MODEL = ROOT / "models/rt_cnn_chrono_global_v3_ch32.urtc"
CONFIG = BENCH / "uno.toml"
SEEDS = (0, 1, 2)


def digest(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def specifications():
    raw_names = {line.split()[0]: Path(line.split()[1]).stem
                 for line in (BENCH / "raw_mode/urls.txt").read_text().splitlines()}
    specs = []
    for ds in PRIMARY:
        for arm in ("biosaur2", "raw"):
            source = (BENCH / f"data/{ds}.features.tsv" if arm == "biosaur2" else
                      Path.home() / f"Data/primary_raw/mzML/{raw_names[ds]}.mzML")
            base = BENCH / "results" / f"uno_{'' if arm == 'biosaur2' else 'raw_'}{ds}.tsv"
            specs.append((ds, "", arm, source, base))
    for stem in ASTRAL_STEMS:
        for arm in ("biosaur2", "raw"):
            source = (ASTRAL / f"features/{stem}.features.tsv" if arm == "biosaur2" else
                      ASTRAL / f"mzML/{stem}.mzML")
            sub = "uno_fixed" if arm == "biosaur2" else "uno_raw_fixed"
            specs.append((ASTRAL_NAME[stem], stem, arm, source, ASTRAL / f"results/{sub}/{stem}.tsv"))
    return specs


def output_path(ds, arm, seed, base):
    if seed == 0:
        return base
    if ds in PRIMARY and arm == "biosaur2":
        return BENCH / f"results/seed_sweep/uno_{ds}_seed{seed}.tsv"
    return OUT / f"{ds}.{arm}.seed{seed}.tsv"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summarize", action="store_true")
    parser.add_argument("--jobs", type=int, default=2)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    specs = specifications()
    for _, _, _, source, base in specs:
        if not source.is_file() or not base.is_file():
            raise SystemExit(f"missing input: {source} or {base}")
    fingerprints = {str(p): digest(p) for p in (UNO, FASTA, MODEL, CONFIG)}
    fingerprints.update({str(p): digest(p) for p in {s[3] for s in specs}})
    jobs = []
    for ds, stem, arm, source, base in specs:
        for seed in SEEDS:
            path = output_path(ds, arm, seed, base)
            if path.parent != OUT:
                if not path.is_file():
                    raise SystemExit(f"missing archived run: {path}")
                continue
            command = [str(UNO), "--features" if arm == "biosaur2" else "--raw", str(source),
                       "--fasta", str(FASTA), "--rt-model", str(MODEL), "--config", str(CONFIG),
                       "--seed", str(seed), "--out-proteins", str(path.with_suffix(".partial"))]
            inputs = {str(p): fingerprints[str(p)] for p in (UNO, FASTA, MODEL, CONFIG, source)}
            jobs.append((path, command, inputs))

    def run(job):
        path, command, inputs = job
        stamp = path.with_suffix(".json")
        if path.exists() and stamp.exists():
            metadata = json.loads(stamp.read_text())
            if metadata["inputs"] == inputs and metadata["output_sha256"] == digest(path):
                return f"reused {path.name}"
            raise RuntimeError(f"provenance changed for {path}; use a new output directory")
        if args.summarize:
            raise RuntimeError(f"incomplete run: {path}")
        with path.with_suffix(".log").open("w") as log, path.with_suffix(".time").open("w") as timing:
            subprocess.run(["/usr/bin/time", "-v", *command], stdout=log, stderr=timing,
                           cwd=ROOT, check=True)
        path.with_suffix(".partial").replace(path)
        stamp.write_text(json.dumps({"command": command, "inputs": inputs,
                                     "output_sha256": digest(path), "serialized": False}, indent=2) + "\n")
        return f"completed {path.name}"

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(run, job) for job in jobs]
        for i, future in enumerate(as_completed(futures), 1):
            print(f"{i}/{len(jobs)} {future.result()}", flush=True)
    rows, paths = [], []
    sets = {}
    for ds, stem, arm, source, base in specs:
        for seed in SEEDS:
            path = output_path(ds, arm, seed, base)
            depths, _, _ = entrapment.curve(entrapment.load(path))
            names = accepted(path, 0.01)
            assert len(names) == depths[0.01], (path, "accepted-list/depth mismatch")
            sets[ds, arm, seed] = names
            rows.append(dict(dataset=ds, astral_stem=stem, arm=arm, seed=seed,
                             depth_1pct=depths[0.01], depth_5pct=depths[0.05],
                             path=str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)))
            paths.append(path)
    for row in rows:
        key = row["dataset"], row["seed"]
        a, b = sets[key[0], "biosaur2", key[1]], sets[key[0], "raw", key[1]]
        # Two empty lists have no defined shared/union fraction. Retain the
        # search in the depth distribution, but do not invent zero agreement.
        row["jaccard_1pct"] = len(a & b) / len(a | b) if a | b else ""
    summary = OUT / "summary.tsv"
    with summary.with_suffix(".partial").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    summary.with_suffix(".partial").replace(summary)
    (OUT / "provenance.json").write_text(json.dumps({"seeds": SEEDS, "script_sha256": digest(__file__), "inputs": fingerprints,
        "outputs": {str(p): digest(p) for p in paths}, "summary_sha256": digest(summary)}, indent=2) + "\n")
    print(f"wrote {summary}: {len(rows)} rows", flush=True)


if __name__ == "__main__":
    main()
