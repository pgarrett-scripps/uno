#!/usr/bin/env python3
"""Matched frozen-versus-candidate search stability study (standard library only).

By default runs both binaries afresh on the same biosaur2 features, model, FASTA,
config, seeds and thread limits; verified completed baselines can be reused.
Archived variability-extension outputs remain read-only references. Entrapment is
used only after searches finish; every seed, zero and loss remains in summaries.
Concurrent elapsed times are provenance, not performance comparisons.

Example (prespecified low-input panel):
  python benchmark/stability_study.py --baseline /path/frozen/uno \
    --candidate /path/candidate/uno --output /tmp/uno-study-lowinput \
    --datasets lowinput --seeds 0-4 --jobs 2 --threads 2
Use a separate output directory for --datasets primary --seeds 0-2.
--dry-run validates and prints the plan without launching searches or writing.
--baseline-study DIR copies verified completed baseline runs into the new study;
candidate searches remain fresh. Seeds absent from a valid source study run a
fresh baseline search. Copied timings are historical, not measured here.
--summarize-only requires complete, hash-verified outputs. Repeating the exact
command resumes matching successful runs; changed provenance requires a new
output directory. Inherited UNO_* experimental knobs are cleared. Explicit
--env UNO_NAME=value applies identically to both binaries and is fingerprinted.
"""
from __future__ import annotations

import argparse
import csv
import fcntl
import hashlib
import itertools
import json
import math
import os
import shutil
from pathlib import Path
import statistics
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import entrapment

BENCH = Path(__file__).resolve().parent
ROOT = BENCH.parent
INVENTORY = BENCH / "results/variability_extension/summary.tsv"
SCHEMA = 1


def digest(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def atomic_json(path, value):
    partial = path.with_suffix(path.suffix + ".partial")
    partial.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    partial.replace(path)


def write_tsv(path, rows):
    if not rows:
        return
    partial = path.with_suffix(path.suffix + ".partial")
    with partial.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    partial.replace(path)


def numbers(values):
    result = set()
    for token in values:
        for part in token.split(","):
            if "-" in part:
                lo, hi = map(int, part.split("-", 1))
                if not 0 <= lo <= hi:
                    raise ValueError(f"invalid seed range: {part}")
                result.update(range(lo, hi + 1))
            else:
                result.add(int(part))
    if not result or min(result) < 0 or max(result) >= 2**64:
        raise ValueError("seeds must be unsigned 64-bit integers")
    return sorted(result)


def inventory():
    datasets, archives = {}, {}
    with INVENTORY.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["arm"] != "biosaur2":
                continue
            ds, stem = row["dataset"], row["astral_stem"]
            # Same feature-path conventions as variability_extension.specifications.
            feature = (Path.home() / f"Data/PXD066701_Astral/features/{stem}.features.tsv"
                       if stem else BENCH / f"data/{ds}.features.tsv")
            panel = ("lowinput" if ds.startswith(("astral_200pg_", "astral_1ng_"))
                     else "astral_highinput" if stem else "primary")
            spec = dict(dataset=ds, astral_stem=stem, panel=panel, features=str(feature.resolve()))
            if ds in datasets and datasets[ds] != spec:
                raise ValueError(f"inconsistent inventory for {ds}")
            datasets[ds] = spec
            key = ds, int(row["seed"])
            if key in archives:
                raise ValueError(f"duplicate inventory entry: {key}")
            path = Path(row["path"])
            archives[key] = str((path if path.is_absolute() else ROOT / path).resolve())
    return datasets, archives


def select_datasets(specs, requested):
    selected = set()
    for token in requested:
        for name in token.split(","):
            if name == "all":
                selected.update(specs)
            elif name in {"lowinput", "primary", "astral_highinput"}:
                selected.update(ds for ds, spec in specs.items() if spec["panel"] == name)
            elif name in specs:
                selected.add(name)
            else:
                raise ValueError(f"unknown dataset/panel {name!r}; choices: {', '.join(sorted(specs))}")
    if not selected:
        raise ValueError("no datasets selected")
    return [specs[ds] for ds in sorted(selected)]


def fasta_counts(path):
    human = yeast = 0
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                accession = line.split()[0]
                human += "_HUMAN" in accession
                yeast += "_YEAST" in accession
    if human == 0 or yeast == 0:
        raise ValueError("evaluation requires target human and yeast entries in the FASTA")
    return human, yeast


def evaluate(path):
    """Use the shared scorer; separately retain internal-FDR calibration data."""
    depths, _, _ = entrapment.curve(entrapment.load(path))
    human = yeast = other = decoy = row_count = 0
    ranked = []
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"accession", "decoy", "binomial_score", "q_value"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"invalid UNO output schema: {path}")
        for row in reader:
            row_count += 1
            score, q = float(row["binomial_score"]), float(row["q_value"])
            if not math.isfinite(score) or not math.isfinite(q) or not 0 <= q <= 1:
                raise ValueError(f"invalid score/q value in {path}")
            is_decoy = row["decoy"] not in ("0", "False", "false")
            accession = row["accession"]
            is_h, is_y = "_HUMAN" in accession, "_YEAST" in accession
            if q <= 0.01:
                if is_decoy:
                    decoy += 1
                elif is_y:
                    yeast += 1
                elif is_h:
                    human += 1
                else:
                    other += 1
            if not is_decoy and (is_h or is_y):
                ranked.append((score, is_y, accession))
    # The shared entrapment scorer orders yeast first within tied scores.
    # Adding accession as the last key only canonicalizes within-organism ties.
    ranked.sort(key=lambda r: (-r[0], not r[1], r[2]))
    h = y = 0
    fdps = []
    for _, is_y, _ in ranked:
        y += is_y
        h += not is_y
        fdps.append(y * entrapment.SCALE / h if h else math.inf)
    running, last = math.inf, -1
    for i in range(len(fdps) - 1, -1, -1):
        running = min(running, fdps[i])
        if running <= 0.01 and last < 0:
            last = i
    accepted = {acc for _, is_y, acc in ranked[:last + 1] if not is_y}
    if len(accepted) != depths[0.01]:
        raise ValueError(f"accepted-set/depth mismatch (duplicate accession?): {path}")
    metrics = dict(depth_1pct=depths[0.01], depth_5pct=depths[0.05], rows=row_count,
                   internal_q1_targets=human + yeast + other, internal_q1_human=human,
                   internal_q1_yeast=yeast, internal_q1_other=other, internal_q1_decoy=decoy,
                   internal_q1_entrapment_fdp=yeast * entrapment.SCALE / human if human else "")
    return metrics, accepted


def summarize(out, jobs, held_out):
    rows, sets = [], {}
    for job in jobs:
        path = Path(job["output"])
        metrics, accepted = evaluate(path)
        ds, variant, seed = job["dataset"], job["variant"], job["seed"]
        sets[variant, ds, seed] = accepted
        rows.append(dict(variant=variant, dataset=ds, panel=job["panel"], seed=seed,
                         seed_role="held_out" if seed in held_out else "development",
                         **metrics, path=str(path), archived_seed_path=job["archived_seed_path"]))
    rows.sort(key=lambda r: (r["dataset"], r["seed"], r["variant"]))
    write_tsv(out / "runs.tsv", rows)
    grouped = []
    for variant, ds in sorted({(r["variant"], r["dataset"]) for r in rows}):
        rs = [r for r in rows if (r["variant"], r["dataset"]) == (variant, ds)]
        group = dict(variant=variant, dataset=ds, panel=rs[0]["panel"], n=len(rs),
                     empty_output_count=sum(r["rows"] == 0 for r in rs),
                     zero_depth_1pct_count=sum(r["depth_1pct"] == 0 for r in rs),
                     zero_internal_q1_count=sum(r["internal_q1_targets"] == 0 for r in rs))
        for metric in ("depth_1pct", "depth_5pct", "internal_q1_human", "rows"):
            values = [r[metric] for r in rs]
            med, low, high = statistics.median(values), min(values), max(values)
            group.update({f"{metric}_median": med, f"{metric}_min": low,
                          f"{metric}_max": high, f"{metric}_range": high - low,
                          f"{metric}_range_over_median": (high - low) / med if med else ""})
        fdps = [r["internal_q1_entrapment_fdp"] for r in rs if r["internal_q1_entrapment_fdp"] != ""]
        group.update(internal_q1_fdp_defined_n=len(fdps),
                     internal_q1_fdp_median=statistics.median(fdps) if fdps else "",
                     internal_q1_fdp_max=max(fdps) if fdps else "")
        grouped.append(group)
    write_tsv(out / "variability.tsv", grouped)
    pairs, overlaps = [], []
    lookup = {(r["dataset"], r["seed"], r["variant"]): r for r in rows}
    for ds, seed in sorted({(r["dataset"], r["seed"]) for r in rows}):
        baseline, candidate = (lookup[ds, seed, v] for v in ("baseline", "candidate"))
        pair = dict(dataset=ds, panel=baseline["panel"], seed=seed, seed_role=baseline["seed_role"])
        for metric in ("depth_1pct", "depth_5pct", "rows", "internal_q1_targets", "internal_q1_human",
                       "internal_q1_yeast", "internal_q1_entrapment_fdp"):
            a, b = baseline[metric], candidate[metric]
            pair.update({f"baseline_{metric}": a, f"candidate_{metric}": b,
                         f"delta_{metric}": b - a if a != "" and b != "" else ""})
        a, b = (sets[v, ds, seed] for v in ("baseline", "candidate"))
        pair["accepted_1pct_jaccard"] = len(a & b) / len(a | b) if a | b else ""
        pairs.append(pair)
    write_tsv(out / "paired_deltas.tsv", pairs)
    for variant, ds in sorted({(r["variant"], r["dataset"]) for r in rows}):
        seeds = sorted(r["seed"] for r in rows if (r["variant"], r["dataset"]) == (variant, ds))
        for seed_a, seed_b in itertools.combinations(seeds, 2):
            a, b = sets[variant, ds, seed_a], sets[variant, ds, seed_b]
            overlaps.append(dict(variant=variant, dataset=ds, seed_a=seed_a, seed_b=seed_b,
                                 n_a=len(a), n_b=len(b), shared=len(a & b), union=len(a | b),
                                 jaccard_1pct=len(a & b) / len(a | b) if a | b else ""))
    write_tsv(out / "within_variant_overlap.tsv", overlaps)
    summary_paths = [out / name for name in ("runs.tsv", "variability.tsv", "paired_deltas.tsv", "within_variant_overlap.tsv")]
    atomic_json(out / "summary_provenance.json", dict(
        study_sha256=digest(out / "study.json"),
        outputs={job["output"]: digest(job["output"]) for job in jobs},
        summaries={p.name: digest(p) for p in summary_paths if p.exists()},
        interpretation="Ranges and pairwise overlaps are descriptive; seed-pair comparisons are not independent. Empty/empty overlap and FDP with zero accepted human proteins are undefined (blank)."))
    print(f"Summarized {len(rows)} matched searches in {out}", flush=True)


def scientific_recipe(command, fingerprints):
    """Normalize only the search-affecting recipe, allowing relocated inputs."""
    flags = ["--features", "--fasta", "--rt-model", "--config", "--seed", "--out-proteins"]
    if len(command) != 13 or command[1::2] != flags:
        raise ValueError("unsupported source search recipe; baseline reuse requires the same search arguments")
    try:
        recipe = {flag: fingerprints[command[i + 1]] for i, flag in zip(range(1, 9, 2), flags[:4])}
        recipe["binary"] = fingerprints[command[0]]
        recipe["seed"] = int(command[10])
    except (KeyError, ValueError) as error:
        raise ValueError("source study lacks a scientific input fingerprint") from error
    return recipe


def baseline_reuse(source, destination, jobs, protocol):
    """Verify a completed source study and bind selected baseline runs to hashes.

    Source scripts, inventory and unrelated candidate hashes may differ. The
    binary, scientific inputs, seed, applied settings, loader settings and cwd
    must match. Archived outputs never supply candidate results.
    """
    source = source.expanduser().resolve()
    if destination == source or destination.is_relative_to(source):
        raise ValueError("baseline reuse requires a new directory outside the source study")
    if (source / "incomplete.json").exists():
        raise ValueError(f"baseline source study is incomplete: {source}")
    study_path, summary_path = source / "study.json", source / "summary_provenance.json"
    old = json.loads(study_path.read_text())
    summary = json.loads(summary_path.read_text())
    study_sha = digest(study_path)
    if summary.get("study_sha256") != study_sha:
        raise ValueError("source study/summary hash mismatch")
    source_hashes = {str(study_path): study_sha, str(summary_path): digest(summary_path)}
    for name, value in summary.get("summaries", {}).items():
        path = (source / name).resolve()
        if not path.is_relative_to(source) or digest(path) != value:
            raise ValueError(f"source summary hash mismatch: {name}")
        source_hashes[str(path)] = value
    if "runs.tsv" not in summary.get("summaries", {}):
        raise ValueError("source study has no certified complete run summary")
    for field in ("environment", "inherited_runtime_environment", "cwd"):
        if old.get(field) != protocol[field]:
            raise ValueError(f"baseline source has different {field}")
    old_jobs = {}
    for job in old["searches"]:
        if job["variant"] == "baseline":
            key = job["dataset"], job["seed"]
            if key in old_jobs:
                raise ValueError(f"duplicate source baseline: {key}")
            old_jobs[key] = job
    completed = list(csv.DictReader((source / "runs.tsv").open(), delimiter="\t"))
    expected_keys = {(j["dataset"], j["variant"], j["seed"]) for j in old["searches"]}
    declared_keys = {(ds["dataset"], variant, seed) for ds in old["datasets"] for seed in old["seeds"] for variant in ("baseline", "candidate")}
    if expected_keys != declared_keys or len(expected_keys) != len(old["searches"]):
        raise ValueError("source study is missing or duplicates a declared search")
    actual_keys = [(r["dataset"], r["variant"], int(r["seed"])) for r in completed]
    if len(actual_keys) != len(set(actual_keys)) or set(actual_keys) != expected_keys:
        raise ValueError("source run summary does not contain every expected search")
    reused = {}
    for job in jobs:
        if job["variant"] != "baseline":
            continue
        key = job["dataset"], job["seed"]
        if key not in old_jobs:
            if job["dataset"] not in {ds for ds, _ in old_jobs}:
                raise ValueError(f"baseline source lacks requested dataset: {job['dataset']}")
            # A genuinely new seed is executed, never fabricated or silently
            # dropped. Missing runs declared by the source study fail above.
            continue
        prior = old_jobs[key]
        if scientific_recipe(prior["command"], old["fingerprints"]) != scientific_recipe(job["command"], protocol["fingerprints"]):
            raise ValueError(f"baseline scientific inputs/seed differ: {key}")
        path = Path(prior["output"]).resolve()
        if not path.is_relative_to(source):
            raise ValueError("source baseline output is outside the source study")
        stamp = path.with_suffix(".json")
        metadata = json.loads(stamp.read_text())
        if metadata.get("study_sha256") != study_sha or metadata.get("command") != prior["command"] or metadata.get("returncode") != 0:
            raise ValueError(f"source baseline run provenance mismatch: {path}")
        hashes = dict(source_hashes)
        hashes[str(stamp)] = digest(stamp)
        files = {"output": str(path), "stdout": str(path.with_suffix(".stdout.log")), "stderr": str(path.with_suffix(".stderr.log"))}
        for kind, filename in files.items():
            value = metadata.get(f"{kind}_sha256")
            if digest(filename) != value:
                raise ValueError(f"source baseline output/log hash mismatch: {filename}")
            hashes[filename] = value
        if summary.get("outputs", {}).get(str(path)) != hashes[str(path)]:
            raise ValueError(f"source summary/output hash mismatch: {path}")
        reused[key] = dict(study=str(source), files=files, source_hashes=hashes,
                           source_command=prior["command"], source_started_unix=metadata.get("started_unix"),
                           source_elapsed_seconds=metadata.get("elapsed_seconds"))
    return reused


def copy_baseline(reused, path, expected):
    # Recheck immediately before copying; source artifacts must remain immutable.
    if any(digest(name) != value for name, value in reused["source_hashes"].items()):
        raise ValueError("baseline source changed after validation")
    destinations = {"output": path, "stdout": path.with_suffix(".stdout.log"), "stderr": path.with_suffix(".stderr.log")}
    if any(p.exists() for p in destinations.values()):
        raise ValueError(f"will not overwrite baseline copy or log: {path}")
    for kind, target in destinations.items():
        source = reused["files"][kind]
        # Exclusive creation protects both archived data and existing new results.
        with open(source, "rb") as reader, target.open("xb") as writer:
            shutil.copyfileobj(reader, writer)
        if digest(target) != reused["source_hashes"][source]:
            raise ValueError(f"baseline changed while copying: {source}")
    metadata = dict(**expected, returncode=0, executed_command=None, reused_from=reused,
                    started_unix=None, elapsed_seconds=None, timing_measured=False, serialized=False,
                    timing_note="No search executed here; source timings are historical provenance, not measurements for this study.",
                    **{f"{kind}_sha256": digest(target) for kind, target in destinations.items()})
    atomic_json(path.with_suffix(".json"), metadata)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--baseline-study", type=Path, help="reuse verified baseline runs from this complete study; new seeds and all candidates run fresh")
    parser.add_argument("--datasets", nargs="+", default=["lowinput", "primary"], help="panel names or dataset IDs, space/comma separated")
    parser.add_argument("--seeds", nargs="+", default=["0-2"], help="integers, inclusive ranges or comma lists; e.g. 0-4 or 0,1,2")
    parser.add_argument("--held-out-seeds", nargs="+", default=["3-4"])
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--threads", type=int, default=2, help="Rayon threads per search; OpenMP/BLAS limits stay one")
    parser.add_argument("--fasta", type=Path, default=BENCH / "data/entrap_human_yeast.fasta")
    parser.add_argument("--model", type=Path, default=ROOT / "models/rt_cnn_chrono_global_v3_ch32.urtc")
    parser.add_argument("--config", type=Path, default=BENCH / "uno.toml")
    parser.add_argument("--env", action="append", default=[], metavar="UNO_NAME=value")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--summarize-only", action="store_true")
    args = parser.parse_args(argv)
    if args.jobs < 1 or args.threads < 1:
        parser.error("jobs and threads must be positive")
    seeds, held_out = numbers(args.seeds), numbers(args.held_out_seeds)
    specs, archives = inventory()
    selected = select_datasets(specs, args.datasets)
    binaries = {v: getattr(args, v).expanduser().resolve() for v in ("baseline", "candidate")}
    fasta, model, config, out = (getattr(args, key).expanduser().resolve() for key in ("fasta", "model", "config", "output"))
    if out == ROOT or out == BENCH or out == BENCH / "results" or out == INVENTORY.parent:
        raise ValueError("use a dedicated new study output directory")
    settings = {"RAYON_NUM_THREADS": str(args.threads), "OMP_NUM_THREADS": "1",
                "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "RUST_LOG": "info"}
    for item in args.env:
        name, separator, value = item.partition("=")
        if not separator or not name.startswith("UNO_"):
            parser.error("--env requires UNO_NAME=value")
        settings[name] = value
    runtime_env = {k: v for k, v in os.environ.items() if not k.startswith("UNO_")}
    runtime_env.update(settings)
    # Capture inherited loader/cache variables without recording unrelated secrets.
    inherited = {k: runtime_env[k] for k in ("LD_LIBRARY_PATH", "LD_PRELOAD", "XDG_CACHE_HOME") if k in runtime_env}
    input_paths = {Path(__file__).resolve(), Path(entrapment.__file__).resolve(), INVENTORY, fasta, model, config,
                   *binaries.values(), *(Path(spec["features"]) for spec in selected)}
    archive_paths = {Path(archives[spec["dataset"], seed]) for spec in selected for seed in seeds
                     if (spec["dataset"], seed) in archives}
    fingerprints = {str(path): digest(path) for path in sorted(input_paths | archive_paths)}
    for binary in binaries.values():
        if not os.access(binary, os.X_OK):
            raise ValueError(f"binary is not executable: {binary}")
    n_human, n_yeast = fasta_counts(fasta)
    entrapment.SCALE = n_human / n_yeast
    jobs = []
    for spec in selected:
        for seed in seeds:
            for variant in ("baseline", "candidate"):
                path = out / "runs" / f"{spec['dataset']}.seed{seed}.{variant}.tsv"
                command = [str(binaries[variant]), "--features", spec["features"], "--fasta", str(fasta),
                           "--rt-model", str(model), "--config", str(config), "--seed", str(seed),
                           "--out-proteins", str(path.with_suffix(".partial.tsv"))]
                jobs.append(dict(**spec, variant=variant, seed=seed, command=command, output=str(path),
                                 archived_seed_path=archives.get((spec["dataset"], seed), "")))
    protocol = dict(schema=SCHEMA, datasets=selected, seeds=seeds, held_out_seeds=held_out,
                    jobs=args.jobs, threads=args.threads, cwd=str(ROOT), environment=settings,
                    inherited_runtime_environment=inherited, fingerprints=fingerprints,
                    fasta_human_entries=n_human, fasta_yeast_entries=n_yeast, entrapment_scale=entrapment.SCALE,
                    searches=jobs, archived_outputs_are_reference_only=True)
    if args.baseline_study:
        reuse = baseline_reuse(args.baseline_study, out, jobs, protocol)
        protocol["baseline_study"] = str(args.baseline_study.expanduser().resolve())
        for job in jobs:
            key = job["dataset"], job["seed"]
            if job["variant"] == "baseline" and key in reuse:
                job["reused_baseline"] = reuse[key]
                fingerprints.update(job["reused_baseline"]["source_hashes"])
        # Bind the valid source even if every requested seed is new.
        source = args.baseline_study.expanduser().resolve()
        for name in ("study.json", "summary_provenance.json"):
            fingerprints[str(source / name)] = digest(source / name)
        print(f"Baseline reuse: {len(reuse)} verified copies; {len(jobs) // 2 - len(reuse)} fresh baseline searches; {len(jobs) // 2} fresh candidate searches", flush=True)
    result_label = "matched results" if args.baseline_study else "searches"
    print(f"Plan: {len(selected)} datasets × {len(seeds)} seeds × 2 binaries = {len(jobs)} {result_label}; "
          f"{args.jobs} jobs × {args.threads} Rayon threads", flush=True)
    if args.dry_run:
        for job in jobs:
            print(json.dumps(dict(variant=job["variant"], dataset=job["dataset"], seed=job["seed"], command=job["command"])))
        return 0
    out.mkdir(parents=True, exist_ok=True)
    with (out / ".study.lock").open("a+") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        study = out / "study.json"
        if study.exists():
            if json.loads(study.read_text()) != protocol:
                raise ValueError("study provenance changed; use a new output directory")
        else:
            if any(p.name != ".study.lock" for p in out.iterdir()):
                raise ValueError("unrecognized nonempty output directory; use a new study directory")
            if args.summarize_only:
                raise ValueError("no completed study to summarize")
            atomic_json(study, protocol)
        (out / "runs").mkdir(exist_ok=True)
        protocol_hash = digest(study)

        def run(job):
            path = Path(job["output"])
            stamp = path.with_suffix(".json")
            stdout, stderr = path.with_suffix(".stdout.log"), path.with_suffix(".stderr.log")
            expected = dict(study_sha256=protocol_hash, command=job["command"])
            if path.exists() or stamp.exists():
                if not (path.is_file() and stamp.is_file()):
                    raise ValueError(f"incomplete output/stamp pair: {path}; inspect before resuming")
                metadata = json.loads(stamp.read_text())
                if any(metadata.get(k) != v for k, v in expected.items()) or metadata.get("returncode") != 0:
                    raise ValueError(f"mismatched completed-run provenance: {path}")
                if "reused_baseline" in job and (metadata.get("reused_from") != job["reused_baseline"] or metadata.get("executed_command", "missing") is not None or metadata.get("timing_measured") is not False):
                    raise ValueError(f"baseline reuse metadata mismatch: {path}")
                for p, key in ((path, "output_sha256"), (stdout, "stdout_sha256"), (stderr, "stderr_sha256")):
                    if not p.is_file() or digest(p) != metadata.get(key):
                        raise ValueError(f"completed output/log changed: {p}")
                return f"reused {path.name}"
            if args.summarize_only:
                raise ValueError(f"missing completed search: {path}")
            if "reused_baseline" in job:
                copy_baseline(job["reused_baseline"], path, expected)
                return f"copied verified baseline {path.name}"
            # Preserve failed/interrupted logs before retrying a run with no final output.
            for old in (stdout, stderr, path.with_suffix(".failure.json")):
                if old.exists():
                    index = 1
                    while old.with_name(old.name + f".attempt{index}").exists():
                        index += 1
                    old.rename(old.with_name(old.name + f".attempt{index}"))
            partial = Path(job["command"][-1])
            if partial.exists():
                partial.unlink()
            started = time.time()
            with stdout.open("w") as so, stderr.open("w") as se:
                completed = subprocess.run(job["command"], cwd=ROOT, env=runtime_env, stdout=so, stderr=se, check=False)
            metadata = dict(**expected, returncode=completed.returncode, started_unix=started,
                            elapsed_seconds=time.time() - started, serialized=False,
                            stdout_sha256=digest(stdout), stderr_sha256=digest(stderr))
            if completed.returncode != 0 or not partial.is_file():
                atomic_json(path.with_suffix(".failure.json"), metadata)
                raise RuntimeError(f"search failed ({completed.returncode}): {path}; see {stderr}")
            evaluate(partial)  # Empty valid tables are outcomes; malformed tables are failures.
            partial.replace(path)
            metadata["output_sha256"] = digest(path)
            atomic_json(stamp, metadata)
            return f"completed {path.name}"

        errors = []
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = {pool.submit(run, job): job for job in jobs}
            for i, future in enumerate(as_completed(futures), 1):
                try:
                    print(f"{i}/{len(jobs)} {future.result()}", flush=True)
                except Exception as error:
                    errors.append(str(error))
                    print(f"{i}/{len(jobs)} ERROR: {error}", file=sys.stderr, flush=True)
        if errors:
            atomic_json(out / "incomplete.json", dict(errors=errors, expected_searches=len(jobs)))
            raise RuntimeError(f"{len(errors)} searches incomplete; no complete-study summaries written")
        # Fail closed if an input or binary changed while the study was running.
        if any(digest(path) != sha for path, sha in fingerprints.items()):
            raise ValueError("an input/binary changed during execution; results cannot be certified")
        summarize(out, jobs, held_out)
        (out / "incomplete.json").unlink(missing_ok=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, RuntimeError, OSError) as error:
        raise SystemExit(str(error))
