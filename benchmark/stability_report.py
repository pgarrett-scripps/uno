#!/usr/bin/env python3
"""Report completed stability_study.py panels without launching any searches.

Usage: python benchmark/stability_report.py --studies RESULTS/lowinput RESULTS/primary
       --output RESULTS/report
Writes report.json and report.md. Verifies study, summary and protein-file hashes.
Use --lowinput-seed-count 7 --protect-discoveries --protocol V2_PLAN for the
prespecified v2 panel. Defaults retain the phase-1 five-seed acceptance rule.
Blank/undefined FDP and relative ranges remain null; no zero imputation. Pooled
FDP is descriptive because repeated files/seeds are correlated. Seed roles refer
to new random draws on previously inspected files, not held-out biological samples.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import statistics

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "docs/stability-study.md"


def sha(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def read_tsv(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def median(values):
    return statistics.median(values) if values else None


def load_study(directory):
    directory = directory.resolve()
    if (directory / "incomplete.json").exists():
        raise ValueError(f"incomplete study: {directory}")
    study = json.loads((directory / "study.json").read_text())
    provenance = json.loads((directory / "summary_provenance.json").read_text())
    paths = {str(directory / "study.json"): provenance["study_sha256"]}
    paths.update({str(directory / name): digest for name, digest in provenance["summaries"].items()})
    paths.update(provenance["outputs"])
    if any(sha(path) != digest for path, digest in paths.items()):
        raise ValueError(f"study/summary/output hash mismatch: {directory}")
    rows = read_tsv(directory / "runs.tsv")
    expected = {(j["dataset"], j["variant"], j["seed"]) for j in study["searches"]}
    for row in rows:
        for field in ("seed", "rows", "depth_1pct", "depth_5pct", "internal_q1_targets", "internal_q1_human", "internal_q1_yeast"):
            row[field] = int(row[field])
    keys = [(r["dataset"], r["variant"], r["seed"]) for r in rows]
    if len(keys) != len(set(keys)) or set(keys) != expected:
        raise ValueError(f"missing, duplicated or unexpected runs: {directory}")
    return study, rows, paths


def aggregates(rows, scale):
    summaries, comparisons = [], []
    for subset in ("all", "development", "held_out"):
        selected = [r for r in rows if subset == "all" or r["seed_role"] == subset]
        if not selected:
            continue
        ranges = {}
        for variant in ("baseline", "candidate"):
            rs = [r for r in selected if r["variant"] == variant]
            by_file = {ds: [r for r in rs if r["dataset"] == ds] for ds in sorted({r["dataset"] for r in rs})}
            defined, excluded = {}, {}
            for ds, data in by_file.items():
                depths = [r["depth_1pct"] for r in data]
                med = median(depths)
                if len(depths) < 2 or not med:
                    excluded[ds] = "fewer than two seeds" if len(depths) < 2 else "zero median"
                else:
                    defined[ds] = (max(depths) - min(depths)) / med
            ranges[variant] = defined, excluded
            h, y = sum(r["internal_q1_human"] for r in rs), sum(r["internal_q1_yeast"] for r in rs)
            summaries.append(dict(subset=subset, variant=variant, searches=len(rs), files=len(by_file),
                empty_output_count=sum(r["rows"] == 0 for r in rs),
                zero_depth_1pct_count=sum(r["depth_1pct"] == 0 for r in rs),
                zero_internal_q1_count=sum(r["internal_q1_targets"] == 0 for r in rs),
                output_rows=sum(r["rows"] for r in rs), internal_q1_targets=sum(r["internal_q1_targets"] for r in rs),
                internal_q1_human=h, internal_q1_yeast=y, pooled_internal_q1_fdp=y * scale / h if h else None,
                median_file_depth_1pct=median([median([r["depth_1pct"] for r in data]) for data in by_file.values()]),
                median_file_depth_5pct=median([median([r["depth_5pct"] for r in data]) for data in by_file.values()]),
                median_relative_1pct_range=median(list(defined.values())), defined_range_files=len(defined), excluded_range_files=excluded))
        common = sorted(set(ranges["baseline"][0]) & set(ranges["candidate"][0]))
        comparisons.append(dict(subset=subset, comparable_files=common,
            baseline_median_relative_range=median([ranges["baseline"][0][ds] for ds in common]),
            candidate_median_relative_range=median([ranges["candidate"][0][ds] for ds in common]),
            excluded_files={v: ranges[v][1] for v in ranges}))
    return summaries, comparisons


def file_medians(rows):
    output = []
    for ds in sorted({r["dataset"] for r in rows}):
        by_variant = {v: [r for r in rows if r["dataset"] == ds and r["variant"] == v] for v in ("baseline", "candidate")}
        if {r["seed"] for r in by_variant["baseline"]} != {r["seed"] for r in by_variant["candidate"]}:
            raise ValueError(f"unmatched seeds for {ds}")
        for metric in ("depth_1pct", "depth_5pct"):
            a, b = (median([r[metric] for r in by_variant[v]]) for v in ("baseline", "candidate"))
            output.append(dict(dataset=ds, metric=metric, n_seeds=len(by_variant["baseline"]),
                               baseline_median=a, candidate_median=b, delta=b - a,
                               relative_delta=(b - a) / a if a else None,
                               loses_more_than_5pct=b < 0.95 * a))
    return output


def paired_outputs(rows):
    output = []
    non_group = lambda r: tuple(sorted((k, v) for k, v in r.items() if k not in {"n_group_members", "group_members"}))
    for ds, seed in sorted({(r["dataset"], r["seed"]) for r in rows}):
        a, b = (next(r for r in rows if (r["dataset"], r["seed"], r["variant"]) == (ds, seed, v)) for v in ("baseline", "candidate"))
        ar, br = read_tsv(a["path"]), read_tsv(b["path"])
        scores_equal = sorted(map(non_group, ar)) == sorted(map(non_group, br))
        full_rows_equal = sorted(tuple(sorted(r.items())) for r in ar) == sorted(tuple(sorted(r.items())) for r in br)
        losses = {metric: dict(baseline=a[metric], candidate=b[metric], delta=b[metric] - a[metric],
                              relative_delta=(b[metric] - a[metric]) / a[metric] if a[metric] else None,
                              loss_above_10pct=b[metric] < 0.9 * a[metric]) for metric in ("depth_1pct", "depth_5pct")}
        output.append(dict(dataset=ds, seed=seed, seed_role=a["seed_role"],
                           baseline_internal_q1_targets=a["internal_q1_targets"],
                           candidate_internal_q1_targets=b["internal_q1_targets"],
                           newly_empty_reported_1pct=a["internal_q1_targets"] > 0 and b["internal_q1_targets"] == 0,
                           score_rows_equal=scores_equal,
                           group_only_changes=scores_equal and not full_rows_equal,
                           byte_equal=sha(a["path"]) == sha(b["path"]), depth_changes=losses))
    return output


def acceptance(panels, lowinput_seed_count=5, protect_discoveries=False):
    checks = []
    def add(name, passed, detail):
        checks.append(dict(criterion=name, passed=passed, detail=detail))
    for name in ("lowinput", "primary"):
        if name not in panels:
            add(f"{name}_panel_complete", None, "Required panel absent")
            continue
        panel = panels[name]
        if protect_discoveries:
            collapses = [dict(dataset=p["dataset"], seed=p["seed"], baseline=p["baseline_internal_q1_targets"],
                              candidate=p["candidate_internal_q1_targets"])
                         for p in panel["paired_outputs"] if p["newly_empty_reported_1pct"]]
            add(f"{name}_no_newly_empty_reported_1pct", not collapses, collapses)
        expected_seeds = list(range(lowinput_seed_count if name == "lowinput" else 3))
        design = panel["seeds_by_dataset"]
        add(f"{name}_panel_complete", len(design) == 6 and all(seeds == expected_seeds for seeds in design.values()),
            dict(expected_files=6, expected_seeds=expected_seeds, observed=design))
        whole = {r["variant"]: r for r in panel["aggregates"] if r["subset"] == "all"}
        a, b = whole["baseline"], whole["candidate"]
        fa, fb = a["pooled_internal_q1_fdp"], b["pooled_internal_q1_fdp"]
        add(f"{name}_pooled_fdp", fb <= fa + 0.005 if fa is not None and fb is not None else None,
            dict(baseline=fa, candidate=fb, maximum_increase=0.005))
        if name == "lowinput":
            az, bz = a["zero_depth_1pct_count"], b["zero_depth_1pct_count"]
            add("lowinput_zero_depth_halved", bz * 2 <= az if az else None, dict(baseline=az, candidate=bz))
            relative = next(r for r in panel["paired_relative_ranges"] if r["subset"] == "all")
            x, y = relative["baseline_median_relative_range"], relative["candidate_median_relative_range"]
            add("lowinput_relative_range_lower", y < x if x is not None and y is not None else None, relative)
            for metric in ("median_file_depth_1pct", "median_file_depth_5pct"):
                add(f"lowinput_{metric}_nondecreasing", b[metric] >= a[metric], dict(baseline=a[metric], candidate=b[metric]))
        else:
            losses = [r for r in panel["file_medians"] if r["loses_more_than_5pct"]]
            add("primary_file_medians_within_5pct", not losses, losses)
    return dict(promote=bool(checks) and all(r["passed"] is True for r in checks), checks=checks,
                indeterminate=any(r["passed"] is None for r in checks),
                depth_interpretation="Median-depth guards are evaluated at both entrapment 1% and 5%; relative range and collapse gates use 1%.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--studies", nargs="+", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True, help="output filename prefix, without .json/.md")
    parser.add_argument("--protocol", type=Path, default=PROTOCOL, help="frozen acceptance-plan snapshot")
    parser.add_argument("--lowinput-seed-count", type=int, default=5, help="required contiguous low-input seeds beginning at 0; use 7 for the prespecified v2 panel")
    parser.add_argument("--protect-discoveries", action="store_true", help="v2: reject newly empty reported-1%% lists when the matched baseline has targets")
    args = parser.parse_args()
    if args.lowinput_seed_count < 1:
        parser.error("--lowinput-seed-count must be positive")
    protocol = args.protocol.expanduser().resolve()
    panels, hashes, inputs = {}, {str(protocol): sha(protocol), str(Path(__file__).resolve()): sha(__file__)}, {}
    for directory in args.studies:
        study, rows, verified = load_study(directory)
        hashes.update(verified)
        if any(path in inputs and inputs[path] != digest for path, digest in study["fingerprints"].items()):
            raise ValueError("studies disagree on the hash of a shared input")
        inputs.update(study["fingerprints"])
        for name in sorted({r["panel"] for r in rows}):
            if name in panels:
                raise ValueError(f"duplicate panel {name}; supply one complete study per panel")
            selected = [r for r in rows if r["panel"] == name]
            summary, relative = aggregates(selected, study["entrapment_scale"])
            pairs = paired_outputs(selected)
            seeds_by_dataset = {ds: sorted({r["seed"] for r in selected if r["dataset"] == ds}) for ds in sorted({r["dataset"] for r in selected})}
            panels[name] = dict(study=str(directory.resolve()), seeds_by_dataset=seeds_by_dataset,
                                aggregates=summary, paired_relative_ranges=relative,
                                file_medians=file_medians(selected), paired_outputs=pairs,
                                losses_above_10pct=[p for p in pairs if any(v["loss_above_10pct"] for v in p["depth_changes"].values())])
    report = dict(acceptance=acceptance(panels, args.lowinput_seed_count, args.protect_discoveries),
                  protocol=str(protocol), lowinput_seed_count=args.lowinput_seed_count,
                  protect_discoveries=args.protect_discoveries, panels=panels, verified_file_hashes=hashes,
                  study_input_hashes=inputs, caveats=[
                      "Repeated seeds/files are correlated: pooled FDP and ranges are descriptive, not independent-sample confidence estimates.",
                      "Undefined FDP and zero-median/fewer-than-two-seed ranges remain null; relative-range acceptance uses only the same defined files in both variants.",
                      "Held-out seeds are new random draws on inspected files, not new biological samples; all individual losses remain reported."])
    lines = ["# UNO stability study", "", f"Acceptance rule: {'PASS' if report['acceptance']['promote'] else 'NOT PASSED'}.", "",
             f"Protocol: `{protocol}`. Required low-input seeds: 0–{args.lowinput_seed_count - 1}. Protect reported discoveries: {args.protect_discoveries}.", ""]
    for name, panel in panels.items():
        lines += [f"## {name}", "", "| Subset | Variant | Runs | Empty outputs | Zero 1% depth | Zero internal 1% | Pooled FDP |", "|---|---|---:|---:|---:|---:|---:|"]
        for row in panel["aggregates"]:
            fdp = row["pooled_internal_q1_fdp"]
            lines.append(f"| {row['subset']} | {row['variant']} | {row['searches']} | {row['empty_output_count']} | {row['zero_depth_1pct_count']} | {row['zero_internal_q1_count']} | {100 * fdp:.3f}% |" if fdp is not None else f"| {row['subset']} | {row['variant']} | {row['searches']} | {row['empty_output_count']} | {row['zero_depth_1pct_count']} | {row['zero_internal_q1_count']} | NA |")
        pairs = panel["paired_outputs"]
        lines += ["", f"Score rows equal in {sum(p['score_rows_equal'] for p in pairs)}/{len(pairs)} pairs; group-only changes in {sum(p['group_only_changes'] for p in pairs)}. Individual pairs with >10% depth loss: {len(panel['losses_above_10pct'])}.", ""]
    lines += ["## Acceptance checks", ""]
    for check in report["acceptance"]["checks"]:
        status = "PASS" if check["passed"] is True else "FAIL" if check["passed"] is False else "UNDEFINED"
        lines += [f"- {status}: {check['criterion']} — `{json.dumps(check['detail'], sort_keys=True)}`"]
    lines += ["", *report["caveats"], ""]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.with_suffix(".json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    args.output.with_suffix(".md").write_text("\n".join(lines))
    print(f"Wrote {args.output.with_suffix('.json')} and {args.output.with_suffix('.md')}")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError) as error:
        raise SystemExit(str(error))
