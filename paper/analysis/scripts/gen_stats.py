#!/usr/bin/env python3
"""Derive every analysis-owned number stated in the Technical Note prose."""
from __future__ import annotations

import csv
import statistics
from pathlib import Path

from _stats import Stats
import seed_data
import koth_matrix_stats

HERE = Path(__file__).resolve().parent
PAPER = HERE.parent.parent
BENCH = PAPER.parent / "benchmark"
SUMMARY = BENCH / "results" / "summary.tsv"
ENTRAPMENT = PAPER / "figdata" / "entrapment_curves.tsv"
OVERLAP = PAPER / "figdata" / "overlap_analysis.tsv"
RAW_MODE = BENCH / "raw_mode" / "summary.tsv"
SPREAD = BENCH / "results" / "spread.tsv"
EXTENSION = BENCH / "results/variability_extension/summary.tsv"
ASTRAL_SUMMARY = Path.home() / "Data/PXD066701_Astral/results/input_series_summary.tsv"
ENTRAP_FASTA = BENCH / "data" / "entrap_human_yeast.fasta"
RERUN = BENCH / "results/publication_rerun/cases.tsv"

# The generating Astral run is external to the repository. These eight values
# remain explicit hand entries with provenance notes in stats.json (plus the
# external-RT-calibration arm added 2026-09-03); every other
# hand-entered value belonged to a manuscript branch removed in the Technical
# Note reduction and is intentionally retired on the next write.
HAND_IDS = {
    "astral.files",
    "astral.levels",
    "astral.memory_max",
    "astral.memory_min",
    "astral.ms1_median_1pct_100ng",
    "astral.ms1_median_1pct_10ng",
    "astral.ms1_median_1pct_1ng",
    "astral.ms1_median_1pct_200pg",
    "astral.ms1_median_5pct_100ng",
    "astral.ms1_median_5pct_10ng",
    "astral.ms1_median_5pct_1ng",
    "astral.ms1_median_5pct_200pg",
    "astral.ms1_stock_rt_sigma_max",
    "astral.ms1_stock_rt_sigma_min",
    "astral.ms1pl_median_1pct_100ng",
    "astral.ms1pl_median_1pct_10ng",
    "astral.ms1pl_median_1pct_1ng",
    "astral.ms1pl_median_1pct_200pg",
    "astral.ms1pl_median_5pct_100ng",
    "astral.ms1pl_median_5pct_10ng",
    "astral.ms1pl_median_5pct_1ng",
    "astral.ms1pl_median_5pct_200pg",
    "astral.ms1pl_rt_sigma_max",
    "astral.ms1pl_rt_sigma_min",
    "astral.pl_list_max",
    "astral.pl_list_min",
    "astral.speed_max",
    "astral.speed_min",
    "astral.uno_median_1pct_100ng",
    "astral.uno_median_1pct_10ng",
    "astral.uno_median_1pct_1ng",
    "astral.uno_median_1pct_200pg",
    "astral.uno_median_5pct_100ng",
    "astral.uno_median_5pct_10ng",
    "astral.uno_median_5pct_1ng",
    "astral.uno_median_5pct_200pg",
}


def table_rows(path: Path) -> list[dict[str, str]]:
    lines = [line for line in path.read_text().splitlines()
             if line and not line.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def fasta_count(path: Path, suffix: str) -> int:
    return sum(line.startswith(">") and suffix in line
               for line in path.read_text().splitlines())


def main() -> int:
    rows = table_rows(SUMMARY)
    spread_rows = seed_data.rows()
    datasets = list(dict.fromkeys(row["dataset"] for row in rows))

    uno_wall = [float(row["wall_s"]) for row in rows if row["tool"] == "uno"]
    ms1_wall = [float(row["wall_s"]) for row in rows if row["tool"] == "ms1"]
    factors: list[float] = []
    margins: list[float] = []
    spreads_abs: list[int] = []
    spreads_pct: list[float] = []
    uno_deeper_5pct = 0
    uno_tied_5pct = 0
    uno_beats_best_1pct = 0

    for dataset in datasets:
        uno = next(row for row in rows
                   if row["dataset"] == dataset and row["tool"] == "uno")
        ms1 = [row for row in rows
               if row["dataset"] == dataset and row["tool"] == "ms1"]
        ms1_wall_dataset = [float(row["wall_s"]) for row in ms1]
        draws = [r for r in spread_rows if r["dataset"] == dataset]
        ms1_1pct = [int(r["depth_1pct"]) for r in draws if r["engine"] == "ms1searchpy"]
        ms1_5pct = [int(r["depth_5pct"]) for r in draws if r["engine"] == "ms1searchpy"]
        uno_1pct = statistics.median(int(r["depth_1pct"]) for r in draws if r["engine"] == "uno")
        uno_5pct = statistics.median(int(r["depth_5pct"]) for r in draws if r["engine"] == "uno")
        median_1pct = statistics.median(ms1_1pct)
        median_5pct = statistics.median(ms1_5pct)

        factors.append(statistics.median(ms1_wall_dataset) / float(uno["wall_s"]))
        margins.append(100 * (uno_1pct - median_1pct) / median_1pct)
        spread = max(ms1_1pct) - min(ms1_1pct)
        spreads_abs.append(spread)
        spreads_pct.append(100 * spread / median_1pct)
        uno_deeper_5pct += uno_5pct > median_5pct
        uno_tied_5pct += uno_5pct == median_5pct
        uno_beats_best_1pct += uno_1pct > max(ms1_1pct)

    st = Stats()
    st.add("bench.dataset_count", len(datasets), fmt="d",
           desc="Primary benchmark datasets", between=(1, 20))
    st.add("bench.ms1_run_count", sum(row["tool"] == "ms1" for row in rows),
           fmt="d", desc="Total primary ms1searchpy benchmark searches",
           sign="+", between=(1, 100))
    st.add("bench.uno_wall_min", min(uno_wall), fmt=".0f", unit="s",
           desc="Minimum uno wall time", sign="+", between=(1, 600))
    st.add("bench.uno_wall_max", max(uno_wall), fmt=".0f", unit="s",
           desc="Maximum uno wall time", sign="+", between=(1, 600))
    st.add("bench.ms1_wall_min", min(ms1_wall), fmt=".0f", unit="s",
           desc="Minimum ms1searchpy wall time", sign="+", between=(1, 3600))
    st.add("bench.ms1_wall_max", max(ms1_wall), fmt=".0f", unit="s",
           desc="Maximum ms1searchpy wall time", sign="+", between=(1, 3600))
    st.add("bench.speed_factor_min", min(factors), fmt=".1f", unit="x",
           desc="Minimum ms1searchpy-median/uno wall-time factor",
           sign="+", between=(1, 50))
    st.add("bench.speed_factor_max", max(factors), fmt=".1f", unit="x",
           desc="Maximum ms1searchpy-median/uno wall-time factor",
           sign="+", between=(1, 50))
    st.add("bench.median_margin_1pct_min", min(margins), fmt=".1f", unit="%",
           desc="Minimum uno versus ms1searchpy median depth margin at 1% entrapment FDR",
           between=(-100, 200))
    st.add("bench.median_margin_1pct_max", max(margins), fmt=".1f", unit="%",
           desc="Maximum uno versus ms1searchpy median depth margin at 1% entrapment FDR",
           between=(-100, 200))
    st.add("bench.uno_median_wins_5pct", uno_deeper_5pct, fmt="d",
           desc="Datasets where uno is deeper at 5% entrapment FDR",
           sign="+", between=(0, len(datasets)))
    st.add("bench.uno_median_ties_5pct", uno_tied_5pct, fmt="d",
           desc="Datasets tied at 5% entrapment FDR",
           between=(0, len(datasets)))
    st.add("bench.uno_median_beats_best_1pct", uno_beats_best_1pct, fmt="d",
           desc="Datasets where uno beats every ms1searchpy replicate at 1%",
           sign="+", between=(0, len(datasets)))
    st.add("bench.ms1_spread_1pct_abs_min", min(spreads_abs), fmt="d",
           desc="Minimum ms1searchpy replicate range at 1%",
           sign="+", between=(1, 5000))
    st.add("bench.ms1_spread_1pct_abs_max", max(spreads_abs), fmt="d",
           desc="Maximum ms1searchpy replicate range at 1%",
           sign="+", between=(1, 5000))
    st.add("bench.ms1_spread_1pct_pct_min", min(spreads_pct), fmt=".1f", unit="%",
           desc="Minimum replicate range as percent of median",
           sign="+", between=(0, 100))
    st.add("bench.ms1_spread_1pct_pct_max", max(spreads_pct), fmt=".1f", unit="%",
           desc="Maximum replicate range as percent of median",
           sign="+", between=(0, 100))

    # -- randomization spread, measured the same way for both engines ---------
    # uno's byte-identical output came from hardcoded constants, not from an
    # algorithm with less randomness in it, so a locked uno against an unlocked
    # ms1searchpy compared the locking rather than the engines. `uno --seed`
    # perturbs the decoy shuffle, RANSAC sampling, fold assignment and booster
    # subsampling together; these values summarize uno's depth over independent
    # seeds beside ms1searchpy's over independent reruns.
    seed_ranges_abs: list[int] = []
    seed_ranges_pct: list[float] = []
    seed_percentiles: list[float] = []
    seed_counts: list[int] = []
    uno_median_wins = 0
    matched_datasets = 0
    for dataset in dict.fromkeys(row["dataset"] for row in spread_rows):
        uno_draws = sorted(int(row["depth_1pct"]) for row in spread_rows
                           if row["dataset"] == dataset and row["engine"] == "uno")
        ms1_draws = sorted(int(row["depth_1pct"]) for row in spread_rows
                           if row["dataset"] == dataset and row["engine"] == "ms1searchpy")
        if len(uno_draws) < 3:
            continue
        seed_counts.append(len(uno_draws))
        span = uno_draws[-1] - uno_draws[0]
        median = statistics.median(uno_draws)
        seed_ranges_abs.append(span)
        seed_ranges_pct.append(100 * span / median)
        published = next((int(row["depth_1pct"]) for row in spread_rows
                          if row["dataset"] == dataset
                          and row["engine"] == "uno" and row["draw"] == "0"), None)
        if published is not None:
            below = sum(value < published for value in uno_draws)
            seed_percentiles.append(100 * below / (len(uno_draws) - 1))
        if ms1_draws:
            matched_datasets += 1
            uno_median_wins += median > statistics.median(ms1_draws)

    ms1_counts = [sum(r["dataset"] == d and r["engine"] == "ms1searchpy" for r in spread_rows)
                  for d in datasets]
    st.add("bench.ms1_draw_count", sum(ms1_counts), fmt="d", desc="Primary comparator draws", between=(1, 1000))
    st.add("bench.ms1_draw_min", min(ms1_counts), fmt="d", desc="Minimum comparator draws per primary file", between=(3, 100))
    st.add("bench.ms1_draw_max", max(ms1_counts), fmt="d", desc="Maximum comparator draws per primary file", between=(3, 100))
    st.add("bench.uno_draw_count", sum(seed_counts), fmt="d", desc="Primary uno seed draws", between=(1, 1000))
    st.add("bench.uno_seed_count", min(seed_counts), fmt="d",
           desc="Independent uno seeds per dataset in the spread comparison",
           sign="+", between=(3, 100))
    st.add("bench.uno_seed_spread_1pct_abs_min", min(seed_ranges_abs), fmt="d",
           desc="Minimum uno seed-to-seed range at 1% entrapment FDR",
           sign="+", between=(1, 5000))
    st.add("bench.uno_seed_spread_1pct_abs_max", max(seed_ranges_abs), fmt="d",
           desc="Maximum uno seed-to-seed range at 1% entrapment FDR",
           sign="+", between=(1, 5000))
    st.add("bench.uno_seed_spread_1pct_pct_min", min(seed_ranges_pct), fmt=".1f",
           unit="%", desc="Minimum uno seed range as percent of its median",
           sign="+", between=(0, 100))
    st.add("bench.uno_seed_spread_1pct_pct_max", max(seed_ranges_pct), fmt=".1f",
           unit="%", desc="Maximum uno seed range as percent of its median",
           sign="+", between=(0, 100))
    st.add("bench.uno_seed0_percentile_max", max(seed_percentiles), fmt=".0f",
           desc="Highest percentile of a published seed-0 run within its own seed draws",
           sign="+", between=(0, 100))
    st.add("bench.uno_seed0_high_datasets",
           sum(value >= 90 for value in seed_percentiles), fmt="d",
           desc="Datasets where the published uno run is at or above the 90th "
                "percentile of its own seed draws",
           between=(0, 20))
    st.add("bench.uno_median_wins", uno_median_wins, fmt="d",
           desc="Datasets where uno's seed median exceeds ms1searchpy's rerun median",
           between=(0, 20))
    st.add("bench.spread_datasets", matched_datasets, fmt="d",
           desc="Datasets with both uno seed draws and ms1searchpy reruns",
           sign="+", between=(1, 20))

    entrapment = table_rows(ENTRAPMENT)
    reported_1pct = [row for row in entrapment
                     if row["cutoff_type"] == "reported_q"
                     and float(row["cutoff"]) == 0.01]
    uno_1pct = [row for row in reported_1pct if row["engine"] == "uno"]
    ms1_1pct = [row for row in reported_1pct if row["engine"] == "ms1searchpy"]
    uno_fdp = [100 * float(row["entrapment_fdp"]) for row in uno_1pct]

    st.add("cal.uno_fdp_1pct_min", min(uno_fdp), fmt=".2f", unit="%",
           desc="Minimum uno empirical FDP at reported 1%", between=(0, 100))
    st.add("cal.uno_fdp_1pct_max", max(uno_fdp), fmt=".2f", unit="%",
           desc="Maximum uno empirical FDP at reported 1%", between=(0, 100))
    st.add("cal.uno_fdp_1pct_median", statistics.median(uno_fdp), fmt=".2f", unit="%",
           desc="Median uno empirical FDP at reported 1%", between=(0, 100))
    st.add("cal.uno_conservative_runs_1pct", sum(value < 1 for value in uno_fdp), fmt="d",
           desc="Searches where uno empirical FDP is below 1%",
           between=(0, len(uno_1pct)))
    st.add("cal.ms1_optimistic_runs_1pct", sum(
        float(row["entrapment_fdp"]) > 0.01 for row in ms1_1pct), fmt="d",
        desc="ms1searchpy searches above 1% empirical FDP",
        sign="+", between=(0, len(ms1_1pct)))

    human_db = fasta_count(ENTRAP_FASTA, "_HUMAN")
    yeast_db = fasta_count(ENTRAP_FASTA, "_YEAST")
    st.add("bench.human_db", human_db, fmt=",d",
           desc="Human target entries in benchmark FASTA",
           sign="+", between=(1, 100000))
    st.add("bench.yeast_db", yeast_db, fmt=",d",
           desc="Yeast target entries in benchmark FASTA",
           sign="+", between=(1, 100000))
    st.add("cal.db_scale", human_db / yeast_db, fmt=".2f",
           desc="Human/yeast database-size correction factor",
           sign="+", between=(1, 10))

    overlap = table_rows(OVERLAP)

    def frac(analysis: str) -> list[float]:
        return [100 * float(r["fraction"]) for r in overlap if r["analysis"] == analysis]

    for key, analysis, desc in (
        ("jaccard_cross", "jaccard_uno_vs_ms1",
         "shared/union of uno and one ms1searchpy 1% list"),
        ("jaccard_within_uno", "jaccard_within_uno",
         "shared/union of two uno seed 1% lists"),
        ("jaccard_within", "jaccard_within_ms1",
         "shared/union of two ms1searchpy replicate 1% lists"),
        ("core_recovery", "uno_recovers_ms1_core",
         "proteins accepted by all available ms1searchpy replicates that uno also accepted"),
        ("uno_confirmed", "uno_confirmed_by_any_ms1",
         "uno 1% proteins accepted by at least one ms1searchpy replicate"),
        ("uno_novel", "uno_in_no_ms1",
         "uno 1% proteins accepted by no ms1searchpy replicate"),
        ("ms1_confirmed", "ms1_confirmed_by_other_ms1",
         "ms1searchpy 1% proteins accepted by at least one other replicate"),
    ):
        vals = frac(analysis)
        st.add(f"overlap.{key}_min", min(vals), fmt=".0f", unit="%",
               desc=f"Minimum {desc}", sign="+", between=(0, 100))
        st.add(f"overlap.{key}_max", max(vals), fmt=".0f", unit="%",
               desc=f"Maximum {desc}", sign="+", between=(0, 100))

    # Built-in feature detection (Figure 6). koth partitions MS1 detection by
    # FAIMS compensation voltage, so fixed-CV, multi-CV, and non-FAIMS files are
    # members of the same comparison population.
    extension = table_rows(EXTENSION)
    st.add("bench.extension_seed_count", len({r["seed"] for r in extension}), fmt="d",
           desc="Seeds per file and arm in the Astral/raw extension", between=(3, 100))
    # Acquisition files are equally weighted: median across seeds/reruns per
    # file, then median of the three file medians per input level.
    astral_rows = [r for r in table_rows(ASTRAL_SUMMARY) if r["engine"] in ("ms1searchpy", "ms1searchpy-pl")]
    for r in extension:
        if r["astral_stem"] and r["arm"] == "biosaur2":
            astral_rows.append(dict(r, dataset=r["astral_stem"], engine="uno",
                                   level=r["dataset"].split("_")[1]))
    generated_astral = set()
    astral_level_medians = {}
    for engine, tag in (("uno", "uno"), ("ms1searchpy", "ms1"), ("ms1searchpy-pl", "ms1pl")):
        for level in ("200pg", "1ng", "10ng", "100ng"):
            rr = [r for r in astral_rows if r["engine"] == engine and r["level"] == level]
            for threshold in ("1pct", "5pct"):
                medians = [statistics.median(int(r[f"depth_{threshold}"]) for r in rr if r["dataset"] == ds)
                           for ds in sorted({r["dataset"] for r in rr})]
                assert len(medians) == 3, (engine, level, medians)
                key = f"astral.{tag}_median_{threshold}_{level}"
                st.add(key, statistics.median(medians), fmt=".0f",
                       desc=f"Median of per-file median depths: {engine}, {level}, {threshold}", between=(0, 10000))
                astral_level_medians.setdefault((tag, threshold), []).append(statistics.median(medians))
                generated_astral.add(key)
    for (tag, threshold), medians in astral_level_medians.items():
        for endpoint, value in (("min", min(medians)), ("max", max(medians))):
            st.add(f"astral.{tag}_level_median_{threshold}_{endpoint}", value, fmt=".0f",
                   desc=f"{endpoint} across four input-level medians: {tag} at {threshold}", between=(0, 10000))
    astral_uno = [r for r in astral_rows if r["engine"] == "uno"]
    st.add("astral.uno_seed_zero_depth", sum(int(r["depth_1pct"]) == 0 for r in astral_uno),
           fmt="d", desc="Astral uno searches with zero depth at 1% entrapment FDR", between=(0, len(astral_uno)))
    st.add("astral.uno_seed_draws", len(astral_uno), fmt="d", desc="Astral uno seed draws", between=(1, 1000))
    st.add("astral.uno_seed_zero_files", len({r["dataset"] for r in astral_uno if int(r["depth_1pct"]) == 0}),
           fmt="d", desc="Astral files with at least one zero-depth uno seed", between=(0, 12))
    # Seed sensitivity at the two lowest inputs, as the full within-file range
    # over the matched seeds divided by that file's median depth. Recovery
    # removes empty lists; it does not make one low-input search sufficient.
    low_ranges = []
    for ds in sorted({r["dataset"] for r in astral_uno if r["level"] in ("200pg", "1ng")}):
        v = [int(r["depth_1pct"]) for r in astral_uno if r["dataset"] == ds]
        median = statistics.median(v)
        if median <= 0:
            raise ValueError(f"undefined low-input range on {ds}; revise the prose")
        low_ranges.append(100 * (max(v) - min(v)) / median)
    for tag, value in (("min", min(low_ranges)), ("max", max(low_ranges))):
        st.add(f"astral.low_range_pct_{tag}", value, fmt=".0f", unit="%",
               desc=f"{tag} within-file seed range at 200 pg or 1 ng, as percent of that file's median depth",
               between=(0, 500))
    level_medians = [statistics.median(statistics.median(int(r["depth_1pct"]) for r in astral_uno
                     if r["dataset"] == ds) for ds in sorted({r["dataset"] for r in astral_uno if r["level"] == lv}))
                     for lv in ("200pg", "1ng", "10ng", "100ng")]
    for tag, value in (("min", min(level_medians)), ("max", max(level_medians))):
        st.add(f"astral.uno_seed_level_median_{tag}", value, fmt=".0f",
               desc=f"{tag} Astral input-level median of per-file seed medians", between=(0, 10000))
    raw_rows = table_rows(RAW_MODE)
    raw_datasets = list(dict.fromkeys(row["dataset"] for row in raw_rows))
    low_names = {d for d in raw_datasets if d.startswith(("astral_200pg", "astral_1ng"))}
    rest_names = set(raw_datasets) - low_names
    # Paired-seed ratios are reported per input group. A ratio needs a nonzero
    # biosaur2 depth in its own seed; an empty built-in list is a defined ratio
    # of zero and is kept. The zero-depth searches are also counted below.
    groups = {"rest": (rest_names, []), "low": (low_names, [])}
    agreement = []
    for dataset in raw_datasets:
        bb = {int(r["seed"]): r for r in extension if r["dataset"] == dataset and r["arm"] == "biosaur2"}
        rr = {int(r["seed"]): r for r in extension if r["dataset"] == dataset and r["arm"] == "raw"}
        assert set(bb) == set(rr) == {0, 1, 2}, dataset
        if any(int(bb[k]["depth_1pct"]) == 0 for k in bb):
            raise ValueError(f"undefined depth ratio on {dataset}; revise the prose")
        group = "rest" if dataset in rest_names else "low"
        groups[group][1].append(
            statistics.median(int(rr[k]["depth_1pct"]) / int(bb[k]["depth_1pct"]) for k in bb))
        if dataset in rest_names:
            if any(rr[k]["jaccard_1pct"] == "" for k in bb):
                raise ValueError(f"undefined higher-input list agreement on {dataset}; revise the prose")
            agreement.append(statistics.median(100 * float(rr[k]["jaccard_1pct"]) for k in bb))
    for key, value, desc in (
        ("dataset_count", len(raw_datasets), "Files with both feature sources"),
        ("low_count", len(low_names), "Files at 200 pg or 1 ng"),
        ("rest_count", len(rest_names), "Files at 10 ng or more"),
    ):
        st.add(f"raw.{key}", value, fmt="d", desc=desc, between=(1, 20))
    for group, label in (("rest", ">=10 ng"), ("low", "200 pg or 1 ng")):
        ratios = groups[group][1]
        for tag, value in (("min", min(ratios)), ("max", max(ratios)),
                           ("median", statistics.median(ratios))):
            st.add(f"raw.paired_ratio_{tag}_{group}", value, fmt=".2f",
                   desc=f"{tag} per-file median paired-seed built-in/biosaur2 depth ratio at {label}",
                   between=(0, 10))
    for tag, value in (("min", min(agreement)), ("max", max(agreement))):
        st.add(f"raw.paired_jaccard_{tag}_rest", value, fmt=".0f", unit="%",
               desc=f"{tag} per-file median paired-seed list agreement at >=10 ng", between=(0, 100))
    for arm in ("biosaur2", "raw"):
        low = [r for r in extension if r["dataset"] in low_names and r["arm"] == arm]
        st.add(f"raw.low_zero_{arm}", sum(int(r["depth_1pct"]) == 0 for r in low), fmt="d",
               desc=f"Zero-depth searches at low input, {arm}", between=(0, len(low)))
    st.add("raw.low_draws_per_arm", len(low), fmt="d", desc="Low-input searches per feature source", between=(1, 1000))
    st.add("raw.wall_max", max(float(r["wall_s"]) for r in raw_rows if r["arm"] == "raw"),
           fmt=".0f", unit="s", desc="Maximum original serialized raw-mode wall time", between=(1, 600))

    # Provenance of the main single-run comparisons. These searches were
    # re-executed on one pinned build (benchmark/publication_rerun.py); the
    # case table records, per search, whether the protein output is
    # byte-identical to the superseded build's and what its depth was before.
    rerun = table_rows(RERUN)
    st.add("astral.uno_seed_recovered", sum(
        bool(r["astral_stem"]) and r["arm"] == "biosaur2" and r["had_prev"] == "1"
        and r["prev_depth_1pct"] == "0" and int(r["depth_1pct"]) > 0 for r in rerun),
        fmt="d", desc="Astral biosaur2 seed searches recovered from zero by the pinned rerun",
        between=(0, 1000))
    prev = [(r, r["prev_depth_1pct"], int(r["depth_1pct"])) for r in rerun if r["had_prev"] == "1"]
    st.add("rerun.searches", len(rerun), fmt="d",
           desc="uno searches re-executed on the pinned build", between=(1, 1000))
    st.add("rerun.unchanged", sum(int(r["identical_to_prev"]) for r in rerun), fmt="d",
           desc="Re-run uno outputs byte-identical to the superseded build", between=(0, 1000))
    st.add("rerun.changed", sum(1 for r, _, _ in prev if r["identical_to_prev"] == "0"), fmt="d",
           desc="Re-run uno outputs that differ from the superseded build", between=(0, 1000))
    st.add("rerun.recovered", sum(1 for _, old, new in prev if old == "0" and new > 0), fmt="d",
           desc="Searches that returned nothing on the superseded build and identify proteins on the pinned build",
           between=(0, 1000))
    st.add("rerun.empty", sum(1 for r in rerun if int(r["depth_1pct"]) == 0), fmt="d",
           desc="Re-run uno searches that still return no proteins at 1% entrapment FDR",
           between=(0, 1000))

    koth_matrix_stats.add(st)

    return st.write(
        inputs=[
            "figdata/koth_matrix_study.json",
            "../benchmark/results/summary.tsv",
            "../benchmark/results/spread.tsv",
            "../benchmark/results/variability_extension/summary.tsv",
            str(ASTRAL_SUMMARY),
            "figdata/entrapment_curves.tsv",
            "figdata/overlap_analysis.tsv",
            "../benchmark/raw_mode/summary.tsv",
            "../benchmark/results/publication_rerun/cases.tsv",
            "../benchmark/data/entrap_human_yeast.fasta",
        ],
        keep_hand_ids=HAND_IDS - generated_astral,
    )


if __name__ == "__main__":
    raise SystemExit(main())
