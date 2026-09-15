#!/usr/bin/env python3
"""Post hoc, matrix-only seed ensemble diagnostic; no new searches or paper edits.

Primary aggregation: minimum saved internal q per accession, for both routes.
Sensitivity: maximum saved binomial score, for both routes. No score averaging,
seed selection, or cross-arm rescue. Reevaluate entrapment after aggregation.
"""
import csv
import hashlib
import itertools
import json
from collections import defaultdict
from decimal import Decimal
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parents[1]
STUDY = ROOT / "benchmark/results/matrix_20260913_koth050"
OUT = ROOT / "benchmark/results/matrix_seed_ensemble_20260913"
SEEDS = (0, 1, 2)
SPECIES = ("HUMAN", "ECOLI", "YEAST")
THRESHOLDS = (1, 3, 5)


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def kind(acc):
    matches = [sp for sp in SPECIES if "_" + sp in acc]
    return matches[0] if len(matches) == 1 else "AMBIGUOUS"


def merge(tables):
    result = {}
    for table in tables:
        for acc, value in table.items():
            result[acc] = max(result.get(acc, Decimal("-Infinity")), value)
    return result


def evaluate(scores, counts):
    blocks = defaultdict(list)
    for acc, score in scores.items():
        blocks[score].append(acc)
    present = yeast = 0
    best = {t: (0, 0, None) for t in THRESHOLDS}
    for score in sorted(blocks, reverse=True):
        present += sum(kind(a) != "YEAST" for a in blocks[score])
        yeast += sum(kind(a) == "YEAST" for a in blocks[score])
        for threshold in THRESHOLDS:
            if (present > best[threshold][0]
                    and yeast * (counts["HUMAN"] + counts["ECOLI"]) * 100
                    <= threshold * present * counts["YEAST"]):
                best[threshold] = present, yeast, score
    result = {}
    for threshold, (n, y, score) in best.items():
        accepted = sorted(a for a, s in scores.items()
                          if score is not None and s >= score)
        assert n == sum(kind(a) != "YEAST" for a in accepted)
        assert y == sum(kind(a) == "YEAST" for a in accepted)
        result[str(threshold / 100)] = dict(
            present=n, yeast=y, rank_threshold=str(score) if score is not None else None,
            empirical_fdp=y * (counts["HUMAN"] + counts["ECOLI"]) / counts["YEAST"] / n if n else None,
            accessions=[a for a in accepted if kind(a) != "YEAST"],
            yeast_accessions=[a for a in accepted if kind(a) == "YEAST"])
    return result


def main():
    report = json.loads((STUDY / "report.json").read_text())
    hashes = {str(STUDY / "report.json"): sha(STUDY / "report.json")}
    cache = {}

    def load(path):
        if path in cache:
            return cache[path]
        receipt_path = path.parent / "job.json"
        receipt = json.loads(receipt_path.read_text())
        assert receipt["returncode"] == 0
        hashes[str(path)] = sha(path)
        hashes[str(receipt_path)] = sha(receipt_path)
        assert hashes[str(path)] == receipt["outputs"][str(path)]
        assert hashes[str(path)] == report["hashes"][str(path)]
        tables = {"min_q": {}, "max_score": {}}
        for row in csv.DictReader(path.open(), delimiter="\t"):
            if row["decoy"] != "0":
                continue
            members = row["group_members"].split(";")
            classes = {kind(a) for a in members}
            if len(classes) != 1 or "AMBIGUOUS" in classes:
                continue
            for method, value in (("min_q", -Decimal(row["q_value"])),
                                  ("max_score", Decimal(row["binomial_score"]))):
                assert value.is_finite()
                for acc in members:
                    tables[method][acc] = max(tables[method].get(acc, Decimal("-Infinity")), value)
        cache[path] = tables
        return tables

    results = {}
    for dataset, dd in report["datasets"].items():
        counts = {sp: 0 for sp in SPECIES}
        fasta = STUDY / dataset / "human_ecoli_yeast.fasta"
        hashes[str(fasta)] = sha(fasta)
        for line in fasta.open():
            if line.startswith(">"):
                sp = kind(line.split()[0])
                assert sp in counts
                counts[sp] += 1
        results[dataset] = {}
        for cohort, cd in dd["cohorts"].items():
            result = dict(n_files=len(cd["runs"]), fasta_counts=counts, arms={})
            results[dataset][cohort] = result
            for arm in ("union", "joint", "joint_cell"):
                paths = {seed: [STUDY / dataset / "searches" / name / f"seed{seed}" / "proteins.tsv"
                                for name in (["single_" + run for run in cd["runs"]] if arm == "union"
                                             else [arm + "_" + cohort])] for seed in SEEDS}
                result["arms"][arm] = {}
                for method in ("min_q", "max_score"):
                    by_seed = {seed: merge([load(p)[method] for p in paths[seed]]) for seed in SEEDS}
                    subsets = {}
                    for size in (1, 2, 3):
                        for seeds in itertools.combinations(SEEDS, size):
                            combined = merge([by_seed[s] for s in seeds])
                            evaluated = evaluate(combined, counts)
                            if size == 1 and method == ("min_q" if arm == "union" else "max_score"):
                                for threshold, point in evaluated.items():
                                    old = cd["seeds"][str(seeds[0])][arm][threshold]
                                    assert (point["present"], point["yeast"], point["accessions"]) == (old["present"], old["yeast"], old["accessions"])
                            internal = [a for a, q in merge([load(p)["min_q"] for s in seeds for p in paths[s]]).items()
                                        if q >= Decimal("-0.01")]
                            n = sum(kind(a) != "YEAST" for a in internal)
                            y = len(internal) - n
                            subsets[",".join(map(str, seeds))] = dict(seeds=seeds,
                                search_invocations=sum(len(paths[s]) for s in seeds),
                                thresholds=evaluated, internal_1pct=dict(present=n, yeast=y,
                                    empirical_fdp=y * (counts["HUMAN"] + counts["ECOLI"]) / counts["YEAST"] / n if n else None))
                    sizes = {str(size): {t: dict(median=median(vals), min=min(vals), max=max(vals))
                             for t in ("0.01", "0.03", "0.05")
                             for vals in [[x["thresholds"][t]["present"] for x in subsets.values() if len(x["seeds"]) == size]]}
                             for size in (1, 2, 3)}
                    singles = [set(subsets[str(s)]["thresholds"]["0.01"]["accessions"]) for s in SEEDS]
                    ensemble = set(subsets["0,1,2"]["thresholds"]["0.01"]["accessions"])
                    naïve = set.union(*singles)
                    naive_yeast = set.union(*(set(subsets[str(s)]["thresholds"]["0.01"]["yeast_accessions"]) for s in SEEDS))
                    result["arms"][arm][method] = dict(subsets=subsets, by_size=sizes,
                        overlap=dict(independently_thresholded_seed_union=len(naïve),
                            independently_thresholded_seed_union_yeast=len(naive_yeast),
                            independently_thresholded_seed_union_fdp=len(naive_yeast) * (counts["HUMAN"] + counts["ECOLI"]) / counts["YEAST"] / len(naïve) if naïve else None,
                            shared_all_seeds=len(set.intersection(*singles)),
                            ensemble_retained_from_seed_union=len(ensemble & naïve),
                            ensemble_lost_from_seed_union=len(naïve - ensemble),
                            ensemble_added_outside_seed_union=len(ensemble - naïve)))
            print(dataset, cohort, {a: {m: [v["by_size"][str(s)]["0.01"]["median"] for s in (1,2,3)]
                                        for m,v in av.items()} for a,av in result["arms"].items()}, flush=True)
    assert len(cache) == 144
    OUT.mkdir(parents=True, exist_ok=True)
    output = dict(status="complete", exploratory=True, new_searches=0,
        method="All nonempty subsets of saved seeds 0,1,2; minimum internal q (primary) and maximum binomial score (sensitivity), applied symmetrically to both routes. Whole score ties and exact integer entrapment FDP checks at 1%, 3%, 5%.",
        limitations=["Post hoc evaluation; entrapment-selected cutoffs are not deployable cohort q-values.",
                     "Matrix ensembles contain matrix search outputs only; no separate-file results enter them.",
                     "Seeds change RT calibration, CV folds and classifier randomness together; effects are not isolated.",
                     "Search counts differ: k matrix searches versus k times N separate-file searches; no matched compute trial.",
                     "Only three seeds; no confidence intervals or significance claims. Main paper remains unchanged."],
        validation="All 144 protein tables match successful original receipts and audited report hashes; 72 baseline curves/accession sets reproduced at each of three thresholds.",
        source_sha256=sha(Path(__file__)), input_hashes=hashes, datasets=results)
    (OUT / "report.json").write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
