#!/usr/bin/env python3
"""Independent integer-score verification and fixed-rule checks on saved cohorts."""
import csv
import itertools
import json
from collections import defaultdict
from decimal import Decimal
from pathlib import Path

from joint_evidence_audit import BASE, OLD, OUT, sha


def evaluate(files, counts):
    scores, qs = {}, {}
    for run, path in enumerate(files):
        with path.open() as f:
            for row in csv.DictReader(f, delimiter='\t'):
                if int(row['decoy']):
                    continue
                members = row['group_members'].split(';')
                classes = [{s for s in ('HUMAN', 'ECOLI', 'YEAST') if '_' + s in a} for a in members]
                if any(len(c) != 1 for c in classes) or len(set.union(*classes)) != 1:
                    continue
                for a in members:
                    s = int(Decimal(row['binomial_score']) * 1000)
                    q = int(Decimal(row['q_value']) * 100000)
                    sv = scores.setdefault(a, [0] * len(files))
                    qv = qs.setdefault(a, [100000] * len(files))
                    sv[run] = max(sv[run], s)
                    qv[run] = min(qv[run], q)
    methods = {
        'separate_union_min_q': {a: -min(v) for a, v in qs.items()},
        'max_protein_score': {a: max(v) for a, v in scores.items()},
        'sum_top_two_protein_scores': {a: sum(sorted(v)[-2:]) for a, v in scores.items()},
        'sum_all_protein_scores': {a: sum(v) for a, v in scores.items()},
    }
    result = {}
    for method, ranking in methods.items():
        blocks = defaultdict(list)
        for a, s in ranking.items():
            blocks[s].append(a)
        n = y = 0
        best = {'.01': (0, 0), '.03': (0, 0), '.05': (0, 0)}
        for s in sorted(blocks, reverse=True):
            for a in blocks[s]:
                y += '_YEAST' in a
                n += '_HUMAN' in a or '_ECOLI' in a
            for threshold in best:
                # Exact rational comparison avoids floating-point threshold ties.
                if n and y * (counts['HUMAN'] + counts['ECOLI']) * 100 <= int(Decimal(threshold) * 100) * n * counts['YEAST']:
                    best[threshold] = max(best[threshold], (n, y))
        result[method] = {str(float(t)): dict(present=n, yeast=y) for t, (n, y) in best.items()}
    return result


def main():
    counts = json.loads((OLD / 'frozen.json').read_text())['fasta_counts']
    groups = defaultdict(list)
    with (OLD / 'inputs.tsv').open() as f:
        for r in csv.DictReader(f, delimiter='\t'):
            groups[int(r['replicate_order'])].append(OLD / 'searches' / ('full_single_' + r['run']) / 'seed0/proteins.tsv')
    results = {str(rep): evaluate(sorted(files), counts) for rep, files in sorted(groups.items())}
    pilot = json.loads((OUT / 'report.json').read_text())['protein_pooling']
    for method, thresholds in pilot.items():
        for threshold, point in thresholds.items():
            assert results['1'][method][threshold] == {k: point[k] for k in ('present', 'yeast')}
    report = dict(pilot_integer_verification='passed for all four rankings and three thresholds',
        scope='Four disjoint five-run IonStar subsets defined by frozen replicate_order; one run per condition each. Rules fixed before checking subsets 2-4. Same dataset, not external validation.',
        cohorts=results,
        hashes={str(p): sha(p) for p in [Path(__file__), OLD / 'frozen.json', OLD / 'inputs.tsv',
                                       *itertools.chain.from_iterable(groups.values())]})
    (OUT / 'validation.json').write_text(json.dumps(report, indent=2) + '\n')
    print(report['pilot_integer_verification'])
    for rep, methods in results.items():
        print(rep, {name: vals['0.01'] for name, vals in methods.items()})


if __name__ == '__main__':
    main()
