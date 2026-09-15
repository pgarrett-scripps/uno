#!/usr/bin/env python3
"""Independent float evaluator cross-check against the exact-integer ensemble audit."""
import csv
import hashlib
import json
from pathlib import Path
from joint_evidence_audit import curve

ROOT = Path(__file__).resolve().parents[1]
STUDY = ROOT / 'benchmark/results/matrix_20260913_koth050'
OUT = ROOT / 'benchmark/results/matrix_seed_ensemble_20260913'
report = json.loads((OUT / 'report.json').read_text())
original = json.loads((STUDY / 'report.json').read_text())
cache = {}
checked = 0
for dataset, cohorts in report['datasets'].items():
    for cohort, cd in cohorts.items():
        scale = (cd['fasta_counts']['HUMAN'] + cd['fasta_counts']['ECOLI']) / cd['fasta_counts']['YEAST']
        runs = original['datasets'][dataset]['cohorts'][cohort]['runs']
        for arm, methods in cd['arms'].items():
            for method, results in methods.items():
                for subset in results['subsets'].values():
                    names = ['single_' + r for r in runs] if arm == 'union' else [arm + '_' + cohort]
                    scores = {}
                    for seed in subset['seeds']:
                        for name in names:
                            p = STUDY / dataset / 'searches' / name / f'seed{seed}' / 'proteins.tsv'
                            if p not in cache:
                                cache[p] = list(csv.DictReader(p.open(), delimiter='\t'))
                            for row in cache[p]:
                                if row['decoy'] != '0':
                                    continue
                                members = row['group_members'].split(';')
                                species = [{s for s in ('HUMAN', 'ECOLI', 'YEAST') if '_' + s in a} for a in members]
                                if any(len(s) != 1 for s in species) or len(set.union(*species)) != 1:
                                    continue
                                score = -float(row['q_value']) if method == 'min_q' else float(row['binomial_score'])
                                for acc in members:
                                    scores[acc] = max(scores.get(acc, float('-inf')), score)
                    independent = curve(scores, scale)
                    for threshold, expected in subset['thresholds'].items():
                        observed = independent[threshold]
                        for field in ('present', 'yeast', 'accessions'):
                            assert observed[field] == expected[field], (dataset, cohort, arm, method, subset['seeds'], threshold, field)
                        checked += 1
assert checked == 1008

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

# Explicitly test whole ties and exact 1% boundary using the production evaluator.
from decimal import Decimal
from matrix_seed_ensemble import evaluate
scores = {f'p{i}_HUMAN': Decimal(1) for i in range(99)}
scores['y_YEAST'] = Decimal(1)
assert evaluate(scores, dict(HUMAN=1, ECOLI=0, YEAST=1))['0.01']['present'] == 0
scores['p99_HUMAN'] = Decimal(1)
assert evaluate(scores, dict(HUMAN=1, ECOLI=0, YEAST=1))['0.01']['present'] == 100
result = dict(status='passed', checked_thresholds=checked, unique_search_tables=len(cache),
              checks=['Independent original float curve matches Decimal/integer ensemble counts and present accession sets at all 1%, 3%, 5% cutoffs.',
                      'Whole tied-score block rejection and exact 1% acceptance boundary passed.'],
              report_sha256=sha(OUT / 'report.json'), checker_sha256=sha(Path(__file__)))
(OUT / 'validation.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
