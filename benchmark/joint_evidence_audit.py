#!/usr/bin/env python3
"""Read-only five-run feasibility audit. Candidate recurrence is NOT identification.

No searches, fitted models, or production changes. Protein score pooling is an
exploratory accession-level diagnostic on already inferred per-run protein groups.
"""
import csv
import hashlib
import itertools
import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'benchmark/results/lfq_matrix_fixes'
OUT = BASE / 'joint_evidence_audit'
OLD = Path('/home/patrick-garrett/Repos/uno/benchmark/results/matrix_20260912')


def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def species(proteins, decoy=False):
    if decoy:
        return 'DECOY'
    kinds = {s for p in proteins.split(';') for s in ('HUMAN', 'ECOLI', 'YEAST') if '_' + s in p}
    if kinds and kinds <= {'HUMAN', 'ECOLI'}:
        return 'PRESENT'
    return next(iter(kinds)) if len(kinds) == 1 else 'AMBIGUOUS'


def key(row, original=False):
    if original:
        return (row['feature_charge'], format(float(row['feature_mass']), '.5f'),
                format(float(row['feature_rt']), '.4f'))
    return row['charge'], row['mass'], row['rt']


def curve(scores, scale):
    """Whole tied-score blocks; threshold selected using entrapment for evaluation.

    These are empirical development curves, not deployable cohort q-values.
    """
    points = []
    present = yeast = 0
    ranked = sorted(scores.items(), key=lambda x: -x[1])
    for score, block in itertools.groupby(ranked, key=lambda x: x[1]):
        for acc, _ in block:
            kind = species(acc)
            if kind == 'PRESENT':
                present += 1
            elif kind == 'YEAST':
                yeast += 1
        if present:
            points.append((present, yeast, yeast * scale / present, score))
    result = {}
    for threshold in (.01, .03, .05):
        eligible = [p for p in points if p[2] <= threshold]
        best = max(eligible, default=(0, 0, None, None), key=lambda p: p[0])
        keep = {a for a, score in scores.items() if best[3] is not None and score >= best[3]
                and species(a) == 'PRESENT'}
        result[str(threshold)] = dict(present=best[0], yeast=best[1], empirical_fdp=best[2],
                                      score_threshold=best[3], accessions=sorted(keep))
    return result


def protein_audit(scale):
    scores, qs = {}, {}
    files = []
    for run in range(5):
        path = BASE / f'E_batch/run{run:03}/proteins.tsv'
        files.append(path)
        with path.open() as f:
            for row in csv.DictReader(f, delimiter='\t'):
                if row['decoy'] != '0':
                    continue
                members = row['group_members'].split(';')
                # Same-species groups only, following the frozen study convention.
                kinds = {next((s for s in ('HUMAN', 'ECOLI', 'YEAST') if '_' + s in p), '') for p in members}
                if len(kinds) != 1 or '' in kinds:
                    continue
                for acc in members:
                    scores.setdefault(acc, [0.] * 5)[run] = max(
                        scores.get(acc, [0.] * 5)[run], float(row['binomial_score']))
                    qs.setdefault(acc, [1.] * 5)[run] = min(
                        qs.get(acc, [1.] * 5)[run], float(row['q_value']))
    methods = {
        'separate_union_min_q': {a: -min(q) for a, q in qs.items()},
        'max_protein_score': {a: max(s) for a, s in scores.items()},
        # Stored scores have three decimals. Preserve ties at that precision.
        'sum_top_two_protein_scores': {a: round(sum(sorted(s)[-2:]), 3) for a, s in scores.items()},
        'sum_all_protein_scores': {a: round(sum(s), 3) for a, s in scores.items()},
    }
    results = {name: curve(s, scale) for name, s in methods.items()}
    original = {a for a, q in qs.items() if min(q) <= .01 and species(a) == 'PRESENT'}
    for name, thresholds in results.items():
        for point in thresholds.values():
            accepted = set(point.pop('accessions'))
            point['new_vs_original_1pct_per_run_union'] = len(accepted - original)
            point['lost_vs_original_1pct_per_run_union'] = len(original - accepted)
            if name != 'separate_union_min_q':
                point['rescued_without_any_run_q_le_1pct'] = len(accepted - original)
    return results, files


def candidate_audit():
    # Peptide TSVs are unthresholded quantification candidate maps, without ML
    # scores or peptide error probabilities. Never interpret these as peptide IDs.
    support = {}
    coherent = defaultdict(dict)
    runs = []
    sources = [BASE / 'aligned5/feature_observations.tsv']
    with sources[0].open() as f:
        reader = csv.DictReader(f, delimiter='\t')
        for run, rows in itertools.groupby(reader, key=lambda r: int(r['run_id'])):
            lookup = {}
            collisions = 0
            for row in rows:
                k = key(row, True)
                if k in lookup:
                    collisions += 1
                    lookup[k] = None  # Rounded-coordinate collision: do not guess identity.
                else:
                    lookup[k] = int(row['consensus_id']) if row['consensus_id'] else -1
            stats = Counter(rounded_coordinate_collisions=collisions)
            path = BASE / f'E_batch/run{run:03}/peptides.tsv'
            sources.append(path)
            seen = set()
            with path.open() as pf:
                for row in csv.DictReader(pf, delimiter='\t'):
                    ident = (row['decoy'], row['peptide'])
                    if ident in seen:
                        raise ValueError('Duplicate sequence/decoy in run: modform-aware export required')
                    seen.add(ident)
                    kind = species(row['proteins'], row['decoy'] == '1')
                    prior = support.get(ident, (kind, 0))
                    assert prior[0] == kind
                    support[ident] = kind, prior[1] | (1 << run)
                    stats[kind + '_candidates'] += 1
                    k = key(row)
                    if k not in lookup:
                        stats['unmatched_coordinates'] += 1
                        continue
                    cid = lookup[k]
                    if cid is None:
                        stats['ambiguous_coordinates'] += 1
                    elif cid < 0:
                        stats[kind + '_not_in_any_consensus'] += 1
                    else:
                        cm = coherent[ident]
                        cm[cid] = cm.get(cid, 0) | (1 << run)
            runs.append(dict(run=run, **stats))
            print('candidate audit run', run, dict(stats), flush=True)
    summary = {}
    for kind in ('PRESENT', 'YEAST', 'DECOY', 'AMBIGUOUS'):
        hist, same = Counter(), Counter()
        for ident, (s, mask) in support.items():
            if s != kind:
                continue
            hist[mask.bit_count()] += 1
            cm = coherent.get(ident, {})
            same[max((m.bit_count() for m in cm.values()), default=0)] += 1
        summary[kind] = dict(distinct_candidates=sum(hist.values()),
            number_of_runs_histogram=dict(sorted(hist.items())),
            largest_same_consensus_run_count_histogram=dict(sorted(same.items())))
    return dict(runs=runs, summary=summary), sources


def main():
    OUT.mkdir(exist_ok=True)
    frozen = OLD / 'frozen.json'
    counts = json.loads(frozen.read_text())['fasta_counts']
    scale = (counts['HUMAN'] + counts['ECOLI']) / counts['YEAST']
    proteins, pfiles = protein_audit(scale)
    (OUT / 'protein-pooling.json').write_text(json.dumps(proteins, indent=2) + '\n')
    print('protein pooling', json.dumps(proteins), flush=True)
    candidates, cfiles = candidate_audit()
    report = dict(scope='Five-run IonStar, seed zero, offline development audit',
                  candidate_tables_are_unthresholded=True, empirical_scale=scale,
                  protein_pooling=proteins, candidate_recurrence=candidates,
                  limitations=[
                      'No peptide scores/probabilities in saved candidate tables.',
                      'Protein scores already underwent per-run inference and competition.',
                      'Protein groups expanded to accessions; not independent protein groups.',
                      'Pooled score is a ranking statistic, not a calibrated probability.',
                      'Entrapment-selected thresholds are exploratory evaluation only.',
                      'Same aligned feature can support multiple wrong peptide sequences.',
                      'Candidate identity is sequence/decoy; script rejects within-run duplicates.',
                  ],
                  hashes={str(p): sha(p) for p in [Path(__file__), frozen, *pfiles, *cfiles]})
    (OUT / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print('candidate summary', json.dumps(candidates['summary']), flush=True)


if __name__ == '__main__':
    main()
