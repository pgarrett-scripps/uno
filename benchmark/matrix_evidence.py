#!/usr/bin/env python3
"""Validate KOTH LFQ sidecars and generate explicit matrix evidence ablations."""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import warnings
import numpy as np
from matrix_study import STUDY, sha, write_json

QC_COLUMNS = [
    'matrix_confident_fraction', 'matrix_confident_direct_fraction',
    'matrix_confident_mbr_fraction', 'matrix_median_log_confidence',
    'matrix_median_isotope_pattern', 'matrix_median_coelution',
    'matrix_median_isotopes', 'matrix_median_abs_ppm', 'matrix_median_abs_rt',
]


def load_evidence(directory, threshold=.01):
    directory = Path(directory)
    with open(directory/'intensity_matrix.tsv') as f:
        reader = csv.DictReader(f, delimiter='\t')
        header = reader.fieldnames
        rows = list(reader)
    first = header.index('n_contributing_runs') + 1
    runs = header[first:]
    if len(set(runs)) != len(runs): raise ValueError('Duplicate run names')
    index = {r: i for i, r in enumerate(runs)}
    n, k = len(rows), len(runs)
    intensity = np.array([[float(r[c]) for c in runs] for r in rows])
    if not np.all(np.isfinite(intensity) & (intensity >= 0)):
        raise ValueError('Invalid intensity cells')
    qvalue = np.full((n,k), np.nan)
    mbr = np.zeros((n,k), dtype=bool)
    seen = np.zeros((n,k), dtype=bool)
    values = np.full((n,k,5), np.nan, dtype=np.float32)
    # Validate identity and coverage even when no row passes the q-value gate.
    with open(directory/'lfq_details.tsv') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            if r['is_decoy'] not in ('true', 'false'):
                raise ValueError('Invalid decoy flag')
            if r['is_decoy'] == 'true': continue
            i = int(r['feature_idx']); j = index[r['run_name']]
            if i < 0 or i >= n or seen[i,j]: raise ValueError('Duplicate/invalid feature-run cell')
            ref = rows[i]
            for a,b in [('massCalib','massCalib'),('mz','mz'),('charge','charge'),
                        ('rtApex','refRtApex'),('im','refIm'),('combined_score','seed_combined_score'),
                        ('seed_run','seed_run'),('n_contributing_runs','n_contributing_runs')]:
                if ref[a] != r[b]: raise ValueError(f'Metadata mismatch at feature {i}: {a}')
            if not np.isclose(float(r['intensity']), intensity[i,j], rtol=1e-6, atol=1e-12):
                raise ValueError(f'Intensity mismatch at feature {i}, run {j}')
            q = float(r['q_value'])
            if not np.isfinite(q) or not 0 <= q <= 1: raise ValueError('Invalid q-value')
            if r['is_mbr'] not in ('true','false'): raise ValueError('Invalid MBR flag')
            qvalue[i,j], mbr[i,j], seen[i,j] = q, r['is_mbr']=='true', True
            values[i,j] = [float(r[c]) for c in ['spectral_bhattacharyya','coelution',
                                                'n_isotopes_found','ppm_error','rt_diff']]
    if not seen.all(): raise ValueError('Incomplete target cell coverage')
    expected = np.array([int(r['n_contributing_runs']) for r in rows])
    if not np.array_equal((~mbr).sum(axis=1), expected):
        raise ValueError('Contributing-run counts disagree with per-cell MBR flags')
    # Matrix q-values are printed to four decimals; details retain six.
    with open(directory/'qvalue_matrix.tsv') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if reader.fieldnames != header: raise ValueError('Q-value matrix header mismatch')
        count = 0
        for i,r in enumerate(reader):
            if i >= n: raise ValueError('Extra q-value row')
            if any(r[c] != rows[i][c] for c in header[:first]):
                raise ValueError('Q-value matrix row identity mismatch')
            q = np.array([float(r[c]) for c in runs])
            if not np.all(np.abs(q-qvalue[i]) <= .000051):
                raise ValueError('Q-value matrix/details disagreement')
            count += 1
        if count != n: raise ValueError('Q-value matrix row count mismatch')
    positive = intensity > 0
    confident = positive & (qvalue <= threshold)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)  # all-missing row => explicit zero + support fraction
        med = np.nanmedian(np.where(confident[:,:,None], np.abs(values), np.nan), axis=1)
    med = np.nan_to_num(med, nan=0, posinf=0, neginf=0)
    evidence = np.column_stack([confident.mean(axis=1), (confident & ~mbr).mean(axis=1),
        (confident & mbr).mean(axis=1), np.median(-np.log10(np.maximum(qvalue,1e-6)),axis=1), med])
    return header, rows, intensity, confident, evidence, {
        'rows': n, 'runs': k, 'positive_cells': int(positive.sum()),
        'direct_cells': int((~mbr).sum()), 'confident_cells': int(confident.sum()),
        'confident_direct_cells': int((confident & ~mbr).sum()),
        'confident_mbr_cells': int((confident & mbr).sum()),
        'fraction_positive_confident': float(confident.sum()/max(positive.sum(),1)),
        'per_run': {run: {'original_contribution_rows':int((~mbr[:,j]).sum()),
                         'positive_cells':int(positive[:,j].sum()),
                         'zero_cells':int((~positive[:,j]).sum()),
                         'confident_direct_cells':int((confident[:,j]&~mbr[:,j]).sum()),
                         'confident_mbr_cells':int((confident[:,j]&mbr[:,j]).sum())}
                    for j,run in enumerate(runs)},
        'q_threshold': threshold, 'q_source': 'lfq_details.tsv (six-decimal precision)',
        'sources': {name:sha(directory/name) for name in
                    ['intensity_matrix.tsv','lfq_details.tsv','qvalue_matrix.tsv']},
    }, mbr


def prepare(cohort):
    directory = STUDY/'aligned'/cohort
    header, rows, intensity, confident, evidence, summary, mbr = load_evidence(directory)
    out = STUDY/'evidence'/cohort
    out.mkdir(parents=True,exist_ok=True)
    pivot = header.index('n_contributing_runs')
    np.savez_compressed(out/'cell_evidence.npz', confident=confident, direct=~mbr,
        mass=np.array([float(r['massCalib']) for r in rows]),
        rt=np.array([float(r['rtApex']) for r in rows]),
        charge=np.array([int(r['charge']) for r in rows]),
        run_names=np.array(header[pivot+1:]))
    for mode in ['masked', 'quality', 'quality_null', 'run_evidence_null']:
        enriched = mode in ('quality', 'quality_null')
        columns = header[:pivot] + (QC_COLUMNS if enriched else []) + header[pivot:]
        with open(out/(mode+'.tsv'),'w') as f:
            writer = csv.DictWriter(f,fieldnames=columns,delimiter='\t')
            writer.writeheader()
            for i,r in enumerate(rows):
                r = dict(r)
                if mode == 'run_evidence_null':
                    # Keep precisely the baseline floor-2 rows; then suppress
                    # recurrence and intensity variation. Search with floor 1.
                    # This ablation is ID-only: its quantities are artificial.
                    if int(r['n_contributing_runs']) < 2: continue
                    r['n_contributing_runs'] = '1'
                    for c in header[pivot+1:]: r[c] = '1'
                if enriched:
                    r.update({c:format(float(v),'.8g') for c,v in zip(QC_COLUMNS,
                        evidence[i] if mode=='quality' else np.zeros(len(QC_COLUMNS)))})
                elif mode == 'masked':
                    for j,c in enumerate(header[pivot+1:]):
                        if not confident[i,j]: r[c] = '0'
                writer.writerow(r)
    summary['outputs'] = {p.name:sha(p) for p in out.glob('*.tsv')}
    summary['cell_evidence_sha256'] = sha(out/'cell_evidence.npz')
    summary['run_evidence_null'] = 'floor-2 membership retained; use --matrix-min-runs 1; ID-only, also changes intensity tie-breaking'
    write_json(out/'summary.json',summary)
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cohort',default='pilot')
    prepare(p.parse_args().cohort)
