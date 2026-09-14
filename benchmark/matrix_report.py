#!/usr/bin/env python3
"""Summarize only the fresh matrix study, including honest cohort unions."""
from __future__ import annotations
import csv
import json
import math
from pathlib import Path
from statistics import median
import re
import numpy as np

from matrix_study import STUDY, cohort_rows, sha, write_json


def species(accession):
    found = [s for s in ['HUMAN', 'ECOLI', 'YEAST'] if '_' + s in accession]
    if len(found) != 1:
        return 'AMBIGUOUS'
    return found[0]


def load_proteins(path):
    rows = []
    for r in csv.DictReader(open(path), delimiter='\t'):
        if r['decoy'] != '0': continue
        members = r.get('group_members', r['accession']).split(';')
        classes = {species(p) for p in members}
        # Mixed present/absent groups cannot be unambiguously scored as false.
        r['species'] = next(iter(classes)) if len(classes) == 1 else 'AMBIGUOUS'
        r['q'] = float(r['q_value'])
        r['score'] = float(r['binomial_score'])
        r['members'] = members
        rows.append(r)
    return rows


def list_counts(rows, scale):
    counts = {s: sum(r['species'] == s for r in rows)
              for s in ['HUMAN', 'ECOLI', 'YEAST', 'AMBIGUOUS']}
    present = counts['HUMAN'] + counts['ECOLI']
    return dict(counts, present=present,
                entrapment_fdp=counts['YEAST'] * scale / present if present else None)


def evaluate(rows, scale):
    result = {'scored': list_counts(rows, scale),
              'internal_1pct': list_counts([r for r in rows if r['q'] <= .01], scale)}
    ranked = sorted((r for r in rows if r['species'] != 'AMBIGUOUS'),
                    key=lambda r: (-r['score'], r['species'] != 'YEAST', r['accession']))
    present = absent = 0
    points = []
    for r in ranked:
        absent += r['species'] == 'YEAST'
        present += r['species'] in ('HUMAN', 'ECOLI')
        if present: points.append((present, absent * scale / present))
    q = float('inf')
    depths = {'0.01': 0, '0.05': 0}
    for n, fdp in reversed(points):
        q = min(q, fdp)
        for t in depths:
            if q <= float(t): depths[t] = max(depths[t], n)
    result['entrapment_depth'] = depths
    return result


def accession_union(tables, threshold=.01):
    # Expand groups to accessions, then deduplicate. These counts deliberately
    # differ from group counts and are labeled at accession level in the report.
    union = {}
    for rows in tables:
        for r in rows:
            if r['q'] > threshold or r['species'] == 'AMBIGUOUS': continue
            for acc in r['members']:
                union[acc] = {'accession': acc, 'species': species(acc)}
    return list(union.values())


def accession_curve(tables, scale):
    return [dict(threshold=t, **list_counts(accession_union(tables,t),scale))
            for t in [.001,.002,.005,.01,.02,.05,.1]]


def compare_cohort(joint, singles):
    """Compare present-species accessions, preserving same-species group expansion."""
    def present(rows):
        return {r['accession'] for r in accession_union([rows]) if r['species'] in ('HUMAN','ECOLI')}
    a=present(joint);runs={run:present(rows) for run,rows in singles.items()}
    b=set().union(*runs.values())
    return {'joint_present':len(a),'separate_union_present':len(b),'intersection':len(a&b),
            'jaccard':len(a&b)/len(a|b) if a|b else None,
            'joint_only_accessions':sorted(a-b),'separate_only_accessions':sorted(b-a),
            'per_run':{run:{'separate_present':len(r),'joint_overlap':len(a&r),
                            'separate_only':len(r-a),'fraction_recovered':len(a&r)/len(r) if r else None}
                       for run,r in runs.items()}}


def resources(name):
    p = STUDY / 'logs' / (name + '.json')
    if not p.exists(): return {}
    j = json.loads(p.read_text())
    out = {'wall_seconds': j['elapsed'], 'exit_code': j['returncode']}
    timefile = p.with_suffix('.time')
    if timefile.exists():
        for line in timefile.read_text().splitlines():
            for key, target in [('Maximum resident set size (kbytes)', 'peak_rss_kb'),
                                ('User time (seconds)', 'user_seconds'),
                                ('System time (seconds)', 'system_seconds')]:
                if key + ':' in line: out[target] = float(line.rsplit(':', 1)[1])
    return out


def quant_metrics(path, fixed_factors=None, common_groups=None):
    """No imputation; label-independent median-ratio normalization, then truth check."""
    with open(path) as f:
        reader=csv.DictReader(f,delimiter='\t');header=reader.fieldnames;rows=list(reader)
    if common_groups is not None:
        rows=[r for r in rows if tuple(sorted(r['group_members'].split(';'))) in common_groups]
    runs=header[4:]
    if len(runs)<2 or not rows: return None
    conditions=[re.search(r'human_ecoli_([A-E])_',r).group(1) for r in runs]
    values=np.array([[float(r[c]) for c in runs] for r in rows])
    positive=values>0
    logged=np.full_like(values,np.nan)
    np.log2(values,out=logged,where=positive)
    complete=positive.all(axis=1)
    if fixed_factors is not None:
        factors=np.array([fixed_factors[r] for r in runs])
    elif complete.any():
        residual=logged[complete]-logged[complete].mean(axis=1,keepdims=True)
        factors=np.median(residual,axis=0);factors-=np.mean(factors)
    else: factors=np.zeros(len(runs))
    normalized=values/np.exp2(factors)
    classes=[]
    for r in rows:
        tags={species(p) for p in r['group_members'].split(';')}
        classes.append(next(iter(tags)) if len(tags)==1 else 'AMBIGUOUS')
    result={'normalization':'baseline median-ratio factors; no imputation' if fixed_factors is not None else 'median ratios on complete protein groups; no imputation',
            'size_factors_log2':dict(zip(runs,factors.tolist())),
            'complete_groups_for_normalization':int(complete.sum()),'species':{}}
    for sp in ['HUMAN','ECOLI','YEAST']:
        sel=np.array([c==sp for c in classes]);v=normalized[sel];pos=positive[sel]
        cvs=[]
        for condition in sorted(set(conditions)):
            cols=np.array([c==condition for c in conditions])
            if cols.sum()<2: continue
            for line in v[:,cols]:
                observed=line[line>0]
                if len(observed)>=2: cvs.append(float(np.std(observed,ddof=1)/np.mean(observed)))
        errors=[]
        if 'A' in conditions and 'E' in conditions:
            acols=np.array([c=='A' for c in conditions]);ecols=np.array([c=='E' for c in conditions])
            for line in v:
                a=line[acols];a=a[a>0];e=line[ecols];e=e[e>0]
                if len(a)>=2 and len(e)>=2:
                    errors.append(float(np.log2(np.median(e)/np.median(a))))
        expected=np.log2(3) if sp=='ECOLI' else 0
        curve=[]
        if 'A' in conditions:
            acols=np.array([c=='A' for c in conditions])
            for condition in sorted(set(conditions)):
                cols=np.array([c==condition for c in conditions]);ratios=[]
                for line in v:
                    a=line[acols];a=a[a>0];b=line[cols];b=b[b>0]
                    if len(a)>=2 and len(b)>=2:
                        ratios.append(float(np.log2(np.median(b)/np.median(a))))
                truth=float(np.log2(dict(zip('ABCDE',[1,1.5,2,2.5,3]))[condition])) if sp=='ECOLI' else 0.
                curve.append({'condition':condition,'groups':len(ratios),
                              'expected_log2_ratio':truth,'median_log2_ratio':median(ratios) if ratios else None})
        result['species'][sp]={'groups':int(sel.sum()),
            'condition_ratio_curve':curve,
            'positive_cell_fraction':float(pos.mean()) if pos.size else None,
            'median_replicate_cv':median(cvs) if cvs else None,
            'replicate_cv_group_condition_pairs':len(cvs),
            'extreme_ratio_groups':len(errors),
            'median_log2_E_over_A':median(errors) if errors else None,
            'median_abs_log2_ratio_error':median([abs(x-expected) for x in errors]) if errors else None}
    return result


def report():
    frozen = json.loads((STUDY / 'frozen.json').read_text())
    n = frozen['fasta_counts']
    scale = (n['HUMAN'] + n['ECOLI']) / n['YEAST']
    searches = {}
    tables = {}
    for path in sorted((STUDY / 'searches').glob('*/seed*/proteins.tsv')):
        name, seed = path.parent.parent.name, path.parent.name
        resource = resources(f'search_{name}_{seed}')
        if resource.get('exit_code') != 0:
            continue
        tables[name, seed] = load_proteins(path)
        searches[f'{name}/{seed}'] = dict(evaluate(tables[name, seed], scale),
            resources=resource, file_sha256=sha(path),
            accession_level_1pct=list_counts(accession_union([tables[name, seed]]), scale),
            accession_calibration_curve=accession_curve([tables[name, seed]],scale))
        if ('matrix' in name or name.startswith(('pilot_', 'full_quality', 'full_masked', 'full_profile'))) and 'run_evidence_null' not in name:
            abundance=path.with_name('abundance.tsv')
            if abundance.exists(): searches[f'{name}/{seed}']['quantification']=quant_metrics(abundance)
    unions = {}
    for cohort in ['pilot', '5', '10', '20']:
        for seed in ['seed0', 'seed1', 'seed2']:
            prefix='full_single_' if all(('full_single_'+r['run'],seed) in tables for r in cohort_rows(cohort)) else 'single_'
            keys = [(prefix + r['run'], seed) for r in cohort_rows(cohort)]
            if all(k in tables for k in keys):
                u = accession_union([tables[k] for k in keys])
                unions[f'{cohort}/{seed}'] = dict(list_counts(u, scale), source_prefix=prefix,
                    accession_calibration_curve=accession_curve([tables[k] for k in keys],scale))
    for seed in ['seed0','seed1','seed2']:
        names=['full_matrix_20','full_masked_20','full_quality_null_20','full_quality_20','full_profile_reassign_20']
        if not all(n+'/'+seed in searches for n in names):continue
        paths={n:STUDY/'searches'/n/seed/'abundance.tsv' for n in names}
        keys=[]
        for path in paths.values():
            with open(path) as f:
                keys.append({tuple(sorted(r['group_members'].split(';'))) for r in csv.DictReader(f,delimiter='\t')})
        common=set.intersection(*keys)
        base=quant_metrics(paths['full_matrix_20'])
        if base is None:continue
        factors=base['size_factors_log2']
        for name,path in paths.items():
            searches[name+'/'+seed]['quantification']=quant_metrics(path,fixed_factors=factors)
            searches[name+'/'+seed]['quantification_common_groups']=quant_metrics(path,fixed_factors=factors,common_groups=common)
    comparisons={}
    for (name,seed),rows in tables.items():
        if not name.startswith('full_') or name.startswith('full_single_'):continue
        cohort=next((n for n in ['5','10'] if name=='full_matrix_'+n),'20')
        singles={r['run']:tables.get(('full_single_'+r['run'],seed)) for r in cohort_rows(cohort)}
        if all(r is not None for r in singles.values()):
            comparisons[name+'/'+seed]=compare_cohort(rows,singles)
    equivalence = {}
    for (name, seed), raw in tables.items():
        if not name.startswith('raw_'): continue
        single = tables.get(('single_' + name[4:], seed))
        if single is None: continue
        a = {r['accession'] for r in raw if r['q'] <= .01}
        b = {r['accession'] for r in single if r['q'] <= .01}
        equivalence[name + '/' + seed] = {'raw': len(a), 'tsv': len(b),
            'intersection': len(a & b), 'raw_only': sorted(a-b), 'tsv_only': sorted(b-a),
            'jaccard': len(a & b)/len(a | b) if a | b else None,
            'byte_identical': sha(STUDY/'searches'/name/seed/'proteins.tsv') ==
                              sha(STUDY/'searches'/('single_' + name[4:])/seed/'proteins.tsv')}
    summary = {'status': 'preliminary; only completed searches', 'scale': scale,
               'searches': searches, 'accession_level_unions_at_internal_1pct': unions,
               'route_equivalence': equivalence,
               'joint_vs_separate':comparisons,
               'warning': 'Union lists are empirically rescored; not assumed to control cohort FDR at 1%.'}
    write_json(STUDY / 'summary.json', summary)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__': report()
