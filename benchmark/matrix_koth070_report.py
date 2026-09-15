#!/usr/bin/env python3
"""Evaluate the fixed two-cohort consensus-row study for the manuscript."""
import collections
import csv
import json
import math
from pathlib import Path
import shutil
from statistics import median

from matrix_report import load_proteins, quant_metrics
from matrix_release_report import count, independent_check
from joint_evidence_audit import sha

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'benchmark/results/matrix_20260915_koth070'


def read(path):
    return json.loads(path.read_text())


def cell_summary(path):
    manifest = read(path/'matrix_manifest.json')
    assert manifest['schema_version'] == 3
    for key in ['cells', 'observations']:
        assert sha(path/manifest[key]['path']) == manifest[key]['sha256']
    total = collections.Counter()
    runs = collections.defaultdict(collections.Counter)
    with (path/'lfq_matrix.long.tsv').open() as f:
        reader = csv.reader(f, delimiter='\t')
        columns = {name: i for i, name in enumerate(next(reader))}
        for row in reader:
            if row[columns['extraction_kind']] != 'target':
                continue
            c = runs[row[columns['run_id']]]
            c['cells'] += 1
            intensity = float(row[columns['lfq_intensity']])
            ownership = row[columns['lfq_ownership_status']]
            c[ownership] += 1
            if ownership in ['shared_signal_only', 'ambiguous_residual']:
                assert intensity == 0
            if intensity > 0:
                c['positive'] += 1
                kind = 'mbr' if row[columns['lfq_is_mbr']] == 'true' else 'direct'
                c[kind+'_positive'] += 1
                q = row[columns['lfq_q_value']]
                if q and float(q) <= .01:
                    c['confident'] += 1
                    c[kind+'_confident'] += 1
    for c in runs.values():
        total.update(c)
    assert total['cells'] == manifest['n_consensus'] * len(manifest['runs'])
    return dict(groups=manifest['n_consensus'], totals=dict(total),
                runs={r['name']: dict(runs[str(r['run_id'])]) for r in manifest['runs']})


def main():
    audit = read(OUT/'execution-audit.json')
    expected = dict(extraction=23, alignment=2, search=75, timing=4,
                    validation=3, pilot_alignment=1, pilot_search=1)
    assert audit['status'] == 'execution_verified' and audit['counts'] == expected
    plan = read(OUT/'plan.json')
    assert plan['arms'] == ['separate', 'consensus']
    sources = {str(OUT/n): sha(OUT/n) for n in ['execution-audit.json', 'frozen.json',
               'plan.json', 'complete.json', 'timing-audit.json']}
    for entry in audit['jobs']:
        p = Path(entry['receipt'])
        assert sha(p) == entry['sha256']
        r = read(p)
        assert r['returncode'] == 0 and '--matrix-joint' not in r['signature']['argv']
        # Recheck all scientific outputs used by the replacement analysis.
        for file, digest in r['outputs'].items():
            p = Path(file)
            if p.name in ['proteins.tsv', 'abundance.tsv']:
                assert sha(p) == digest
                sources[str(p)] = digest
    frozen = read(OUT/'frozen.json')
    result = dict(status='complete', study='KOTH 0.7.0 consensus rows',
                  publication_seeds=[0, 1, 2], arms=plan['arms'], job_counts=expected,
                  provenance=dict(uno_commit=frozen['sources']['uno']['head'],
                      koth_commit=frozen['sources']['koth_rust']['head'],
                      binaries={n: frozen['runtime_inputs'][str(OUT/'bin'/n)] for n in ['uno','koth_ff','koth_align']},
                      runtime_library=frozen['runtime_inputs'][str(OUT/'lib/lib_lightgbm.so')],
                      driver_revisions=frozen['driver_revisions']),
                  timing=dict(usable_for_speed_claim=False, audit=read(OUT/'timing-audit.json')),
                  datasets={})
    for cohort in plan['cohorts']:
        ds, name = cohort['dataset'], cohort['cohort']
        d = OUT/ds
        fc = collections.Counter({sp: 0 for sp in ['HUMAN','ECOLI','YEAST']})
        for line in (d/'database.fasta').open():
            if line.startswith('>'):
                species = line.split()[0].rsplit('_', 1)[-1]
                assert species in fc
                fc[species] += 1
        scale = (fc['HUMAN']+fc['ECOLI'])/fc['YEAST']
        cd = dict(runs=[r['run'] for r in cohort['inputs']], seeds={}, cells=cell_summary(d/'aligned/publication'))
        result['datasets'][ds] = dict(fasta_counts=dict(fc), scale=scale, cohorts={name: cd})
        sources[str(d/'database.fasta')] = sha(d/'database.fasta')
        sources[str(d/'aligned/publication/matrix_manifest.json')] = sha(d/'aligned/publication/matrix_manifest.json')
        for seed in plan['seeds']:
            paths = [d/'search'/f"{ds}_search_separate_{r['run']}_seed{seed}"/'proteins.tsv' for r in cohort['inputs']]
            consensus = d/'search'/f'{ds}_search_consensus_publication_seed{seed}'
            sd = dict(union=count([load_proteins(p) for p in paths], scale, True),
                      consensus=count([load_proteins(consensus/'proteins.tsv')], scale))
            independent_check(paths, fc, sd['union'], True)
            independent_check([consensus/'proteins.tsv'], fc, sd['consensus'])
            base = set(sd['union']['0.01']['accessions'])
            keep = set(sd['consensus']['0.01']['accessions'])
            sd['consensus']['comparison'] = dict(added=sorted(keep-base), lost=sorted(base-keep),
                 intersection=len(keep & base), jaccard=len(keep & base)/len(keep | base) if keep | base else None)
            if ds == 'ionstar':
                sd['consensus']['quant'] = quant_metrics(consensus/'abundance.tsv')
            cd['seeds'][str(seed)] = sd
        print(ds, {arm: [s[arm]['0.01']['present'] for s in cd['seeds'].values()] for arm in ['union','consensus']}, flush=True)
    for filename in ['matrix_koth070_report.py','matrix_report.py','matrix_release_report.py','joint_evidence_audit.py']:
        sources[str(ROOT/'benchmark'/filename)] = sha(ROOT/'benchmark'/filename)
    result['source_hashes'] = sources
    result['limitations'] = ['Development cohorts; no unseen validation or low-input/scaling claim.',
        'Entrapment thresholds are retrospective accession-level evaluations, not deployable q-values.',
        'Quantification is conditional on internal-q accepted protein groups and observed positive cells; no imputation.',
        'Timing controls overlapped other work and are excluded from speed claims.',
        'Other raw-mode and IM manuscript evidence retains its previous version; those reruns are separate.']
    payload = json.dumps(result, indent=2, allow_nan=False)+'\n'
    (OUT/'paper-report.json').write_text(payload)
    dest = ROOT/'paper/figdata/koth_matrix_study.json'
    archive = OUT/'previous-paper-evidence'
    archive.mkdir(exist_ok=True)
    if not (archive/'koth_matrix_study_060.json').exists():
        shutil.copy2(dest, archive/'koth_matrix_study_060.json')
    dest.write_text(payload)
    print('Wrote audited replacement evidence:', dest, flush=True)


if __name__ == '__main__':
    main()
