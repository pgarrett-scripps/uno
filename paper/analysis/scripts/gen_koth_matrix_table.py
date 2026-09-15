#!/usr/bin/env python3
"""Three compact SI tables for the audited 0.7.0 consensus-row comparison."""
import json
from statistics import median
import style
from _assets import record


def spread(values, digits=0, scale=1):
    vals = [v * scale for v in values if v is not None]
    if not vals:
        return 'NA'
    fmt = lambda v: f'{v:,.{digits}f}'
    return fmt(median(vals)) if min(vals) == max(vals) else f'{fmt(median(vals))} ({fmt(min(vals))}–{fmt(max(vals))})'


def table(name, headers, rows, desc):
    text = ['// Generated; do not edit by hand.', '#table(', f'  columns: {len(headers)},',
            '  align: (' + ', '.join(['left'] + ['right'] * (len(headers)-1)) + '),',
            '  table.header(repeat: true, ' + ', '.join('['+h+']' for h in headers) + '),']
    text += ['  '+', '.join('['+str(v)+']' for v in row)+',' for row in rows]
    text.append(')')
    path = 'si/koth_matrix_'+name+'.typ'
    (style.PAPER/path).write_text('\n'.join(text)+'\n')
    record('tbl.koth-matrix-'+name, path, kind='table', inputs=['figdata/koth_matrix_study.json'], desc=desc)


def main():
    d = json.loads((style.FIGDATA/'koth_matrix_study.json').read_text())
    assert d['status'] == 'complete' and d['arms'] == ['separate', 'consensus']
    # Retire the old generated joint/scaling tables. Their exact versions remain
    # in Git and the study's previous-paper-evidence archive.
    manifest = style.PAPER/'assets.json'
    assets = json.loads(manifest.read_text())
    for name in ['astral', 'additional', 'resources', 'sample-support', 'routes']:
        old = assets['values'].pop('tbl.koth-matrix-'+name, None)
        if old:
            (style.PAPER/old['path']).unlink(missing_ok=True)
    manifest.write_text(json.dumps(assets, indent=2)+'\n')
    rows = []
    for dataset, cohort, label in [('astral', '100ng', 'Astral 100 ng'), ('ionstar', '20', 'IonStar 20 runs')]:
        seeds = d['datasets'][dataset]['cohorts'][cohort]['seeds'].values()
        for arm, title in [('union', 'Separate union'), ('consensus', 'Consensus rows')]:
            vals = [s[arm] for s in seeds]
            rows.append([label+' / '+title, *[spread([v[t]['present'] for v in vals]) for t in ['0.01','0.03','0.05']],
                         spread([v['internal_1pct']['entrapment_fdp'] for v in vals], 2, 100)])
    table('identification', ['Cohort / input', '1% FDP', '3% FDP', '5% FDP', 'Internal-list FDP (%)'], rows,
          'Present-species accession depth at empirical cutoffs and FDP of internal 1% lists; median and full range over three seeds')
    seeds = d['datasets']['ionstar']['cohorts']['20']['seeds'].values()
    rows = []
    for sp, label in [('HUMAN','Human'),('ECOLI','_E. coli_')]:
        vals = [s['consensus']['quant']['species'][sp] for s in seeds]
        rows.append([label, spread([v['extreme_ratio_groups'] for v in vals]),
                     spread([v['median_abs_log2_ratio_error'] for v in vals], 2),
                     spread([v['median_replicate_cv'] for v in vals], 1, 100),
                     spread([v['positive_cell_fraction'] for v in vals], 1, 100)])
    table('quantification', ['Species', 'Ratio groups', 'Abs. log2 error', 'CV (%)', 'Positive cells (%)'], rows,
          'IonStar consensus-matrix quantities; internal-q accepted protein groups, median-ratio normalization, no imputation, three-seed medians and ranges')
    p = d['provenance']
    rows = [[name.replace('_','\\_'), p['uno_commit' if name == 'uno' else 'koth_commit'][:12], p['binaries'][name][:12]]
            for name in ['uno','koth_ff','koth_align']]
    rows.append(['LightGBM runtime','4.6.0',p['runtime_library'][:12]])
    table('versions', ['Program', 'Source revision / version', 'Binary SHA-256 prefix'], rows,
          'Frozen UNO and released KOTH 0.7.0 build identities; complete working-tree snapshots and driver amendments archived')


if __name__ == '__main__':
    main()
