#!/usr/bin/env python3
"""Bounded genuine matrix-search pilot; never merges independent protein scores."""
import argparse
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

from joint_evidence_audit import curve, sha, species
from matrix_report import load_proteins, evaluate, accession_union, list_counts

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'benchmark/results/lfq_matrix_fixes'
OUT = BASE / 'joint_search_v1'
OLD = Path('/home/patrick-garrett/Repos/uno/benchmark/results/matrix_20260912')


def run(name, args):
    directory = OUT / name
    directory.mkdir(parents=True, exist_ok=True)
    binary = OUT / 'uno'
    argv = [str(binary), *map(str, args), '--fasta', str(OLD / 'human_ecoli_yeast.fasta'),
            '--rt-model', str(OLD / 'rt.urtc'), '--seed', '0',
            '--out-proteins', str(directory / 'proteins.tsv'),
            '--out-lfq', str(directory / 'lfq.tsv'), '--out-abundance', str(directory / 'abundance.tsv')]
    signature = dict(argv=argv, binary=sha(binary), fasta=sha(OLD / 'human_ecoli_yeast.fasta'),
                     rt_model=sha(OLD / 'rt.urtc'), config=sha(Path(argv[argv.index('--config') + 1])),
                     matrix_manifest=sha(BASE / 'aligned5/matrix_manifest.json'))
    if '--matrix-long' in argv:
        signature['matrix_manifest'] = sha(Path(argv[argv.index('--matrix-long') + 1]) / 'matrix_manifest.json')
    if '--features' in argv:
        signature['features'] = sha(Path(argv[argv.index('--features') + 1]))
    stamp = directory / 'job.json'
    if stamp.exists():
        prior = json.loads(stamp.read_text())
        assert prior['signature'] == signature and prior['returncode'] == 0, name
        assert all(sha(Path(p)) == h for p, h in prior['outputs'].items()), name
        return
    env = dict(os.environ, RAYON_NUM_THREADS='8', OMP_NUM_THREADS='8', OPENBLAS_NUM_THREADS='1',
               RUST_LOG='info', UNO_RT_CACHE='1', UNO_RT_CACHE_DIR=str(BASE / 'rt-cache'))
    # Ambient diagnostic flags must not silently change this pilot.
    env = {k: v for k, v in env.items() if not k.startswith('UNO_') or k in ('UNO_RT_CACHE', 'UNO_RT_CACHE_DIR')}
    started = time.time()
    print('Starting', name, flush=True)
    with (directory / 'run.log').open('w') as log:
        p = subprocess.run(['/usr/bin/time', '-v', '-o', str(directory / 'time.txt'), *argv],
                           env=env, stdout=log, stderr=subprocess.STDOUT)
    outputs = {str(directory / n): sha(directory / n) for n in ['proteins.tsv', 'lfq.tsv', 'abundance.tsv']
               if (directory / n).is_file()}
    stamp.write_text(json.dumps(dict(signature=signature, seconds=time.time() - started,
        returncode=p.returncode, outputs=outputs, started=started), indent=2) + '\n')
    assert p.returncode == 0, name
    print('Finished', name, round(time.time() - started, 1), flush=True)


def pilot():
    OUT.mkdir(parents=True, exist_ok=True)
    build = ROOT.parent / 'target/release/uno'
    if not (OUT / 'uno').exists():
        shutil.copy2(build, OUT / 'uno')
        (OUT / 'source-manifest.json').write_text(json.dumps({str(p.relative_to(ROOT)): sha(p)
            for p in sorted((ROOT / 'crates').rglob('*.rs'))}, indent=2) + '\n')
    config = BASE / 'uno.toml'
    ablation = OUT / 'no-agreement.toml'
    ablation.write_text(config.read_text() + '\n[scoring]\nmatrix_joint_agreement = false\n')
    common = ['--matrix-long', BASE / 'aligned5']
    run('joint', [*common, '--matrix-joint', '--config', config])
    run('no_agreement', [*common, '--matrix-joint', '--config', ablation])
    run('representative_control', [*common, '--matrix-min-runs', '2', '--matrix-prefer-isotope-support',
                                   '--matrix-legacy-values', '--config', config])
    report()


def guarded():
    OUT.mkdir(parents=True, exist_ok=True)
    if not (OUT / 'uno').exists():
        shutil.copy2(ROOT.parent / 'target/release/uno', OUT / 'uno')
        (OUT / 'source-manifest.json').write_text(json.dumps({str(p.relative_to(ROOT)): sha(p)
            for p in sorted((ROOT / 'crates').rglob('*.rs'))}, indent=2) + '\n')
    run('replicated_joint', ['--matrix-long', BASE / 'aligned5-replicated', '--matrix-joint', '--config', BASE / 'uno.toml'])
    run('joint', ['--matrix-long', BASE / 'aligned5', '--matrix-joint', '--config', BASE / 'uno.toml'])
    parity = {name: sha(OUT / 'joint' / name) == sha(BASE / 'joint_search_v2/joint' / name)
              for name in ['proteins.tsv', 'lfq.tsv', 'abundance.tsv']}
    assert all(parity.values()), parity
    (OUT / 'strong-group-parity.json').write_text(json.dumps(parity, indent=2) + '\n')
    report()


def report():
    counts = json.loads((OLD / 'frozen.json').read_text())['fasta_counts']
    scale = (counts['HUMAN'] + counts['ECOLI']) / counts['YEAST']
    results = {}
    for directory in sorted(OUT.iterdir()):
        if not (directory / 'job.json').is_file():
            continue
        # Alignment jobs may share the output tree with search jobs.
        if not (directory / 'proteins.tsv').is_file():
            continue
        rows = load_proteins(directory / 'proteins.tsv')
        scores = {}
        for row in rows:
            if row['species'] == 'AMBIGUOUS': continue
            for a in row['members']: scores[a] = max(scores.get(a, 0.0), row['score'])
        accession_curve = curve(scores, scale)
        for v in accession_curve.values(): v.pop('accessions')
        results[directory.name] = dict(groups=evaluate(rows, scale),
            internal_1pct_accessions=list_counts(accession_union([rows]), scale),
            empirical_accession_curve=accession_curve,
            seconds=json.loads((directory / 'job.json').read_text())['seconds'])
    control = OUT / 'representative_control'
    if (control / 'proteins.tsv').exists():
        parity = {name: sha(control / name) == sha(BASE / 'F_final_validation' / name)
                  for name in ['proteins.tsv', 'lfq.tsv', 'abundance.tsv']}
        assert all(parity.values()), parity
        results['legacy_parity'] = parity
    reference = json.loads((BASE / 'joint_evidence_audit/report.json').read_text())
    results['separate_union_reference'] = reference['protein_pooling']['separate_union_min_q']['0.01']
    results['scope'] = 'Actual long-matrix searches, one cohort result per arm; no individual-search results used as input.'
    singles = [r for name, r in results.items() if name.startswith('single_run')]
    if len(singles) == 5:
        results['timing_comparison'] = dict(individual_sum_seconds=sum(r['seconds'] for r in singles),
            joint_seconds=results['joint']['seconds'],
            search_speed_ratio=sum(r['seconds'] for r in singles) / results['joint']['seconds'],
            conditions='Serial runs, eight threads, same FASTA/model/config/output types, warmed prediction cache. Matrix alignment/feature finding excluded.')
    (OUT / 'report.json').write_text(json.dumps(results, indent=2) + '\n')
    print(json.dumps(results, indent=2), flush=True)


def timing():
    bundle = json.loads((BASE / 'aligned5/matrix_manifest.json').read_text())
    for r in bundle['runs']:
        name = f"single_run{r['run_id']:03}"
        run(name, ['--features', OLD / 'batches/5' / r['name'] / 'features.tsv', '--config', BASE / 'uno.toml'])
        assert sha(OUT / name / 'proteins.tsv') == sha(BASE / 'E_batch' / f"run{r['run_id']:03}" / 'proteins.tsv')
    report()


def expanded():
    aligned = BASE / 'aligned5-replicated'
    aligned.mkdir(exist_ok=True)
    config = BASE / 'align5-replicated.toml'
    binary = OUT / 'koth_align'
    if not binary.exists():
        shutil.copy2(ROOT.parent / 'target/release/koth_align', binary)
        source = ROOT.parent / 'koth_rust'
        (OUT / 'koth-source-manifest.json').write_text(json.dumps({str(p.relative_to(source)): sha(p)
            for p in sorted((source / 'koth_ff/src').rglob('*.rs'))}, indent=2) + '\n')
    signature = dict(binary=sha(binary), config=sha(config), input_manifest=sha(BASE / 'aligned5/matrix_manifest.json'))
    stamp = aligned / 'job.json'
    if stamp.exists():
        prior = json.loads(stamp.read_text())
        assert prior['signature'] == signature and prior['returncode'] == 0
        assert all(sha(Path(p)) == h for p, h in prior['outputs'].items())
    else:
        started = time.time()
        argv = [str(binary), str(OLD / 'batches/5'), '--config', str(config), '--output', str(aligned)]
        print('Starting replicated-group alignment', flush=True)
        with (aligned / 'run.log').open('w') as log:
            result = subprocess.run(['/usr/bin/time', '-v', '-o', str(aligned / 'time.txt'), *argv],
                env=dict(os.environ, RAYON_NUM_THREADS='8'), stdout=log, stderr=subprocess.STDOUT)
        outputs = {str(p): sha(p) for p in aligned.glob('*.tsv')}
        if (aligned / 'matrix_manifest.json').exists():
            outputs[str(aligned / 'matrix_manifest.json')] = sha(aligned / 'matrix_manifest.json')
        stamp.write_text(json.dumps(dict(signature=signature, argv=argv, returncode=result.returncode,
            seconds=time.time() - started, outputs=outputs), indent=2) + '\n')
        assert result.returncode == 0
        print('Finished replicated-group alignment', round(time.time() - started, 1), flush=True)
    run('replicated_joint', ['--matrix-long', aligned, '--matrix-joint', '--config', BASE / 'uno.toml'])
    run('replicated_no_agreement', ['--matrix-long', aligned, '--matrix-joint', '--config', OUT / 'no-agreement.toml'])
    report()


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('action', choices=['pilot', 'report', 'timing', 'expanded', 'guarded'])
    p.add_argument('--version', choices=['v1', 'v2', 'v3'], default='v3')
    args = p.parse_args()
    OUT = BASE / ('joint_search_' + args.version)
    globals()[args.action]()
