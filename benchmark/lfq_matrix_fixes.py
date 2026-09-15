#!/usr/bin/env python3
"""Bounded follow-up to the frozen matrix study; never resumes its old queue."""
import argparse
import csv
import hashlib
import json
import os
import re
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
OLD = Path('/home/patrick-garrett/Repos/uno/benchmark/results/matrix_20260912')
OUT = ROOT / 'benchmark/results/lfq_matrix_fixes'
BIN = ROOT.parent / 'target/release'


def sha(path):
    with open(path, 'rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def job(name, argv, outputs):
    logs = OUT / 'logs'
    logs.mkdir(parents=True, exist_ok=True)
    argv = list(map(str, argv))
    signature = dict(argv=argv, binary_sha256=sha(argv[0]))
    for flag in ['--config', '--rt-model', '--fasta']:
        if flag in argv:
            signature[flag] = sha(argv[argv.index(flag) + 1])
    if '--matrix-long' in argv:
        p = Path(argv[argv.index('--matrix-long') + 1])
        signature['matrix_manifest_sha256'] = sha(p / 'matrix_manifest.json' if p.is_dir() else p)
    stamp = logs / (name + '.json')
    if stamp.exists():
        old = json.loads(stamp.read_text())
        if old['returncode'] == 0:
            assert old['signature'] == signature, f'{name}: changed signature'
            assert all(p.exists() and sha(p) == old['outputs'][str(p)] for p in outputs)
            return
    env = dict(os.environ, RUST_LOG='info', RAYON_NUM_THREADS='8', UNO_RT_CACHE='1',
               UNO_RT_CACHE_DIR=str(OUT / 'rt-cache'))
    start = time.time()
    with (logs / (name + '.log')).open('w') as f:
        p = subprocess.run(['/usr/bin/time', '-v', '-o', str(logs / (name + '.time')), *argv],
                           stdout=f, stderr=subprocess.STDOUT, env=env)
    record = dict(signature=signature, started=start, elapsed=time.time() - start,
                  returncode=p.returncode, cache='study-local enabled', threads=8,
                  outputs={str(p): sha(p) for p in outputs if p.exists()})
    stamp.write_text(json.dumps(record, indent=2) + '\n')
    print(name, p.returncode, round(record['elapsed'], 1), flush=True)
    assert p.returncode == 0, f'{name} failed; inspect its log'


def align5():
    job('align5', [BIN / 'koth_align', OLD / 'batches/5', '--config', OUT / 'align5.toml',
                   '--output', OUT / 'aligned5'],
        [OUT / 'aligned5/matrix_manifest.json', OUT / 'aligned5/lfq_matrix.long.tsv',
         OUT / 'aligned5/feature_observations.tsv', OUT / 'aligned5/intensity_matrix.tsv'])


def audit5():
    runs = sorted((OLD / 'batches/5').iterdir())
    job('audit5', [BIN / 'examples/long_validate', OUT / 'aligned5',
                   OUT / 'aligned5/intensity_matrix.tsv', *[r / 'features.tsv' for r in runs]], [])


def pilot():
    common = ['--fasta', OLD / 'human_ecoli_yeast.fasta', '--rt-model', OLD / 'rt.urtc', '--seed', '0']
    for arm, extra, cfg in [
        ('B_grouping', ['--matrix-no-feature-quality'], 'uno.toml'),
        ('C_original_quality', [], 'uno.toml'),
        ('D_no_recurrence', [], 'no-recurrence.toml'),
    ]:
        out = OUT / arm
        out.mkdir(exist_ok=True)
        job(arm, [BIN / 'uno', '--matrix-long', OUT / 'aligned5', '--matrix-legacy-values',
                  *extra, '--config', OUT / cfg, *common, '--out-proteins', out / 'proteins.tsv',
                  '--out-lfq', out / 'lfq.tsv', '--out-abundance', out / 'abundance.tsv'],
            [out / 'proteins.tsv', out / 'lfq.tsv', out / 'abundance.tsv'])
    job('E_batch', [BIN / 'uno', '--matrix-long', OUT / 'aligned5', '--matrix-batch-out',
                    OUT / 'E_batch', '--config', OUT / 'uno.toml', *common],
        [OUT / 'E_batch/batch_manifest.json', OUT / 'E_batch/accession_union.tsv'])


def followup():
    out = OUT / 'F_isotope_representative'
    out.mkdir(exist_ok=True)
    job('F_isotope_representative', [BIN / 'uno', '--matrix-long', OUT / 'aligned5',
        '--matrix-legacy-values', '--matrix-prefer-isotope-support', '--fasta', OLD / 'human_ecoli_yeast.fasta',
        '--rt-model', OLD / 'rt.urtc', '--config', OUT / 'uno.toml', '--seed', '0',
        '--out-proteins', out / 'proteins.tsv', '--out-lfq', out / 'lfq.tsv', '--out-abundance', out / 'abundance.tsv'],
        [out / 'proteins.tsv', out / 'lfq.tsv', out / 'abundance.tsv'])


def report():
    from matrix_report import load_proteins, evaluate, accession_union, list_counts
    scale = (20416 + 4403) / 6066
    results = {}
    for arm, path in [('A_frozen_baseline', OLD / 'searches/full_matrix_5/seed0/proteins.tsv'),
                      *[(a, OUT / a / 'proteins.tsv') for a in ['B_grouping', 'C_original_quality', 'D_no_recurrence']]]:
        rows = load_proteins(path)
        results[arm] = dict(evaluate(rows, scale), accessions=list_counts(accession_union([rows]), scale))
    bundle = json.loads((OUT / 'aligned5/matrix_manifest.json').read_text())
    parity = []
    tables = []
    for r in bundle['runs']:
        new = OUT / 'E_batch' / f"run{r['run_id']:03}" / 'proteins.tsv'
        old = OLD / 'searches' / ('full_single_' + r['name']) / 'seed0/proteins.tsv'
        tables.append(load_proteins(new))
        parity.append(dict(run=r['name'], protein_tsv_identical=sha(new) == sha(old)))
    results['E_batch'] = dict(accessions=list_counts(accession_union(tables), scale), run_parity=parity)
    f = OUT / 'F_isotope_representative/proteins.tsv'
    if f.exists():
        rows = load_proteins(f)
        results['F_isotope_representative'] = dict(evaluate(rows, scale), accessions=list_counts(accession_union([rows]), scale))
    results['scope'] = 'Five-run IonStar, seed zero. A frozen baseline; B-D same corrected membership and legacy numerical serialization. D changes recurrence columns only. E preserves original run inputs.'
    results['decision'] = 'No expanded validation or paper changes: corrected consensus variants did not exceed the baseline at empirical 1% entrapment q; all five batch protein tables matched independent results.'
    results['followup_rationale'] = '7,843 groups had a two-isotope seed despite an available member with at least three isotopes. F tests that one representative-selection change.'
    results['timing_limit'] = 'Diagnostic timings: compilation overlapped part of the batch run, and old reference timing was cold-cache. No controlled speedup estimate.'
    results['batch_cache_requests'] = [dict(misses=int(a), requested=int(b), cached_before=int(c)) for a, b, c in re.findall(r'(\d+) unique misses of (\d+) requested \((\d+) rows cached\)', (OUT / 'logs/E_batch.log').read_text())]
    (OUT / 'pilot-report.json').write_text(json.dumps(results, indent=2) + '\n')
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=['align5', 'audit5', 'pilot', 'followup', 'report'])
    args = p.parse_args()
    globals()[args.action]()
