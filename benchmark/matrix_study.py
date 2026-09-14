#!/usr/bin/env python3
"""Fresh, restartable KOTH/UNO matrix study. Never reads historical results."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor

ROOT = Path(__file__).resolve().parents[1]
STUDY = ROOT / 'benchmark/results/matrix_20260912'
DATA = Path('/home/patrick-garrett/Data/PXD003881')


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(8 << 20), b''):
            h.update(b)
    return h.hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def fasta_records(path):
    header, seq = None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            if header is not None:
                yield header, ''.join(seq)
            header, seq = line[1:], []
        elif line.strip():
            if header is None:
                raise ValueError(f'Sequence before header: {path}')
            seq.append(line.strip())
    if header is not None:
        yield header, ''.join(seq)


def freeze():
    marker = STUDY / 'frozen.json'
    if marker.exists():
        print('Already frozen:', marker, flush=True)
        return
    STUDY.mkdir(parents=True, exist_ok=True)
    repos = {}
    for name, dirs in [('uno', ['crates', 'models']), ('koth_rust', ['koth_ff']),
                       ('sage', ['crates']), ('Sequoia', ['crates'])]:
        src = ROOT.parent / name
        dst = STUDY / 'source' / name
        if not dst.exists():
            dst.mkdir(parents=True)
            for file in ['Cargo.toml', 'Cargo.lock', 'build.rs']:
                if (src / file).exists():
                    shutil.copy2(src / file, dst / file)
            for folder in dirs:
                shutil.copytree(src / folder, dst / folder,
                                ignore=shutil.ignore_patterns('target', '__pycache__', '.venv', 'build'))
            patch = subprocess.check_output(['git', '-C', str(src), 'diff', 'HEAD', '--binary'])
            (dst / 'working-tree.patch').write_bytes(patch)
        if not (dst / 'working-tree.patch').exists():
            raise RuntimeError(f'Incomplete source snapshot: {dst}')
        commit = subprocess.check_output(['git', '-C', str(src), 'rev-parse', 'HEAD'], text=True).strip()
        repos[name] = {'commit': commit, 'patch_sha256': sha(dst / 'working-tree.patch'),
                       'files': {str(p.relative_to(dst)): sha(p)
                                 for p in sorted(dst.rglob('*')) if p.is_file()}}
    cfg = STUDY / 'config'
    cfg.mkdir(exist_ok=True)
    koth = (STUDY / 'source/uno/crates/uno/src/koth_default.toml').read_text()
    (cfg / 'koth.toml').write_text(koth)
    # A single explicit config is shared by external and embedded KOTH.
    uno = '[digest]\nmin_len = 7\nmax_len = 30\nmissed_cleavages = 0\n\n'
    nested = re.sub(r'^\[([^\]]+)\]', r'[koth.\1]', koth, flags=re.M)
    (cfg / 'uno.toml').write_text(uno + nested)
    # Keep singletons upstream so the prespecified UNO floors 1/2/3 are real.
    # All other alignment settings are the frozen current KOTH defaults.
    (cfg / 'align.toml').write_text(
        '[lfq]\nrun_tdc = true\nnormalize = "none"\n'
        'quant_estimator = "sum"\ndetected_use_grid = true\n'
        'averagine_projection = false\n\n[lfq.consensus]\nmin_group_size = 1\n\n'
        '[output]\nformat = "tsv"\nexport_details = true\nexport_decoys = true\n')
    inputs = []
    seen = {}
    for p in sorted((DATA / 'mzml').glob('*.mzML.gz')):
        condition = re.search(r'human_ecoli_([A-E])_', p.name).group(1)
        seen[condition] = seen.get(condition, 0) + 1
        inputs.append({'run': p.stem, 'path': str(p), 'condition': condition,
                       'replicate_order': seen[condition],
                       'expected_fold': dict(zip('ABCDE', [1, 1.5, 2, 2.5, 3]))[condition],
                       'pilot': condition in 'AE' and seen[condition] <= 2,
                       'bytes': p.stat().st_size, 'sha256': sha(p)})
        print('Hashed', p.name, flush=True)
    assert seen == dict.fromkeys('ABCDE', 4), seen
    with open(STUDY / 'inputs.tsv', 'w') as f:
        w = csv.DictWriter(f, fieldnames=list(inputs[0]), delimiter='\t')
        w.writeheader(); w.writerows(inputs)
    # Both sample species are targets; yeast is the absent entrapment species.
    fasta = STUDY / 'human_ecoli_yeast.fasta'
    entries = {}
    for p in [DATA / 'fasta/human_ecoli.fasta', ROOT / 'benchmark/data/entrap_human_yeast.fasta']:
        for header, sequence in fasta_records(p):
            name = header.split()[0]
            if p.name == 'entrap_human_yeast.fasta' and '_YEAST' not in name:
                continue
            if 'DECOY' in name:
                continue
            if name in entries and entries[name][1] != sequence:
                raise ValueError(f'Conflicting FASTA accession: {name}')
            entries[name] = (header, sequence)
    fasta.write_text(''.join(f'>{h}\n{s}\n' for h, s in entries.values()))
    species = {sp: sum('_' + sp in name for name in entries) for sp in ['HUMAN', 'ECOLI', 'YEAST']}
    assert sum(species.values()) == len(entries), species
    shutil.copy2(ROOT / 'models/rt_cnn_chrono_global_v3_ch32.urtc', STUDY / 'rt.urtc')
    write_json(marker, {'created_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        'repos': repos, 'fasta_counts': species, 'fasta_sha256': sha(fasta),
        'model_sha256': sha(STUDY / 'rt.urtc'), 'seeds': [0, 1, 2], 'threads': 8,
        'config_sha256': {p.name: sha(p) for p in cfg.iterdir()},
        'inputs_sha256': sha(STUDY / 'inputs.tsv'),
        'protocol': str(ROOT / 'docs/koth-matrix-benchmark-plan.md'),
        'note': 'Fresh current defaults; no historical matrix, fitted values or old results reused.'})
    print('Frozen:', marker, species, flush=True)


def run_job(name, argv, outputs=(), env_extra=None):
    logs = STUDY / 'logs'
    logs.mkdir(exist_ok=True)
    stamp = logs / (name + '.json')
    argv = list(map(str, argv))
    outputs = list(map(Path, outputs))
    signature = {'argv': argv, 'executable_sha256': sha(Path(argv[0])),
                 'config_sha256': {p.name: sha(p) for p in (STUDY / 'config').iterdir()}}
    if stamp.exists():
        old = json.loads(stamp.read_text())
        if old.get('returncode') == 0 and old.get('signature') == signature and all(
                p.exists() and sha(p) == old['outputs'].get(str(p)) for p in outputs):
            print('Cached', name, flush=True)
            return
        if old.get('returncode') == 0:
            raise RuntimeError(f'Stale job {name}; inspect before replacing outputs')
        suffix = '.failed-' + str(time.time_ns())
        for p in [stamp, stamp.with_suffix('.time'), stamp.with_suffix('.log')]:
            if p.exists(): p.rename(p.with_name(p.name + suffix))
    env = dict(os.environ, RAYON_NUM_THREADS='8', OMP_NUM_THREADS='8', OPENBLAS_NUM_THREADS='1', RUST_LOG='info')
    env.update(env_extra or {})
    started = time.time()
    print('Starting', name, flush=True)
    with open(logs / (name + '.log'), 'w') as log:
        result = subprocess.run(['/usr/bin/time', '-v', '-o', str(logs / (name + '.time')),
                                 *argv], stdout=log, stderr=subprocess.STDOUT, env=env)
    record = {'signature': signature, 'started': started, 'elapsed': time.time()-started,
              'returncode': result.returncode,
              'outputs': {str(p): sha(p) for p in outputs if p.exists()}}
    write_json(stamp, record)
    if result.returncode or any(not p.exists() for p in outputs):
        raise RuntimeError(f'{name} failed; inspect {logs / (name + ".log")}')
    print('Finished', name, round(record['elapsed'], 1), 's', flush=True)


def extract(pilot=True):
    rows = list(csv.DictReader(open(STUDY / 'inputs.tsv'), delimiter='\t'))
    for row in rows:
        if pilot and row['pilot'] != 'True':
            continue
        run = row['run']
        out = STUDY / 'koth' / run
        run_job('extract_' + run, [STUDY / 'bin/koth_ff', row['path'], '--config',
                STUDY / 'config/koth.toml', '--output', STUDY / 'koth', '--threads', '8'],
                [out / 'features.tsv', out / 'hills.tsv', out / 'config.toml'])


def cohort_rows(cohort):
    rows = list(csv.DictReader(open(STUDY / 'inputs.tsv'), delimiter='\t'))
    if cohort == 'pilot': return [r for r in rows if r['pilot'] == 'True']
    n = int(cohort)
    assert n in (5, 10, 20)
    return [r for r in rows if int(r['replicate_order']) <= n // 5]


def align(cohort):
    batch = STUDY / 'batches' / cohort
    batch.mkdir(parents=True, exist_ok=True)
    rows = cohort_rows(cohort)
    for row in rows:
        src = STUDY / 'koth' / row['run']
        assert (src / 'features.tsv').exists() and (src / 'hills.tsv').exists(), src
        dst = batch / row['run']
        # KOTH discovers real subdirectories (DirEntry metadata does not follow
        # directory symlinks). File symlinks preserve the exact extracted input.
        if dst.is_symlink(): dst.unlink()
        dst.mkdir(exist_ok=True)
        for name in ['features.tsv', 'hills.tsv', 'config.toml']:
            link = dst / name
            if not link.exists(): link.symlink_to(src / name)
    assert set(p.name for p in batch.iterdir()) == {r['run'] for r in rows}
    out = STUDY / 'aligned' / cohort
    run_job('align_' + cohort, [STUDY / 'bin/koth_align', batch, '--config',
            STUDY / 'config/align.toml', '--output', out],
            [out / n for n in ['intensity_matrix.tsv', 'qvalue_matrix.tsv',
                              'lfq_details.tsv', 'consensus_features.tsv']])


def search(name, input_flag, source, seed=0, floor=2, extra=(), cache=None, binary='uno'):
    out = STUDY / 'searches' / name / f'seed{seed}'
    out.mkdir(parents=True, exist_ok=True)
    argv = [STUDY / 'bin' / binary, input_flag, source,
            '--fasta', STUDY / 'human_ecoli_yeast.fasta', '--rt-model', STUDY / 'rt.urtc',
            '--config', STUDY / 'config/uno.toml', '--seed', str(seed),
            '--out-proteins', out / 'proteins.tsv', '--out-lfq', out / 'lfq.tsv',
            '--out-abundance', out / 'abundance.tsv']
    if input_flag == '--matrix': argv += ['--matrix-min-runs', str(floor)]
    run_job('search_' + name + f'_seed{seed}', argv + list(extra),
            [out / n for n in ['proteins.tsv', 'lfq.tsv', 'abundance.tsv']],
            env_extra={'UNO_RT_CACHE': cache} if cache is not None else None)


def pilot_search():
    rows = cohort_rows('pilot')
    for row in rows:
        search('single_' + row['run'], '--features', STUDY / 'koth' / row['run'] / 'features.tsv')
    # Exactly the first pilot file in the frozen manifest: no outcome-based selection.
    row = rows[0]
    search('raw_' + row['run'], '--raw', row['path'])
    for floor in [2, 1, 3]:
        search(f'matrix_pilot_floor{floor}', '--matrix',
               STUDY / 'aligned/pilot/intensity_matrix.tsv', floor=floor)


def full_baseline():
    # Fresh uncached predictions for every timed search. The pilot inherited a
    # pre-existing cache and is explicitly excluded from publication timing.
    extract(pilot=False)
    for cohort in ['5', '10', '20']: align(cohort)
    def single_job(item):
        row, seed = item
        search('full_single_' + row['run'], '--features',
               STUDY/'koth'/row['run']/'features.tsv', seed=seed, cache='0')
    for row in cohort_rows('20'): single_job((row,0))
    for cohort in ['5','10','20']:
        search('full_matrix_' + cohort, '--matrix',
               STUDY/'aligned'/cohort/'intensity_matrix.tsv', seed=0, cache='0')
    # Only seed 0 contributes serialized timing; extra seeds assess ID variation.
    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(single_job, [(r,s) for s in [1,2] for r in cohort_rows('20')]))
    for seed in [1,2]:
        for cohort in ['5','10','20']:
            search('full_matrix_' + cohort, '--matrix',
                   STUDY/'aligned'/cohort/'intensity_matrix.tsv', seed=seed, cache='0')
    for floor in [1, 3]:
        for seed in [0, 1, 2]:
            search(f'full_matrix_20_floor{floor}', '--matrix',
                   STUDY/'aligned/20/intensity_matrix.tsv', floor=floor, seed=seed, cache='0')


def pilot_evidence():
    for arm in ['masked', 'quality_null', 'quality']:
        search('pilot_' + arm, '--matrix', STUDY/'evidence/pilot'/(arm+'.tsv'),
               cache='0', binary='uno_evidence')
    # One unchanged matrix on the new binary verifies the absent-column route.
    search('pilot_baseline_new_binary', '--matrix', STUDY/'aligned/pilot/intensity_matrix.tsv',
           cache='0', binary='uno_evidence')


def full_evidence():
    from matrix_evidence import prepare
    prepare('20')
    for seed in [0,1,2]:
        for arm in ['masked','quality_null','quality']:
            search('full_' + arm + '_20', '--matrix', STUDY/'evidence/20'/(arm+'.tsv'),
                   seed=seed, cache='0', binary='uno_evidence')
        search('full_run_evidence_null_20','--matrix',STUDY/'evidence/20/run_evidence_null.tsv',
               seed=seed,floor=1,cache='0')
        search('full_profile_reassign_20','--matrix',STUDY/'aligned/20/intensity_matrix.tsv',
               seed=seed,cache='0',extra=['--profile-reassign'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['freeze', 'pilot-extract', 'extract', 'align', 'pilot-search', 'full-baseline', 'pilot-evidence', 'full-evidence'])
    parser.add_argument('--cohort', default='pilot', choices=['pilot', '5', '10', '20'])
    args = parser.parse_args()
    if args.action == 'freeze': freeze()
    elif args.action == 'align': align(args.cohort)
    elif args.action == 'pilot-search': pilot_search()
    elif args.action == 'full-baseline': full_baseline()
    elif args.action == 'pilot-evidence': pilot_evidence()
    elif args.action == 'full-evidence': full_evidence()
    else: extract(pilot=args.action == 'pilot-extract')
