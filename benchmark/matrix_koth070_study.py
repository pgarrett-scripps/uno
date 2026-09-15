#!/usr/bin/env python3
"""Frozen, resumable 23-file study: separate features vs consensus rows.

No joint-member or LFQ-cell scoring arm. Results never overwrite older studies.
"""
import argparse
import collections
from concurrent.futures import ThreadPoolExecutor
import csv
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import threading
import time
import traceback

ROOT = Path(__file__).resolve().parents[1]
RELEASE = '88f9e7c54c9ad3d6ca5a56a906703a3d46982f43'


def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def write(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(obj, indent=2) + '\n')
    tmp.replace(path)


def files(path):
    return sorted(p for p in path.rglob('*') if p.is_file())


def hashes(paths):
    return {str(p): sha(p) for p in paths}


class Study:
    def __init__(self, out, threads=4, workers=2):
        self.out = out.resolve()
        self.threads, self.workers = threads, workers
        self.mutex = threading.Lock()
        self.active = set()

    def freeze(self):
        marker = self.out / 'frozen.json'
        if marker.exists():
            self.frozen = json.loads(marker.read_text())
            assert self.frozen['driver_sha256'] == sha(__file__), 'Driver changed since freeze'
            assert (self.threads, self.workers) == (self.frozen['threads'], self.frozen['workers'])
            self.plan = json.loads((self.out / 'plan.json').read_text())
            for p, digest in self.frozen['runtime_inputs'].items():
                assert sha(p) == digest, ('Frozen input changed', p)
            return
        plan_path = ROOT / 'benchmark/matrix_publication_koth070.json'
        self.plan = json.loads(plan_path.read_text())
        assert self.plan['arms'] == ['separate', 'consensus']
        assert self.plan['seeds'] == [0, 1, 2]
        assert self.plan['koth_commit'] == RELEASE
        rows = [r for c in self.plan['cohorts'] for r in c['inputs']]
        assert len(rows) == len({r['path'] for r in rows}) == 23
        assert collections.Counter(r['condition'] for r in self.plan['cohorts'][1]['inputs']) == dict.fromkeys('ABCDE', 4)
        assert shutil.disk_usage(self.out).free > 15 * 2**30, 'Need at least 15 GiB free'
        self.status('freezing')
        for row in rows:
            p = Path(row['path'])
            assert p.stat().st_size == int(row['bytes']) and sha(p) == row['sha256'], p
            print('Verified raw input', row['run'], flush=True)
        sources = {}
        for name in ['uno', 'koth_rust', 'sage', 'Sequoia']:
            source = ROOT.parent / name
            head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=source, text=True).strip()
            patch = subprocess.check_output(['git', 'diff', 'HEAD', '--binary'], cwd=source)
            if name == 'koth_rust':
                assert head == RELEASE and not patch, 'KOTH source differs from release'
            dest = self.out / 'source' / name
            tracked = subprocess.check_output(['git', 'ls-files', '-z'], cwd=source).decode().split('\0')
            for rel in filter(None, tracked):
                src = source / rel
                if src.is_file():
                    dst = dest / rel
                    dst.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(src, dst)
            (dest / 'working-tree.patch').write_bytes(patch)
            sources[name] = {'head': head, 'files': hashes(files(dest))}
        # Runtime resources may be ignored by Git. Freeze and load the native
        # library from this study, rather than the mutable build-tree rpath.
        for rel in ['models/im_cnn_public_20260914.uimc', 'models/rt_cnn_chrono_global_v3_ch32.urtc',
                    'crates/uno-lightgbm-sys/vendor/lib_lightgbm.so']:
            dst = self.out / 'source/uno' / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(ROOT / rel, dst)
            sources['uno']['files'][str(dst)] = sha(dst)
        bin_dir = self.out / 'bin'
        bin_dir.mkdir(exist_ok=True)
        for name in ['uno', 'koth_ff', 'koth_align', 'examples/long_validate']:
            shutil.copy2(ROOT / 'target/release' / name, bin_dir / Path(name).name)
        (self.out / 'lib').mkdir(exist_ok=True)
        shutil.copy2(ROOT / 'crates/uno-lightgbm-sys/vendor/lib_lightgbm.so', self.out / 'lib/lib_lightgbm.so')
        cfg = self.out / 'config'
        cfg.mkdir(exist_ok=True)
        for key, rel in self.plan['configurations'].items():
            shutil.copy2(ROOT / rel, cfg / (key + '.toml'))
        for cohort in self.plan['cohorts']:
            assert sha(ROOT / cohort['source_manifest']) == cohort['source_manifest_sha256']
            d = self.out / cohort['dataset']
            d.mkdir(exist_ok=True)
            shutil.copy2(ROOT / cohort['fasta'], d / 'database.fasta')
            shutil.copy2(ROOT / cohort['rt_model'], d / 'rt.urtc')
        write(self.out / 'plan.json', self.plan)
        shutil.copy2(__file__, self.out / 'driver.py')
        versions = {n: subprocess.check_output([str(bin_dir / n), '--version'], text=True).strip()
                    for n in ['uno', 'koth_ff', 'koth_align']}
        assert all('0.7.0' in versions[n] for n in ['koth_ff', 'koth_align'])
        resources = [*files(bin_dir), *files(cfg), *files(self.out / 'lib'), self.out / 'plan.json']
        resources += [self.out / c['dataset'] / n for c in self.plan['cohorts'] for n in ['database.fasta', 'rt.urtc']]
        self.frozen = dict(created_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                           sources=sources, runtime_inputs=hashes(resources), versions=versions,
                           driver_sha256=sha(__file__), threads=self.threads, workers=self.workers,
                           compiler=subprocess.check_output(['rustc', '--version'], text=True).strip(),
                           timing='Parallel jobs are throughput records; four warm-cache controls run serially at the end.',
                           caches='Study-owned RT/IM directories start empty and populate automatically.',
                           scope='75 publication searches plus 4 timing controls; 4-file IonStar pilot is separate.')
        write(marker, self.frozen)
        print('Frozen 23 inputs, 2 publication matrices, 75 searches.', flush=True)

    def environment(self):
        env = {k: v for k, v in os.environ.items() if not k.startswith('UNO_')}
        env.update(RAYON_NUM_THREADS=str(self.threads), OMP_NUM_THREADS=str(self.threads),
                   OPENBLAS_NUM_THREADS='1', RUST_LOG='info', LD_LIBRARY_PATH=str(self.out / 'lib'),
                   UNO_RT_CACHE='1', UNO_IM_CACHE='1', UNO_RT_CACHE_DIR=str(self.out / 'rt-cache'),
                   UNO_IM_CACHE_DIR=str(self.out / 'im-cache'))
        return env

    def status(self, state='running', error=None):
        with self.mutex:
            counts = collections.Counter()
            for p in self.out.glob('jobs/*/receipt.json'):
                receipt = json.loads(p.read_text())
                if receipt['returncode'] == 0:
                    counts[receipt['signature']['kind']] += 1
            write(self.out / 'status.json', dict(state=state, pid=os.getpid(), updated=time.time(),
                                                completed=dict(counts), active=sorted(self.active), error=error))

    def job(self, name, argv, inputs, outputs, kind):
        d = self.out / 'jobs' / name
        d.mkdir(parents=True, exist_ok=True)
        argv = list(map(str, argv))
        assert '--matrix-joint' not in argv
        env = self.environment()
        signature = dict(argv=argv, kind=kind, inputs=hashes(inputs), binary=sha(argv[0]),
                         library=sha(self.out / 'lib/lib_lightgbm.so'),
                         environment={k: env[k] for k in ['RAYON_NUM_THREADS', 'OMP_NUM_THREADS',
                             'OPENBLAS_NUM_THREADS', 'LD_LIBRARY_PATH', 'UNO_RT_CACHE', 'UNO_IM_CACHE',
                             'UNO_RT_CACHE_DIR', 'UNO_IM_CACHE_DIR']})
        stamp = d / 'receipt.json'
        if stamp.exists():
            old = json.loads(stamp.read_text())
            if old['returncode'] == 0:
                assert old['signature'] == signature and old['outputs'] == hashes(outputs), ('Stale job', name)
                return
        # Keep every failed/interrupted attempt, including logs and partial output hashes.
        attempt = d / ('attempt-' + str(time.time_ns()))
        attempt.mkdir()
        started = time.time()
        with self.mutex:
            self.active.add(name)
        self.status()
        print('Starting', name, flush=True)
        with (attempt / 'run.log').open('w') as log:
            result = subprocess.run(['/usr/bin/time', '-v', '-o', str(attempt / 'time.txt'), *argv],
                                    cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
        present = [p for p in outputs if p.is_file()]
        code = result.returncode if result.returncode else (0 if len(present) == len(outputs) else 99)
        record = dict(signature=signature, started=started, seconds=time.time()-started,
                      returncode=code, outputs=hashes(present), attempt=str(attempt))
        write(attempt / 'receipt.json', record)
        write(stamp, record)
        with self.mutex:
            self.active.remove(name)
        self.status()
        assert code == 0, ('Failed job', name, str(attempt / 'run.log'))
        print('Finished', name, round(record['seconds'], 1), 'seconds', flush=True)

    def parallel(self, fn, items):
        with ThreadPoolExecutor(max_workers=self.workers) as pool:
            list(pool.map(fn, items))

    def extract(self, item):
        cohort, row = item
        d = self.out / cohort['dataset'] / 'features'
        raw = Path(row['path'])
        assert sha(raw) == row['sha256'], raw
        cfg = self.out / 'config/detection.toml'
        self.job('extract_' + row['run'], [self.out/'bin/koth_ff', raw, '--config', cfg,
                 '--output', d, '--threads', self.threads], [raw, cfg],
                 [d/row['run']/n for n in ['features.tsv', 'hills.tsv', 'config.toml']], 'extraction')

    def align(self, cohort, rows, pilot=False):
        d = self.out / cohort['dataset']
        name = 'pilot' if pilot else 'publication'
        batch = d / 'batches' / name
        batch.mkdir(parents=True, exist_ok=True)
        inputs = [self.out/'config/alignment.toml']
        for row in rows:
            src = d/'features'/row['run']
            link = batch/row['run']
            if not link.exists():
                link.symlink_to(src, target_is_directory=True)
            assert link.resolve() == src
            inputs += files(src)
        assert len(list(batch.iterdir())) == len(rows)
        out = d/'aligned'/name
        outputs = [out/n for n in ['matrix_manifest.json', 'feature_observations.tsv',
            'lfq_matrix.long.tsv', 'intensity_matrix.tsv', 'qvalue_matrix.tsv',
            'decoy_intensity_matrix.tsv', 'align_config.toml', 'align_report.json',
            'consensus_features.tsv', 'lfq_details.tsv']]
        self.job(cohort['dataset']+'_'+name+'_align', [self.out/'bin/koth_align', batch,
                 '--config', inputs[0], '--output', out], inputs, outputs,
                 'pilot_alignment' if pilot else 'alignment')
        manifest = json.loads((out/'matrix_manifest.json').read_text())
        assert manifest['schema_version'] == 3 and manifest['n_consensus'] > 0
        source_tables = [d/'features'/r['name']/'features.tsv' for r in manifest['runs']]
        self.job(cohort['dataset']+'_'+name+'_validate', [self.out/'bin/long_validate', out,
                 out/'intensity_matrix.tsv', *source_tables], [*outputs, *source_tables], [], 'validation')
        return out

    def search(self, item):
        cohort, arm, source, seed, kind = item
        name = f"{cohort['dataset']}_{kind}_{arm}_{source.parent.name if arm == 'separate' else source.name}_seed{seed}"
        d = self.out/cohort['dataset']
        out = d/kind/name
        out.mkdir(parents=True, exist_ok=True)
        inputs = [self.out/'config/search.toml', d/'database.fasta', d/'rt.urtc']
        argv = [self.out/'bin/uno', '--config', inputs[0], '--fasta', inputs[1],
                '--rt-model', inputs[2], '--seed', seed, '--out-proteins', out/'proteins.tsv']
        outputs = [out/'proteins.tsv']
        if arm == 'separate':
            argv += ['--features', source]
            inputs += [source]
        else:
            assert arm == 'consensus'
            argv += ['--matrix-long', source, '--matrix-min-runs', '2',
                     '--out-lfq', out/'lfq.tsv', '--out-abundance', out/'abundance.tsv']
            inputs += files(source)
            outputs += [out/'lfq.tsv', out/'abundance.tsv']
        self.job(name, argv, inputs, outputs, kind)

    def run(self):
        self.freeze()
        astral, ionstar = self.plan['cohorts']
        pilot_rows = [r for r in ionstar['inputs'] if r['pilot'] == 'True']
        assert len(pilot_rows) == 4
        self.parallel(self.extract, [(ionstar, r) for r in pilot_rows])
        pilot = self.align(ionstar, pilot_rows, True)
        self.search((ionstar, 'consensus', pilot, 0, 'pilot_search'))
        for cohort in self.plan['cohorts']:
            self.parallel(self.extract, [(cohort, r) for r in cohort['inputs']])
            matrix = self.align(cohort, cohort['inputs'])
            tasks = []
            for seed in self.plan['seeds']:
                tasks += [(cohort, 'separate', self.out/cohort['dataset']/'features'/r['run']/'features.tsv',
                           seed, 'search') for r in cohort['inputs']]
                tasks += [(cohort, 'consensus', matrix, seed, 'search')]
            self.parallel(self.search, tasks)
        for row in astral['inputs']:
            self.search((astral, 'separate', self.out/'astral/features'/row['run']/'features.tsv', 0, 'timing'))
        self.search((astral, 'consensus', self.out/'astral/aligned/publication', 0, 'timing'))
        self.status()
        counts = json.loads((self.out/'status.json').read_text())['completed']
        for key, n in dict(extraction=23, alignment=2, search=75, timing=4, pilot_search=1).items():
            assert counts[key] == n, (key, counts)
        write(self.out/'complete.json', dict(completed=time.time(), counts=counts,
                                           note='Execution complete; scientific evaluation and manuscript updates remain.'))
        self.status('complete')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=['freeze', 'run'])
    p.add_argument('--out', required=True, type=Path)
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--workers', type=int, default=2)
    args = p.parse_args()
    assert args.threads > 0 and args.workers > 0
    args.out.mkdir(parents=True, exist_ok=True)
    study = Study(args.out, args.threads, args.workers)
    with (study.out/'driver.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        try:
            getattr(study, args.action)()
        except BaseException:
            study.status('failed', traceback.format_exc())
            raise


if __name__ == '__main__':
    main()
