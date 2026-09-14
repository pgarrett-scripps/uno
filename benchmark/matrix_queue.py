#!/usr/bin/env python3
"""Resume the fixed study in bounded batches after interrupted background drivers."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import sys
import time
import matrix_study as study
from matrix_finalize import expected_cases


def directory(dataset):return study.ROOT/'benchmark/results/matrix_20260912'/('astral' if dataset=='astral' else '')


def execute(dataset,name,seed):
    study.STUDY=directory(dataset);root=study.STUDY;seed=int(seed)
    single_prefix='full_single_' if dataset=='ionstar' else 'single_'
    if name.startswith(single_prefix):
        study.search(name,'--features',root/'koth'/name[len(single_prefix):]/'features.tsv',seed=seed,cache='0')
        return
    short=name.removeprefix('full_')
    floor=2;binary='uno';extra=[]
    if short.startswith('matrix_'):
        group=short.split('_')[1]
        if '_floor' in short:floor=int(short.rsplit('floor',1)[1])
        source=root/'aligned'/group/'intensity_matrix.tsv'
    elif short.startswith('profile_reassign_'):
        source=root/'aligned/20/intensity_matrix.tsv';extra=['--profile-reassign']
    elif short.startswith('run_evidence_null_'):
        source=root/'evidence/20/run_evidence_null.tsv';floor=1
    else:
        arm,group=short.rsplit('_',1)
        if arm not in ['quality','quality_null','masked']:raise ValueError('Unknown arm: '+name)
        source=root/'evidence'/group/(arm+'.tsv');binary='uno_evidence'
    study.search(name,'--matrix',source,seed=seed,floor=floor,binary=binary,extra=extra,cache='0')


def pending():
    jobs=[]
    for dataset,cases in expected_cases().items():
        for name,seed in cases:
            stamp=directory(dataset)/'logs'/f'search_{name}_seed{seed}.json'
            if stamp.exists() and json.loads(stamp.read_text()).get('returncode')==0:continue
            jobs.append((dataset,name,seed))
    priority=[('ionstar','full_profile_reassign_20',0),('astral','matrix_all',0),
              ('astral','matrix_200pg',0),('astral','quality_all',0),
              ('astral','quality_null_all',0),('astral','masked_all',0)]
    def order(job):
        if job in priority:return (-1,priority.index(job),'',0)
        dataset,name,seed=job
        # Finish a complete seed-zero comparison before expanding repeated searches.
        if seed==0:return (0,0,dataset,name)
        single=name.startswith(('full_single_','single_'))
        return (1+int(single),seed,dataset,name)
    return sorted(jobs,key=order)


def preserve_interrupted(job):
    dataset,name,seed=job;root=directory(dataset)
    stem=f'search_{name}_seed{seed}';stamp=root/'logs'/(stem+'.json')
    log=stamp.with_suffix('.log')
    if not log.exists() or stamp.exists():return
    archive=root/'interrupted_attempts'/stem/str(time.time_ns());archive.mkdir(parents=True)
    moved={}
    for p in [log,stamp.with_suffix('.time'),*(root/'searches'/name/f'seed{seed}').glob('*')]:
        if p.is_file():
            moved[str(p)]=study.sha(p);p.rename(archive/p.name)
    study.write_json(archive/'interruption.json',{'reason':'Prior background exec session ended with SIGTERM (143); no successful completion stamp',
                                                'preserved_files_sha256':moved})


def batch(count,workers):
    selected=pending()[:count]
    if not selected:print('All prespecified searches have completion stamps.',flush=True);return 0
    for job in selected:preserve_interrupted(job)
    print(f'Running {len(selected)} pending cases with {workers} workers; {len(pending())} cases remain.',flush=True)
    def run(job):
        dataset,name,seed=job
        print('Queued',dataset,name,seed,flush=True)
        result=subprocess.run([sys.executable,__file__,'--one',dataset,name,str(seed)],cwd=study.ROOT)
        print('Completed' if result.returncode==0 else 'Failed',dataset,name,seed,flush=True)
        return result.returncode
    with ThreadPoolExecutor(max_workers=workers) as pool:
        results=list(pool.map(run,selected))
    print('Pending after batch:',len(pending()),flush=True)
    return int(any(results))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--one',nargs=3,metavar=('DATASET','ARM','SEED'))
    p.add_argument('--batch-size',type=int,default=3);p.add_argument('--workers',type=int,default=3)
    args=p.parse_args()
    if args.one:execute(*args.one)
    else:raise SystemExit(batch(args.batch_size,args.workers))
