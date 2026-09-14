#!/usr/bin/env python3
"""Execute distant queued repeats early, after the serialized timing cohort.

The original drivers will checksum-reuse these completed outputs. The primary
matrix tail is scheduled first, well before its driver finishes the remaining
individual searches; the Astral matrix tail follows while validation is still
in its individual-search phases. Parameters and case inventory do not change.
"""
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
from matrix_astral import ASTUDY
import matrix_study as study


def phase(tasks):
    def run(task):
        name,source,seed,floor,binary=task
        log=study.STUDY/'logs'/f'search_{name}_seed{seed}.log'
        stamp=log.with_suffix('.json')
        if log.exists() and not stamp.exists():
            raise RuntimeError('Another driver has started '+name+'; do not duplicate it')
        study.search(name,'--matrix',source,seed=seed,floor=floor,binary=binary,cache='0')
    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(run,tasks))


def main():
    root=study.STUDY
    if not (root/'serialized-cohort-complete.json').exists():raise RuntimeError('Timing phase is not complete')
    primary=[('full_matrix_'+g,root/'aligned'/g/'intensity_matrix.tsv',s,2,'uno')
             for s in [1,2] for g in ['5','10','20']]
    primary += [(f'full_matrix_20_floor{f}',root/'aligned/20/intensity_matrix.tsv',s,f,'uno')
                for f in [1,3] for s in [0,1,2]]
    phase(primary)
    study.STUDY=ASTUDY
    validation=[]
    for seed in [1,2]:
        for group in ['200pg','1ng','10ng','100ng','all']:
            validation.append(('matrix_'+group,ASTUDY/'aligned'/group/'intensity_matrix.tsv',seed,2,'uno'))
            if group=='all':
                validation += [(arm+'_all',ASTUDY/'evidence/all'/(arm+'.tsv'),seed,2,'uno_evidence')
                               for arm in ['quality_null','quality','masked']]
    phase(validation)
    study.write_json(root/'tail-complete.json',{'primary_cases':len(primary),'astral_cases':len(validation),
        'purpose':'Prespecified repeats run ahead of original sequential tails; no resource claims from these overlapping searches'})


if __name__=='__main__':main()
