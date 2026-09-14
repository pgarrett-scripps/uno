#!/usr/bin/env python3
"""Run fixed evidence/validation arms after the serialized resource cohort."""
import json
from pathlib import Path
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from matrix_study import ROOT, STUDY, cohort_rows, write_json


def main():
    names=['search_full_single_'+r['run']+'_seed0' for r in cohort_rows('20')]
    names+=['search_full_matrix_'+n+'_seed0' for n in ['5','10','20']]
    deadline=time.monotonic()+4*3600
    print('Waiting for all 23 serialized baseline searches before extensions.',flush=True)
    while True:
        complete=True
        for name in names:
            p=STUDY/'logs'/(name+'.json')
            if not p.exists():complete=False;continue
            try:j=json.loads(p.read_text())
            except json.JSONDecodeError:complete=False;continue
            if j['returncode']!=0:raise RuntimeError('Baseline failed: '+name)
        if complete:break
        if time.monotonic()>deadline:raise TimeoutError('Serialized cohort did not finish in four hours')
        time.sleep(15)
    write_json(STUDY/'serialized-cohort-complete.json',{'searches':names,'completed_at':time.time()})
    print('Serialized measurements complete; remaining experiments may overlap.',flush=True)
    def run(item):
        name,argv=item
        with open(STUDY/(name+'-driver.log'),'w') as log:
            result=subprocess.run([sys.executable,*argv],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
        if result.returncode:raise RuntimeError(name+' failed; inspect its driver log')
        return name
    jobs=[('full-evidence',['benchmark/matrix_study.py','full-evidence']),
          ('astral',['benchmark/matrix_astral.py','run'])]
    with ThreadPoolExecutor(max_workers=2) as pool:
        done=list(pool.map(run,jobs))
    write_json(STUDY/'extensions-complete.json',{'completed':done,'completed_at':time.time()})


if __name__=='__main__':main()
