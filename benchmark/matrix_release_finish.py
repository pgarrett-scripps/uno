#!/usr/bin/env python3
"""Finish this active study after both queues drain, then audit its evidence."""
import json,os,subprocess,time
from pathlib import Path
from matrix_release_study import OUT,ROOT
if __name__=='__main__':
    while True:
        astral=list((OUT/'astral/searches').glob('*/seed*/job.json'))
        if len(astral)==66 and all(json.loads(p.read_text())['returncode']==0 for p in astral) and (OUT/'ionstar-queue-complete.json').exists():break
        time.sleep(10)
    print('Both search queues complete; verifying receipts and running serialized controls.',flush=True)
    subprocess.run(['python','benchmark/matrix_release_study.py','run'],cwd=ROOT,check=True)
    env=dict(os.environ,UV_CACHE_DIR='/tmp/uno-paper-uv-cache',UV_NO_SYNC='1')
    subprocess.run(['uv','run','--project','paper/analysis','--quiet','python','benchmark/matrix_release_report.py'],cwd=ROOT,env=env,check=True)
