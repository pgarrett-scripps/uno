#!/usr/bin/env python3
"""Check the experimental branch's IM-off output against the frozen main binary."""
import json,subprocess,time
from im_timstof_study import OUT,HOME_REPO,env,read,write,sha
reference=OUT/'searches/A2/seed0/off'
r=read(reference/'receipt.json')
exe=HOME_REPO/'experiments/consolidation-20260913/uno-paper'
out=OUT/'main_control';out.mkdir(exist_ok=True)
argv=r['argv'].copy();argv[0]=str(exe);argv[-1]=str(out/'proteins.tsv')
start=time.time()
with (out/'run.log').open('w') as f:result=subprocess.run(argv,env=env(),stdout=f,stderr=subprocess.STDOUT)
assert result.returncode==0
identical=sha(out/'proteins.tsv')==sha(reference/'proteins.tsv')
write(out/'receipt.json',dict(status='passed' if identical else 'failed',argv=argv,seconds=time.time()-start,binary_sha256=sha(exe),experimental_output_sha256=sha(reference/'proteins.tsv'),main_output_sha256=sha(out/'proteins.tsv')))
assert identical
print('Frozen main and experimental IM-off protein tables are byte-identical.')
