import csv,hashlib,json
from pathlib import Path
src=Path('experiments/im-timstof/uno/benchmark/results/im_timstof_20260913/model');out=Path('benchmark/results/im_cnn_20260914')
counts={}
for old,new,limit in [('train','train',100000),('validation','validation',10000)]:
 rows=list(csv.DictReader((src/f'{old}.csv').open()))
 seqs={r['modified_seq'].split('/')[0] for r in rows}
 selected=set(sorted(seqs,key=lambda s:hashlib.sha256(('uno-im-cnn-v1:'+s).encode()).digest())[:limit])
 with (out/f'{new}.tsv').open('w') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['sequence','charge','im'])
  n=0
  for r in rows:
   s=r['modified_seq'].split('/')[0]
   if s in selected:w.writerow([s,r['charge'],r['predicted IM']]);n+=1
 counts[new]={'sequences':len(selected),'observations':n}
 if old=='validation':
  with (out/'test.tsv').open('w') as f:
   w=csv.writer(f,delimiter='\t');w.writerow(['sequence','charge','im'])
   for r in rows:
    s=r['modified_seq'].split('/')[0]
    if s not in selected:w.writerow([s,r['charge'],r['predicted IM']])
(out/'data-plan.json').write_text(json.dumps({'selection':'Deterministic SHA256 uno-im-cnn-v1:sequence; first 100000 public training sequences and 10000 public validation sequences; remaining original validation sequences reserved for independent test. All charges grouped.', 'counts':counts,'source':json.loads((src/'training-plan.json').read_text())},indent=2))
(out/'train-config.json').write_text(json.dumps({'threads':4,'batch_size':512,'patience':5,'lr':.001}))
print(counts)
