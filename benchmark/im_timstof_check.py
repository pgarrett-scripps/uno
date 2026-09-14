#!/usr/bin/env python3
"""Independently validate floating-point ranked curves and the frozen input archive."""
import csv,json,hashlib,itertools
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'benchmark/results/im_timstof_20260913'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
r=json.loads((OUT/'report.json').read_text());counts=r['plan']['fasta_counts'];scale=(counts['HUMAN']+counts['ECOLI']+counts['YEAST'])/counts['ARATH'];checked=0
for key,seeds in r['results'].items():
 for seed,arms in seeds.items():
  for arm,expected in arms.items():
   path=OUT/'searches'/key/f'seed{seed}'/arm/'proteins.tsv';scores={}
   for row in csv.DictReader(path.open(),delimiter='\t'):
    if row['decoy']!='0':continue
    aa=row['group_members'].split(';');classes=[{s for s in counts if '_'+s in a} for a in aa]
    if any(len(c)!=1 for c in classes) or len(set.union(*classes))!=1:continue
    score=float(row['binomial_score'])
    for a in aa:scores[a]=max(scores.get(a,float('-inf')),score)
   n=y=0;best={t:(0,0,None) for t in [.01,.03,.05]}
   for score,items in itertools.groupby(sorted(scores.items(),key=lambda x:-x[1]),key=lambda x:x[1]):
    for a,_ in items:
     if '_ARATH' in a:y+=1
     else:n+=1
    for t in best:
     if n>best[t][0] and y*scale/n<=t:best[t]=(n,y,score)
   for t,(n,y,cut) in best.items():
    e=expected['thresholds'][str(t)];acc=sorted(a for a,s in scores.items() if cut is not None and s>=cut and '_ARATH' not in a)
    assert (n,y,acc)==(e['present'],e['entrapment'],e['accessions']),(key,seed,arm,t)
    checked+=1
assert checked==81
# Hash the snapshot inputs used by each job, not a silently replaced current file.
for name,h in r['verified_inputs'].items():assert sha(Path(name))==h,name
source=json.loads((OUT/'source.json').read_text())
for name,h in source['files'].items():assert sha(ROOT/name)==h,name
assert json.loads((OUT/'model/batch_parity.json').read_text())['status']=='passed'
result=dict(status='passed',independent_cutoffs=checked,report_sha256=sha(OUT/'report.json'),source_sha256=sha(OUT/'source.json'),checker_sha256=sha(Path(__file__)))
(OUT/'validation.json').write_text(json.dumps(result,indent=2)+'\n');print(result)
