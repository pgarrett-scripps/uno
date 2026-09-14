#!/usr/bin/env python3
"""Prepare a fixed, sequence-disjoint public CCS training/validation split.

CCS * sqrt(reduced molecular mass) / charge is proportional to 1/K0 by
Mason-Schamp. Divide by 1000 to keep labels near unity; the search calibrates
this arbitrary global constant independently by charge. No absolute CCS or
mobility accuracy is inferred from that constant.
"""
import csv,hashlib,io,json,math,zipfile
from collections import defaultdict
from pathlib import Path
from statistics import median
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'benchmark/results/im_timstof_20260913'
(OUT/'model').mkdir(parents=True,exist_ok=True)
SOURCE=Path('/home/patrick-garrett/Repos/pepboost/benchmarks/data/ccs_features/combined_sm.zip')
AA=dict(zip('ACDEFGHIKLMNPQRSTVWY',[71.037113805,103.009184505,115.026943065,129.042593135,147.068413945,57.021463735,137.058911875,113.084063975,128.094963015,113.084063975,131.040484645,114.042927470,97.052763875,128.058577540,156.101111050,87.032028435,101.047678505,99.068413945,186.079312980,163.063328575]))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(SOURCE)=='097f84e3ffc4edd9a88d5dc0fc858069a90a72c050a59bc6bfbc74258ac7dd91'
values=defaultdict(list);counts=defaultdict(int)
with zipfile.ZipFile(SOURCE) as z:
 with z.open('combined_sm.csv') as stream:
  for row in csv.DictReader(io.TextIOWrapper(stream)):
   counts['raw']+=1;s=row['Modified sequence'].strip('_')
   if not 7<=len(s)<=30 or not set(s)<=set(AA):continue
   charge=float(row['Charge']);mass=float(row['Mass']);ccs=float(row['CCS'])
   if charge not in (2,3,4) or not (mass>0 and ccs>0 and math.isfinite(ccs)):continue
   expected=18.010564684+sum(AA[a] for a in s)+57.021463735*s.count('C')
   if abs(mass-expected)>.02:counts['mass_mismatch']+=1;continue
   mu=mass*28.0134/(mass+28.0134)
   values[(s,int(charge))].append(ccs*math.sqrt(mu)/charge/1000)
train=[];val=[]
for (s,z),v in sorted(values.items()):
 split=int.from_bytes(hashlib.sha256(('uno-im-public-v1:'+s).encode()).digest()[:8],'big')%10
 (val if split==0 else train).append((s,z,median(v)))
assert not {r[0] for r in train}&{r[0] for r in val}
for name,rows in [('train',train),('validation',val)]:
 with (OUT/'model'/f'{name}.csv').open('w') as f:
  w=csv.writer(f);w.writerow(['modified_seq','charge','predicted IM'])
  for s,z,v in rows:w.writerow([s+'/'+str(z),z,v])
(OUT/'model/validation_predict.tsv').write_text(''.join(f'{s}\t{z}\n' for s,z,v in val))
result=dict(source=str(SOURCE),source_sha256=sha(SOURCE),counts=dict(counts),training_rows=len(train),validation_rows=len(val),training_sequences=len({r[0] for r in train}),validation_sequences=len({r[0] for r in val}),
 split='SHA256 uno-im-public-v1:sequence modulo 10; bucket 0 external validation; group by sequence, all charge states together.',
 label='CCS * sqrt(mass * 28.0134 / (mass + 28.0134)) / charge / 1000; proportional inverse mobility; dataset-specific multiplicative calibration is required.',
 chemistry='Unmodified sequences with fixed carbamidomethyl C; measured mass within 0.02 Da of that chemistry; lengths 7-30, charges 2-4; replicate medians.',
 parameters=dict(backend='sequoia',max_depth=5,iterations=300,shrinkage=.1,seed=42),
 hashes={str(p.name):sha(p) for p in (OUT/'model').glob('*') if p.is_file()})
(OUT/'model/training-plan.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
