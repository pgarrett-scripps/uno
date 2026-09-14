#!/usr/bin/env python3
"""Read-only matched evaluation of every completed IM arm and seed."""
import csv,json,hashlib,re
from collections import defaultdict
from decimal import Decimal
from pathlib import Path
from statistics import median
from im_timstof_study import OUT,EVAL,SEEDS,ARMS,ROOT,read,write,sha
SPECIES=['HUMAN','ECOLI','YEAST','ARATH']
def kind(acc):
 hits=[s for s in SPECIES if '_'+s in acc]
 return hits[0] if len(hits)==1 else 'AMBIGUOUS'
def curve(scores,counts):
 blocks=defaultdict(list)
 for acc,score in scores.items():blocks[score].append(acc)
 best={t:(0,0,None) for t in [1,3,5]};n=y=0
 for score in sorted(blocks,reverse=True):
  n+=sum(kind(a)!='ARATH' for a in blocks[score]);y+=sum(kind(a)=='ARATH' for a in blocks[score])
  for t in best:
   if n>best[t][0] and y*sum(counts[s] for s in SPECIES[:3])*100<=t*n*counts['ARATH']:best[t]=(n,y,score)
 out={}
 for t,(n,y,cut) in best.items():
  keep=sorted(a for a,s in scores.items() if cut is not None and s>=cut)
  assert sum(kind(a)=='ARATH' for a in keep)==y and len(keep)-y==n
  out[str(t/100)]=dict(present=n,entrapment=y,fdp=y*sum(counts[s] for s in SPECIES[:3])/counts['ARATH']/n if n else None,
    threshold=str(cut) if cut is not None else None,accessions=[a for a in keep if kind(a)!='ARATH'],entrapment_accessions=[a for a in keep if kind(a)=='ARATH'])
 return out

def main():
 assert read(OUT/'complete.json')==dict(searches=27)
 plan=read(OUT/'plan.json');counts=plan['fasta_counts'];support=read(OUT/'inputs/sage_support.json');result={};hashes={}
 for key in EVAL:
  result[key]={}
  for seed in SEEDS:
   result[key][str(seed)]={}
   for arm in ARMS:
    p=OUT/'searches'/key/f'seed{seed}'/arm;receipt=read(p/'receipt.json');assert receipt['returncode']==0
    for name,h in receipt['inputs'].items():
     if name not in hashes:hashes[name]=sha(Path(name))
     assert hashes[name]==h,name
    for name,h in receipt['outputs'].items():assert sha(p/name)==h,(p,name)
    scores={};qs={};amb=0
    for r in csv.DictReader((p/'proteins.tsv').open(),delimiter='\t'):
     if r['decoy']!='0':continue
     aa=r['group_members'].split(';');sp={kind(a) for a in aa}
     if len(sp)!=1 or 'AMBIGUOUS' in sp:amb+=1;continue
     score=Decimal(r['binomial_score']);q=Decimal(r['q_value']);assert score.is_finite() and q.is_finite()
     for a in aa:scores[a]=max(scores.get(a,Decimal('-Infinity')),score);qs[a]=min(qs.get(a,Decimal('Infinity')),q)
    values=curve(scores,counts)
    for v in values.values():
     acc=v['accessions'];v['sage_supported']=sum(a.split('|')[1] in support[key] for a in acc)
    accepted=[a for a,q in qs.items() if q<=Decimal('.01')];n=sum(kind(a)!='ARATH' for a in accepted);y=len(accepted)-n
    log=(p/'run.log').read_text();gate=re.search(r'IM gate:.*kept (\d+)/(\d+) PFMs',log)
    assert (gate is not None)==(arm!='off'),(key,seed,arm)
    result[key][str(seed)][arm]=dict(thresholds=values,ambiguous_groups=amb,
      internal_1pct=dict(present=n,entrapment=y,fdp=y*sum(counts[s] for s in SPECIES[:3])/counts['ARATH']/n if n else None),
      gate=dict(kept=int(gate[1]),before=int(gate[2]),rejected_fraction=1-int(gate[1])/int(gate[2])) if gate else None,
      seconds=receipt['seconds'])
 summaries={};deltas={arm:[] for arm in ARMS[1:]}
 for key,sd in result.items():
  summaries[key]={}
  for arm in ARMS:
   ns=[sd[str(s)][arm]['thresholds']['0.01']['present'] for s in SEEDS]
   summaries[key][arm]=dict(values=ns,median=median(ns),min=min(ns),max=max(ns))
   if arm!='off':deltas[arm].extend(ns[i]-sd[str(s)]['off']['thresholds']['0.01']['present'] for i,s in enumerate(SEEDS))
 outcomes={arm:dict(paired_deltas=v,improved=sum(x>0 for x in v),unchanged=sum(x==0 for x in v),worse=sum(x<0 for x in v),
   median_delta=median(v),all_files_median_improved=all(summaries[k][arm]['median']>summaries[k]['off']['median'] for k in EVAL)) for arm,v in deltas.items()}
 for arm,v in outcomes.items():v['passes_prespecified_gain_rule']=v['improved']>=7 and v['all_files_median_improved']
 report=dict(status='complete',plan=plan,calibration=read(OUT/'calibration.json'),public_model=read(OUT/'model/training-plan.json'),
  results=result,summary=summaries,outcomes=outcomes,verified_inputs=hashes,producer_sha256=sha(Path(__file__)))
 write(OUT/'report.json',report)
 print(json.dumps(dict(summary=summaries,outcomes=outcomes),indent=2))
if __name__=='__main__':main()
