#!/usr/bin/env python3
"""Independent verification of the published IM CNN study numbers and provenance."""
import csv,hashlib,json,re
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'benchmark/results/im_cnn_20260914'
OLD=ROOT/'experiments/im-timstof/uno/benchmark/results/im_timstof_20260913'
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
 report=json.loads((OUT/'report.json').read_text());counts=json.loads((OLD/'plan.json').read_text())['fasta_counts']
 archive=json.loads((OUT/'binary-archive.json').read_text())
 for p,h in report['plan']['source_hashes'].items():
  actual=archive['archived_path'] if p==archive['original_path'] else p
  assert sha(actual)==h,p
 assert sha(OUT/'im_public.uimc')==sha(ROOT/'models/im_cnn_public_20260914.uimc')
 splits={}
 for part in ['train','validation','test']:
  splits[part]={r['sequence'] for r in csv.DictReader((OUT/f'{part}.tsv').open(),delimiter='\t')}
 assert not splits['train']&splits['validation'] and not splits['train']&splits['test'] and not splits['validation']&splits['test']
 checked=0
 for run,seeds in report['results'].items():
  for seed,arms in seeds.items():
   for arm,result in arms.items():
    folder=OUT/'searches'/run/f'seed{seed}'/arm;receipt=json.loads((folder/'receipt.json').read_text())
    assert receipt['returncode']==0 and receipt['output_sha256']==sha(folder/'proteins.tsv')
    scores={}
    for row in csv.DictReader((folder/'proteins.tsv').open(),delimiter='\t'):
     if int(row['decoy']):continue
     names=row['group_members'].split(';')
     species=[{s for s in counts if '_'+s in a} for a in names]
     if any(len(s)!=1 for s in species) or len(set.union(*species))!=1:continue
     for name in names:scores[name]=max(scores.get(name,Fraction('-10000000')),Fraction(row['binomial_score']))
    ordered=sorted(scores,key=lambda a:(-scores[a],a))
    h=a=0;best={1:0,3:0,5:0}
    for i,name in enumerate(ordered):
     if '_ARATH' in name:a+=1
     else:h+=1
     if i+1<len(ordered) and scores[ordered[i+1]]==scores[name]:continue
     if not h:continue
     fdp=Fraction(a*sum(counts[s] for s in ['HUMAN','ECOLI','YEAST']),h*counts['ARATH'])
     for percent in best:
      if fdp<=Fraction(percent,100):best[percent]=max(best[percent],h)
    for percent,value in best.items():
     assert value==result['thresholds'][str(percent/100)]['present'],(run,seed,arm,percent,value)
     checked+=1
    log=(folder/'run.log').read_text()
    assert 'IM gate:' not in log
    if arm=='off':assert not result['im_stages']
    elif result['im_stages']:
     baseline=(OUT/'searches'/run/f'seed{seed}'/'off/run.log').read_text()
     initial=re.search(r'\(conf=(\d+) dim=(\d+) decoy=(\d+)\)',baseline)
     assert initial is not None,(run,seed,arm)
     assert sum(map(int,initial.groups()))==result['im_stages'][0][1],(run,seed,arm,'initial RT shortlist changed')
 result={'status':'passed','independently_verified_thresholds':checked,'sequence_disjoint_splits':{k:len(v) for k,v in splits.items()},'report_sha256':sha(OUT/'report.json'),'model_sha256':sha(OUT/'im_public.uimc'),'checker_sha256':sha(__file__)}
 (OUT/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
