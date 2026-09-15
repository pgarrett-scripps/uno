#!/usr/bin/env python3
"""Reproducible RT-first IM-soft comparison; existing frozen experiment is read-only."""
import argparse,csv,hashlib,json,os,re,subprocess,time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal
from pathlib import Path
from statistics import median
ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT/'experiments/im-timstof/uno/benchmark/results/im_timstof_20260913'
OUT=ROOT/'benchmark/results/im_cnn_20260914'
ARMS=['off','tree_soft','cnn_soft'];RUNS=['A2','B2','C2'];SEEDS=[0,1,2]
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write(p,x):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def read(p):return json.loads(p.read_text())
def prepare():
 assert not (OUT/'plan.json').exists()
 paths=[ROOT/'target/release/uno',OUT/'im_public.uimc',OLD/'model/im_public.sequoia',OLD/'inputs/hye_arabidopsis.fasta',OLD/'inputs/uno.toml',ROOT/'benchmark/results/matrix_20260913_koth050/ionstar/rt.urtc']
 for key in RUNS:
  receipt=read(OLD/'searches'/key/'seed0/off/receipt.json'); paths.append(Path(receipt['argv'][2]))
 write(OUT/'plan.json',dict(arms=ARMS,runs=RUNS,seeds=SEEDS,source_hashes={str(p):sha(p) for p in paths},
  design='Paired current-source runs. RT and mass finish shortlisting before IM prediction. No IM cutoff. Only a calibrated relative IM residual and availability flag enter existing semi-supervised scoring. Same targets/decoys, RT model, features, config and seed in each triple. No tuning on search results.',
  calibration='Unsupervised per-charge median observed/predicted on all RT-shortlisted target and decoy candidates equally; minimum 20 candidates. Independent A1 predictions are diagnostic only, not used for search tuning.',
  model='Same emb32/ch32 residual dilated RT CNN with an added charge/4 global input; trained on a fixed public subset with all charges of each sequence in one split. UIMC artifacts use separate magic and schema from RT.',
  decision='Primary cnn_soft vs off: improved median 1% empirical entrapment FDP depth in all three files and at least seven of nine paired searches. tree_soft is a composition-model control. Report all arms and seeds regardless of outcome.',
  limitations='Retrospective empirical accession FDP thresholds; one HYE timsTOF dataset. Public pretrained sequences can occur in search acquisitions. Current source differs from frozen previous IM study; only within-study pairs support comparisons.',
  data=read(OUT/'data-plan.json'),driver_sha256=sha(__file__)))
 subprocess.run(['git','diff','--','crates'],cwd=ROOT,stdout=(OUT/'source.patch').open('w'),check=True)
 print('Frozen study plan',flush=True)
def run():
 plan=read(OUT/'plan.json')
 archive=read(OUT/'binary-archive.json') if (OUT/'binary-archive.json').exists() else None
 for p,h in plan['source_hashes'].items():
  actual=archive['archived_path'] if archive and p==archive['original_path'] else p
  assert sha(actual)==h,p
 exe=Path(archive['archived_path']) if archive else ROOT/'target/release/uno'
 jobs=[(key,seed,arm) for seed in SEEDS for key in RUNS for arm in ARMS]
 def one(job):
  key,seed,arm=job;folder=OUT/'searches'/key/f'seed{seed}'/arm;folder.mkdir(parents=True,exist_ok=True)
  if (folder/'receipt.json').exists():
   r=read(folder/'receipt.json');assert r['returncode']==0 and sha(folder/'proteins.tsv')==r['output_sha256'];return
  old=read(OLD/'searches'/key/f'seed{seed}'/'off/receipt.json');argv=old['argv'];argv[0]=str(exe);argv[-1]=str(folder/'proteins.tsv')
  if arm!='off':argv+=['--im-model',str(OUT/'im_public.uimc' if arm=='cnn_soft' else OLD/'model/im_public.sequoia')]
  env={k:v for k,v in os.environ.items() if not k.startswith('UNO_')}
  env.update(RAYON_NUM_THREADS='4',OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='1',RUST_LOG='info',UNO_RT_CACHE='1',UNO_RT_CACHE_DIR=str(OUT/'rt-cache'))
  start=time.time();print('Starting',*job,flush=True)
  with (folder/'run.log').open('w') as f:r=subprocess.run(['/usr/bin/time','-v','-o',str(folder/'time.txt'),*argv],env=env,stdout=f,stderr=subprocess.STDOUT)
  receipt=dict(argv=argv,environment={k:v for k,v in env.items() if k.startswith(('UNO_','RAYON','OMP','OPENBLAS'))},seconds=time.time()-start,returncode=r.returncode,output_sha256=sha(folder/'proteins.tsv') if (folder/'proteins.tsv').exists() else None)
  write(folder/'receipt.json',receipt);assert r.returncode==0,job
  print('Finished',*job,round(receipt['seconds'],1),flush=True)
 with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(one,jobs))
 write(OUT/'complete.json',dict(searches=len(jobs)))
def report():
 assert read(OUT/'complete.json')==dict(searches=27)
 counts=read(OLD/'plan.json')['fasta_counts'];ratio=sum(counts[s] for s in ['HUMAN','ECOLI','YEAST'])/counts['ARATH']
 result={};summaries={};deltas={a:[] for a in ARMS[1:]}
 for key in RUNS:
  result[key]={}
  for seed in SEEDS:
   result[key][seed]={}
   for arm in ARMS:
    folder=OUT/'searches'/key/f'seed{seed}'/arm;receipt=read(folder/'receipt.json')
    assert receipt['returncode']==0 and sha(folder/'proteins.tsv')==receipt['output_sha256']
    blocks=defaultdict(set);accepted=set();seen={}
    for r in csv.DictReader((folder/'proteins.tsv').open(),delimiter='\t'):
     if r['decoy']!='0':continue
     accs=r['group_members'].split(';');species={s for a in accs for s in counts if '_'+s in a}
     if len(species)!=1:continue
     for a in accs:
      seen[a]=max(seen.get(a,Decimal('-Infinity')),Decimal(r['binomial_score']))
      if Decimal(r['q_value'])<=Decimal('.01'):accepted.add(a)
    for a,score in seen.items():blocks[score].add(a)
    thresholds={str(t):dict(present=0,entrapment=0,fdp=None) for t in [.01,.03,.05]};n=y=0
    for score,accs in sorted(blocks.items(),reverse=True):
     y+=sum('_ARATH' in a for a in accs);n+=sum('_ARATH' not in a for a in accs)
     for t in [.01,.03,.05]:
      # Integer arithmetic preserves whole-score blocks at the exact boundary.
      if n>thresholds[str(t)]['present'] and y*sum(counts[s] for s in ['HUMAN','ECOLI','YEAST'])*100<=round(t*100)*n*counts['ARATH']:
       thresholds[str(t)]=dict(present=n,entrapment=y,fdp=y*ratio/n)
    log=(folder/'run.log').read_text();im=[list(map(int,m)) for m in re.findall(r'IM soft feature: (\d+) unique sequence/charge predictions for (\d+) RT-shortlisted candidates; (\d+) calibrated',log)]
    result[key][seed][arm]=dict(thresholds=thresholds,seconds=receipt['seconds'],im_stages=im,internal_1pct_present=sum('_ARATH' not in a for a in accepted),internal_1pct_entrapment=sum('_ARATH' in a for a in accepted))
  summaries[key]={a:dict(values=[result[key][s][a]['thresholds']['0.01']['present'] for s in SEEDS],median=median(result[key][s][a]['thresholds']['0.01']['present'] for s in SEEDS)) for a in ARMS}
  for a in ARMS[1:]:deltas[a].extend(result[key][s][a]['thresholds']['0.01']['present']-result[key][s]['off']['thresholds']['0.01']['present'] for s in SEEDS)
 outcomes={a:dict(deltas=d,improved=sum(x>0 for x in d),tied=sum(x==0 for x in d),worse=sum(x<0 for x in d),consistent_gain=sum(x>0 for x in d)>=7 and all(summaries[k][a]['median']>summaries[k]['off']['median'] for k in RUNS)) for a,d in deltas.items()}
 write(OUT/'report.json',dict(summary=summaries,outcomes=outcomes,results=result,plan=read(OUT/'plan.json')))
 print(json.dumps(dict(summary=summaries,outcomes=outcomes),indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run','report']);globals()[p.parse_args().action]()
