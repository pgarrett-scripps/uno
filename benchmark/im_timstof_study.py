#!/usr/bin/env python3
"""Fixed paired IM-gate test; model unchanged, calibration run excluded from searches."""
import argparse,csv,hashlib,json,math,os,re,shutil,subprocess,time
from collections import Counter,defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from statistics import median
ROOT=Path(__file__).resolve().parents[1]
HOME_REPO=Path('/home/patrick-garrett/Repos/uno')
OUT=ROOT/'benchmark/results/im_timstof_20260913'
STEMS={x:f'LFQ_Ultra2_PASEF_15min_50ng_Condition_{x[0]}_REP{x[1]}' for x in ['A1','A2','B2','C2']}
EVAL=['A2','B2','C2'];SEEDS=[0,1,2];ARMS=['off','im3','im2']
RAW=Path('/home/patrick-garrett/Data/dnoise/data/dda_15min/raw')
SAGE=Path('/home/patrick-garrett/Data/dnoise/results/dda_15min/original/results.sage.tsv')
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write(p,v):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def read(p):return json.loads(p.read_text())
def env():
 e={k:v for k,v in os.environ.items() if not k.startswith('UNO_')};e.update(RAYON_NUM_THREADS='8',OMP_NUM_THREADS='8',OPENBLAS_NUM_THREADS='1',RUST_LOG='info',UNO_RT_CACHE='1',UNO_RT_CACHE_DIR=str(OUT/'rt-cache'));return e

def prepare():
 OUT.mkdir(parents=True,exist_ok=True)
 assert not (OUT/'plan.json').exists(),'Already frozen; use the existing plan.'
 inputs=OUT/'inputs';inputs.mkdir(exist_ok=True)
 targets=HOME_REPO/'benchmark/results/matrix_20260913_koth050/ionstar/human_ecoli_yeast.fasta'
 plant=Path('/home/patrick-garrett/Data/fasta/arabidopsis_UP000006548_reviewed.fasta')
 fasta=inputs/'hye_arabidopsis.fasta';fasta.write_bytes(targets.read_bytes()+b'\n'+plant.read_bytes())
 counts=Counter();accessions=set()
 for l in fasta.open():
  if l.startswith('>'):
   a=l.split()[0][1:];assert a not in accessions;accessions.add(a);counts[a.split('_')[-1]]+=1
 assert set(counts)=={'HUMAN','YEAST','ECOLI','ARATH'}
 shutil.copy2(ROOT/'crates/uno/src/koth_bruker.toml',inputs/'koth.toml')
 (inputs/'uno.toml').write_text('[digest]\nmin_len = 7\nmax_len = 30\nmissed_cleavages = 0\n\n[matching]\nmin_charge = 2\nppm_tolerance = 20.0\n\n[scoring]\nmass_anchor_ppm = 20.0\n')
 sources={str(p):sha(p) for p in [targets,plant,SAGE]}
 for key in EVAL:
  for p in sorted((RAW/(STEMS[key]+'.d')).rglob('*')):
   if p.is_file():sources[str(p)]=sha(p)
 obs=defaultdict(list);support=defaultdict(set)
 with SAGE.open() as f:
  for r in csv.DictReader(f,delimiter='\t'):
   key=next((k for k,s in STEMS.items() if Path(r['filename']).stem==s),None)
   if key is None or r['label']!='1' or float(r['peptide_q'])>.01 or float(r['protein_q'])>.01:continue
   for acc in r['proteins'].split(';'):
    if acc.count('|')>=2:support[key].add(acc.split('|')[1])
   if key!='A1':continue
   mods=re.findall(r'\[([^]]+)\]',r['peptide'])
   if any(abs(float(m)-57.021464)>.001 for m in mods):continue
   seq=re.sub(r'\[[^]]+\]','',r['peptide']);z=int(r['charge']);im=float(r['ion_mobility'])
   if 7<=len(seq)<=30 and set(seq)<=set('ACDEFGHIKLMNPQRSTVWY') and z>=2 and im>0 and math.isfinite(im):obs[(seq,z)].append(im)
 with (inputs/'calibration.tsv').open('w') as f:
  wr=csv.writer(f,delimiter='\t');wr.writerow(['sequence','charge','im'])
  for (s,z),values in sorted(obs.items()):wr.writerow([s,z,median(values)])
 (inputs/'predict.tsv').write_text(''.join(f'{s}\t{z}\n' for s,z in sorted(obs)))
 write(inputs/'sage_support.json',{k:sorted(v) for k,v in support.items()})
 for p in inputs.iterdir():sources[str(p)]=sha(p)
 write(OUT/'plan.json',dict(calibration_run='A1',evaluation_runs=EVAL,seeds=SEEDS,arms=ARMS,
  primary='im3 versus off; im2 is prespecified sensitivity. All 27 searches retained.',
  model='Compatible 100-feature Sequoia model trained on public measured CCS transformed to inverse mobility; no evaluation-run labels used. Per-charge widths estimated independently on A1.',
  calibration='A1 used only for error-width estimation; each searched run gets the original gate\'s unsupervised per-charge median obs/pred scale from its candidate pool.',
  evaluation='Whole tied-binomial-score blocks, expanded same-species accessions; absent Arabidopsis count / present HYE count times FASTA HYE/Arabidopsis size, at 1%, 3%, 5%. Also internal q<=1% lists and per-file Sage corroboration.',
  decision='Call a consistent gain only if the primary arm improves median 1% depth in each of the three files and improves at least seven of nine paired searches; otherwise report mixed/no benefit. No parameter tuning after results.',
  limitations='One timsTOF Ultra 2 HYE dataset, three evaluation acquisitions and three seeds. MS/MS corroboration is not a complete truth set. Entrapment thresholds are retrospective, not deployable q-values. No paper update.',
  fasta_counts=dict(counts),source_hashes=sources,driver_sha256=sha(Path(__file__)),paper_tag='paper-2026-09-13'))
 print('Frozen inputs and independent calibration observations:',len(obs),dict(counts),flush=True)

def extract():
 koth=HOME_REPO/'benchmark/results/matrix_20260913_koth050/bin/koth_ff'
 def one(key):
  out=OUT/'features';out.mkdir(exist_ok=True);job=OUT/'extraction'/key;job.mkdir(parents=True,exist_ok=True)
  argv=[str(koth),str(RAW/(STEMS[key]+'.d')),'--config',str(OUT/'inputs/koth.toml'),'--output',str(out),'--threads','8']
  start=time.time()
  with (job/'run.log').open('w') as f:result=subprocess.run(argv,env=env(),stdout=f,stderr=subprocess.STDOUT)
  assert result.returncode==0,key
  p=out/STEMS[key]/'features.tsv'
  write(job/'receipt.json',dict(argv=argv,seconds=time.time()-start,binary_sha256=sha(koth),features_sha256=sha(p),features=str(p)))
  print('Extracted',key,flush=True)
 with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(one,EVAL))

def calibrate():
 trainer=OUT/'bin/uno-rt-train'
 argv=[str(trainer),'--labels',str(OUT/'inputs/calibration.tsv'),'--out',str(OUT/'model/im_public.sequoia'),'--predict-im',str(OUT/'inputs/predict.tsv')]
 with (OUT/'prediction.tsv').open('w') as f:subprocess.run(argv,stdout=f,check=True,env=env())
 predictions={}
 for r in csv.reader((OUT/'prediction.tsv').open(),delimiter='\t'):predictions[(r[0],int(r[1]))]=float(r[2])
 byz=defaultdict(list)
 for r in csv.DictReader((OUT/'inputs/calibration.tsv').open(),delimiter='\t'):
  z=int(r['charge']);observed=float(r['im']);pred=predictions[(r['sequence'],z)];assert pred>0 and math.isfinite(pred)
  byz[z].append((observed,pred))
 summaries={};residuals=defaultdict(list)
 for z,values in sorted(byz.items()):
  sc=median(o/p for o,p in values);errs=[abs((o-p*sc)/o) for o,p in values]
  residuals[2 if z==2 else 3].extend(errs)
  summaries[z]=dict(n=len(values),scale=sc,median_absolute_relative_error=median(errs))
 sigmas={}
 for z,errs in residuals.items():
  assert len(errs)>=100,z
  errs.sort();sigmas[z]=errs[int(.6827*len(errs))];assert 0<sigmas[z]<.3,(z,sigmas[z])
 result=dict(per_charge=summaries,relative_sigma=sigmas,model_sha256=sha(OUT/'model/im_public.sequoia'),prediction_sha256=sha(OUT/'prediction.tsv'))
 write(OUT/'calibration.json',result);print('Independent calibration',result,flush=True)

def run():
 cal=read(OUT/'calibration.json');sig=cal['relative_sigma'];exe=OUT/'bin/uno';fasta=OUT/'inputs/hye_arabidopsis.fasta';cfg=OUT/'inputs/uno.toml'
 rt=HOME_REPO/'benchmark/results/matrix_20260913_koth050/ionstar/rt.urtc'
 jobs=[(key,seed,arm) for seed in SEEDS for key in EVAL for arm in ARMS]
 def one(item):
  key,seed,arm=item;out=OUT/'searches'/key/f'seed{seed}'/arm;out.mkdir(parents=True,exist_ok=True)
  if (out/'receipt.json').exists():
   old=read(out/'receipt.json');assert old['returncode']==0 and sha(out/'proteins.tsv')==old['outputs']['proteins.tsv'];return
  feature=OUT/'features'/STEMS[key]/'features.tsv'
  argv=list(map(str,[exe,'--features',feature,'--fasta',fasta,'--rt-model',rt,'--config',cfg,'--seed',seed,'--out-proteins',out/'proteins.tsv']))
  e=env()
  if arm!='off':
   argv+=['--im-model',str(OUT/'model/im_public.sequoia')]
   e.update(UNO_IM_GATE=arm[-1],UNO_IM_SIG2=str(sig['2']),UNO_IM_SIG3=str(sig['3']))
  start=time.time();print('Starting',key,seed,arm,flush=True)
  with (out/'run.log').open('w') as f:res=subprocess.run(['/usr/bin/time','-v','-o',str(out/'time.txt'),*argv],env=e,stdout=f,stderr=subprocess.STDOUT)
  result=dict(argv=argv,env={k:v for k,v in e.items() if k.startswith('UNO_')},returncode=res.returncode,seconds=time.time()-start,
    inputs={str(p):sha(p) for p in [exe,fasta,cfg,rt,feature,OUT/'model/im_public.sequoia',OUT/'calibration.json']},
    outputs={n:sha(out/n) for n in ['proteins.tsv','run.log','time.txt'] if (out/n).exists()})
  write(out/'receipt.json',result);assert res.returncode==0,item
  if arm!='off':assert 'IM gate:' in (out/'run.log').read_text(),item
  print('Finished',key,seed,arm,round(result['seconds'],1),flush=True)
 with ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(one,jobs))
 write(OUT/'complete.json',dict(searches=27))

if __name__=='__main__':
 action=argparse.ArgumentParser();action.add_argument('action',choices=['prepare','extract','calibrate','run']);args=action.parse_args();globals()[args.action]()
