#!/usr/bin/env python3
"""Evaluate the frozen single-arm confidence-weighted bootstrap plan, never launching searches.

Stage lowinput-screen evaluates 30 pairs; screening adds primary/raw controls;
confirmation requires a passing screening report for the identical candidate.
Pooled FDP is descriptive and undefined values never become zero.
"""
from pathlib import Path
from fractions import Fraction
import argparse,csv,hashlib,itertools,json,math,statistics
import stability_study as helper

ROOT=Path(__file__).resolve().parent.parent
PLAN=ROOT/'benchmark/results/soft_bootstrap/plan-manifest.json'
REGISTRY=ROOT/'benchmark/results/soft_bootstrap/candidate-registry.json'
ADDITIONAL=ROOT/'benchmark/results/soft_bootstrap/additional-inputs.json'
V3_SHA='50b4227bca68a58283368d7189dd65dbbd54854f340d05d2ca355824c4ad1738'
DATASETS={
 'lowinput':{'astral_1ng_rep1','astral_1ng_rep2','astral_1ng_rep4','astral_200pg_rep1','astral_200pg_rep2','astral_200pg_rep3'},
 'primary':{'hela5min_09746','hela5min_09766','hela5min_09734','hela5min_09778','hela5min_09950','faims_LUM2_05433_VG'},
 'raw':{'astral_200pg_rep1','astral_200pg_rep2','astral_200pg_rep3'},
 'additional':{'astral_10ng_rep1','astral_10ng_rep4','astral_10ng_rep6'}}
SEEDS={'screening':{'lowinput':[0,1,2,5,8],'primary':[0,1,2],'raw':[0,1,2]},'confirmation':{'lowinput':[9,10],'primary':[11,12],'raw':[5,6],'additional':[3,4]}}
METRICS=('depth_1pct','depth_5pct','rows','internal_q1_targets','internal_q1_human','internal_q1_yeast','internal_q1_other')

def sha(path):return helper.digest(path)
def read(path):return json.loads(Path(path).read_text())
def tsv(path):
 with Path(path).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def median(values):return statistics.median(values) if values else None
def verify(path,expected,hashes):
 path=Path(path).resolve();actual=sha(path)
 if actual!=expected:raise ValueError(f'hash mismatch: {path}')
 hashes[str(path)]=actual

def verify_candidate_registry(path,arm,candidate_sha,hashes):
 registry=read(path);hashes[str(Path(path).resolve())]=sha(path)
 entry=registry['arms'][arm]
 if entry['binary_sha256']!=candidate_sha:raise ValueError('candidate does not match frozen arm registry')
 verify(entry['manifest_path'],entry['manifest_sha256'],hashes)
 verify(entry['binary_path'],entry['binary_sha256'],hashes)
 return entry

def verify_additional_input(job,fingerprints,registered):
 entry=registered[job['dataset']];actual=Path(job['command'][2]).resolve()
 if actual!=Path(entry['features']).resolve() or fingerprints.get(job['command'][2])!=entry['sha256']:
  raise ValueError('additional acquisition does not match frozen dataset mapping')

def verify_reuse(job,stamp):
 reused=stamp.get('reused_from') or stamp.get('reused_baseline')
 if reused!=job.get('reused_baseline'):raise ValueError('copied baseline metadata differs from declared reuse')
 if reused:
  for kind in ('output','stdout','stderr'):
   if stamp.get(kind+'_sha256')!=reused['source_hashes'][reused['files'][kind]]:
    raise ValueError('copied artifact differs from verified source')
 return reused

def exact_pooled_fdp_delta(baseline,candidate,scale):
 if not baseline['human'] or not candidate['human']:return None
 # The production FASTA scale is a ratio of integer entry counts. Preserve that
 # ratio exactly rather than turning its rounded float into a new rational.
 factor=Fraction(20647,6066) if scale==20647/6066 else Fraction(str(scale))
 return factor*(Fraction(candidate['yeast'],candidate['human'])-Fraction(baseline['yeast'],baseline['human']))

def recipe(command,fingerprints):
 if len(command)!=13 or command[1] not in ('--features','--raw') or command[3::2]!=['--fasta','--rt-model','--config','--seed','--out-proteins']:
  raise ValueError('unexpected scientific search recipe')
 return dict(mode=command[1],seed=int(command[10]),inputs=tuple(fingerprints[command[i]] for i in (2,4,6,8)))

def load_panel(directory,panel,seeds,candidate_sha,hashes,additional_inputs=None):
 directory=directory.resolve()
 if (directory/'incomplete.json').exists():raise ValueError(f'incomplete panel: {directory}')
 study=read(directory/'study.json');cert=read(directory/'summary_provenance.json')
 verify(directory/'study.json',cert['study_sha256'],hashes);hashes[str(directory/'summary_provenance.json')]=sha(directory/'summary_provenance.json')
 for name,h in cert['summaries'].items():verify(directory/name,h,hashes)
 expected={(ds,s,v) for ds in DATASETS[panel] for s in seeds for v in ('baseline','candidate')}
 jobs={(j['dataset'],int(j['seed']),j['variant']):j for j in study['searches']}
 summary=tsv(directory/'runs.tsv');keys=[(r['dataset'],int(r['seed']),r['variant']) for r in summary]
 if set(jobs)!=expected or len(jobs)!=len(study['searches']) or set(keys)!=expected or len(keys)!=len(expected):raise ValueError(f'wrong, missing, or duplicate panel coverage: {panel}')
 wanted_env={'RAYON_NUM_THREADS':'4','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','RUST_LOG':'info','UNO_RT_CACHE_DIR':'/tmp/uno-stability/rt-cache'}
 if study['environment']!=wanted_env or study['cwd']!=str(ROOT):raise ValueError('unexpected search environment/cwd')
 if study['entrapment_scale']!=20647/6066:raise ValueError('unexpected FASTA entrapment scale')
 rows=[];tables={};scale=study['entrapment_scale'];helper.entrapment.SCALE=scale
 for r in summary:
  key=r['dataset'],int(r['seed']),r['variant'];job=jobs[key];command=job['command'];binary_hash=V3_SHA if key[2]=='baseline' else candidate_sha
  if study['fingerprints'].get(command[0])!=binary_hash:raise ValueError(f'wrong engine fingerprint: {key}')
  for i in (0,2,4,6,8):
   p=command[i]
   if p not in hashes:verify(p,study['fingerprints'][p],hashes)
   elif hashes[p]!=study['fingerprints'][p]:raise ValueError('inconsistent scientific input hash')
  other=jobs[key[0],key[1],'candidate' if key[2]=='baseline' else 'baseline']
  if recipe(command,study['fingerprints'])!=recipe(other['command'],study['fingerprints']):raise ValueError('unmatched scientific inputs/seed')
  if recipe(command,study['fingerprints'])['seed']!=key[1]:raise ValueError('command and record seed disagree')
  if panel=='additional':
   if additional_inputs is None:raise ValueError('missing frozen additional acquisition mapping')
   verify_additional_input(job,study['fingerprints'],additional_inputs)
  path=Path(job['output']).resolve()
  if not path.is_relative_to(directory):raise ValueError('output outside panel directory')
  verify(path,cert['outputs'][str(path)],hashes);stamp_path=path.with_suffix('.json');stamp=read(stamp_path);hashes[str(stamp_path)]=sha(stamp_path)
  if stamp.get('study_sha256')!=cert['study_sha256'] or stamp.get('command')!=command or stamp.get('returncode')!=0 or stamp.get('output_sha256')!=cert['outputs'][str(path)]:raise ValueError('run provenance mismatch')
  for stream in ('stdout','stderr'):verify(path.with_suffix('.'+stream+'.log'),stamp[stream+'_sha256'],hashes)
  reused=verify_reuse(job,stamp)
  if key[2]=='candidate' and (reused or job.get('reused_baseline')):raise ValueError('candidate may not be copied from another study')
  if reused:
   if key[2]!='baseline' or stamp.get('executed_command','missing') is not None or stamp.get('elapsed_seconds','missing') is not None or stamp.get('started_unix','missing') is not None or stamp.get('timing_measured') is not False:raise ValueError('copied baseline has fabricated execution/timing')
   for name,h in reused['source_hashes'].items():verify(name,h,hashes)
  if command[1]!=('--raw' if panel=='raw' else '--features'):raise ValueError('wrong input mode for panel')
  got,accepted=helper.evaluate(path)
  for metric in METRICS:
   if got[metric]!=int(r[metric]):raise ValueError(f'summary metric mismatch: {key} {metric}')
  fdp=got['internal_q1_entrapment_fdp'];fdp=None if fdp=='' else fdp
  recorded=r.get('internal_q1_entrapment_fdp',r.get('internal_q1_fdp',''))
  if fdp!=(None if recorded=='' else float(recorded)):raise ValueError('reported FDP mismatch')
  table=tsv(path)
  reported={r['accession'] for r in table if r['decoy'] in ('0','False','false') and float(r['q_value'])<=.01}
  rows.append(dict(dataset=key[0],seed=key[1],variant=key[2],**{k:got[k] for k in METRICS},fdp=fdp,entrapment1_human_set=accepted,reported_q1_target_set=reported,path=str(path)))
  signature=lambda data:hashlib.sha256(json.dumps(sorted(data),separators=(',',':')).encode()).hexdigest()
  tables[key]=dict(scores=signature(tuple(sorted((k,v) for k,v in row.items() if k not in ('n_group_members','group_members'))) for row in table),all_columns=signature(tuple(sorted(row.items())) for row in table),sha=cert['outputs'][str(path)])
 return rows,scale,tables,study.get('inherited_runtime_environment',{})

def aggregate(rows,scale):
 h=sum(r['internal_q1_human'] for r in rows);y=sum(r['internal_q1_yeast'] for r in rows);defined=[r for r in rows if r['fdp'] is not None]
 return dict(n=len(rows),empty_q1=sum(not r['internal_q1_targets'] for r in rows),zero_depth1=sum(not r['depth_1pct'] for r in rows),empty_outputs=sum(not r['rows'] for r in rows),output_rows=sum(r['rows'] for r in rows),human=h,yeast=y,other=sum(r['internal_q1_other'] for r in rows),targets=sum(r['internal_q1_targets'] for r in rows),pooled_fdp=y*scale/h if h else None,individual_fdp_defined=len(defined),individual_fdp_above1=sum(r['fdp']>.01 for r in defined),maximum_individual_fdp=max((r['fdp'] for r in defined),default=None))

def reported_variability(rows):
 result={}
 for metric in ('internal_q1_targets','internal_q1_human'):
  files=[]
  for ds in sorted({r['dataset'] for r in rows}):
   item=dict(dataset=ds)
   for variant in ('baseline','candidate'):
    values=[r[metric] for r in rows if r['dataset']==ds and r['variant']==variant];med=median(values);mean=statistics.mean(values)
    item[variant]=dict(n=len(values),minimum=min(values),maximum=max(values),median=med,mean=mean,relative_range=(max(values)-min(values))/med if len(values)>=2 and med else None,sample_cv=statistics.stdev(values)/mean if len(values)>=2 and mean else None)
   files.append(item)
  common=[r for r in files if all(r[v]['relative_range'] is not None for v in ('baseline','candidate'))]
  result[metric]=dict(files=files,common_range_files=[r['dataset'] for r in common],median_relative_range={v:median([r[v]['relative_range'] for r in common]) for v in ('baseline','candidate')},excluded_range_files=[r['dataset'] for r in files if r not in common],interpretation='Target-count relative range supplies the low-input improvement criterion; human-count median supplies a depth guard. Relative range uses the same positive-median files in both variants. Sample CV is descriptive, uses sample SD (ddof=1)/mean, and is null with fewer than two seeds or zero mean.')
 return result

def summarize_panel(rows,scale,panel,stage,tables=None):
 subsets={'all':rows}
 if stage!='confirmation' and panel=='lowinput':subsets.update(core=[r for r in rows if r['seed'] in (0,1,2)],stress=[r for r in rows if r['seed'] in (5,8)])
 stats={};ranges={};medians={};overlaps=[]
 for subset,rr in subsets.items():
  stats[subset]={v:aggregate([r for r in rr if r['variant']==v],scale) for v in ('baseline','candidate')}
  ranges_by_variant={};excluded={};medians[subset]=[]
  for v in ('baseline','candidate'):
   rv={};ex={}
   for ds in sorted({r['dataset'] for r in rr}):
    values=[r['depth_1pct'] for r in rr if r['dataset']==ds and r['variant']==v];med=median(values)
    if len(values)>=2 and med:rv[ds]=(max(values)-min(values))/med
    else:ex[ds]='zero median' if not med else 'fewer than two seeds'
   ranges_by_variant[v]=rv;excluded[v]=ex
  common=sorted(ranges_by_variant['baseline'].keys()&ranges_by_variant['candidate'].keys())
  ranges[subset]=dict(common_files=common,excluded=excluded,**{v:median([ranges_by_variant[v][ds] for ds in common]) for v in ('baseline','candidate')})
  for ds in sorted({r['dataset'] for r in rr}):
   for metric in ('depth_1pct','depth_5pct'):
    a,b=[median([r[metric] for r in rr if r['dataset']==ds and r['variant']==v]) for v in ('baseline','candidate')]
    medians[subset].append(dict(dataset=ds,metric=metric,baseline=a,candidate=b,relative_delta=(b-a)/a if a else None,loss_over5pct=b<.95*a))
  for v in ('baseline','candidate'):
   for metric in ('depth_1pct','depth_5pct'):stats[subset][v]['median_file_'+metric]=median([r[v] for r in medians[subset] if r['metric']==metric])
 pairs=[];rescued=[];new_empty=[];losses=[];large_losses=[]
 lookup={(r['dataset'],r['seed'],r['variant']):r for r in rows}
 for ds,seed in sorted({(r['dataset'],r['seed']) for r in rows}):
  a,b=(lookup[ds,seed,v] for v in ('baseline','candidate'))
  p=dict(dataset=ds,seed=seed,baseline_q1=a['internal_q1_targets'],candidate_q1=b['internal_q1_targets'])
  for label,key in [('entrapment1_human','entrapment1_human_set'),('reported_q1_target','reported_q1_target_set')]:
   union=a[key]|b[key];p['paired_'+label+'_jaccard']=len(a[key]&b[key])/len(union) if union else None
  if tables:
   ta,tb=(tables[ds,seed,v] for v in ('baseline','candidate'));p.update(score_rows_equal=ta['scores']==tb['scores'],group_only_change=ta['scores']==tb['scores'] and ta['all_columns']!=tb['all_columns'],byte_equal=ta['sha']==tb['sha'])
  for metric in METRICS:p[metric]=dict(baseline=a[metric],candidate=b[metric],delta=b[metric]-a[metric])
  pairs.append(p)
  if a['internal_q1_targets'] and not b['internal_q1_targets']:new_empty.append(dict(dataset=ds,seed=seed))
  if not a['internal_q1_targets'] and b['internal_q1_targets']:rescued.append(b)
  for metric in ('depth_1pct','depth_5pct'):
   if b[metric]<a[metric]:
    loss=dict(dataset=ds,seed=seed,metric=metric,baseline=a[metric],candidate=b[metric],relative_loss=1-b[metric]/a[metric]);losses.append(loss)
    if b[metric]<.9*a[metric]:large_losses.append(loss)
 for v in ('baseline','candidate'):
  for ds in sorted({r['dataset'] for r in rows}):
   rr=sorted([r for r in rows if r['dataset']==ds and r['variant']==v],key=lambda r:r['seed'])
   for a,b in itertools.combinations(rr,2):
    overlap=dict(dataset=ds,variant=v,seeds=[a['seed'],b['seed']])
    for label,key in [('entrapment1_human','entrapment1_human_set'),('reported_q1_target','reported_q1_target_set')]:
     union=a[key]|b[key];overlap[label+'_jaccard']=len(a[key]&b[key])/len(union) if union else None
    overlaps.append(overlap)
 rescued_stats={subset:aggregate([r for r in rescued if subset=='all' or r['seed'] in ({0,1,2} if subset=='core' else {5,8})],scale) for subset in subsets}
 a,b=stats['all']['baseline'],stats['all']['candidate'];fdp_delta=None if a['pooled_fdp'] is None or b['pooled_fdp'] is None else b['pooled_fdp']-a['pooled_fdp']
 exact_delta=exact_pooled_fdp_delta(a,b,scale)
 checks=[dict(criterion=f'{panel}_no_new_empty_q1',passed=not new_empty,detail=new_empty),dict(criterion=f'{panel}_no_file_median_loss_over5pct',passed=not any(r['loss_over5pct'] for r in medians['all']),detail=[r for r in medians['all'] if r['loss_over5pct']]),dict(criterion=f'{panel}_pooled_fdp_delta_at_most0.5pp',passed=exact_delta<=Fraction(1,200) if exact_delta is not None else None,detail=dict(baseline=a['pooled_fdp'],candidate=b['pooled_fdp'],delta=fdp_delta,exact_delta=str(exact_delta) if exact_delta is not None else None,maximum_delta=.005))]
 reported = reported_variability(rows)
 human_losses=[dict(dataset=f['dataset'],baseline=f['baseline']['median'],candidate=f['candidate']['median']) for f in reported['internal_q1_human']['files'] if f['candidate']['median'] < .95*f['baseline']['median']]
 checks.append(dict(criterion=f'{panel}_no_file_median_reported_human_loss_over5pct',passed=not human_losses,detail=human_losses))
 if panel=='lowinput':
  span=reported['internal_q1_targets']['median_relative_range'];fewer=b['empty_q1']<a['empty_q1']
  less=span['baseline'] is not None and span['candidate'] is not None and span['baseline']>0 and span['candidate']<=.9*span['baseline']
  checks.append(dict(criterion='lowinput_improves_failures_or_reported_count_variability',passed=fewer or less,detail=dict(fewer_actual_empty_q1=fewer,reported_target_relative_range_reduced_at_least10pct=less,ranges=span)))
 for ds,seed in sorted({(r['dataset'],r['seed']) for r in rows}):
  aa,bb=(lookup[ds,seed,v] for v in ('baseline','candidate'))
  metric='internal_q1_human'
  if bb[metric]<aa[metric]:
   loss=dict(dataset=ds,seed=seed,metric=metric,baseline=aa[metric],candidate=bb[metric],relative_loss=1-bb[metric]/aa[metric]);losses.append(loss)
   if bb[metric]<.9*aa[metric]:large_losses.append(loss)

 return dict(descriptive_reported_q1_variability=reported_variability(rows),aggregates=stats,relative_ranges=ranges,file_medians=medians,pairs=pairs,newly_empty=new_empty,all_losses=losses,losses_over10pct=large_losses,rescued_only=rescued_stats,within_variant_overlap=overlaps,individual_results=[{k:v for k,v in r.items() if k not in ('entrapment1_human_set','reported_q1_target_set')} for r in rows],checks=checks)

def rejection_label(checks):
 # A definite guard failure takes precedence over a different undefined metric.
 return 'FAIL: KEEP V3' if any(c['passed'] is False for c in checks) else 'INCONCLUSIVE: KEEP V3'

def main(argv=None):
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--stage',choices=('lowinput-screen','screening','confirmation'),required=True);parser.add_argument('--arm',choices=('soft',),required=True)
 parser.add_argument('--candidate-sha256',required=True);parser.add_argument('--lowinput',type=Path,required=True);parser.add_argument('--primary',type=Path);parser.add_argument('--raw',type=Path);parser.add_argument('--additional',type=Path)
 parser.add_argument('--screening-report',type=Path);parser.add_argument('--plan-manifest',type=Path,default=PLAN);parser.add_argument('--output',type=Path,required=True)
 args=parser.parse_args(argv)
 if args.stage!='lowinput-screen' and (not args.primary or not args.raw):parser.error('screening and confirmation require primary and raw panels')
 if args.stage=='lowinput-screen' and (args.primary or args.raw or args.additional):parser.error('lowinput-screen accepts only its lowinput panel')
 if args.stage=='confirmation' and not args.additional:parser.error('confirmation requires additional acquisition files')
 if args.stage=='confirmation' and not args.screening_report:parser.error('confirmation requires a passing screening report')
 hashes={str(Path(__file__).resolve()):sha(__file__),str(Path(helper.__file__).resolve()):sha(helper.__file__),str(Path(helper.entrapment.__file__).resolve()):sha(helper.entrapment.__file__)}
 manifest=read(args.plan_manifest);hashes[str(args.plan_manifest.resolve())]=sha(args.plan_manifest)
 for p,h in manifest['files'].items():verify(p,h,hashes)
 verify_candidate_registry(REGISTRY,args.arm,args.candidate_sha256,hashes)
 additional_inputs=read(ADDITIONAL) if args.stage=='confirmation' else None
 if args.stage=='confirmation':
  prior=read(args.screening_report)
  if prior['stage']!='screening' or not prior['gate_pass'] or prior['arm']!=args.arm or prior['candidate_sha256']!=args.candidate_sha256 or prior['plan_manifest_sha256']!=sha(args.plan_manifest):raise ValueError('confirmation requires passing matching frozen screening evidence')
  for p,h in prior['source_hashes'].items():verify(p,h,hashes)
  hashes[str(args.screening_report.resolve())]=sha(args.screening_report)
 panels={};checks=[];loader=None
 for name in ('lowinput','primary','raw','additional'):
  directory=getattr(args,name)
  if directory is None:continue
  seeds=SEEDS['confirmation' if args.stage=='confirmation' else 'screening'][name]
  rows,scale,tables,current_loader=load_panel(directory,name,seeds,args.candidate_sha256,hashes,additional_inputs)
  if loader is not None and loader!=current_loader:raise ValueError('panels have different inherited loader/cache settings')
  loader=current_loader;panels[name]=summarize_panel(rows,scale,name,args.stage,tables);checks+=panels[name]['checks']
 passed=all(c['passed'] is True for c in checks);decision=('ADVANCE TO CONTROL SCREEN' if args.stage=='lowinput-screen' else 'ADVANCE TO UNSEEN CONFIRMATION' if args.stage=='screening' else 'PASS: ELIGIBLE FOR PROMOTION') if passed else rejection_label(checks)
 result=dict(stage=args.stage,arm=args.arm,candidate_sha256=args.candidate_sha256,baseline_sha256=V3_SHA,gate_pass=passed,eligible_for_promotion=passed and args.stage=='confirmation',decision=decision,plan_manifest_sha256=sha(args.plan_manifest),checks=checks,panels=panels,source_hashes=hashes,caveats=['FDP and ranges are descriptive; seeds/files are correlated.','Unseen random draws on known files are not biological holdouts.','Individual losses above10% in entrapment depth or reported human count are flagged, not an additional rejection threshold.','Group and score changes are descriptive; V3 fallback-only byte-equality criteria do not apply to this arm.'])
 args.output.parent.mkdir(parents=True,exist_ok=True);args.output.with_suffix('.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
 def pct(v):return 'undefined' if v is None else f'{100*v:.3f}%'
 lines=[f'# {args.arm}: {args.stage}','',decision,'',f'Candidate `{args.candidate_sha256}` versus accepted V3 `{V3_SHA}`.','','| Panel/subset | Paired runs | Actual q1 empty V3→arm | Honest1% zero V3→arm | Reported H/Y/other V3→arm | Pooled FDP V3→arm | Individual FDP>1% V3→arm |','|---|---:|---:|---:|---|---|---|']
 for name,p in panels.items():
  for subset,ab in p['aggregates'].items():
   a,b=ab['baseline'],ab['candidate'];lines.append(f"| {name}/{subset} | {a['n']} | {a['empty_q1']}→{b['empty_q1']} | {a['zero_depth1']}→{b['zero_depth1']} | {a['human']}/{a['yeast']}/{a['other']}→{b['human']}/{b['yeast']}/{b['other']} | {pct(a['pooled_fdp'])}→{pct(b['pooled_fdp'])} | {a['individual_fdp_above1']}/{a['individual_fdp_defined']}→{b['individual_fdp_above1']}/{b['individual_fdp_defined']} |")
 for name,p in panels.items():
  span=p['relative_ranges']['all'];lines+=['',f"{name} common-file median relative1% range: {pct(span['baseline'])}→{pct(span['candidate'])}; {len(span['common_files'])} comparable files. Individual >10% losses: {len(p['losses_over10pct'])}."]
  for metric,desc in p['descriptive_reported_q1_variability'].items():lines.append(f"Reported-count {name}/{metric} relative range: {pct(desc['median_relative_range']['baseline'])}→{pct(desc['median_relative_range']['candidate'])}; {len(desc['common_range_files'])} same defined files. Per-file ranges/sample CV are in JSON; target range is an improvement gate and human median is a depth guard.")
  for subset,r in p['rescued_only'].items():lines.append(f"Rescued {name}/{subset}: {r['n']} runs, {r['human']}H/{r['yeast']}Y, pooled FDP {pct(r['pooled_fdp'])}, individual >1% {r['individual_fdp_above1']}/{r['individual_fdp_defined']}.")
 lines+=['','## Gates','']+[f"- {'PASS' if c['passed'] is True else 'UNDEFINED' if c['passed'] is None else 'FAIL'}: {c['criterion']}." for c in checks]
 lines+=['','Per-file medians, all paired changes, >10% loss flags, individual FDP, undefined-range exclusions, and verified source hashes are retained in the JSON. Paired and within-variant overlaps separately label actual reported-q1 nondecoy targets and entrapment-1%-depth human proteins; empty/empty is undefined for both.','']+result['caveats']
 args.output.with_suffix('.md').write_text('\n'.join(lines)+'\n');print(decision)
if __name__=='__main__':main()
