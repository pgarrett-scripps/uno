#!/usr/bin/env python3
"""Audit the LFQ-only 0.6.0 update and emit manuscript evidence, without selection.
Run with paper/analysis's locked environment (numpy required). All tied score
blocks enter together. Entrapment cutoffs are retrospective evaluation only.
"""
import csv,json,math,shutil,re
from collections import Counter,defaultdict
from pathlib import Path
from statistics import median
from decimal import Decimal
from matrix_koth060_study import OUT,ROOT,OLD,cohorts,rows
from joint_evidence_audit import sha,curve
from matrix_report import load_proteins,list_counts,accession_union,quant_metrics

def read(p):return json.loads(p.read_text())
def write(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def scores(tables,union=False):
    out={}
    for tt in tables:
        for r in tt:
            if r['species']=='AMBIGUOUS':continue
            for a in r['members']:out[a]=max(out.get(a,-math.inf),-r['q'] if union else r['score'])
    return out

def count(tables,scale,union=False):
    result=curve(scores(tables,union),scale)
    result['internal_1pct']=list_counts(accession_union(tables),scale)
    result['ambiguous_groups']=sum(r['species']=='AMBIGUOUS' for tt in tables for r in tt)
    return result

def independent_check(files,counts,expected,union=False):
    """Independent Decimal ranking and integer FDP arithmetic, including sets."""
    blocks=defaultdict(set);ranks={}
    for path in files:
        for row in csv.DictReader(path.open(),delimiter='\t'):
            if row['decoy']!='0':continue
            members=row['group_members'].split(';')
            kinds=[{sp for sp in counts if '_'+sp in acc} for acc in members]
            if any(len(k)!=1 for k in kinds) or len(set.union(*kinds))!=1:continue
            rank=-Decimal(row['q_value']) if union else Decimal(row['binomial_score'])
            for acc in members:ranks[acc]=max(ranks.get(acc,Decimal('-Infinity')),rank)
    for acc,rank in ranks.items():blocks[rank].add(acc)
    n=y=0;accepted=set();best={t:(0,0,set()) for t in [1,3,5]}
    for rank in sorted(blocks,reverse=True):
        for acc in blocks[rank]:
            if '_YEAST' in acc:y+=1
            else:n+=1;accepted.add(acc)
        for t in best:
            if n>best[t][0] and y*(counts['HUMAN']+counts['ECOLI'])*100<=t*n*counts['YEAST']:
                best[t]=(n,y,accepted.copy())
    for t,(n,y,acc) in best.items():
        e=expected[str(t/100)]
        assert (n,y,acc)==(e['present'],e['yeast'],set(e['accessions'])),(files,t,n,y,e['present'],e['yeast'])

def resources(p):
    j=read(p/'job.json');r=dict(wall_seconds=j['seconds'])
    for line in (p/'time.txt').read_text().splitlines():
        for key,dst in [('Maximum resident set size (kbytes)','peak_rss_kb'),('User time (seconds)','user_seconds'),('System time (seconds)','system_seconds')]:
            if key+':' in line:r[dst]=float(line.rsplit(':',1)[1])
    return r

def cells(p):
    manifest=read(p/'matrix_manifest.json')
    source_hashes={name:sha(p/name) for name in ['matrix_manifest.json','lfq_matrix.long.tsv','feature_observations.tsv']}
    source_hashes['producer']=sha(Path(__file__))
    for key in ['cells','observations']:
        assert source_hashes[manifest[key]['path']]==manifest[key]['sha256']
    cache=p/'extraction-summary-cache.json'
    if cache.exists():
        saved=read(cache)
        if saved['source_hashes']==source_hashes:return saved['result']
    rr=defaultdict(Counter);groups=set();total=Counter()
    with (p/'lfq_matrix.long.tsv').open() as f:
        reader=csv.reader(f,delimiter='\t');header=next(reader)
        columns=[header.index(k) for k in ['consensus_id','run_id','extraction_kind','lfq_intensity','lfq_is_mbr','lfq_q_value']]
        group_i,run_i,kind_i,intensity_i,mbr_i,q_i=columns
        for row in reader:
            if row[kind_i]!='target':continue
            groups.add(row[group_i]);c=rr[row[run_i]];c['cells']+=1
            if float(row[intensity_i] or 0)<=0:continue
            kind='mbr' if row[mbr_i]=='true' else 'direct'
            c['positive']+=1;c[kind+'_positive']+=1
            q=row[q_i]
            if q and math.isfinite(float(q)):
                c['q_available']+=1
                if float(q)<=.01:c['confident']+=1;c[kind+'_confident']+=1
    for c in rr.values():total.update(c)
    assert total['cells']==manifest['n_consensus']*len(manifest['runs'])
    assert len(groups)==manifest['n_consensus']
    result=dict(groups=len(groups),totals=dict(total),runs={r['name']:dict(rr[str(r['run_id'])]) for r in manifest['runs']})
    write(cache,dict(source_hashes=source_hashes,result=result))
    return result

def main():
    assert read(OUT/'complete.json')==dict(new_searches=48,reused_searches=96,serialized_searches=4,reused_feature_tables=32,matrices=8)
    print('Verifying complete job inventory, sources, and input/output hashes.',flush=True)
    frozen=read(OUT/'frozen.json');assert frozen['seeds']==[0,1,2]
    hashes={};kinds=Counter()
    for p in sorted(OUT.rglob('job.json')):
        j=read(p);assert j['returncode']==0,p;kinds[j['signature']['kind']]+=1
        assert sha(Path(j['signature']['argv'][0]))==j['signature']['binary'],p
        for name,h in {**j['signature']['inputs'],**j['outputs']}.items():
            if name not in hashes:hashes[name]=sha(Path(name))
            assert hashes[name]==h,name
        hashes[str(p)]=sha(p)
    assert kinds==dict(alignment=8,search=48,**{'serial-timing':4}),kinds
    reuse=read(OUT/'reuse.json');reused=Counter()
    for item in reuse:
        p=Path(item['job']);assert sha(p)==item['sha256'],p
        j=read(p);assert j['returncode']==0,p;reused[item['kind']]+=1
        assert sha(Path(j['signature']['argv'][0]))==j['signature']['binary'],p
        # Reuse the already audited raw-input receipts; verify every reused output.
        checks=j['outputs'] if item['kind']=='extraction' else {**j['signature']['inputs'],**j['outputs']}
        for name,h in checks.items():
            if name not in hashes:hashes[name]=sha(Path(name))
            assert hashes[name]==h,name
        hashes[str(p)]=sha(p)
    assert reused==dict(extraction=32,**{'single-search':96}),reused
    assert sha(OLD/'report.json')==frozen['parent_report_sha256']
    for repo,manifest in frozen['source'].items():
        for name,h in manifest['files'].items():assert sha(OUT/'source'/repo/name)==h
    for name,h in frozen['binary'].items():assert sha(OUT/'bin'/name)==h
    supplement=read(OUT/'source-supplement.json')
    for repo,files in supplement['files'].items():
        for name,h in files.items():assert sha(OUT/'source'/repo/name)==h
    for name,h in supplement['runtime'].items():assert sha(Path(name))==h
    oldms=Path('/home/patrick-garrett/Repos/uno/benchmark/results/matrix_20260912/astral/msms_evidence.json')
    ms=read(oldms);assert sha(Path(ms['source']))==ms['source_sha256']
    shutil.copy2(oldms,OUT/'astral/msms_evidence.json');hashes[str(oldms)]=sha(oldms)
    result=dict(status='complete',publication_seeds=[0,1,2],provenance=frozen,execution_amendment=dict(workers=4,scope='LFQ-only rerun; unchanged single searches reused'),reuse_counts=dict(reused),source_supplement=supplement,job_counts=dict(kinds),datasets={})
    assert sha(OUT/'matrix_koth060_study.py')==frozen['driver_sha256']
    for name in ['frozen.json','reuse.json','matrix_koth060_study.py','koth-release.json','source-supplement.json']:
        hashes[str(OUT/name)]=sha(OUT/name)
    for dataset in ['astral','ionstar']:
        d=OUT/dataset;fc=read(d/'previous-frozen.json')['fasta_counts']
        # Independently count the exact searched FASTA, not merely inherited metadata.
        observed=Counter()
        for line in (d/'human_ecoli_yeast.fasta').open():
            if line.startswith('>'):
                for sp in ['HUMAN','ECOLI','YEAST']:
                    if '_'+sp in line.split()[0]:observed[sp]+=1
        assert all(observed[sp]==fc[sp] for sp in fc),(observed,fc)
        scale=(fc['HUMAN']+fc['ECOLI'])/fc['YEAST'];dd=dict(scale=scale,cohorts={});result['datasets'][dataset]=dd
        for name,rr in cohorts(dataset).items():
            print('Evaluating',dataset,name,'all seeds, extraction cells, and independent checks.',flush=True)
            cd=dict(runs=[r['run'] for r in rr],seeds={},cells=cells(d/'aligned'/name));dd['cohorts'][name]=cd
            supported={p for r in rr for p in ms['runs'][r['run']]['proteins']} if dataset=='astral' else set()
            for seed in [0,1,2]:
                single_paths=[d/'searches'/('single_'+r['run'])/f'seed{seed}/proteins.tsv' for r in rr]
                singles=[load_proteins(p) for p in single_paths]
                sd={'union':count(singles,scale,True)};cd['seeds'][str(seed)]=sd
                independent_check(single_paths,fc,sd['union'],True)
                for arm in ['joint','joint_cell']:
                    path=d/'searches'/(arm+'_'+name)/f'seed{seed}'
                    sd[arm]=count([load_proteins(path/'proteins.tsv')],scale)
                    independent_check([path/'proteins.tsv'],fc,sd[arm])
                for arm,v in sd.items():
                    accepted=set(v['0.01']['accessions']);base=set(sd['union']['0.01']['accessions'])
                    added=accepted-base;lost=base-accepted
                    v['comparison']=dict(added=sorted(added),lost=sorted(lost),intersection=len(accepted&base),jaccard=len(accepted&base)/len(accepted|base) if accepted|base else None)
                    if dataset=='astral':
                        access=lambda x:x.split('|')[1] if x.count('|')>=2 else x
                        v['corroboration']=dict(accepted=len(accepted),supported=sum(access(x) in supported for x in accepted),added=len(added),added_supported=sum(access(x) in supported for x in added))
                if dataset=='ionstar' and name in ['10','20']:
                    paths={a:d/'searches'/(a+'_'+name)/f'seed{seed}/abundance.tsv' for a in ['joint','joint_cell']}
                    def groups(p):return {tuple(sorted(r['group_members'].split(';'))) for r in csv.DictReader(p.open(),delimiter='\t')}
                    common=set.intersection(*(groups(p) for p in paths.values()))
                    baseline=quant_metrics(paths['joint']);factors=baseline['size_factors_log2'] if baseline else None
                    for a,p in paths.items():
                        sd[a]['quant_all']=quant_metrics(p,fixed_factors=factors)
                        sd[a]['quant_common']=quant_metrics(p,fixed_factors=factors,common_groups=common)
    timing=OUT/'astral/timing';r=cohorts('astral')['100ng'];single=[resources(timing/('single_'+x['run'])/'seed0') for x in r]
    combined={k:(max(x[k] for x in single) if k=='peak_rss_kb' else sum(x[k] for x in single)) for k in single[0]}
    joint=resources(timing/'joint_cell_100ng/seed0')
    for p in timing.glob('*/seed0'):
        # Timing repetition must reproduce the prespecified seed-zero identification.
        assert sha(p/'proteins.tsv')==sha(OUT/'astral/searches'/p.parent.name/'seed0/proteins.tsv'),p
        misses=re.findall(r'RT cache .*: (\d+) unique misses of', (p/'run.log').read_text())
        assert misses and all(int(n)==0 for n in misses),('timing cache not fully warm',p,misses)
    parity=[]
    for dataset in ['astral','ionstar']:
        for name in cohorts(dataset):
            file='feature_observations.tsv'
            assert sha(OUT/dataset/'aligned'/name/file)==sha(OLD/dataset/'aligned'/name/file)
            for seed in [0,1,2]:
                rel=Path(dataset)/'searches'/('joint_'+name)/f'seed{seed}'/'proteins.tsv'
                parity.append(dict(dataset=dataset,cohort=name,seed=seed,identical=sha(OUT/rel)==sha(OLD/rel)))
    result['coordinate_parity']='All eight original-observation tables are byte-identical to 0.5.0.'
    result['joint_identification_parity']=parity
    result['koth_release']=read(OUT/'koth-release.json')
    assert not result['koth_release']['isDraft'] and result['koth_release']['tagName']=='v0.6.0'
    result['timing']=dict(cohort='Astral 100 ng triplicate',cache='warm shared RT cache',repeats=1,single=combined,joint_cell=joint,speedup=combined['wall_seconds']/joint['wall_seconds'])
    result['independent_validation']='Decimal rankings, exact integer FDP comparisons, and accession-set checks passed for every cohort/arm/seed at 1%, 3%, and 5%';result['hashes']=hashes;result['report_source_sha256']=sha(Path(__file__))
    write(OUT/'report.json',result)
    write(ROOT/'paper/figdata/koth_matrix_study.json',result)
    for dataset,dd in result['datasets'].items():
        for name,cd in dd['cohorts'].items():
            print(dataset,name,{a:[v[a]['0.01']['present'] for v in cd['seeds'].values()] for a in ['union','joint','joint_cell']})
    print('Serialized speedup',result['timing']['speedup'])
if __name__=='__main__':main()
