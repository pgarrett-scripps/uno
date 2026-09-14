#!/usr/bin/env python3
"""Release-pinned matrix study: 32 fresh feature tables, 8 matrices, 144 searches.
Astral is scheduled first; all cohorts and three seeds are fixed before execution.
"""
import argparse,csv,json,os,shutil,subprocess,time
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
from joint_evidence_audit import sha
ROOT=Path(__file__).resolve().parents[1]
OLD=Path('/home/patrick-garrett/Repos/uno/benchmark/results/matrix_20260912')
OUT=ROOT/'benchmark/results/matrix_20260913_koth050'

def write(path,obj):path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(obj,indent=2)+'\n')
def rows(dataset):return list(csv.DictReader((OUT/dataset/'inputs.tsv').open(),delimiter='\t'))
def cohorts(dataset):
    rr=rows(dataset)
    if dataset=='astral':return {g:[r for r in rr if g=='all' or r['condition']==g] for g in ['200pg','1ng','10ng','100ng','all']}
    return {str(n):[r for r in rr if int(r['replicate_order'])<=n//5] for n in [5,10,20]}
def freeze():
    OUT.mkdir(parents=True,exist_ok=True)
    if (OUT/'frozen.json').exists():return
    (OUT/'bin').mkdir(exist_ok=True)
    for name in ['uno','koth_ff','koth_align']:shutil.copy2(ROOT.parent/'target/release'/name,OUT/'bin'/name)
    version=subprocess.check_output([str(OUT/'bin/koth_ff'),'--version'],text=True).strip()
    assert '0.5.0' in version,version
    source={}
    for repo in ['uno','koth_rust','sage','Sequoia']:
        root=ROOT if repo=='uno' else ROOT.parent/repo
        files=[p for p in root.rglob('*.rs') if 'target' not in p.parts and 'benchmark' not in p.parts and 'experiments' not in p.relative_to(root).parts]
        files+=[p for p in root.rglob('Cargo.toml') if 'target' not in p.parts and 'experiments' not in p.relative_to(root).parts]
        if (root/'Cargo.lock').exists():files.append(root/'Cargo.lock')
        checks={}
        for p in sorted(set(files)):
            rel=p.relative_to(root);dst=OUT/'source'/repo/rel;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dst);checks[str(rel)]=sha(dst)
        source[repo]=dict(head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),files=checks)
    for dataset in ['astral','ionstar']:
        old=OLD/'astral' if dataset=='astral' else OLD
        d=OUT/dataset;cfg=d/'config';cfg.mkdir(parents=True,exist_ok=True)
        for name in ['inputs.tsv','human_ecoli_yeast.fasta','rt.urtc','frozen.json']:shutil.copy2(old/name,d/name)
        (d/'frozen.json').rename(d/'previous-frozen.json')
        shutil.copy2(OLD/'config/koth.toml',cfg/'koth.toml')
        shutil.copy2(ROOT/'benchmark/results/lfq_matrix_fixes/uno.toml',cfg/'uno.toml')
        (cfg/'cell.toml').write_text((cfg/'uno.toml').read_text()+'\n[scoring]\nmatrix_joint_cell_quality = true\n')
        (cfg/'align.toml').write_text('[lfq]\nrun_tdc = true\nnormalize = "none"\nquant_estimator = "sum"\ndetected_use_grid = true\naveragine_projection = false\n\n[lfq.consensus]\nmin_group_size = 1\nallow_replicated_weak_seeds = true\n\n[output]\nexport_long = true\nformat = "tsv"\nexport_details = true\nexport_decoys = true\n')
        assert len(rows(dataset))==(12 if dataset=='astral' else 20)
    write(OUT/'frozen.json',dict(koth_version=version,seeds=[0,1,2],threads_per_job=8,workers=2,planned_searches=144,
        datasets=['astral','ionstar'],matrix_arms=['joint','joint_cell'],
        extraction='Fresh KOTH 0.5.0 detection from all 32 raw inputs using the existing UNO feature-finding profile.',
        grouping='Fixed member floor 0.5, seed floor 0.75, replicated weak seeds enabled, min group 1, bounded default spans, automatic reference.',
        evaluation='Whole tied-score empirical entrapment FDP<=1% accession curves; independent union ranked by minimum per-run internal q.',
        timing='Main jobs overlap. Dedicated serialized warm-cache controls are reported separately.',
        validation='All four Astral input levels plus all-file matrix; all nested IonStar cohorts. No selection of successful seeds.',
        binary={n:sha(OUT/'bin'/n) for n in ['uno','koth_ff','koth_align']},source=source,
        compiler=subprocess.check_output(['rustc','--version'],text=True).strip(),driver_sha256=sha(Path(__file__))))

def job(directory,argv,outputs,inputs,kind):
    directory.mkdir(parents=True,exist_ok=True)
    argv=list(map(str,argv));outputs=list(map(Path,outputs));inputs=list(map(Path,inputs))
    signature=dict(argv=argv,binary=sha(Path(argv[0])),inputs={str(p):sha(p) for p in inputs},kind=kind)
    stamp=directory/'job.json'
    if stamp.exists():
        old=json.loads(stamp.read_text())
        if old['returncode']==0:
            assert old['signature']==signature,('changed successful job',directory)
            assert all(sha(p)==old['outputs'][str(p)] for p in outputs),directory
            return
        stamp.rename(directory/f'failed-{time.time_ns()}.json')
    env={k:v for k,v in os.environ.items() if not k.startswith('UNO_')}
    env.update(RAYON_NUM_THREADS='8',OMP_NUM_THREADS='8',OPENBLAS_NUM_THREADS='1',RUST_LOG='info',
               UNO_RT_CACHE='1',UNO_RT_CACHE_DIR=str(OUT/'rt-cache'))
    start=time.time();print('Starting',directory.relative_to(OUT),flush=True)
    with (directory/'run.log').open('w') as log:
        p=subprocess.run(['/usr/bin/time','-v','-o',str(directory/'time.txt'),*argv],env=env,stdout=log,stderr=subprocess.STDOUT)
    record=dict(signature=signature,started=start,seconds=time.time()-start,returncode=p.returncode,
                outputs={str(p):sha(p) for p in outputs if p.is_file()})
    write(stamp,record)
    assert p.returncode==0 and len(record['outputs'])==len(outputs),directory
    print('Finished',directory.relative_to(OUT),round(record['seconds'],1),flush=True)

def extract(dataset,row):
    d=OUT/dataset;out=d/'koth'/row['run'];source=Path(row['path'])
    assert sha(source)==row['sha256'],source
    job(d/'jobs'/('extract_'+row['run']),[OUT/'bin/koth_ff',source,'--config',d/'config/koth.toml','--output',d/'koth','--threads','8'],
        [out/n for n in ['features.tsv','hills.tsv','config.toml']],[source,d/'config/koth.toml'],'extraction')
def align(dataset,name,rr):
    d=OUT/dataset;batch=d/'batches'/name;batch.mkdir(parents=True,exist_ok=True)
    for r in rr:
        dst=batch/r['run'];dst.mkdir(exist_ok=True)
        for n in ['features.tsv','hills.tsv','config.toml']:
            p=dst/n
            if not p.exists():p.symlink_to(d/'koth'/r['run']/n)
    assert sorted(p.name for p in batch.iterdir())==sorted(r['run'] for r in rr)
    out=d/'aligned'/name
    job(d/'jobs'/('align_'+name),[OUT/'bin/koth_align',batch,'--config',d/'config/align.toml','--output',out],
        [out/n for n in ['matrix_manifest.json','feature_observations.tsv','lfq_matrix.long.tsv','intensity_matrix.tsv','qvalue_matrix.tsv']],
        [d/'config/align.toml',*[d/'koth'/r['run']/n for r in rr for n in ['features.tsv','hills.tsv']]],'alignment')
def search(dataset,name,seed,arm,source,timing=False):
    d=OUT/dataset;out=d/('timing' if timing else 'searches')/name/f'seed{seed}'
    cfg=d/'config'/('cell.toml' if arm=='joint_cell' else 'uno.toml')
    argv=[OUT/'bin/uno','--fasta',d/'human_ecoli_yeast.fasta','--rt-model',d/'rt.urtc','--seed',str(seed),'--config',cfg,
        '--out-proteins',out/'proteins.tsv','--out-lfq',out/'lfq.tsv','--out-abundance',out/'abundance.tsv']
    inputs=[cfg,d/'human_ecoli_yeast.fasta',d/'rt.urtc']
    if arm=='single':argv+=['--features',source];inputs+=[source]
    else:argv+=['--matrix-long',source,'--matrix-joint'];inputs+=[source/'matrix_manifest.json']
    job(out,argv,[out/n for n in ['proteins.tsv','lfq.tsv','abundance.tsv']],inputs,'serial-timing' if timing else 'search')
def parallel(tasks):
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures=[pool.submit(fn,*args) for fn,args in tasks]
        errors=[]
        for f in as_completed(futures):
            try:f.result()
            except Exception as e:errors.append(repr(e));print('FAILED',repr(e),flush=True)
        if errors:raise RuntimeError(errors)
def run():
    freeze()
    for dataset in ['astral','ionstar']:
        parallel([(extract,(dataset,r)) for r in rows(dataset)])
        parallel([(align,(dataset,n,rr)) for n,rr in cohorts(dataset).items()])
        tasks=[]
        for seed in [0,1,2]:
            tasks += [(search,(dataset,'single_'+r['run'],seed,'single',OUT/dataset/'koth'/r['run']/'features.tsv')) for r in rows(dataset)]
            tasks += [(search,(dataset,arm+'_'+n,seed,arm,OUT/dataset/'aligned'/n)) for n in cohorts(dataset) for arm in ['joint','joint_cell']]
        parallel(tasks)
    # One explicitly serialized, warmed comparison; ordinary multiworker timings
    # are never presented as sequential resource measurements.
    for r in cohorts('astral')['100ng']:
        search('astral','single_'+r['run'],0,'single',OUT/'astral/koth'/r['run']/'features.tsv',True)
    search('astral','joint_cell_100ng',0,'joint_cell',OUT/'astral/aligned/100ng',True)
    write(OUT/'complete.json',dict(searches=144,serialized_searches=4,feature_tables=32,matrices=8))
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['freeze','run']);a=parser.parse_args();globals()[a.action]()
