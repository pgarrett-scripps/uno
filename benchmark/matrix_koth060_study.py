#!/usr/bin/env python3
"""LFQ-only KOTH 0.6.0 update; retain every cohort, arm and seed from 0.5.0."""
import argparse,json,shutil,subprocess,tarfile,io
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import matrix_release_study as base
from joint_evidence_audit import sha
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'benchmark/results/matrix_20260914_koth060'
OLD=Path('/home/patrick-garrett/Repos/uno/benchmark/results/matrix_20260913_koth050').resolve()
RELEASE=ROOT.parent/'release'
base.ROOT=ROOT;base.OUT=OUT
rows=base.rows;cohorts=base.cohorts
write=base.write

def freeze():
    assert not (OUT/'frozen.json').exists(), 'Use the frozen study; do not replace its inputs.'
    OUT.mkdir(parents=True,exist_ok=True);(OUT/'bin').mkdir()
    old=json.loads((OLD/'frozen.json').read_text());release=json.loads((RELEASE/'koth-release.json').read_text())
    assert release['tagName']=='v0.6.0' and not release['isDraft']
    archive=next(RELEASE.glob('*.tar.gz'))
    asset=next(a for a in release['assets'] if a['name']==archive.name)
    assert 'sha256:'+sha(archive)==asset['digest']
    for name in ['uno','koth_ff']:shutil.copy2(OLD/'bin'/name,OUT/'bin'/name)
    shutil.copy2(next((RELEASE/'unpacked').glob('*/koth_align')),OUT/'bin/koth_align')
    version=subprocess.check_output([str(OUT/'bin/koth_align'),'--version'],text=True).strip();assert '0.6.0' in version
    repo=Path('/home/patrick-garrett/Repos/koth_rust');tag=subprocess.check_output(['git','rev-parse','v0.6.0^{commit}'],cwd=repo,text=True).strip()
    assert tag=='7cc80fb846254f958fa3df6f2950cf83aaa460d5'
    for name in ['uno','sage','Sequoia']:shutil.copytree(OLD/'source'/name,OUT/'source'/name)
    shutil.copytree(OLD/'source/koth_rust',OUT/'source/koth_detection')
    with tarfile.open(fileobj=io.BytesIO(subprocess.check_output(['git','archive','v0.6.0'],cwd=repo))) as t:
        t.extractall(OUT/'source/koth_rust',filter='data')
    source={k:v for k,v in old['source'].items() if k!='koth_rust'}
    source['koth_detection']=old['source']['koth_rust']
    kr=OUT/'source/koth_rust';source['koth_rust']=dict(head=tag,files={str(p.relative_to(kr)):sha(p) for p in kr.rglob('*') if p.is_file()})
    reuse=[]
    for dataset in ['astral','ionstar']:
        d=OUT/dataset;d.mkdir();shutil.copytree(OLD/dataset/'config',d/'config')
        for name in ['inputs.tsv','human_ecoli_yeast.fasta','rt.urtc','previous-frozen.json']:shutil.copy2(OLD/dataset/name,d/name)
        # Explicitly pin the release defaults so future versions cannot change this study.
        cfg=d/'config/align.toml';cfg.write_text(cfg.read_text().replace('[lfq]\n','[lfq]\nrt_window_pct = 0.005\nim_tolerance = 0.015\n'))
        (d/'koth').symlink_to(OLD/dataset/'koth');(d/'searches').mkdir()
        for r in rows(dataset):
            p=OLD/dataset/'jobs'/('extract_'+r['run'])/'job.json';reuse.append(dict(kind='extraction',job=str(p),sha256=sha(p)))
            dest=d/'searches'/('single_'+r['run']);dest.symlink_to(OLD/dataset/'searches'/dest.name)
            for seed in [0,1,2]:
                p=dest/f'seed{seed}/job.json';reuse.append(dict(kind='single-search',job=str(p.resolve()),sha256=sha(p)))
    supplement=json.loads((OLD/'source-supplement.json').read_text())
    supplement['files']['koth_detection']=supplement['files'].pop('koth_rust')
    supplement['note']='Inherited unchanged UNO/runtime and feature-detection resources; complete released 0.6.0 LFQ source is archived separately.'
    write(OUT/'source-supplement.json',supplement);write(OUT/'reuse.json',reuse);write(OUT/'koth-release.json',release)
    frozen=dict(koth_version=version,detection_version='0.5.0',seeds=[0,1,2],threads_per_job=8,workers=4,
        planned_searches=48,reused_single_searches=96,reused_feature_tables=32,planned_matrices=8,serialized_searches=4,
        source=source,binary={n:sha(OUT/'bin'/n) for n in ['uno','koth_ff','koth_align']},
        parent_archive=str(OLD),parent_report_sha256=sha(OLD/'report.json'),
        changes='LFQ RT half-window 0.01 -> 0.005 of observed run span; IM half-window 0.05 -> 0.015 absolute 1/K0. IM window is inert on these Orbitrap inputs.',
        reuse_reason='Release diff changes only LFQ defaults, version and documentation/tests. Detection, alignment, grouping, schemas, UNO executable, FASTA, RT model and scoring settings unchanged.',
        protocol='All eight cohorts, both matrix arms and all three seeds retained. Regenerate LFQ and all 48 matrix searches, including abundance outputs. Reuse 96 unchanged separate searches. Verify original-observation bundles and non-cell protein outputs for exact parity. No outcome-based parameter changes.',
        evaluation=old['evaluation'],grouping=old['grouping'],timing='Four serialized warm-cache controls after all parallel jobs end; concurrent timings are not speed evidence.',
        driver_sha256=sha(Path(__file__)),release_asset_sha256=sha(archive))
    write(OUT/'frozen.json',frozen);shutil.copy2(__file__,OUT/Path(__file__).name)
    # Reuse predictions read-only at the input stage, then let each study own its cache.
    shutil.copytree(OLD/'rt-cache',OUT/'rt-cache')
    print('Frozen 8 matrices, 48 new searches, 96 reused separate searches.',flush=True)

def parallel(tasks,workers=4):
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures=[pool.submit(fn,*args) for fn,args in tasks]
        for f in futures:f.result()

def run():
    if not (OUT/'frozen.json').exists():freeze()
    for dataset in ['astral','ionstar']:
        parallel([(base.align,(dataset,n,rr)) for n,rr in cohorts(dataset).items()],2)
        # Grouping/coordinate evidence must be unchanged before searching.
        for n in cohorts(dataset):
            a=OUT/dataset/'aligned'/n/'feature_observations.tsv';b=OLD/dataset/'aligned'/n/'feature_observations.tsv'
            assert sha(a)==sha(b), ('unexpected observation change',dataset,n)
        parallel([(base.search,(dataset,arm+'_'+n,seed,arm,OUT/dataset/'aligned'/n))
                  for seed in [0,1,2] for n in cohorts(dataset) for arm in ['joint','joint_cell']])
    for r in cohorts('astral')['100ng']:
        base.search('astral','single_'+r['run'],0,'single',OUT/'astral/koth'/r['run']/'features.tsv',True)
    base.search('astral','joint_cell_100ng',0,'joint_cell',OUT/'astral/aligned/100ng',True)
    write(OUT/'complete.json',dict(new_searches=48,reused_searches=96,serialized_searches=4,reused_feature_tables=32,matrices=8))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['freeze','run']);a=p.parse_args();globals()[a.action]()
