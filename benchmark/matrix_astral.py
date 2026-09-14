#!/usr/bin/env python3
"""Frozen-settings validation on exactly the manuscript's twelve Astral files."""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor

import matrix_study as study
PARENT = study.STUDY
ASTUDY = PARENT/'astral'
sys.path.insert(0,str(study.ROOT/'benchmark/raw_mode'))
from summarize import ASTRAL, ASTRAL_STEMS


def freeze():
    if (ASTUDY/'frozen.json').exists(): return
    ASTUDY.mkdir(parents=True,exist_ok=True)
    for name in ['config','bin']:
        dst=ASTUDY/name
        if not dst.exists(): dst.symlink_to(PARENT/name,target_is_directory=True)
    shutil.copy2(PARENT/'rt.urtc',ASTUDY/'rt.urtc')
    # search() uses this fixed filename; here it deliberately contains NO E. coli.
    fasta=ASTUDY/'human_ecoli_yeast.fasta'
    shutil.copy2(study.ROOT/'benchmark/data/entrap_human_yeast.fasta',fasta)
    counts={s:0 for s in ['HUMAN','ECOLI','YEAST']}
    for h,_ in study.fasta_records(fasta):
        for sp in counts: counts[sp]+=('_'+sp in h.split()[0])
    assert counts['ECOLI']==0 and counts['HUMAN']>0 and counts['YEAST']>0
    rows=[]
    for stem in ASTRAL_STEMS:
        p=ASTRAL/'mzML'/(stem+'.mzML')
        assert p.is_file(),p
        level=stem.split('OldFAIMS_')[1].split('__')[0]
        rows.append({'run':stem,'path':str(p),'condition':level,'pilot':False,
                     'bytes':p.stat().st_size,'sha256':study.sha(p)})
        print('Frozen Astral input',stem,flush=True)
    with open(ASTUDY/'inputs.tsv','w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    study.write_json(ASTUDY/'frozen.json',{'parent':str(PARENT/'frozen.json'),
        'settings_frozen_before_validation':True,'fasta_counts':counts,
        'fasta_sha256':study.sha(fasta),'inputs_sha256':study.sha(ASTUDY/'inputs.tsv'),
        'seeds':[0,1,2], 'threads':8,
        'groups':['200pg','1ng','10ng','100ng','all'],
        'cohort_source':'benchmark/raw_mode/summarize.py ASTRAL_STEMS',
        'cohort_source_sha256':study.sha(study.ROOT/'benchmark/raw_mode/summarize.py')})


def run():
    freeze()
    study.STUDY=ASTUDY
    from matrix_evidence import prepare
    rows=list(csv.DictReader(open(ASTUDY/'inputs.tsv'),delimiter='\t'))
    study.extract(pilot=False)
    for group in ['200pg','1ng','10ng','100ng','all']:
        chosen=[r for r in rows if group=='all' or r['condition']==group]
        batch=ASTUDY/'batches'/group;batch.mkdir(parents=True,exist_ok=True)
        for r in chosen:
            sub=batch/r['run'];sub.mkdir(exist_ok=True)
            for name in ['features.tsv','hills.tsv','config.toml']:
                if not (sub/name).exists(): (sub/name).symlink_to(ASTUDY/'koth'/r['run']/name)
        out=ASTUDY/'aligned'/group
        study.run_job('align_'+group,[ASTUDY/'bin/koth_align',batch,'--config',ASTUDY/'config/align.toml',
                                     '--output',out],[out/n for n in
                    ['intensity_matrix.tsv','qvalue_matrix.tsv','lfq_details.tsv','consensus_features.tsv']])
        prepare(group)
    for r in rows:
        study.search('single_'+r['run'],'--features',ASTUDY/'koth'/r['run']/'features.tsv',cache='0')
    for group in ['200pg','1ng','10ng','100ng','all']:
        study.search('matrix_'+group,'--matrix',ASTUDY/'aligned'/group/'intensity_matrix.tsv',cache='0')
        for arm in (['quality_null','quality','masked'] if group=='all' else []):
            study.search(arm+'_'+group,'--matrix',ASTUDY/'evidence'/group/(arm+'.tsv'),
                         cache='0',binary='uno_evidence')
    def extra(item):
        row,seed=item
        study.search('single_'+row['run'],'--features',ASTUDY/'koth'/row['run']/'features.tsv',seed=seed,cache='0')
    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(extra,[(r,s) for s in [1,2] for r in rows]))
    for seed in [1,2]:
        for group in ['200pg','1ng','10ng','100ng','all']:
            study.search('matrix_'+group,'--matrix',ASTUDY/'aligned'/group/'intensity_matrix.tsv',seed=seed,cache='0')
            for arm in (['quality_null','quality','masked'] if group=='all' else []):
                study.search(arm+'_'+group,'--matrix',ASTUDY/'evidence'/group/(arm+'.tsv'),
                             seed=seed,cache='0',binary='uno_evidence')
    study.write_json(ASTUDY/'completed.json',{'searches':60,'status':'all prespecified searches completed'})


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['freeze','run'])
    if p.parse_args().action=='freeze':freeze()
    else:run()
