#!/usr/bin/env python3
"""Verify the entire fresh study before publishing its compact evidence bundle."""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys
import time
from matrix_study import ROOT, STUDY, cohort_rows, sha, write_json


def expected_cases():
    ionstar=['full_single_'+r['run'] for r in cohort_rows('20')]
    ionstar+=['full_matrix_'+n for n in ['5','10','20']]
    ionstar+=['full_matrix_20_floor'+str(n) for n in [1,3]]
    ionstar+=['full_'+arm+'_20' for arm in
              ['masked','quality_null','quality','run_evidence_null','profile_reassign']]
    astral=['single_'+r['run'] for r in csv.DictReader(open(STUDY/'astral/inputs.tsv'),delimiter='\t')]
    astral+=['matrix_'+g for g in ['200pg','1ng','10ng','100ng','all']]
    astral+=[a+'_all' for a in ['quality_null','quality','masked']]
    return {'ionstar':[(n,s) for n in ionstar for s in range(3)],
            'astral':[(n,s) for n in astral for s in range(3)]}


def validate_job(directory, name):
    p=directory/'logs'/(name+'.json')
    if not p.exists():raise ValueError('Missing completion record: '+str(p))
    record=json.loads(p.read_text())
    if record['returncode']!=0:raise ValueError('Failed job: '+name)
    if not record['outputs']:raise ValueError('No output checksums: '+name)
    for filename,digest in record['outputs'].items():
        path=Path(filename)
        if not path.is_file() or sha(path)!=digest:raise ValueError('Changed output: '+filename)
    signature=record['signature']
    if sha(Path(signature['argv'][0]))!=signature['executable_sha256']:
        raise ValueError('Changed executable: '+name)
    for filename,digest in signature['config_sha256'].items():
        if sha(directory/'config'/filename)!=digest:raise ValueError('Changed configuration: '+filename)
    return record


def verify_sources():
    for filename,source in [('frozen.json','source'),('enhanced-source.json','source_enhanced')]:
        doc=json.loads((STUDY/filename).read_text())
        if 'repos' in doc:
            for repo,item in doc['repos'].items():
                for name,digest in item['files'].items():
                    p=STUDY/source/repo/name
                    if sha(p)!=digest:raise ValueError('Frozen source changed: '+str(p))
        else:
            for name,digest in doc['files'].items():
                p=STUDY/source/name
                if sha(p)!=digest:raise ValueError('Frozen enhanced source changed: '+str(p))


def resource_totals(records):
    required={'wall_seconds','user_seconds','system_seconds','peak_rss_kb'}
    if any(not required.issubset(r) or r.get('exit_code')!=0 for r in records):
        raise ValueError('Incomplete successful resource measurements')
    return {'wall_seconds':sum(r['wall_seconds'] for r in records),
            'cpu_seconds':sum(r.get('user_seconds',0)+r.get('system_seconds',0) for r in records),
            'peak_rss_kb':max((r.get('peak_rss_kb',0) for r in records),default=0)}


def finalize(seeds=(0,1,2)):
    planned=expected_cases()
    cases={d:[(n,s) for n,s in rows if s in seeds] for d,rows in planned.items()}
    inventory={};sources={}
    for dataset,expected in cases.items():
        directory=STUDY if dataset=='ionstar' else STUDY/'astral'
        inventory[dataset]={}
        for name,seed in expected:
            key=name+'/seed'+str(seed)
            record=validate_job(directory,'search_'+name+'_seed'+str(seed))
            inventory[dataset][key]=record
        inputs=list(csv.DictReader(open(directory/'inputs.tsv'),delimiter='\t'))
        frozen=json.loads((directory/'frozen.json').read_text())
        if sha(directory/'human_ecoli_yeast.fasta')!=frozen['fasta_sha256']:
            raise ValueError('Changed FASTA')
        model_hash=frozen.get('model_sha256',json.loads((STUDY/'frozen.json').read_text())['model_sha256'])
        if sha(directory/'rt.urtc')!=model_hash:raise ValueError('Changed RT model')
        if sha(directory/'inputs.tsv')!=frozen.get('inputs_sha256',sha(directory/'inputs.tsv')):
            raise ValueError('Changed input manifest')
        for row in inputs:
            # Input content was hashed at freeze; check size now without rereading all raw data.
            if Path(row['path']).stat().st_size!=int(row['bytes']):raise ValueError('Changed raw input size')
            validate_job(directory,'extract_'+row['run'])
        groups=['5','10','20'] if dataset=='ionstar' else ['200pg','1ng','10ng','100ng','all']
        for group in groups:validate_job(directory,'align_'+group)
        for p in (directory/'evidence').glob('*/summary.json'):
            if dataset=='ionstar' and p.parent.name!='20':continue
            doc=json.loads(p.read_text())
            for name,digest in doc['sources'].items():
                if sha(directory/'aligned'/p.parent.name/name)!=digest:raise ValueError('Changed LFQ source')
            for name,digest in doc['outputs'].items():
                if sha(p.parent/name)!=digest:raise ValueError('Changed evidence input')
            if sha(p.parent/'cell_evidence.npz')!=doc['cell_evidence_sha256']:raise ValueError('Changed cell evidence')
    verify_sources()
    completion={}
    for dataset,expected in planned.items():
        directory=STUDY if dataset=='ionstar' else STUDY/'astral'
        completed=[];pending=[]
        for name,seed in expected:
            key=name+'/seed'+str(seed)
            stamp=directory/'logs'/f'search_{name}_seed{seed}.json'
            if stamp.exists() and json.loads(stamp.read_text()).get('returncode')==0:
                if key not in inventory[dataset]:validate_job(directory,f'search_{name}_seed{seed}')
                completed.append(key)
            else:pending.append(key)
        completion[dataset]={'completed':completed,'pending':pending}
    serial=[inventory['ionstar']['full_single_'+r['run']+'/seed0'] for r in cohort_rows('20')]
    serial += [inventory['ionstar']['full_matrix_'+n+'/seed0'] for n in ['5','10','20']]
    serial.sort(key=lambda r:r['started'])
    if any(a['started']+a['elapsed']>b['started']+.001 for a,b in zip(serial,serial[1:])):
        raise ValueError('Timed searches overlap')
    old=STUDY/'searches/matrix_pilot_floor2/seed0/proteins.tsv'
    new=STUDY/'searches/pilot_baseline_new_binary/seed0/proteins.tsv'
    if sha(old)!=sha(new):raise ValueError('Ordinary matrix regression failed')
    for seed in seeds:
        a=STUDY/'searches/full_matrix_20'/f'seed{seed}'/'proteins.tsv'
        b=STUDY/'searches/full_profile_reassign_20'/f'seed{seed}'/'proteins.tsv'
        if sha(a)!=sha(b):raise ValueError('Profile reassignment changed protein IDs')
    for script in ['matrix_report.py','matrix_astral_report.py']:
        subprocess.run([sys.executable,str(ROOT/'benchmark'/script)],check=True,stdout=subprocess.DEVNULL)
    primary=json.loads((STUDY/'summary.json').read_text())
    astral=json.loads((STUDY/'astral/summary.json').read_text())
    if not all(n+'/seed'+str(s) in primary['searches'] for n,s in cases['ionstar']):
        raise ValueError('Primary summary omitted cases')
    if not set(inventory['astral']).issubset(astral['searches']):raise ValueError('Astral summary case mismatch')
    for result in astral['searches'].values():
        if 'per_sample_support' in result:
            support=result['per_sample_support']
            if support['unmatched_coordinate_rows']:
                raise ValueError('Unmapped validation LFQ rows need investigation')
    scope='full three-seed study' if len(seeds)==3 else 'seed-zero core comparison; repeated-search extension incomplete'
    primary['status']='verified publication scope: '+scope
    astral['status']='verified publication scope: '+scope
    write_json(STUDY/'summary.json',primary)
    write_json(STUDY/'astral/summary.json',astral)
    import matrix_report
    timings={}
    for group in ['5','10','20']:
        single=[matrix_report.resources('search_full_single_'+r['run']+'_seed0') for r in cohort_rows(group)]
        extraction=[matrix_report.resources('extract_'+r['run']) for r in cohort_rows(group)]
        joint=matrix_report.resources('search_full_matrix_'+group+'_seed0')
        alignment=matrix_report.resources('align_'+group)
        timings[group]={'single_searches':resource_totals(single),'joint_search':resource_totals([joint]),
            'feature_extraction':resource_totals(extraction),'alignment':resource_totals([alignment]),
            'single_end_to_end':resource_totals(extraction+single),
            'joint_end_to_end':resource_totals(extraction+[alignment,joint])}
    for name in ['frozen.json','enhanced-source.json','binaries.json','protocol-decisions.json',
                 'inputs.tsv','astral/frozen.json','astral/inputs.tsv','astral/msms_evidence.json',
                 'summary.json','astral/summary.json','route_audit/feature_comparison.json']:
        p=STUDY/name
        if p.exists():sources[name]=sha(p)
    scripts={str(p.relative_to(ROOT)):sha(p) for p in sorted((ROOT/'benchmark').glob('matrix_*.py'))}
    bundle={'status':'complete' if len(seeds)==3 else 'complete-core',
        'scope':scope,'publication_seeds':list(seeds),
        'execution_completion':completion,
        'planned_case_counts':{d:len(v) for d,v in planned.items()},
        'study':'matrix_20260912','case_counts':{k:len(v) for k,v in inventory.items()},
        'sources':sources,'analysis_scripts_sha256':scripts,'ionstar':primary,'astral':astral,
        'serialized_resources':timings,
        'evidence':{'ionstar_20':json.loads((STUDY/'evidence/20/summary.json').read_text()),
                    'astral':{g:json.loads((STUDY/'astral/evidence'/g/'summary.json').read_text())
                              for g in ['200pg','1ng','10ng','100ng','all']}},
        'verification':{'profile_ids_unchanged_publication_seeds':True,'output_checksums_verified':True,
            'ordinary_matrix_binary_regression_byte_identical':True,
            'raw_input_checks':'Content SHA256 at freeze; file size at completion',
            'timing_scope':'IonStar baseline seed 0 only; extraction, alignment, and searches measured separately; RT cache disabled for timed searches'},
        'provenance':{'baseline':json.loads((STUDY/'frozen.json').read_text()),
                      'enhanced':json.loads((STUDY/'enhanced-source.json').read_text()),
                      'binaries':json.loads((STUDY/'binaries.json').read_text())}}
    bundle['route_feature_comparison']=json.loads((STUDY/'route_audit/feature_comparison.json').read_text())
    inventory_name='case-inventory.json' if len(seeds)==3 else 'core-case-inventory.json'
    write_json(STUDY/inventory_name,inventory)
    bundle['case_inventory_sha256']=sha(STUDY/inventory_name)
    write_json(STUDY/('final.json' if len(seeds)==3 else 'core.json'),bundle)
    target=ROOT/'paper/figdata/koth_matrix_study.json'
    write_json(target,bundle)
    print('Verified and published',len(inventory['ionstar']),'IonStar and',len(inventory['astral']),'Astral searches')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--status',action='store_true');p.add_argument('--core',action='store_true',help='Verify and publish the complete default-seed comparison, explicitly labeled as a core subset');args=p.parse_args()
    if args.status:
        expected=expected_cases()
        for dataset,cases in expected.items():
            directory=STUDY if dataset=='ionstar' else STUDY/'astral'
            stamps=[directory/'logs'/f'search_{n}_seed{s}.json' for n,s in cases]
            completed=sum(p.exists() and json.loads(p.read_text()).get('returncode')==0 for p in stamps)
            print(dataset,completed,'/',len(cases))
    else:finalize((0,) if args.core else (0,1,2))
