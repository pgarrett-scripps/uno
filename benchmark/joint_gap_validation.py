#!/usr/bin/env python3
"""Fixed cross-dataset validation on the preselected Astral 100ng triplicate."""
import json,os,shutil,subprocess,time
from pathlib import Path
import joint_matrix_pilot as p
from matrix_report import load_proteins
from joint_evidence_audit import curve
p.OUT=p.BASE/'joint_gap_v1/astral_validation'
p.OLD=p.OLD/'astral'
frozen=json.loads((p.OUT/'frozen.json').read_text())
if not (p.OUT/'uno').exists():
    shutil.copy2(p.BASE/'joint_gap_alternatives_v1/uno',p.OUT/'uno')
    shutil.copy2(p.BASE/'joint_gap_alternatives_v1/source-manifest.json',p.OUT/'source-manifest.json')
if not (p.OUT/'koth_align').exists():shutil.copy2(p.BASE/'joint_search_v2/koth_align',p.OUT/'koth_align')
aligned=p.OUT/'aligned';aligned.mkdir(exist_ok=True)
stamp=aligned/'job.json'
av=[str(p.OUT/'koth_align'),str(p.OLD/'batches/100ng'),'--config',str(p.OUT/'align.toml'),'--output',str(aligned)]
sig=dict(argv=av,binary=p.sha(p.OUT/'koth_align'),config=p.sha(p.OUT/'align.toml'),
         sources={str(p.OLD/'koth'/r/f):p.sha(p.OLD/'koth'/r/f) for r in frozen['runs'] for f in ['features.tsv','hills.tsv']})
if stamp.exists():
    old=json.loads(stamp.read_text());assert old['signature']==sig and old['returncode']==0
    assert all(p.sha(Path(f))==h for f,h in old['outputs'].items())
else:
    t=time.time();print('Starting Astral matrix',flush=True)
    with (aligned/'run.log').open('w') as log:
        proc=subprocess.run(av,stdout=log,stderr=subprocess.STDOUT,env=dict(os.environ,RAYON_NUM_THREADS='8'))
    outputs={str(f):p.sha(f) for f in aligned.glob('*.tsv')}
    if (aligned/'matrix_manifest.json').exists():outputs[str(aligned/'matrix_manifest.json')]=p.sha(aligned/'matrix_manifest.json')
    stamp.write_text(json.dumps(dict(signature=sig,returncode=proc.returncode,seconds=time.time()-t,outputs=outputs),indent=2)+'\n')
    assert proc.returncode==0
    print('Finished Astral matrix',round(time.time()-t,1),flush=True)
report=json.loads((p.BASE/'joint_gap_v1/report.json').read_text())
# Select only on development IonStar, before any new Astral search output.
variants={'full':(True,True),'counts':(True,False),'q':(False,True)}
pilot_names={'full':'expanded_control','counts':'expanded_counts','q':'expanded_q'}
winner=max(variants,key=lambda name:report[pilot_names[name]]['empirical_accession_curve']['0.01']['present'])
alt=json.loads((p.BASE/'joint_gap_alternatives_v1/report.json').read_text())['expanded_alternatives']['empirical_accession_curve']['0.01']['present']
choices=[('no_cell_scores',False,True,True,False),('full',True,True,True,False)]
if winner!='full':choices.append((winner,True,*variants[winner],False))
if alt>report['expanded_control']['empirical_accession_curve']['0.01']['present']:
    choices.append(('alternatives',True,True,True,True))
selection=dict(pilot_winner=winner,comparisons=choices,selection_dataset='IonStar only')
selection_file=p.OUT/'selection.json'
if selection_file.exists():assert json.loads(selection_file.read_text())==json.loads(json.dumps(selection))
else:selection_file.write_text(json.dumps(selection,indent=2)+'\n')
for name,enabled,counts,quality,alternatives in choices:
    cfg=p.OUT/(name+'.toml')
    cfg.write_text((p.BASE/'uno.toml').read_text()+'\n[scoring]\n'+
        '\n'.join(f'{k} = {str(v).lower()}' for k,v in dict(matrix_joint_cell_quality=enabled,
        matrix_joint_cell_counts=counts,matrix_joint_cell_q=quality,matrix_joint_alternative_groups=alternatives).items())+'\n')
    p.run(name,['--matrix-long',aligned,'--matrix-joint','--config',cfg])
c=json.loads((p.OLD/'frozen.json').read_text())['fasta_counts'];scale=c['HUMAN']/c['YEAST']
results={}
for name,*_ in choices:
    scores={}
    for row in load_proteins(p.OUT/name/'proteins.tsv'):
        if row['species']=='AMBIGUOUS':continue
        for a in row['members']:scores[a]=max(scores.get(a,0),row['score'])
    values=curve(scores,scale)
    for v in values.values():v.pop('accessions')
    results[name]=values
qs={}
for run in frozen['runs']:
    for row in load_proteins(p.OLD/'searches'/('single_'+run)/'seed0/proteins.tsv'):
        if row['species']=='AMBIGUOUS':continue
        for a in row['members']:qs[a]=min(qs.get(a,1),row['q'])
values=curve({a:-q for a,q in qs.items()},scale)
for v in values.values():v.pop('accessions')
results['separate_union']=values
(p.OUT/'report.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2),flush=True)
