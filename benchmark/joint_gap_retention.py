#!/usr/bin/env python3
"""One predefined lower member-quality floor, keeping all other gates fixed."""
import json,os,shutil,subprocess,time
from pathlib import Path
import joint_matrix_pilot as p
p.OUT=p.BASE/'joint_gap_retention_v1';p.OUT.mkdir(exist_ok=True)
if not (p.OUT/'uno').exists():shutil.copy2(p.BASE/'joint_gap_alternatives_v1/uno',p.OUT/'uno')
if not (p.OUT/'source-manifest.json').exists():shutil.copy2(p.BASE/'joint_gap_alternatives_v1/source-manifest.json',p.OUT/'source-manifest.json')
if not (p.OUT/'koth_align').exists():shutil.copy2(p.BASE/'joint_search_v2/koth_align',p.OUT/'koth_align')
config=p.OUT/'align.toml'
config.write_text((p.BASE/'align5-replicated.toml').read_text().replace('min_group_size = 1','min_group_size = 1\nmin_member_combined_score = 0.3'))
aligned=p.OUT/'aligned';aligned.mkdir(exist_ok=True)
av=[str(p.OUT/'koth_align'),str(p.OLD/'batches/5'),'--config',str(config),'--output',str(aligned)]
signature=dict(argv=av,binary=p.sha(p.OUT/'koth_align'),config=p.sha(config),
               sources_manifest=p.sha(p.BASE/'aligned5-replicated/matrix_manifest.json'))
stamp=aligned/'job.json'
if stamp.exists():
    old=json.loads(stamp.read_text());assert old['signature']==signature and old['returncode']==0
    assert all(p.sha(Path(f))==h for f,h in old['outputs'].items())
else:
    started=time.time();print('Starting member-floor 0.3 alignment',flush=True)
    with (aligned/'run.log').open('w') as log:
        proc=subprocess.run(av,stdout=log,stderr=subprocess.STDOUT,env=dict(os.environ,RAYON_NUM_THREADS='8'))
    outputs={str(f):p.sha(f) for f in aligned.glob('*.tsv')}
    if (aligned/'matrix_manifest.json').exists():outputs[str(aligned/'matrix_manifest.json')]=p.sha(aligned/'matrix_manifest.json')
    stamp.write_text(json.dumps(dict(signature=signature,returncode=proc.returncode,seconds=time.time()-started,outputs=outputs),indent=2)+'\n')
    assert proc.returncode==0
    print('Finished alignment',round(time.time()-started,1),flush=True)
cfg=p.OUT/'search.toml';cfg.write_text((p.BASE/'uno.toml').read_text()+'\n[scoring]\nmatrix_joint_cell_quality = true\n')
p.run('member_floor_03',['--matrix-long',aligned,'--matrix-joint','--config',cfg])
p.report()
