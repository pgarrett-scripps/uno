#!/usr/bin/env python3
"""Fresh same-binary individual controls for the Astral cross-dataset result."""
import json
from pathlib import Path
import joint_matrix_pilot as p
from matrix_report import load_proteins
from joint_evidence_audit import curve
p.OUT=p.BASE/'joint_gap_v1/astral_validation';p.OLD=p.OLD/'astral'
runs=json.loads((p.OUT/'frozen.json').read_text())['runs'];qs={};parity={}
for i,run in enumerate(runs):
    name=f'single_{i:03}'
    p.run(name,['--features',p.OLD/'koth'/run/'features.tsv','--config',p.BASE/'uno.toml'])
    parity[run]=p.sha(p.OUT/name/'proteins.tsv')==p.sha(p.OLD/'searches'/('single_'+run)/'seed0/proteins.tsv')
    for row in load_proteins(p.OUT/name/'proteins.tsv'):
        if row['species']=='AMBIGUOUS':continue
        for a in row['members']:qs[a]=min(qs.get(a,1),row['q'])
c=json.loads((p.OLD/'frozen.json').read_text())['fasta_counts']
values=curve({a:-q for a,q in qs.items()},c['HUMAN']/c['YEAST'])
for v in values.values():v.pop('accessions')
report=dict(curve=values,protein_parity_with_frozen=parity,same_executable=p.sha(p.OUT/'uno'),
    same_config=p.sha(p.BASE/'uno.toml'),runs=runs)
(p.OUT/'fresh-controls.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2),flush=True)
