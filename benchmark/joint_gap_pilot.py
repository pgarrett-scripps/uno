#!/usr/bin/env python3
"""Predefined count/q ablations and selected-candidate diagnostics; no ID pooling."""
import json
import shutil
from pathlib import Path
import joint_matrix_pilot as p
p.OUT = p.BASE / 'joint_gap_v1'
p.OUT.mkdir(exist_ok=True)
if not (p.OUT/'uno').exists():
    shutil.copy2(p.ROOT.parent/'target/release/uno',p.OUT/'uno')
    (p.OUT/'source-manifest.json').write_text(json.dumps({str(f.relative_to(p.ROOT)):p.sha(f)
        for f in sorted((p.ROOT/'crates').rglob('*.rs'))},indent=2)+'\n')
for name,bundle,counts,quality in [
    ('expanded_control','aligned5-replicated',True,True),
    ('expanded_counts','aligned5-replicated',True,False),
    ('expanded_q','aligned5-replicated',False,True),
    ('standard_counts','aligned5',True,False),
    ('standard_q','aligned5',False,True),
]:
    directory=p.OUT/name;directory.mkdir(exist_ok=True)
    cfg=p.OUT/(name+'.toml')
    cfg.write_text((p.BASE/'uno.toml').read_text()+ '\n[scoring]\nmatrix_joint_cell_quality = true\n'+
        f'matrix_joint_cell_counts = {str(counts).lower()}\nmatrix_joint_cell_q = {str(quality).lower()}\n'+
        f'diagnostic_peptides_path = "{directory / "diagnostic.tsv"}"\n')
    p.run(name,['--matrix-long',p.BASE/bundle,'--matrix-joint','--config',cfg])
    assert (directory/'diagnostic.tsv').stat().st_size>0
    if name=='expanded_control':
        parity={f:p.sha(directory/f)==p.sha(p.BASE/'joint_cell_quality_v1/expanded_cell_quality'/f)
                for f in ['proteins.tsv','lfq.tsv','abundance.tsv']}
        (p.OUT/'control-parity.json').write_text(json.dumps(parity,indent=2)+'\n')
        assert all(parity.values()),parity
    p.report()
(p.OUT/'diagnostic-hashes.json').write_text(json.dumps({str(f):p.sha(f) for f in p.OUT.glob('*/diagnostic.tsv')},indent=2)+'\n')
