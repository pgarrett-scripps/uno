#!/usr/bin/env python3
"""Single fixed test of coherent alternative rows during protein competition."""
import json,shutil
import joint_matrix_pilot as p
p.OUT=p.BASE/'joint_gap_alternatives_v1';p.OUT.mkdir(exist_ok=True)
if not (p.OUT/'uno').exists():
    shutil.copy2(p.ROOT.parent/'target/release/uno',p.OUT/'uno')
    (p.OUT/'source-manifest.json').write_text(json.dumps({str(f.relative_to(p.ROOT)):p.sha(f)
        for f in sorted((p.ROOT/'crates').rglob('*.rs'))},indent=2)+'\n')
cfg=p.OUT/'alternatives.toml'
cfg.write_text((p.BASE/'uno.toml').read_text()+'\n[scoring]\nmatrix_joint_cell_quality = true\nmatrix_joint_alternative_groups = true\n')
p.run('expanded_alternatives',['--matrix-long',p.BASE/'aligned5-replicated','--matrix-joint','--config',cfg])
p.report()
