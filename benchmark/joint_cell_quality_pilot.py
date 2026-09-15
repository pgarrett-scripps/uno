#!/usr/bin/env python3
"""Fixed seed0 pilot of optional MBR and extraction-q row covariates."""
import json
import shutil
from pathlib import Path
import joint_matrix_pilot as pilot

pilot.OUT = pilot.BASE / 'joint_cell_quality_v1'
pilot.OUT.mkdir(exist_ok=True)
binary = pilot.OUT / 'uno'
if not binary.exists():
    shutil.copy2(pilot.ROOT.parent / 'target/release/uno', binary)
    (pilot.OUT / 'source-manifest.json').write_text(json.dumps({str(p.relative_to(pilot.ROOT)): pilot.sha(p)
        for p in sorted((pilot.ROOT / 'crates').rglob('*.rs'))}, indent=2) + '\n')
config = pilot.BASE / 'uno.toml'
enabled = pilot.OUT / 'cell-quality.toml'
enabled.write_text(config.read_text() + '\n[scoring]\nmatrix_joint_cell_quality = true\n')
# Predefined comparison: baseline parity, added descriptors, expanded retention.
for name, bundle, cfg in [('baseline', 'aligned5', config),
                          ('cell_quality', 'aligned5', enabled),
                          ('expanded_cell_quality', 'aligned5-replicated', enabled)]:
    pilot.run(name, ['--matrix-long', pilot.BASE / bundle, '--matrix-joint', '--config', cfg])
    if name == 'baseline':
        parity = {f: pilot.sha(pilot.OUT / name / f) == pilot.sha(pilot.BASE / 'joint_search_v3/joint' / f)
                  for f in ['proteins.tsv', 'lfq.tsv', 'abundance.tsv']}
        (pilot.OUT / 'baseline-parity.json').write_text(json.dumps(parity, indent=2) + '\n')
        assert all(parity.values()), parity
pilot.report()
result = json.loads((pilot.OUT / 'report.json').read_text())
prior = json.loads((pilot.BASE / 'joint_search_v3/report.json').read_text())
summary = dict(
    scope='Five-run IonStar, seed0; present-species accessions at empirical entrapment FDP <=1%; development pilot.',
    arms={name: dict(**result[name]['empirical_accession_curve']['0.01'], seconds=result[name]['seconds'],
                    internal_1pct=result[name]['internal_1pct_accessions'])
          for name in ['baseline', 'cell_quality', 'expanded_cell_quality']},
    expanded_without_cell_scores=prior['replicated_joint']['empirical_accession_curve']['0.01'],
    separate_union=result['separate_union_reference'],
    baseline_parity=json.loads((pilot.OUT / 'baseline-parity.json').read_text()),
    tests=dict(uno_core=63, uno_cli=7),
    config_option='scoring.matrix_joint_cell_quality = true',
    option_default=False,
    binary_sha256=pilot.sha(binary),
    driver_sha256=pilot.sha(Path(__file__)),
)
(pilot.ROOT / 'benchmark/joint_cell_quality_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
