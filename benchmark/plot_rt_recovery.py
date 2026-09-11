#!/usr/bin/env python3
"""Plot every prospective low-input RT-recovery pair, including zero results."""
import hashlib
import json
import os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR', '/tmp/uno-rt-recovery-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'benchmark/results/rt_recovery_investigation'


def main():
    source = BASE / 'screen-report.json'
    report = json.loads(source.read_text())
    rows = report['panels']['lowinput']['individual_results']
    datasets = ['astral_1ng_rep1', 'astral_1ng_rep2', 'astral_1ng_rep4',
                'astral_200pg_rep1', 'astral_200pg_rep2', 'astral_200pg_rep3']
    assert len(rows) == 24
    assert {(r['dataset'], r['seed'], r['variant']) for r in rows} == {
        (ds, seed, variant) for ds in datasets for seed in (9, 10) for variant in ('baseline', 'candidate')}
    fig, axes = plt.subplots(2, 3, figsize=(11, 6.5), constrained_layout=True)
    for ax, dataset in zip(axes.flat, datasets):
        for seed, offset, color in [(9, -.06, '#2878b5'), (10, .06, '#d95f02')]:
            counts = [next(r['internal_q1_targets'] for r in rows if r['dataset'] == dataset and
                           r['seed'] == seed and r['variant'] == variant) for variant in ('baseline', 'candidate')]
            ax.plot([offset, 1 + offset], counts, 'o-', color=color, label=f'Seed {seed}', linewidth=1.4)
        ax.set_title(dataset.replace('astral_', '').replace('_', ' '))
        ax.set_xticks([0, 1], ['V3', 'RT recovery'])
        ax.set_xlim(-.3, 1.3)
        ax.set_ylim(bottom=-20)
        ax.set_ylabel('Reported targets (q ≤ 1%)')
        ax.grid(axis='y', alpha=.2)
        ax.spines[['top', 'right']].set_visible(False)
    axes[0, 0].legend()
    fig.suptitle('Broader RT only after V3 failure: new seed draws', fontsize=15)
    for suffix in ('png', 'pdf'):
        fig.savefig(BASE / f'new-draw-counts.{suffix}', dpi=180, bbox_inches='tight')
    plt.close(fig)
    provenance = dict(report_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                      script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      caption='All 12 paired cases are shown. Lines connect matching parent seeds; zero is an empty target list. Two correlated algorithmic draws per known acquisition file do not establish general biological reproducibility. Panel scales differ.')
    (BASE / 'figure-provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')


if __name__ == '__main__':
    main()
