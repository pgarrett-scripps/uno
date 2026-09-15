#!/usr/bin/env python3
"""Two-cohort consensus-row identification and known-ratio quantities."""
import json
import numpy as np
import style
import matplotlib.pyplot as plt
from _assets import record


def main():
    d = json.loads((style.FIGDATA/'koth_matrix_study.json').read_text())
    assert d['status'] == 'complete' and not d['timing']['usable_for_speed_claim']
    fig, axes = plt.subplots(1, 2, figsize=(style.W2, 2.9))
    ax = axes[0]
    for i, (dataset, cohort) in enumerate([('astral','100ng'),('ionstar','20')]):
        seeds = d['datasets'][dataset]['cohorts'][cohort]['seeds'].values()
        for j, (arm, label, color) in enumerate([('union','Separate union',style.C['accent2']),
                                               ('consensus','Consensus rows',style.C['uno'])]):
            vals = [s[arm]['0.01']['present'] for s in seeds]
            style.spread_marker(ax, i+(j-.5)*.25, vals, color, label=label if i == 0 else None,
                                half_width=.075, jitter=.025, ms=3)
    ax.set_xticks([0,1], ['Astral 100 ng\n3 runs','IonStar\n20 runs'])
    ax.set_ylim(bottom=0)
    ax.set_ylabel('Accessions at ≤1% entrapment FDP')
    ax.legend(frameon=False, fontsize=7, loc='lower right')
    ax = axes[1]
    seeds = d['datasets']['ionstar']['cohorts']['20']['seeds'].values()
    expected = np.array([1,1.5,2,2.5,3])
    for sp, label, color in [('ECOLI','E. coli',style.C['uno']),('HUMAN','Human',style.C['accent2'])]:
        vals = np.array([[2**v['median_log2_ratio'] if v['median_log2_ratio'] is not None else np.nan
                          for v in s['consensus']['quant']['species'][sp]['condition_ratio_curve']] for s in seeds])
        mid = np.nanmedian(vals, axis=0)
        ax.errorbar(expected, mid, yerr=[mid-np.nanmin(vals,axis=0),np.nanmax(vals,axis=0)-mid],
                    color=color, marker='o', ms=3, capsize=2, label=label, lw=1)
    ax.plot(expected, expected, ':', color=style.C['neutral'], lw=.8, label='Expected E. coli')
    ax.axhline(1, color=style.C['accent2'], ls=':', lw=.7)
    ax.set_xlabel('Expected E. coli ratio versus A')
    ax.set_ylabel('Observed protein-group ratio versus A')
    ax.legend(frameon=False, fontsize=7)
    for label, ax in zip('ab',axes):
        style.panel_label(ax,label)
    fig.tight_layout(w_pad=2)
    for ext, suffix in [('png',''),('pdf','-pdf')]:
        path='figures/figS7_koth_matrix.'+ext
        fig.savefig(style.PAPER/path,dpi=300,metadata={'Software':None} if ext=='png' else {'CreationDate':None,'ModDate':None})
        record('fig.figs7-koth-matrix'+suffix,path,kind='figure',inputs=['figdata/koth_matrix_study.json'],
               desc='KOTH 0.7.0: individual-file union versus consensus-row identification and IonStar known-ratio recovery; three seeds; no timing claim')
    plt.close(fig)


if __name__ == '__main__':
    main()
