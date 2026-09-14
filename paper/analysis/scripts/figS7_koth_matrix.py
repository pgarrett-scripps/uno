#!/usr/bin/env python3
"""Astral-first release-pinned matrix study; matched error thresholds, all seeds."""
import json
import numpy as np
import style
import matplotlib.pyplot as plt
from _assets import record
ARMS=[('union','Separate union',style.C['accent2']),('joint','Joint',style.C['uno']),('joint_cell','Joint + cells',style.C['accent'])]
def main():
    d=json.loads((style.FIGDATA/'koth_matrix_study.json').read_text());assert d['status']=='complete'
    fig,axes=plt.subplots(2,2,figsize=(style.W2,5.4))
    for ax,dataset,title in [(axes[0,0],'astral','Astral'),(axes[0,1],'ionstar','IonStar')]:
        cohorts=d['datasets'][dataset]['cohorts']
        for j,(arm,label,color) in enumerate(ARMS):
            mid=[];xx=[]
            for i,(name,c) in enumerate(cohorts.items()):
                vals=[s[arm]['0.01']['present'] for s in c['seeds'].values()];x=i+(j-1)*.19
                style.spread_marker(ax,x,vals,color,label=label if i==0 else None,half_width=.065,jitter=.025,ms=3);mid.append(np.median(vals));xx.append(x)
            ax.plot(xx[:-1] if dataset=='astral' else xx,mid[:-1] if dataset=='astral' else mid,color=color,lw=.7,alpha=.5)
        labels={'200pg':'200 pg','1ng':'1 ng','10ng':'10 ng','100ng':'100 ng','all':'All 12'}
        ax.set_xticks(np.arange(len(cohorts)),[labels.get(n,n) for n in cohorts]);ax.set_ylim(bottom=0)
        if dataset=='astral':ax.axvline(3.5,color=style.C['light'],lw=.8,ls=':')
        ax.set_title(title,fontsize=9);ax.set_xlabel('Input level / pooled files' if dataset=='astral' else 'Runs in cohort')
        ax.set_ylabel('Accessions at ≤1% entrapment FDP')
    axes[0,0].legend(frameon=False,fontsize=7,loc='upper left')
    t=d['timing'];ax=axes[1,0]
    ax.bar([0,1],[t[k]['wall_seconds'] for k in ['single','joint_cell']],color=[style.C['accent2'],style.C['accent']],width=.55)
    ax.set_xticks([0,1],['Three separate','One joint + cells']);ax.set_ylabel('Serialized search wall time (s)');ax.set_title('Astral 100 ng; warm RT cache',fontsize=9)
    ax=axes[1,1];truth=np.log2([1,1.5,2,2.5,3])
    for arm,label,color in ARMS[1:]:
        qq=[s[arm]['quant_common'] for s in d['datasets']['ionstar']['cohorts']['20']['seeds'].values()]
        values=np.array([[x['median_log2_ratio'] for x in q['species']['ECOLI']['condition_ratio_curve']] for q in qq if q],dtype=float)
        if not len(values):continue
        mid=np.nanmedian(values,axis=0);lo=np.nanmin(values,axis=0);hi=np.nanmax(values,axis=0)
        ax.errorbar(truth,mid,yerr=[mid-lo,hi-mid],color=color,marker='o',ms=3,capsize=2,label=label,lw=1)
    ax.plot(truth,truth,color=style.C['neutral'],ls=':',lw=.8,label='Expected')
    ax.set_xlabel('Expected log2 ratio versus A');ax.set_ylabel('Observed E. coli log2 ratio');ax.set_title('IonStar: common protein groups',fontsize=9);ax.legend(frameon=False,fontsize=7)
    for label,ax in zip('abcd',axes.flat):style.panel_label(ax,label)
    fig.tight_layout(h_pad=2,w_pad=2)
    for ext,suffix in [('png',''),('pdf','-pdf')]:
        name='figures/figS7_koth_matrix.'+ext
        fig.savefig(style.PAPER/name,dpi=300,metadata={'Software':None} if ext=='png' else {'CreationDate':None,'ModDate':None})
        record('fig.figs7-koth-matrix'+suffix,name,kind='figure',inputs=['figdata/koth_matrix_study.json'],desc='Astral and IonStar matched-entrapment accession depth, warm-cache serialized timing, and known-ratio recovery; three seeds')
    plt.close(fig)
if __name__=='__main__':main()
