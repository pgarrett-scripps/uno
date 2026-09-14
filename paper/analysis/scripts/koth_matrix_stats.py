"""Traceable prose statistics for the audited KOTH LFQ update."""
import json
from statistics import median
from pathlib import Path
def add(st):
    d=json.loads((Path(__file__).resolve().parents[2]/'figdata/koth_matrix_study.json').read_text());assert d['status']=='complete'
    def stat(key,value,desc,fmt='.0f',bounds=(0,1000000)):
        st.add('matrix.'+key,value,fmt=fmt,desc=desc,between=bounds)
    for dataset,cohort,prefix in [('astral','100ng','astral'),('astral','200pg','astral_low'),('ionstar','20','ionstar')]:
        cc=d['datasets'][dataset]['cohorts'][cohort]['seeds']
        for arm in ['union','joint','joint_cell']:
            v=[s[arm]['0.01']['present'] for s in cc.values()]
            for suffix,fn in [('median',median),('min',min),('max',max)]:
                stat(prefix+'_'+arm+'_'+suffix,fn(v),f'{dataset} {cohort} {arm}: {suffix} present accessions at empirical FDP ≤1%, three seeds')
        gains=[100*(s['joint_cell']['0.01']['present']/s['union']['0.01']['present']-1) for s in cc.values()]
        stat(prefix+'_gain_percent',median(gains),f'{dataset} {cohort}: median paired relative gain, joint + cells versus independent union, three seeds','.1f',(-100,10000))
    c=d['datasets']['ionstar']['cohorts']['20']['cells']['totals']
    stat('positive_cells_millions',c['positive']/1e6,'All positive target IonStar 20-run extraction cells, millions','.2f')
    stat('confident_positive_percent',100*c['confident']/c['positive'],'Positive IonStar target cells passing extraction q ≤0.01, percent','.1f',(0,100))
    stat('search_speedup',d['timing']['speedup'],'Dedicated serialized Astral 100 ng search-only wall-time ratio, warm RT cache','.1f',(0,1000))
    for arm in ['joint','joint_cell']:
        ratios=[s[arm]['quant_common']['species']['ECOLI']['median_log2_E_over_A'] for s in d['datasets']['ionstar']['cohorts']['20']['seeds'].values()]
        assert all(r is not None for r in ratios)
        stat('ecoli_ratio_'+arm,2**median(ratios),f'Median common-group E. coli E/A ratio, {arm}, three seeds','.2f',(0,10))
