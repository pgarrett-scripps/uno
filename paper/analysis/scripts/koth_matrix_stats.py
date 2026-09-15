"""Traceable prose statistics for the KOTH 0.7.0 consensus-row study."""
import json
from statistics import median
from pathlib import Path


def add(st):
    d=json.loads((Path(__file__).resolve().parents[2]/'figdata/koth_matrix_study.json').read_text())
    assert d['status']=='complete' and d['arms']==['separate','consensus']
    def stat(key,value,desc,fmt='.0f',bounds=(0,1000000)):
        st.add('matrix.'+key,value,fmt=fmt,desc=desc,between=bounds)
    for dataset,cohort in [('astral','100ng'),('ionstar','20')]:
        seeds=d['datasets'][dataset]['cohorts'][cohort]['seeds'].values()
        assert all(s['union']['0.01']['present'] > s['consensus']['0.01']['present'] for s in seeds), 'Review the claim of lower consensus depth in both cohorts'
        for arm in ['union','consensus']:
            vals=[s[arm]['0.01']['present'] for s in seeds]
            for suffix,fn in [('median',median),('min',min),('max',max)]:
                stat(dataset+'_'+arm+'_'+suffix,fn(vals),f'KOTH 0.7.0 {dataset} {cohort} {arm}: {suffix} present accessions at empirical FDP ≤1%, three seeds')
        stat(dataset+'_groups',d['datasets'][dataset]['cohorts'][cohort]['cells']['groups'],f'KOTH 0.7.0 {dataset}: consensus feature rows')
    c=d['datasets']['ionstar']['cohorts']['20']['cells']['totals']
    stat('positive_cells_millions',c['positive']/1e6,'KOTH 0.7.0 positive IonStar target extraction cells, millions','.2f')
    stat('confident_positive_percent',100*c['confident']/c['positive'],'KOTH 0.7.0 positive IonStar target cells with extraction q ≤0.01, percent','.1f',(0,100))
    seeds=d['datasets']['ionstar']['cohorts']['20']['seeds'].values()
    for sp in ['HUMAN','ECOLI']:
        vals=[s['consensus']['quant']['species'][sp]['median_log2_E_over_A'] for s in seeds]
        assert all(v is not None for v in vals)
        stat(sp.lower()+'_ratio_consensus',2**median(vals),f'KOTH 0.7.0 consensus-row {sp} median E/A protein-group ratio across three seeds','.2f',(0,10))
