#!/usr/bin/env python3
"""Generate all matrix SI tables from the audited KOTH 0.5.0 three-seed study."""
import json,re
from statistics import median
import style
from _assets import record
SOURCE=style.FIGDATA/'koth_matrix_study.json'
ARMS=[('union','Separate union'),('joint','Joint'),('joint_cell','Joint + cells')]
def spread(values,digits=0,scale=1):
    v=[x*scale for x in values if x is not None]
    if not v:return 'NA'
    fmt=lambda x:f'{x:,.{digits}f}'
    return fmt(median(v)) if min(v)==max(v) else f'{fmt(median(v))} ({fmt(min(v))}–{fmt(max(v))})'
def table(name,headers,rows,desc):
    out=['// Generated; do not edit by hand.','#table(',f'  columns: {len(headers)},',
         '  align: ('+', '.join(['left']+['right']*(len(headers)-1))+'),',
         '  table.header(repeat: true, '+', '.join('['+h+']' for h in headers)+'),']
    out+=['  '+', '.join('['+str(v)+']' for v in row)+',' for row in rows];out.append(')')
    path='si/koth_matrix_'+name+'.typ';(style.PAPER/path).write_text('\n'.join(out)+'\n')
    record('tbl.koth-matrix-'+name,path,kind='table',inputs=['figdata/koth_matrix_study.json'],desc=desc)
def main():
    data=json.loads(SOURCE.read_text());assert data['status']=='complete' and data['publication_seeds']==[0,1,2]
    ds=data['datasets'];out=[]
    for n,c in ds['astral']['cohorts'].items():
        for arm,label in ARMS:
            vv=[s[arm] for s in c['seeds'].values()];points=[v['0.01'] for v in vv]
            out.append([n+' / '+label,spread([v['present'] for v in points]),spread([v['yeast'] for v in points]),
                spread([v['empirical_fdp'] for v in points],2,100),
                spread([v['corroboration']['supported']/v['corroboration']['accepted'] if v['corroboration']['accepted'] else None for v in vv],1,100)])
    table('astral',['Astral cohort / arm','Human accessions','Yeast','FDP (%)','MS/MS (%)'],out,'Astral accession counts at matched empirical 1% entrapment cutoffs; three-seed medians and ranges')
    out=[]
    for n,c in ds['ionstar']['cohorts'].items():
        for arm,label in ARMS:
            vv=[s[arm] for s in c['seeds'].values()]
            out.append([n+' / '+label,spread([v['0.01']['present'] for v in vv]),spread([v['0.01']['yeast'] for v in vv]),
                spread([v['0.01']['empirical_fdp'] for v in vv],2,100),spread([v['internal_1pct']['entrapment_fdp'] for v in vv],1,100)])
    table('identification',['IonStar runs / arm','Present accessions','Yeast','FDP (%)','Internal-list FDP (%)'],out,'IonStar matched-cutoff accession coverage and separately measured FDP of internal 1% lists')
    out=[]
    for dataset in ['astral','ionstar']:
        for n,c in ds[dataset]['cohorts'].items():
            for arm,label in ARMS[1:]:
                vv=[s[arm] for s in c['seeds'].values()]
                out.append([dataset.capitalize()+' '+n+' / '+label,spread([len(v['comparison']['added']) for v in vv]),
                    spread([len(v['comparison']['lost']) for v in vv]),spread([v['comparison']['jaccard'] for v in vv],1,100)])
    table('additional',['Cohort / arm','Joint only','Union only','Jaccard (%)'],out,'Added and lost accessions at matched empirical cutoffs; full three-seed ranges')
    out=[]
    for scope,key in [('All','quant_all'),('Common','quant_common')]:
        for arm,label in ARMS[1:]:
            for sp,labelsp in [('HUMAN','Human'),('ECOLI','_E. coli_')]:
                vals=[s[arm][key]['species'][sp] for s in ds['ionstar']['cohorts']['20']['seeds'].values() if s[arm][key]]
                out.append([scope+' / '+label+' / '+labelsp,spread([v['extreme_ratio_groups'] for v in vals]),
                    spread([v['median_abs_log2_ratio_error'] for v in vals],2),spread([v['median_replicate_cv'] for v in vals],1,100),
                    spread([v['positive_cell_fraction'] for v in vals],1,100)])
    table('quantification',['20-run comparison','Ratio groups','Abs. log2 error','CV (%)','Positive (%)'],out,'IonStar known-ratio measurements from internal-q accepted abundance outputs; shared normalization and no imputation')
    t=data['timing'];out=[]
    for arm,label in [('single','Three separate searches'),('joint_cell','One joint + cells search')]:
        r=t[arm];out.append([label,f"{r['wall_seconds']:.1f}",f"{r['user_seconds']+r['system_seconds']:.1f}",f"{r['peak_rss_kb']/1024**2:.2f}"])
    table('resources',['Astral 100 ng route','Wall (s)','CPU (s)','Peak RSS (GiB)'],out,'One dedicated serialized seed-zero repetition per route with warm shared RT cache; search-only times')
    out=[]
    for run,c in ds['astral']['cohorts']['all']['cells']['runs'].items():
        level=run.split('OldFAIMS_')[1].split('__')[0];rep=re.search(r'_rep(\d+)_',run).group(1)
        out.append([level+' / '+rep,f"{c.get('direct_positive',0):,}",f"{c.get('mbr_positive',0):,}",f"{c.get('direct_confident',0):,}",f"{c.get('mbr_confident',0):,}"])
    table('sample-support',['Input / replicate','Direct positive','MBR positive','Direct q ≤ 1%','MBR q ≤ 1%'],out,'All target extraction cells before UNO identification; positive and extraction-q-qualified cells do not imply identified peptides')
    p=data['provenance'];out=[]
    for binary,repo in [('uno','uno'),('koth_ff','koth_rust'),('koth_align','koth_rust')]:
        out.append([binary.replace('_','\\_'),p['source'][repo]['head'][:12],p['binary'][binary][:12]])
    out.append(['LightGBM runtime','4.6.0',data['source_supplement']['files']['uno']['crates/uno-lightgbm-sys/vendor/lib_lightgbm.so'][:12]])
    table('versions',['Program','Revision / version','Binary SHA-256 prefix'],out,'Frozen KOTH 0.5.0, UNO, and runtime library identities; complete source snapshots and hashes archived')
    table('routes',['Route','Measurements used','Searches per cohort'],[
        ['Separate union','Original per-run koth feature table','One per file'],
        ['Joint','Original observations in long bundle','One'],
        ['Joint + cells','Original observations + extraction summaries','One']],
        'Freshly rerun routes; ordinary --raw mode was not rerun as a separate matrix-study arm')
if __name__=='__main__':main()
