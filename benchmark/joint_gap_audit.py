#!/usr/bin/env python3
"""Trace independent quant assignments through original KOTH observations.
These are competition-resolved quant assignments, not validated peptide IDs.
Independent search output is used only to diagnose losses, never to score a matrix.
"""
import csv,json,itertools
from collections import Counter,defaultdict
from pathlib import Path
from joint_evidence_audit import BASE,OLD,curve,species,sha
from matrix_report import load_proteins
OUT=BASE/'joint_gap_v1';OUT.mkdir(exist_ok=True)
c=json.loads((OLD/'frozen.json').read_text())['fasta_counts'];scale=(c['HUMAN']+c['ECOLI'])/c['YEAST']
qs={};eligible=[]
for run in range(5):
    current={}
    for r in load_proteins(BASE/f'E_batch/run{run:03}/proteins.tsv'):
        if r['species']=='AMBIGUOUS':continue
        for a in r['members']:
            qs[a]=min(qs.get(a,1),r['q']);current[a]=r['q']
    eligible.append(current)
u=curve({a:-q for a,q in qs.items()},scale)['0.01'];union=set(u.pop('accessions'));qcut=-u['score_threshold']
scores={}
matrix_rows=load_proteins(BASE/'joint_cell_quality_v1/expanded_cell_quality/proteins.tsv')
for r in matrix_rows:
    if r['species']=='AMBIGUOUS':continue
    for a in r['members']:scores[a]=max(scores.get(a,0),r['score'])
m=curve(scores,scale)['0.01'];matrix=set(m.pop('accessions'));lost=union-matrix;gained=matrix-union
print('overlap',len(union&matrix),'union only',len(lost),'matrix only',len(gained),flush=True)
# Retain all confident-union assignments for a grouping-conflict diagnostic.
wanted=[defaultdict(list) for _ in range(5)]
for run in range(5):
    with (BASE/f'E_batch/run{run:03}/lfq.tsv').open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            if row['decoy']!='0':continue
            prots=set(row['proteins'].split(';'))
            accepted={a for a in prots if eligible[run].get(a,1)<=qcut and a in union}
            if not accepted:continue
            key=(row['charge'],row['mass'],row['rt'])
            wanted[run][key].append(dict(peptide=row['peptide'],proteins=sorted(accepted),
                unique=len(prots)==1))
records=[];collisions=0;matched=0;counts=Counter();group_sequences=defaultdict(set)
with (BASE/'aligned5-replicated/feature_observations.tsv').open() as f:
    for run,rows in itertools.groupby(csv.DictReader(f,delimiter='\t'),key=lambda r:int(r['run_id'])):
        found={}
        for o in rows:
            key=(o['feature_charge'],format(float(o['feature_mass']),'.5f'),format(float(o['feature_rt']),'.4f'))
            if key not in wanted[run]:continue
            if key in found:found[key]=None;collisions+=1
            else:found[key]=o
        for key,assignments in wanted[run].items():
            if key not in found:counts['unmatched_quant_coordinates']+=len(assignments);continue
            o=found[key]
            if o is None:counts['ambiguous_quant_coordinates']+=len(assignments);continue
            for a in assignments:
                matched+=1
                cid=int(o['consensus_id']) if o['consensus_id'] else None
                r=dict(run=run,observation=int(o['observation_id']),consensus=cid,
                       reason=o['exclusion_reason'],quality=float(o['feature_combined_score']),**a)
                records.append(r)
                if cid is not None and a['unique']:group_sequences[cid].add(a['peptide'])
per={}
for a in sorted(lost):
    rs=[r for r in records if a in r['proteins']]
    seq={r['peptide'] for r in rs};retained={r['peptide'] for r in rs if r['consensus'] is not None}
    groups={r['consensus'] for r in rs if r['consensus'] is not None}
    conflict={g for g in groups if len(group_sequences[g])>1}
    per[a]=dict(quant_peptides=len(seq),retained_quant_peptides=len(retained),groups=len(groups),
        unique_quant_peptides=len({r['peptide'] for r in rs if r['unique']}),
        risky_groups=len(conflict),excluded_reasons=dict(Counter(r['reason'] for r in rs if r['consensus'] is None)),
        matrix_scored=a in scores,matrix_score=scores.get(a),assignments=rs)
    counts['proteins_'+('no_quant_assignment' if not seq else 'all_quant_observations_excluded' if not retained else
                       'some_quant_peptides_excluded' if len(retained)<len(seq) else 'all_quant_peptides_retained')]+=1
    if len(retained)>=2:counts['proteins_with_at_least_two_retained_quant_peptides']+=1
    if conflict:counts['proteins_with_potential_group_conflict']+=1
def retention_summary(accessions):
    import statistics
    by=defaultdict(lambda: [set(),set()])
    for r in records:
        for acc in r['proteins']:
            if acc not in accessions:continue
            by[acc][0].add(r['peptide'])
            if r['consensus'] is not None:by[acc][1].add(r['peptide'])
    return dict(proteins=len(accessions),with_quant_assignments=len(by),
        median_quant_peptides=statistics.median(len(v[0]) for v in by.values()),
        median_retained_quant_peptides=statistics.median(len(v[1]) for v in by.values()),
        median_retained_fraction=statistics.median(len(v[1])/len(v[0]) for v in by.values()))
retention_comparison=dict(union_only=retention_summary(lost),shared=retention_summary(union&matrix))
result=dict(retention_comparison=retention_comparison,scope='Pilot five-run IonStar at empirical entrapment FDP<=1%; accession-level comparison.',
    separate_union=u,matrix=m,overlap=len(union&matrix),union_only=len(lost),matrix_only=len(gained),
    lost=sorted(lost),gained=sorted(gained),counts=dict(counts),per_protein=per,
    limitations=['Quant assignments are not accepted peptide IDs or the exact surviving peptides of every inference pass.',
                 'Distinct peptide assignments in a group flag possible mixing; they do not prove an incorrect alignment.'],
    hashes={str(f):sha(f) for f in [Path(__file__),BASE/'aligned5-replicated/matrix_manifest.json',
        BASE/'joint_cell_quality_v1/expanded_cell_quality/proteins.tsv']})
(OUT/'retention-audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(counts),flush=True)
