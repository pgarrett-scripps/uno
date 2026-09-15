#!/usr/bin/env python3
"""Compare independent quant-associated sequences with selected matrix candidates."""
import csv,json
from collections import Counter,defaultdict
from joint_evidence_audit import BASE
OUT=BASE/'joint_gap_v1'
a=json.loads((OUT/'retention-audit.json').read_text())
needed={r['peptide'] for p in a['per_protein'].values() for r in p['assignments']}
diag={}
with (OUT/'expanded_control/diagnostic.tsv').open() as f:
    for r in csv.DictReader(f,delimiter='\t'):
        if r['decoy']=='0' and r['peptide'] in needed:diag[r['peptide']]=r
anchors={(int(r['group_anchor_run']),int(r['group_anchor_observation'])) for r in diag.values()}
cids={}
with (BASE/'aligned5-replicated/feature_observations.tsv').open() as f:
    for r in csv.DictReader(f,delimiter='\t'):
        key=int(r['run_id']),int(r['observation_id'])
        if key in anchors:cids[key]=int(r['consensus_id'])
counts=Counter();rows={};sequence_counts=Counter()
for acc,p in a['per_protein'].items():
    seqs=defaultdict(set)
    for r in p['assignments']:
        if r['consensus'] is not None:seqs[r['peptide']].add(r['consensus'])
    local=Counter();details=[]
    for seq,groups in seqs.items():
        r=diag.get(seq)
        if r is None:category='not_in_final_matrix_candidate_set';cid=None
        else:
            cid=cids[int(r['group_anchor_run']),int(r['group_anchor_observation'])]
            same=cid in groups
            included=r['in_ensemble']=='1'
            category=('same_group' if same else 'different_group')+('_in_ensemble' if included else '_below_ensemble_rank')
            assert r['ml_selected']=='1','deviation result selected; cannot attribute to ML ensemble'
        local[category]+=1;sequence_counts[category]+=1
        details.append(dict(peptide=seq,category=category,selected_group=cid,
                            source_groups=sorted(groups),rank_fraction=float(r['rank_fraction']) if r else None))
    included=local['same_group_in_ensemble']+local['different_group_in_ensemble']
    if included>=2:category='at_least_two_candidates_enter_protein_inference'
    elif included:category='one_candidate_enters_protein_inference'
    else:category='no_retained_quant_sequence_enters_protein_inference'
    counts[category]+=1
    rows[acc]=dict(category=category,counts=dict(local),details=details)
result=dict(protein_counts=dict(counts),sequence_association_counts=dict(sequence_counts),per_protein=rows,
    limitations=a['limitations']+['Sequence-association counts count a shared sequence once per associated protein. Availability does not establish correct identity.'])
(OUT/'stage-audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='per_protein'}),flush=True)
