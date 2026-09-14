#!/usr/bin/env python3
"""Assess fresh Astral joint searches against deposited DIA-NN evidence."""
from __future__ import annotations
import csv
import json
from pathlib import Path
import re
import pyarrow.parquet as pq
import numpy as np
import matrix_report as report
import matrix_study as study
from matrix_astral import ASTUDY, ASTRAL


def accession(value):
    p=value.split('|')
    return p[1] if len(p)>=3 else value


def deposited_evidence():
    out=ASTUDY/'msms_evidence.json'
    source=ASTRAL/'external/CW0710_DN_all_/report.parquet'
    source_hash=study.sha(source)
    if out.exists():
        cached=json.loads(out.read_text())
        if cached['source_sha256']==source_hash:return cached
        raise ValueError('Deposited MS/MS evidence changed')
    runs={r['run'] for r in csv.DictReader(open(ASTUDY/'inputs.tsv'),delimiter='\t')}
    data=pq.read_table(source,columns=['Run','Decoy','Q.Value','PG.Q.Value','Protein.Group',
                                      'Stripped.Sequence','Modified.Sequence','Precursor.Charge'],
                       filters=[('Run','in',sorted(runs)),('Decoy','=',0),
                                ('Q.Value','<=',.01),('PG.Q.Value','<=',.01)])
    evidence={r:{'proteins':set(),'precursors':set(),'filtered_report_rows':0} for r in runs}
    for r in data.to_pylist():
        e=evidence[r['Run']];e['filtered_report_rows']+=1
        e['proteins'].update(r['Protein.Group'].split(';'))
        modified=str(r['Modified.Sequence']).replace('(UniMod:4)','').replace('UniMod:4','')
        if any(c in modified for c in '()[]'):continue
        e['precursors'].add(str(r['Stripped.Sequence'])+'/'+str(int(r['Precursor.Charge'])))
    for e in evidence.values():
        e['proteins']=sorted(e['proteins']);e['precursors']=sorted(e['precursors'])
    cached={'source':str(source),'source_sha256':source_hash,
            'filter':'target rows, precursor Q.Value <= 0.01 and PG.Q.Value <= 0.01; precursor comparison excludes variable modifications',
            'runs':evidence}
    study.write_json(out,cached)
    return cached


def precursor_support(path, group, accepted, msms, cache):
    """Protein-supported assignments, not peptide-FDR-qualified identifications."""
    if group not in cache:
        directory=ASTUDY/'evidence'/group
        meta=json.loads((directory/'summary.json').read_text())
        p=directory/'cell_evidence.npz'
        if study.sha(p)!=meta['cell_evidence_sha256']:raise ValueError('Cell-evidence checksum mismatch')
        with np.load(p,allow_pickle=False) as z:
            cells={k:z[k] for k in z.files}
        lookup={}
        for i,(m,t,c) in enumerate(zip(cells['mass'],cells['rt'],cells['charge'])):
            key=(format(m,'.5f'),format(t,'.4f'),str(c))
            lookup[key]=i if key not in lookup else None
        cache[group]=(cells,lookup)
    cells,lookup=cache[group]
    runs=cells['run_names'].tolist()
    supported={r:set(msms[r]['precursors']) for r in runs}
    metrics={r:{'positive_assignment_rows':0,'confident_direct_rows':0,'confident_mbr_rows':0,
                'confident_precursors':set()} for r in runs}
    mapped=ambiguous=unmapped=0
    with open(path) as f:
        reader=csv.DictReader(f,delimiter='\t')
        if reader.fieldnames[6:]!=runs:raise ValueError('LFQ run order differs from evidence')
        for r in reader:
            if r['decoy']!='0':continue
            if not ({accession(a) for a in r['proteins'].split(';')}&accepted):continue
            key=(r['mass'],r['rt'],r['charge'])
            if key not in lookup:unmapped+=1;continue
            i=lookup[key]
            if i is None:ambiguous+=1;continue
            mapped+=1
            precursor=r['peptide']+'/'+r['charge']
            for j,run in enumerate(runs):
                if float(r[run])<=0:continue
                metrics[run]['positive_assignment_rows']+=1
                if not cells['confident'][i,j]:continue
                field='confident_direct_rows' if cells['direct'][i,j] else 'confident_mbr_rows'
                metrics[run][field]+=1
                metrics[run]['confident_precursors'].add(precursor)
    for run,m in metrics.items():
        precursors=m.pop('confident_precursors')
        m['unique_confident_precursors']=len(precursors)
        m['msms_supported_confident_precursors']=len(precursors&supported[run])
        m['msms_support_fraction']=len(precursors&supported[run])/len(precursors) if precursors else None
    return {'mapped_assignment_rows':mapped,'ambiguous_coordinate_rows':ambiguous,
            'unmatched_coordinate_rows':unmapped,'runs':metrics,
            'scope':'Assignments belonging to accepted human proteins; no peptide-level FDR claim. MS/MS corroboration is same-run exact sequence and charge.'}


def main():
    report.STUDY=ASTUDY
    frozen=json.loads((ASTUDY/'frozen.json').read_text())
    counts=frozen['fasta_counts'];scale=counts['HUMAN']/counts['YEAST']
    inputs=list(csv.DictReader(open(ASTUDY/'inputs.tsv'),delimiter='\t'))
    msms=deposited_evidence()['runs']
    tables={};summaries={};cell_cache={}
    for path in sorted((ASTUDY/'searches').glob('*/seed*/proteins.tsv')):
        name,seed=path.parent.parent.name,path.parent.name
        resource=report.resources(f'search_{name}_{seed}')
        if resource.get('exit_code')!=0:continue
        rows=report.load_proteins(path);tables[name,seed]=rows
        result=report.evaluate(rows,scale)
        accepted={accession(p) for r in rows if r['q']<=.01 and r['species']=='HUMAN' for p in r['members']}
        if name.startswith('single_'):
            runs=[name[7:]]
        else:
            group=name.rsplit('_',1)[1]
            runs=[r['run'] for r in inputs if group=='all' or r['condition']==group]
        msms_proteins={p for r in runs for p in msms[r]['proteins']}
        result['msms_protein_agreement']={'human_accessions':len(accepted),
            'supported_anywhere_in_same_cohort':len(accepted&msms_proteins),
            'fraction':len(accepted&msms_proteins)/len(accepted) if accepted else None}
        result['resources']=resource
        result['accession_level_1pct']=report.list_counts(report.accession_union([rows]),scale)
        result['accession_calibration_curve']=report.accession_curve([rows],scale)
        result['file_sha256']=study.sha(path)
        if not name.startswith('single_'):
            result['per_sample_support']=precursor_support(path.with_name('lfq.tsv'),group,accepted,msms,cell_cache)
        summaries[name+'/'+seed]=result
    unions={};comparisons={}
    for group in ['200pg','1ng','10ng','100ng','all']:
        runs=[r['run'] for r in inputs if group=='all' or r['condition']==group]
        for seed in ['seed0','seed1','seed2']:
            keys=[('single_'+r,seed) for r in runs]
            if all(k in tables for k in keys):
                union=report.accession_union([tables[k] for k in keys])
                accepted={accession(r['accession']) for r in union if r['species']=='HUMAN'}
                supported={p for r in runs for p in msms[r]['proteins']}
                unions[group+'/'+seed]=dict(report.list_counts(union,scale),
                    accession_calibration_curve=report.accession_curve([tables[k] for k in keys],scale),
                    msms_protein_agreement={'human_accessions':len(accepted),
                        'supported_anywhere_in_same_cohort':len(accepted&supported),
                        'fraction':len(accepted&supported)/len(accepted) if accepted else None})
    for (name,seed),rows in tables.items():
        if name.startswith('single_'):continue
        group=name.rsplit('_',1)[1]
        runs=[r['run'] for r in inputs if group=='all' or r['condition']==group]
        singles={r:tables.get(('single_'+r,seed)) for r in runs}
        if not all(r is not None for r in singles.values()):continue
        result=report.compare_cohort(rows,singles)
        added={accession(a) for a in result['joint_only_accessions']}
        supported={p for r in runs for p in msms[r]['proteins']}
        result['additional_msms_corroboration']={'human_accessions':len(added),
            'supported_anywhere_in_same_cohort':len(added&supported),
            'fraction':len(added&supported)/len(added) if added else None}
        comparisons[name+'/'+seed]=result
    study.write_json(ASTUDY/'summary.json',{'searches':summaries,'unions':unions,
        'joint_vs_separate':comparisons,
        'msms_source_sha256':study.sha(ASTUDY/'msms_evidence.json'),
        'interpretation':'Cohort-level corroboration only; pooled identification does not imply detection in every sample.'})
    print(f'Summarized {len(summaries)} completed Astral searches')


if __name__=='__main__':main()
