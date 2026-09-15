#!/usr/bin/env python3
"""Compact results plus independent integer verification of the Astral win."""
import csv,json
from collections import defaultdict
from decimal import Decimal
from pathlib import Path
from joint_evidence_audit import BASE,OLD,ROOT,sha
G=BASE/'joint_gap_v1';A=G/'astral_validation'
read=lambda p:json.loads(p.read_text())
r=read(G/'report.json');ret=read(G/'retention-audit.json');stage=read(G/'stage-audit.json')
a=read(A/'report.json');fresh=read(A/'fresh-controls.json')
alt=read(BASE/'joint_gap_alternatives_v1/report.json')['expanded_alternatives']['empirical_accession_curve']['0.01']
loose=read(BASE/'joint_gap_retention_v1/report.json')['member_floor_03']['empirical_accession_curve']['0.01']
# Independent exact-integer calculation, avoiding floats in threshold comparisons.
def exact(files,union,counts):
    scores={}
    for f in files:
        for row in csv.DictReader(f.open(),delimiter='\t'):
            if row['decoy']!='0':continue
            members=row['group_members'].split(';')
            kinds={x.rsplit('_',1)[-1] for x in members}
            if len(kinds)!=1 or not kinds<={'HUMAN','YEAST'}:continue
            score=-int(Decimal(row['q_value'])*100000) if union else int(Decimal(row['binomial_score'])*1000)
            for x in members:scores[x]=max(scores.get(x,-10**10),score)
    blocks=defaultdict(list)
    for x,s in scores.items():blocks[s].append(x)
    n=y=best=0
    for s in sorted(blocks,reverse=True):
        for x in blocks[s]:
            n+=x.endswith('_HUMAN');y+=x.endswith('_YEAST')
        if n and y*counts['HUMAN']*100<=n*counts['YEAST']:best=max(best,n)
    return best
counts=read(OLD/'astral/frozen.json')['fasta_counts']
checks={name:exact([A/name/'proteins.tsv'],False,counts) for name in ['no_cell_scores','full']}
checks['separate_union']=exact([A/f'single_{i:03}/proteins.tsv' for i in range(3)],True,counts)
assert checks['no_cell_scores']==a['no_cell_scores']['0.01']['present']
assert checks['full']==a['full']['0.01']['present']
assert checks['separate_union']==fresh['curve']['0.01']['present']
# Validate complete outputs for every new search/alignment job.
verified=0
for tree in [G,BASE/'joint_gap_alternatives_v1',BASE/'joint_gap_retention_v1']:
    for job in tree.rglob('job.json'):
        j=read(job);assert j['returncode']==0,job
        assert all(sha(Path(f))==h for f,h in j['outputs'].items()),job
        verified+=1
assert all(fresh['protein_parity_with_frozen'].values())
assert all(sha(Path(f))==h for f,h in read(G/'diagnostic-hashes.json').items())
manifest=read(BASE/'joint_gap_alternatives_v1/source-manifest.json')
assert all(sha(ROOT/f)==h for f,h in manifest.items())
summary=dict(scope='Development IonStar five-run pilot and preselected Astral 100ng triplicate cross-dataset check; empirical entrapment FDP<=1%, accession counts.',
    ionstar={name:value['empirical_accession_curve']['0.01'] for name,value in r.items()
             if isinstance(value,dict) and 'empirical_accession_curve' in value},
    ionstar_alternative_rows=alt,ionstar_member_floor_03=loose,
    overlap={k:ret[k] for k in ['overlap','union_only','matrix_only']},
    retention_comparison=ret['retention_comparison'],stage_counts=stage['protein_counts'],
    missing_proteins_scored=sum(v['matrix_scored'] for v in ret['per_protein'].values()),
    astral={**{k:a[k]['0.01'] for k in ['no_cell_scores','full']},'fresh_separate_union':fresh['curve']['0.01']},
    validation=dict(exact_integer_counts=checks,fresh_control_parity=fresh['protein_parity_with_frozen'],
                    control_parity=read(G/'control-parity.json'),verified_jobs=verified,source_archive_matches=True,
                    tests={'uno_core':65,'uno_cli':7}),
    decisions=['Keep cell scores optional; q-derived information explains nearly all expanded-pilot improvement.',
               'Do not promote alternative-row competition or the lower KOTH member-quality floor.',
               'Matrix search beats the fresh individual union on this Astral cohort; IonStar gap remains.'])
(ROOT/'benchmark/joint_gap_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
no=a['no_cell_scores']['0.01']['present'];full=a['full']['0.01']['present'];union=fresh['curve']['0.01']['present']
text=f'''# Matrix gap investigation

On the preselected Astral 100ng triplicate, the joint matrix search identified
{full:,} present-species accessions versus {union:,} in the freshly rerun individual-search
union ({(full/union-1)*100:.1f}% higher), at empirical entrapment FDP <=1%. The matrix without
cell scores already identified {no:,}; the new cell scores added only {full-no} here.
All three fresh individual protein tables matched the frozen outputs exactly.
This is one cross-dataset check, not evidence of a universal matrix advantage.

## IonStar: what accounts for the gap?

The 1,533-accession matrix result shared 1,362 accessions with the 1,820 union,
gained 171 and missed 458. Of those 458, 427 received a matrix protein score but
fell below the empirical cutoff. At least two retained peptide-associated
candidates entered protein inference for 443 of them, so simple loss of all
input evidence is uncommon. Partial loss still matters: missing proteins retained
a median 70% of their independent quant-associated sequences, versus 86% for
shared proteins. These quant assignments are diagnostic associations, not
validated peptide IDs or a causal partition of the missing proteins.

## Fixed ablations and targeted tests

| IonStar matrix | No cell scores | MBR counts only | Q-derived scores only | All cell scores |
| --- | ---: | ---: | ---: | ---: |
| Standard | 1,458 | 1,304 | 1,268 | 1,309 |
| Expanded | 1,446 | 1,404 | 1,530 | 1,533 |

Q-derived scores provide nearly all the expanded-matrix gain in this pilot.
The ablation keeps confident counts with q-derived scores, since computing them
requires q-values; the count-only arm uses MBR-filled count and fraction.
Allowing coherent alternative rows during protein competition reduced depth to
1,372. Lowering KOTH's member-quality floor from 0.5 to 0.3 increased groups from
182,889 to 242,623 but reduced depth to 1,096. Both changes remain unpromoted;
more candidate evidence did not translate into more reliable identifications.

## Validation and interpretation

All counts use whole tied-score cutoffs selected for empirical entrapment
evaluation, not a claim that internal q-values were recalibrated. Astral used
one current executable, matching FASTA/model/settings, and fresh individual
controls. An independent integer calculation reproduced the matrix and union
counts. The Astral files were part of the earlier manuscript benchmark; they
were selected for this follow-up before examining new scores, but are not a
pristine held-out dataset. Runtime observations are diagnostic because some
jobs overlapped; no new speed claim is made.

All 65 core and 7 CLI tests passed. The diagnostic control reproduced protein,
LFQ and abundance outputs byte for byte, and all {verified} job output manifests were
verified. The selected-calibration diagnostic export and optional ablation
switches are in the isolated worktree. No paper or production defaults changed.
The matrix searches consume only KOTH bundles; independent-search outputs were
used solely as comparators and to trace losses.

The compact record is `benchmark/joint_gap_summary.json`. Reproduction drivers
are `joint_gap_pilot.py`, `joint_gap_audit.py`, `joint_gap_stage_audit.py`,
`joint_gap_alternatives.py`, `joint_gap_retention.py`, `joint_gap_validation.py`,
and `joint_gap_astral_controls.py` under `benchmark/`. Full receipts are under
`benchmark/results/lfq_matrix_fixes/joint_gap*/`.
'''
(ROOT/'docs/joint-gap-results.md').write_text(text)
print(json.dumps(dict(astral=summary['astral'],validation=summary['validation']),indent=2))
