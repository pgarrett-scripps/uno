#!/usr/bin/env python3
"""Build a concise final summary from verified joint-search pilot artifacts."""
import json
import re
from pathlib import Path
from joint_evidence_audit import BASE, OLD, ROOT, curve, sha
from matrix_report import load_proteins


def main():
    v2 = json.loads((BASE / 'joint_search_v2/report.json').read_text())
    v3 = json.loads((BASE / 'joint_search_v3/report.json').read_text())
    retention = json.loads((BASE / 'aligned5-replicated/retention-audit.json').read_text())
    timing = v3['timing_comparison']
    counts = json.loads((OLD / 'frozen.json').read_text())['fasta_counts']
    scale = (counts['HUMAN'] + counts['ECOLI']) / counts['YEAST']
    scores = {}
    for row in load_proteins(OLD / 'searches/full_matrix_5/seed0/proteins.tsv'):
        if row['species'] == 'AMBIGUOUS': continue
        for a in row['members']: scores[a] = max(scores.get(a, 0.0), row['score'])
    frozen = curve(scores, scale)['0.01']
    frozen.pop('accessions')
    tests = {}
    for name, path in [('uno_core', '/tmp/uno-joint-tests.log'), ('uno_cli', '/tmp/uno-joint-cli-tests.log'),
                       ('koth', '/tmp/koth-joint-tests.log')]:
        text = Path(path).read_text()
        m = re.findall(r'test result: ok\. (\d+) passed; 0 failed', text)
        assert len(m) == 1, path
        tests[name] = int(m[0])
    assert all(v2['legacy_parity'].values())
    assert all(json.loads((BASE / 'joint_search_v3/strong-group-parity.json').read_text()).values())
    for i in range(5):
        path = BASE / f'joint_search_v3/single_run{i:03}/proteins.tsv'
        assert sha(path) == sha(BASE / f'E_batch/run{i:03}/proteins.tsv')
    summary = dict(scope='IonStar five-run pilot, seed zero; actual matrix searches, not pooled independent scores.',
        frozen_matrix=frozen, representative=v2['representative_control'],
        joint=v3['joint'], no_agreement=v2['no_agreement'],
        replicated_without_bootstrap_guard=v2['replicated_joint'],
        replicated_no_agreement=v2['replicated_no_agreement'],
        replicated_with_bootstrap_guard=v3['replicated_joint'],
        separate_union_at_empirical_1pct=v3['separate_union_reference'],
        retention=retention, timing=timing, tests=tests,
        parity='Original representative outputs, unchanged strong-group joint outputs, and all five independent protein tables are byte-identical to their references.',
        decision='Keep both new modes explicit and experimental. Joint scoring improves the matrix pilot but does not beat the separate-search union. Expanded KOTH retention does not improve 1% depth.',
        hashes={str(p): sha(p) for p in [BASE / 'joint_search_v2/report.json', BASE / 'joint_search_v3/report.json',
            BASE / 'joint_search_v3/uno', BASE / 'joint_search_v2/koth_align',
            ROOT / 'crates/uno-core/src/joint.rs', ROOT / 'crates/uno-core/src/score.rs',
            ROOT / 'crates/uno-core/src/long_input.rs']})
    (ROOT / 'benchmark/joint_matrix_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    def depth(result): return result['empirical_accession_curve']['0.01']['present']
    joint, rep, expanded = depth(v3['joint']), depth(v2['representative_control']), depth(v3['replicated_joint'])
    gain = 100 * (joint / frozen['present'] - 1)
    doc = f'''# Joint KOTH matrix search: implementation and pilot

UNO now performs one joint search from the long KOTH matrix, scoring original
observations rather than merging independent-search results. On the five-run
IonStar pilot it recovered {joint:,} present-species accessions at an empirical
1% error threshold, {gain:.1f}% above the frozen matrix baseline. It remains below
the {v3['separate_union_reference']['present']:,}-accession separate-search union
at the same error criterion.

| Search | Present accessions at empirical 1% |
| --- | ---: |
| Frozen original matrix | {frozen['present']:,} |
| Previous isotope-supported representative | {rep:,} |
| Joint original-observation scoring | {joint:,} |
| Joint scoring, agreement descriptors disabled | {depth(v2['no_agreement']):,} |
| Expanded KOTH retention, guarded calibration | {expanded:,} |
| Separate-search union, threshold adjusted | {v3['separate_union_reference']['present']:,} |

Counts expand same-species protein groups into accessions. Whole score ties are
kept together, and thresholds are selected using entrapment for development
evaluation. These are one-cohort, one-seed results with only a few entrapment
accessions at the cutoff; they do not establish general sensitivity or calibrated
cohort error on independent data. Several targeted development variants were
examined and remain recorded in the JSON summary.

## Implemented changes

`--matrix-long BUNDLE --matrix-joint` searches the union of original member mass
windows once per aligned group, shares digestion and RT prediction, and estimates
one cohort protein result. Candidate-specific mass/RT consistency and original
isotope/shape descriptors feed the classifier. Alternatives contribute at most
one observation per source run; missing runs and extraction-only MBR cells do not
create identification evidence. Joint mode retains KOTH singleton groups by
default, and recalibration uses actual compatible original observations.

KOTH's optional `allow_replicated_weak_seeds = true` retains moderate-quality
groups supported by at least two runs under the existing coordinate bounds.
It increased searchable original observations from 524,060 to 853,370 and groups
from 99,325 to 182,889. The audit preserved every original measurement, all prior
group memberships, and all 496,625 prior integrated LFQ cells exactly. Newly
retained groups still require quantitative accuracy validation.

Expanded retention raised the initial random-match background from 0.325 to
0.442. Excluding groups below KOTH's original seed floor from initial calibration
restored 0.325 while keeping them available for identification. This did not
improve 1% depth over the smaller joint matrix, so expanded retention remains
opt-in. Its higher depth at looser error thresholds is recorded separately.

## Runtime and validation

The joint search took {v3['joint']['seconds']:.1f} seconds versus
{timing['individual_sum_seconds']:.1f} seconds for five independent searches
({timing['search_speed_ratio']:.1f} times faster). Both used eight threads, warm
prediction caches, the same FASTA/model/config and protein/LFQ/abundance outputs,
with serial execution. These are single timing measurements with different
identification coverage; they exclude feature finding and matrix alignment.
The expanded matrix separately required
{json.loads((BASE / 'aligned5-replicated/job.json').read_text())['seconds']:.1f}
seconds for alignment/LFQ construction and {v3['replicated_joint']['seconds']:.1f}
seconds for search.

All {sum(tests.values())} tests pass ({tests['koth']} KOTH, {tests['uno_core']} UNO
core, {tests['uno_cli']} CLI). Golden comparisons preserve the previous
representative outputs, the unaffected strong-group joint outputs and all five
independent protein tables byte for byte. The new routes remain experimental in
the paired worktrees; no paper claims or default search behavior were promoted.

Reproduction is in `benchmark/joint_matrix_pilot.py`; full per-run outputs,
commands, hashes and timings are under `benchmark/results/lfq_matrix_fixes/`.
The compact result is `benchmark/joint_matrix_summary.json`. The remaining gap is
joint identification sensitivity at matched error, rather than preservation of
original observations or shared prediction caching.
'''
    (ROOT / 'docs/joint-matrix-search-results.md').write_text(doc)
    print(json.dumps(dict(joint=joint, frozen=frozen['present'], separate_union=v3['separate_union_reference']['present'],
        seconds=v3['joint']['seconds'], independent_seconds=timing['individual_sum_seconds'], tests=sum(tests.values())), indent=2))


if __name__ == '__main__':
    main()
