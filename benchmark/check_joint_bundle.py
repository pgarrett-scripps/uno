#!/usr/bin/env python3
"""Verify expanded KOTH retention preserves original measurements and LFQ cells."""
import csv
import itertools
import json
from pathlib import Path
from joint_evidence_audit import BASE, sha


def cells(path):
    groups, values = {}, {}
    with path.open() as f:
        for r in csv.DictReader(f, delimiter='\t'):
            if r['extraction_kind'] != 'target': continue
            c = r['consensus_id']
            seed = (r['seed_run_id'], r['seed_observation_id'])
            groups[c] = seed
            key = (*seed, r['run_id'])
            assert key not in values
            values[key] = (r['lfq_intensity'], r['observation_id'], r['mass'], r['rt'], r['charge'])
    return groups, values


def main():
    old, new = BASE / 'aligned5', BASE / 'aligned5-replicated'
    og, ov = cells(old / 'lfq_matrix.long.tsv')
    ng, nv = cells(new / 'lfq_matrix.long.tsv')
    assert set(og.values()) <= set(ng.values())
    assert all(nv.get(k) == v for k, v in ov.items()), 'Existing LFQ cells or metadata changed'
    with (old / 'feature_observations.tsv').open() as a, (new / 'feature_observations.tsv').open() as b:
        ar, br = csv.DictReader(a, delimiter='\t'), csv.DictReader(b, delimiter='\t')
        assert ar.fieldnames == br.fieldnames
        fields = [k for k in ar.fieldnames if k not in ('consensus_id', 'is_primary', 'is_seed', 'exclusion_reason')]
        n = members_before = members_after = 0
        for x, y in itertools.zip_longest(ar, br):
            assert x is not None and y is not None
            assert all(x[k] == y[k] for k in fields), (x['run_id'], x['observation_id'])
            if x['consensus_id']:
                assert y['consensus_id'] and og[x['consensus_id']] == ng[y['consensus_id']]
                assert x['is_primary'] == y['is_primary'] and x['is_seed'] == y['is_seed']
                members_before += 1
            members_after += bool(y['consensus_id'])
            n += 1
    report = dict(original_observations=n, original_fields_and_alignment_identical=True,
        prior_group_memberships_preserved=True, prior_integrated_cells_exact=len(ov),
        original_groups=len(og), expanded_groups=len(ng), original_member_observations=members_before,
        expanded_member_observations=members_after,
        note='Existing intensities preserved; extraction q-values may change with the expanded target/decoy population. This is not a quantitative accuracy validation of newly retained groups.',
        hashes={str(p): sha(p) for p in [Path(__file__), old / 'matrix_manifest.json', new / 'matrix_manifest.json']})
    (new / 'retention-audit.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
