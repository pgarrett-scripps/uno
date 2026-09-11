#!/usr/bin/env python3
"""Verify the frozen bounded-RT recovery evaluation without launching searches."""
import argparse
import json
import math
from fractions import Fraction
from pathlib import Path
import stability_soft_report as report
from rt_recovery_report import preservation_checks

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'benchmark/results/bounded_rt_recovery'
SEEDS = {'known': [0, 1, 2, 5, 8], 'fresh': [11, 12, 13, 14, 15],
         'primary': [3, 4], 'raw': [5, 6], 'additional': [5, 6]}


def recovery_checks(name, panel):
    rescue = panel['rescued_only']['all']
    pairs = panel['pairs']
    rescued_keys = {(p['dataset'], p['seed']) for p in pairs if p['baseline_q1'] == 0 and p['candidate_q1'] > 0}
    rescued_rows = [r for r in panel['individual_results'] if r['variant'] == 'candidate' and
                    (r['dataset'], r['seed']) in rescued_keys]
    pooled = Fraction(20647 * rescue['yeast'], 6066 * rescue['human']) if rescue['human'] else None
    checks = preservation_checks(name, pairs)
    if name in ('known', 'fresh'):
        checks.append(dict(criterion=f'{name}_at_least_one_rescue', passed=rescue['n'] >= 1, rescued_runs=rescue['n']))
    checks.append(dict(criterion=f'{name}_rescued_pooled_fdp_at_most1pct',
                       passed=True if not rescue['n'] else pooled <= Fraction(1, 100) if pooled is not None else None,
                       applicable=bool(rescue['n']), detail=rescue))
    too_high = [dict(dataset=r['dataset'], seed=r['seed'], fdp=r['fdp']) for r in rescued_rows
                if not r['internal_q1_human'] or Fraction(20647 * r['internal_q1_yeast'], 6066 * r['internal_q1_human']) > Fraction(1, 20)]
    checks.append(dict(criterion=f'{name}_individual_rescued_fdp_at_most5pct', passed=not too_high,
                       applicable=bool(rescued_rows), detail=too_high))
    return checks


def require_fresh(study, both):
    for job in study['searches']:
        if not both and job['variant'] == 'baseline':
            continue
        stamp = report.read(Path(job['output']).with_suffix('.json'))
        assert not job.get('reused_baseline') and not stamp.get('reused_baseline') and not stamp.get('reused_from'), 'execution must be fresh'
        assert all(isinstance(stamp.get(k), (int, float)) and math.isfinite(stamp[k]) and stamp[k] > 0
                   for k in ('started_unix', 'elapsed_seconds')), 'missing fresh execution timing'


def verify_prior(path, candidate_sha, plan_sha, hashes, expected_stage=None):
    prior = report.read(path)
    hashes[str(path)] = report.sha(path)
    assert prior['candidate_sha256'] == candidate_sha and prior['baseline_sha256'] == report.V3_SHA
    assert prior['plan_manifest_sha256'] == plan_sha and prior['gate_pass']
    if expected_stage:
        assert prior['stage'] == expected_stage
    for name, expected in prior['source_hashes'].items():
        report.verify(name, expected, hashes)
    return prior


def verify_timing(candidate_sha, plan_sha, hashes):
    timing = verify_prior(BASE / 'timing/report.json', candidate_sha, plan_sha, hashes, 'timing')
    assert timing['schema'] == 1 and timing['executions'] == 16
    cases = {'normal_lowinput', 'normal_primary', 'failure_1ng', 'failure_200pg'}
    assert {c['case'] for c in timing['cases']} == cases and len(timing['cases']) == 4
    expected = {(case, round_no, variant) for case in cases for round_no in (1, 2) for variant in ('baseline', 'candidate')}
    actual = [(r['case'], r['round'], r['variant']) for r in timing['executed_runs']]
    assert len(actual) == 16 and set(actual) == expected
    checks = {c['criterion']: c for c in timing['checks']}
    assert len(checks) == len(timing['checks']), 'duplicate timing checks'
    required = {'normal_median_case_speed', 'failure_median_case_speed', 'all_fresh_positive_times',
                'all_inputs_binaries_stable', 'all_warm_zero_miss_cache', 'all_no_visible_competing_workers',
                'all_candidate_outputs_match_efficacy', 'all_required_v3_outputs_preserved'}
    for case in cases:
        required.add(f'{case}_case_speed')
        for variant in ('baseline', 'candidate'):
            required.update({f'{case}_{variant}_exact_repeat', f'{case}_{variant}_spread_at_most1.10'})
        for round_no in (1, 2):
            required.update({f'{case}_r{round_no}_single_fallback_structure', f'{case}_r{round_no}_declared_baseline_regime'})
    assert required.issubset(checks) and all(c['passed'] is True for c in checks.values())
    return timing


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=['known', 'fresh', 'complete'], required=True)
    args = parser.parse_args()
    hashes = {}
    plan_path = BASE / 'plan-manifest.json'
    plan = report.read(plan_path)
    plan_sha = report.sha(plan_path)
    hashes[str(plan_path)] = plan_sha
    for name, expected in plan['files'].items():
        report.verify(name, expected, hashes)
    registry = report.read(BASE / 'candidate-registry.json')
    candidate_sha = registry['arms']['bounded']['binary_sha256']
    report.verify_candidate_registry(BASE / 'candidate-registry.json', 'bounded', candidate_sha, hashes)
    if args.stage != 'known':
        verify_prior(BASE / 'known-report.json', candidate_sha, plan_sha, hashes, 'known')
        verify_timing(candidate_sha, plan_sha, hashes)
    if args.stage == 'complete':
        verify_prior(BASE / 'fresh-report.json', candidate_sha, plan_sha, hashes, 'fresh')
    names = ['known'] if args.stage == 'known' else ['fresh'] if args.stage == 'fresh' else list(SEEDS)
    panels, checks = {}, []
    inherited_reference = None
    for name in names:
        directory = BASE / name
        study = report.read(directory / 'study.json')
        require_fresh(study, both=name != 'known')
        loader_name = 'lowinput' if name in ('known', 'fresh') else name
        rows, scale, tables, inherited = report.load_panel(directory, loader_name, SEEDS[name], candidate_sha,
                                                          hashes, report.read(BASE / 'additional-inputs.json'))
        if inherited_reference is None:
            inherited_reference = inherited
        assert inherited == inherited_reference
        panel = report.summarize_panel(rows, scale, name, 'confirmation', tables)
        panel['checks'] += recovery_checks(name, panel)
        panels[name] = panel
        checks += panel['checks']
    if args.stage == 'complete':
        timing = report.read(BASE / 'timing/report.json')
        checks += timing['checks']
    passed = all(c['passed'] is True for c in checks)
    decision = ('PASS: KNOWN-CASE ACCURACY' if args.stage == 'known' else
                'PASS: FRESH-SEED ACCURACY' if args.stage == 'fresh' else
                'PASS: ELIGIBLE FOR PROMOTION') if passed else report.rejection_label(checks)
    result = dict(stage=args.stage, candidate_sha256=candidate_sha, baseline_sha256=report.V3_SHA,
                  plan_manifest_sha256=plan_sha, gate_pass=passed, decision=decision,
                  eligible_for_promotion=passed and args.stage == 'complete', checks=checks,
                  panels=panels, source_hashes=hashes,
                  caveats=['Known cases are development evidence; fresh seeds are new draws on known acquisitions, not biological holdouts.',
                           'Recovery can improve empty results but does not change successful V3 rankings.',
                           'Every rescued-only error result is evaluated separately from preserved V3 successes.',
                           'Relative ranges use identical positive-median files and explicitly retain exclusions.',
                           'Finite empirical entrapment and timing tolerances do not prove nominal per-run FDR or universal performance.'])
    (BASE / f'{args.stage}-report.json').write_text(json.dumps(result, indent=2) + '\n')
    lines = [f'# Bounded RT recovery: {args.stage}', '', f'**{decision}**', '',
             '| Panel | Pairs | Empty V3→candidate | H/Y V3→candidate | Pooled FDP V3→candidate | Rescues H/Y; FDP |',
             '|---|---:|---:|---|---|---|']
    def fdp(value):
        return 'undefined' if value is None else f'{100 * value:.3f}%'
    for name, panel in panels.items():
        a, b = (panel['aggregates']['all'][v] for v in ('baseline', 'candidate'))
        rr = panel['rescued_only']['all']
        lines.append(f"| {name} | {a['n']} | {a['empty_q1']}→{b['empty_q1']} | {a['human']}/{a['yeast']}→{b['human']}/{b['yeast']} | {fdp(a['pooled_fdp'])}→{fdp(b['pooled_fdp'])} | {rr['human']}/{rr['yeast']}; {fdp(rr['pooled_fdp'])} |")
    lines += ['', '## Checks', '']
    lines += [f"- {'PASS' if c['passed'] is True else 'FAIL' if c['passed'] is False else 'UNDEFINED'}: {c['criterion']}" for c in checks]
    lines += ['', *result['caveats'], '', 'All individual metrics, rescue membership, losses, error estimates and count-range exclusions are retained in the JSON.']
    (BASE / f'{args.stage}-report.md').write_text('\n'.join(lines) + '\n')
    print(decision)


if __name__ == '__main__':
    main()
