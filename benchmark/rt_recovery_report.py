#!/usr/bin/env python3
"""Verify and evaluate the frozen recovery-only RT investigation."""
import argparse
import json
import math
from pathlib import Path
import stability_soft_report as report

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'benchmark/results/rt_recovery_investigation'
SEEDS = {'lowinput': [9, 10], 'primary': [11, 12], 'raw': [5, 6], 'additional': [3, 4]}


def preservation_checks(panel, pairs):
    positive = [p for p in pairs if p['baseline_q1'] > 0]
    failed = [p for p in pairs if p['baseline_q1'] == p['candidate_q1'] == 0]
    return [dict(criterion=f'{panel}_{name}_full_output_preserved',
                 passed=all(p['byte_equal'] for p in subset),
                 checked_pairs=len(subset),
                 detail=[{'dataset': p['dataset'], 'seed': p['seed']} for p in subset if not p['byte_equal']])
            for name, subset in [('v3_positive', positive), ('both_empty', failed)]]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=['screen', 'complete'], required=True)
    args = parser.parse_args()
    hashes = {}
    manifest_path = BASE / 'plan-manifest.json'
    manifest = report.read(manifest_path)
    hashes[str(manifest_path)] = report.sha(manifest_path)
    for path, expected in manifest['files'].items():
        report.verify(path, expected, hashes)
    registry = report.read(BASE / 'candidate-registry.json')
    candidate_sha = registry['arms']['recovery']['binary_sha256']
    report.verify_candidate_registry(BASE / 'candidate-registry.json', 'recovery', candidate_sha, hashes)
    development_path = BASE / 'development-verification.json'
    development = report.read(development_path)
    assert development['passed'] and development['candidate_sha256'] == candidate_sha
    assert development['plan_manifest_sha256'] == hashes[str(manifest_path)]
    hashes[str(development_path)] = report.sha(development_path)
    for path, expected in development['source_hashes'].items():
        report.verify(path, expected, hashes)
    if args.stage == 'complete':
        prior_path = BASE / 'screen-report.json'
        prior = report.read(prior_path)
        hashes[str(prior_path)] = report.sha(prior_path)
        assert prior['stage'] == 'screen' and prior['gate_pass'] and prior['candidate_sha256'] == candidate_sha
        assert prior['plan_manifest_sha256'] == hashes[str(manifest_path)]
        for path, expected in prior['source_hashes'].items():
            report.verify(path, expected, hashes)
    names = ['lowinput'] if args.stage == 'screen' else list(SEEDS)
    panels, checks = {}, []
    runtime = None
    for name in names:
        study = report.read(BASE / name / 'study.json')
        for job in study['searches']:
            stamp = report.read(Path(job['output']).with_suffix('.json'))
            assert not job.get('reused_baseline') and not stamp.get('reused_baseline') and not stamp.get('reused_from'), 'new panels require both engines fresh'
            assert all(isinstance(stamp.get(k), (int, float)) and math.isfinite(stamp[k]) and stamp[k] > 0
                       for k in ('started_unix', 'elapsed_seconds')), 'missing fresh execution timing'
        rows, scale, tables, inherited = report.load_panel(BASE / name, name, SEEDS[name], candidate_sha,
                                                          hashes, report.read(BASE / 'additional-inputs.json'))
        if runtime is None:
            runtime = inherited
        assert runtime == inherited, 'runtime environments differ across panels'
        panel = report.summarize_panel(rows, scale, name, 'confirmation', tables)
        panel['checks'] += preservation_checks(name, panel['pairs'])
        panels[name] = panel
        checks += panel['checks']
    passed = all(c['passed'] is True for c in checks)
    decision = ('PASS: ADVANCE TO CONTROLS' if args.stage == 'screen' else 'PASS: BOUNDED RECOVERY CANDIDATE') if passed else report.rejection_label(checks)
    result = dict(stage=args.stage, candidate_sha256=candidate_sha, baseline_sha256=report.V3_SHA,
                  plan_manifest_sha256=hashes[str(manifest_path)], gate_pass=passed,
                  decision=decision, candidate_promoted=False, checks=checks, panels=panels,
                  source_hashes=hashes,
                  caveats=['New random draws on known files are correlated, not independent biological replication.',
                           'The recovery policy was selected after inspecting prior RT results.',
                           'Preserving successful outputs does not prove error control among selected rescues.',
                           'Relative ranges use exactly the same positive-median files in both arms; excluded files remain explicit.',
                           'This investigation does not promote the default engine.'])
    (BASE / f'{args.stage}-report.json').write_text(json.dumps(result, indent=2) + '\n')
    lines = [f'# RT recovery: {args.stage}', '', f'**{decision}**', '',
             '| Panel | Pairs | Empty V3→recovery | Human / yeast V3→recovery | Pooled FDP V3→recovery |',
             '|---|---:|---:|---|---|']
    def fdp(value):
        return 'undefined' if value is None else f'{100 * value:.3f}%'
    for name, panel in panels.items():
        a, b = (panel['aggregates']['all'][v] for v in ('baseline', 'candidate'))
        lines.append(f"| {name} | {a['n']} | {a['empty_q1']}→{b['empty_q1']} | {a['human']}/{a['yeast']}→{b['human']}/{b['yeast']} | {fdp(a['pooled_fdp'])}→{fdp(b['pooled_fdp'])} |")
    lines += ['', '## Checks', '']
    lines += [f"- {'PASS' if c['passed'] is True else 'FAIL' if c['passed'] is False else 'UNDEFINED'}: {c['criterion']}" for c in checks]
    lines += ['', *result['caveats'], '', 'All paired metrics, individual FDPs, rescued-only error, output preservation and excluded ranges are retained in the JSON.']
    (BASE / f'{args.stage}-report.md').write_text('\n'.join(lines) + '\n')
    print(decision)


if __name__ == '__main__':
    main()
