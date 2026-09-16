"""Generate reproducible native evidence, run unmodified reference and Rust, audit parity.

Run with the locked environment: uv run --project tools/quant-reference
    python tools/quant-reference/check_parity.py --out target/parity-baseline
Outputs are deliberately fresh; previously generated evidence is never replaced.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import random
import subprocess
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
REFERENCE = Path(__file__).parent / 'reference' / 'directms1quant.py'
REFERENCE_SHA = 'cb15b7c43b276eb45e77da2b8cf7cd97cd20765f95264202283212fde85cfc90'


def write_tsv(path, columns, rows):
    with path.open('w', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(columns)
        writer.writerows(rows)


def sequence(i):
    alphabet = 'ACDEFGHIKLMNPQRSTVWY'
    letters = []
    for _ in range(5):
        letters.append(alphabet[i % len(alphabet)])
        i //= len(alphabet)
    return 'PEP' + ''.join(letters) + 'K'


def fixture(path, missing=False, structural=False):
    path.mkdir()
    rng = np.random.default_rng(76109)
    samples = [[] for _ in range(6)]
    ml = []
    for p in range(120):
        decoy = p >= 60
        name = ('DECOY_' if decoy else '') + f'P{p % 60:03}'
        if structural and p % 60 == 5:
            name += ';' + name + '_ALT'
        truth = 2.2 if p < 10 else (-2.0 if p < 20 else 0.)
        for k in range(8):
            i = p * 8 + k
            seq = sequence(i)
            baseline = 2 ** rng.uniform(12, 26)
            effect = truth + rng.normal(0, 0.55 if not decoy else 0.8)
            noise = rng.normal(0, 0.2, 6)
            for sample in range(6):
                if missing and ((i % 17 == 0 and sample >= 3) or (i % 19 == 0 and sample < 3) or (i % 23 == sample)):
                    continue
                value = baseline * 2 ** (noise[sample] + (effect if sample >= 3 else 0))
                samples[sample].append([seq, name, 2, 0.0, value, 12.5 + i / 100])
                if structural and i % 11 == 0:
                    samples[sample].append([seq, name, 2, 0.0, value * 0.5, 13.5 + i / 100])
                if structural and i % 13 == 0:
                    samples[sample].append([seq, name, 3, 0.0, value * (1.01 + sample / 100), 12.5 + i / 100])
            ml.append([seq, 1, len(seq)])
    paths = []
    for sample, rows in enumerate(samples):
        base = path / f'run{sample}'
        proteins = [[f'P{p:03}'] for p in range(60)]
        write_tsv(Path(str(base) + '_proteins.tsv'), ['dbname'], proteins)
        write_tsv(Path(str(base) + '_proteins_full.tsv'), ['dbname'], proteins + [[f'DECOY_P{p:03}'] for p in range(60)])
        write_tsv(Path(str(base) + '_PFMs_ML.tsv'), ['seqs', 'qpreds', 'plen'], ml)
        write_tsv(Path(str(base) + '_PFMs.tsv'), ['sequence', 'proteins', 'charge', 'ion_mobility', 'Intensity', 'RT'], rows)
        paths.append(str(base) + '_proteins_full.tsv')
    return paths


def table(path, keys):
    with path.open() as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    result = {tuple(row[key] for key in keys): row for row in rows}
    assert len(result) == len(rows), f'duplicate comparison keys in {path}'
    return result


def compare(reference, rust, numeric, categorical, abs_tol, rel_tol):
    differences = []
    maxima = {name: 0. for name in numeric}
    if reference.keys() != rust.keys():
        differences.append({'kind': 'keys', 'reference_only': list(reference.keys() - rust.keys())[:10], 'rust_only': list(rust.keys() - reference.keys())[:10]})
    counts = {name: 0 for name in numeric + categorical}
    for key in reference.keys() & rust.keys():
        a, b = reference[key], rust[key]
        for col in numeric:
            x, y = float(a[col]), float(b[col])
            if math.isfinite(x) and math.isfinite(y):
                maxima[col] = max(maxima[col], abs(x - y))
            if not (math.isclose(x, y, abs_tol=abs_tol, rel_tol=rel_tol) or (math.isnan(x) and math.isnan(y))):
                counts[col] += 1
                if len(differences) < 20:
                    differences.append({'key': key, 'column': col, 'reference': x, 'rust': y})
        for col in categorical:
            if a[col] != b[col]:
                counts[col] += 1
                if len(differences) < 20:
                    differences.append({'key': key, 'column': col, 'reference': a[col], 'rust': b[col]})
    return {'passed': not differences, 'n_reference': len(reference), 'n_rust': len(rust), 'mismatch_counts': counts, 'max_absolute_error': maxima, 'examples': differences}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--binary', type=Path, default=ROOT/'target/debug/uno-quant')
    parser.add_argument('--missing', action='store_true')
    parser.add_argument('--structural', action='store_true')
    parser.add_argument('--normalization', type=int, choices=(0, 1, 2), default=2)
    parser.add_argument('--absolute-fold-change', action='store_true')
    args = parser.parse_args()
    assert hashlib.sha256(REFERENCE.read_bytes()).hexdigest() == REFERENCE_SHA
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    paths = fixture(out/'inputs', missing=args.missing, structural=args.structural)
    commands = [
        [sys.executable, str(REFERENCE), '-S1', *paths[:3], '-S2', *paths[3:], '-out', str(out/'reference'), '-intensity_norm', str(args.normalization)],
        [str(args.binary.resolve()), '--s1', *paths[:3], '--s2', *paths[3:], '--out', str(out/'rust'), '--normalization', str(args.normalization)]
    ]
    if args.absolute_fold_change:
        commands[0] += ['-fold_change_abs']
        commands[1] += ['--absolute-fold-change']
    for name, command in zip(('reference', 'rust'), commands):
        with (out/f'{name}.log').open('w') as log:
            subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT, env={**os.environ, "NPY_DISABLE_CPU_FEATURES": "AVX512F,AVX2"})
    pep_a = table(out/'reference_quant_peptides.tsv', ['proteins', 'origseq'])
    pep_b = table(out/'rust/peptides.tsv', ['proteins', 'origseq'])
    prot_a = table(out/'reference_quant_full.tsv', ['dbname'])
    prot_b = table(out/'rust/proteins.tsv', ['dbname'])
    peptide = compare(pep_a, pep_b, ['FC_raw', 'FC', 'FC_corrected', 'p-value', 'FC_gr_mean', 'intensity_median'], ['peptide', 'decoy', 'nummissing', 'nonmissing_S1', 'nonmissing_S2', 'iq', 'sign', 'up', 'down'], 1e-8, 1e-6)
    # Upstream serializes bp as floats, Rust as integers: compare numerically and exactly.
    bp = compare(pep_a, pep_b, ['bp'], [], 0., 0.)
    protein = compare(prot_a, prot_b, ['score', 'p-value', 'log2FoldChange(S2/S1)', 'log2FoldChange(S2/S1) no normalization', 'log2FoldChange(S2/S1) using all peptides', 'log2FoldChange(S2/S1) using all peptides and no normalization'], ['differentially expressed peptides', 'identified peptides', 'FC_pass', 'FDR_pass', 'BH_pass'], 1e-8, 1e-6)
    # Group identifiers may be renumbered; the partition must be identical.
    def partition(t):
        groups = {}
        for key, row in t.items():
            groups.setdefault(row['protein_quant_group'], set()).add(key)
        return {frozenset(members) for members in groups.values()}
    group_parity = partition(prot_a) == partition(prot_b)
    report = {'status': 'pass' if all(x['passed'] for x in (peptide, bp, protein)) and group_parity else 'fail', 'reference_sha256': REFERENCE_SHA, 'reference_environment': {'NPY_DISABLE_CPU_FEATURES': 'AVX512F,AVX2'}, 'fixture_seed': 76109, 'missing': args.missing, 'structural': args.structural, 'normalization': args.normalization, 'absolute_fold_change': args.absolute_fold_change, 'commands': commands, 'peptides': peptide, 'resampling': bp, 'proteins': protein, 'group_partition_equal': group_parity, 'scope': 'Synthetic numerical parity, not biological accuracy or FDR validation'}
    (out/'parity.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    return 0 if report['status'] == 'pass' else 1

if __name__ == '__main__':
    sys.exit(main())
