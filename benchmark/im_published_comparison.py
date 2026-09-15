#!/usr/bin/env python3
"""Matched, peptide-disjoint calibration comparison with published IM2Deep.

Run in an environment containing im2deep==2.0.2 and CPU PyTorch >=2.6.
Inputs and predictions are archived; no test labels select either predictor.
"""
from __future__ import annotations
import csv
import hashlib
import importlib.metadata
import json
import math
import platform
import struct
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'benchmark/results/im_cnn_20260914'
OUT = ROOT / 'benchmark/results/im_published_20260914'
MODEL = ROOT / 'models/im_cnn_public_20260914.uimc'
OLD = ROOT / 'experiments/im-timstof/uno/benchmark/results/im_timstof_20260913'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return list(csv.DictReader(path.open(), delimiter='\t'))


def lineage():
    blob = MODEL.read_bytes()
    assert blob[:4] == b'UIMC'
    nb = struct.unpack_from('<I', blob, 28)[0]
    pos = 32 + 4 * nb + 8
    size = struct.unpack_from('<I', blob, pos)[0]
    return set(json.loads(blob[pos + 4:pos + 4 + size])['seen_peptide_hashes'])


def prepare():
    OUT.mkdir(parents=True, exist_ok=True)
    inputs = [SOURCE / 'test-predictions.tsv', SOURCE / 'A1-predictions.tsv', OLD / 'inputs/calibration.tsv', MODEL]
    plan = dict(
        comparator='im2deep 2.0.2 bundled TIMS IM2DeepUni.ckpt, single-output, CPU; fixed carbamidomethyl C',
        cohorts=['public reserved test', 'A1 MS/MS-assigned observations'],
        exclusions='Canonical sequences, length 7-30, charge 2-4, exclude every uno training/validation lineage sequence from both cohorts.',
        calibration='First floor(N_sequences/5) sequences ordered by SHA256 uno-im-published-v1:sequence; all charges together. Per-charge median observed/predicted ratio for BOTH models fitted only on calibration partition.',
        outcomes='Median and 95th percentile absolute relative error (%) and R2 on remaining observations; all charge states pooled. No model fitting or selection on either cohort.',
        caveat='IM2Deep training membership cannot be certified; public Meier data may overlap its training. This is not a mutually unseen-model ranking.',
        inputs={str(p.relative_to(ROOT)): sha(p) for p in inputs},
    )
    plan_path = OUT / 'plan.json'
    if plan_path.exists():
        assert json.loads(plan_path.read_text()) == plan
    else:
        plan_path.write_text(json.dumps(plan, indent=2) + '\n')
    seen = lineage()
    obs = {(r['sequence'], r['charge']): float(r['im']) for r in read(inputs[2])}
    cohorts = {}
    for name, rows in [('public', read(inputs[0])), ('A1', read(inputs[1]))]:
        values = []
        for r in rows:
            s, z = r['sequence'], int(r['charge'])
            if not (7 <= len(s) <= 30 and z in (2, 3, 4) and set(s) <= set('ACDEFGHIKLMNPQRSTVWY')):
                continue
            if hashlib.sha256(s.encode()).hexdigest() in seen:
                continue
            y = float(r['observed_im']) if name == 'public' else obs.get((s, str(z)))
            if y is None or not math.isfinite(y) or y <= 0:
                continue
            values.append(dict(sequence=s, charge=z, observed=y, uno=float(r['predicted_im'])))
        seqs = sorted({r['sequence'] for r in values}, key=lambda s: hashlib.sha256(('uno-im-published-v1:' + s).encode()).digest())
        cal = set(seqs[:len(seqs)//5])
        for r in values:
            r['partition'] = 'calibration' if r['sequence'] in cal else 'test'
        assert values and len({(r['sequence'], r['charge']) for r in values}) == len(values)
        cohorts[name] = values
    return cohorts, plan


def main():
    cohorts, plan = prepare()
    import numpy as np
    import torch
    import im2deep
    from im2deep import predict
    from psm_utils import PSM, PSMList
    torch.set_num_threads(4)
    summary = dict(plan=plan, plan_sha256=sha(OUT/'plan.json'), versions={p: importlib.metadata.version(p) for p in ('im2deep', 'torch', 'psm-utils', 'deeplc', 'numpy')}, python=platform.python_version(), cohorts={})
    weights = Path(im2deep.__file__).parent / 'models/TIMS/IM2DeepUni.ckpt'
    summary['published_model_sha256'] = sha(weights)
    summary['uno_model_sha256'] = sha(MODEL)
    # Same physical conversion used for the public labels, followed by identical
    # charge-wise calibration. Fixed CAM is supplied explicitly to IM2Deep.
    mass = dict(zip('ARNDCEQGHILKMFPSTWYV', [71.037113805,156.10111105,114.04292747,115.026943065,103.009184505,129.042593135,128.05857754,57.021463735,137.058911875,113.084063975,113.084063975,128.094963015,131.040484645,147.068413945,97.052763875,87.032028435,101.047678505,186.07931298,163.063328575,99.068413945]))
    for name, rows in cohorts.items():
        cache = OUT / f'{name}-native-ccs.json'
        keys = [(r['sequence'], r['charge']) for r in rows]
        query_sha = hashlib.sha256(json.dumps(keys).encode()).hexdigest()
        if cache.exists():
            saved = json.loads(cache.read_text())
            assert saved['query_sha256'] == query_sha and saved['model_sha256'] == summary['published_model_sha256']
            ccs = saved['predictions']
        else:
            start = time.monotonic()
            psms = PSMList(psm_list=[PSM(peptidoform=s.replace('C', 'C[UNIMOD:4]') + '/' + str(z), spectrum_id=str(i)) for i, (s,z) in enumerate(keys)])
            ccs = np.asarray(predict(psms, multi=False, predict_kwargs=dict(device='cpu', batch_size=512, num_workers=0, num_threads=4))).reshape(-1).tolist()
            cache.write_text(json.dumps(dict(query_sha256=query_sha, model_sha256=summary['published_model_sha256'], seconds=time.monotonic()-start, predictions=ccs)) + '\n')
        assert len(ccs) == len(rows)
        for r, pred in zip(rows, ccs):
            neutral = 18.010564684 + sum(mass[a] + (57.021463735 if a == 'C' else 0) for a in r['sequence'])
            r['im2deep_ccs'] = pred
            r['im2deep'] = pred * math.sqrt(neutral * 28.0134 / (neutral + 28.0134)) / r['charge'] / 1000
        counts = {part: dict(observations=sum(r['partition']==part for r in rows), sequences=len({r['sequence'] for r in rows if r['partition']==part})) for part in ('calibration', 'test')}
        result = dict(counts=counts, models={})
        for model in ('uno', 'im2deep'):
            scales = {z: float(np.median([r['observed']/r[model] for r in rows if r['partition']=='calibration' and r['charge']==z])) for z in (2,3,4)}
            assert all(math.isfinite(x) and x>0 for x in scales.values())
            for r in rows:
                r[model+'_calibrated'] = r[model]*scales[r['charge']]
            test = [r for r in rows if r['partition']=='test']
            y = np.array([r['observed'] for r in test]); p = np.array([r[model+'_calibrated'] for r in test])
            assert np.isfinite(p).all() and (p>0).all()
            err = 100*np.abs(p-y)/y
            result['models'][model] = dict(calibration_scales=scales, median_are_percent=float(np.median(err)), p95_are_percent=float(np.quantile(err,.95)), r2=float(1-np.sum((p-y)**2)/np.sum((y-y.mean())**2)))
        path = OUT / f'{name}-predictions.tsv'
        with path.open('w') as f:
            w=csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t');w.writeheader();w.writerows(rows)
        result['predictions_sha256'] = sha(path)
        summary['cohorts'][name]=result
        print(name, json.dumps(result), flush=True)
    summary['script_sha256'] = sha(Path(__file__))
    (ROOT/'benchmark/im_published_20260914.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    (OUT/'environment.txt').write_text(subprocess.check_output([__import__('sys').executable,'-m','pip','freeze'], text=True))

if __name__ == '__main__':
    main()
