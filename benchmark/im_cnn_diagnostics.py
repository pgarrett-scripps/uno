#!/usr/bin/env python3
"""Prediction diagnostics independent of search outcome selection."""
import csv,json,statistics,sys,hashlib
from collections import defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'benchmark/results/im_cnn_20260914';OLD=ROOT/'experiments/im-timstof/uno/benchmark/results/im_timstof_20260913'
sys.path.insert(0,str(ROOT/'crates/uno-rt-train/python'))
import numpy as np
import torch
from uno_rt.artifact import load_model
from uno_rt.training import predict,metrics

def main():
 torch.set_num_threads(2)
 model,info,meta=load_model(OUT/'im_public.uimc',target='im')
 rows=list(csv.DictReader((OUT/'test-predictions.tsv').open(),delimiter='\t'))
 # Parity on length/charge-diverse real observations, including all charge states.
 selected=sorted(rows,key=lambda r:hashlib.sha256((r['sequence']+r['charge']).encode()).digest())[:1000]
 oracle=predict(model,[(r['sequence'],int(r['charge'])) for r in selected],info)
 rust=np.array([float(r['predicted_im']) for r in selected]);delta=float(np.max(np.abs(oracle-rust)))
 assert delta<2e-5,delta
 obs={(r['sequence'],int(r['charge'])):float(r['im']) for r in csv.DictReader((OLD/'inputs/calibration.tsv').open(),delimiter='\t')}
 cnn={(r['sequence'],int(r['charge'])):float(r['predicted_im']) for r in csv.DictReader((OUT/'A1-predictions.tsv').open(),delimiter='\t')}
 tree={(r[0],int(r[1])):float(r[2]) for r in csv.reader((OLD/'prediction.tsv').open(),delimiter='\t')}
 diag={}
 for name,pred in [('cnn',cnn),('tree',tree)]:
  byz=defaultdict(list)
  for k,o in obs.items():
   if k in pred:byz[k[1]].append((o,pred[k]))
  diag[name]={}
  for z,v in byz.items():
   sc=statistics.median(o/p for o,p in v);err=[abs(o-p*sc)/o for o,p in v]
   diag[name][z]=dict(n=len(v),scale=sc,median_absolute_relative_error=statistics.median(err))
 # Order ablation with exactly preserved composition, length and terminal residues.
 queries=[(r['sequence'],int(r['charge'])) for r in selected]
 reversed_internal=[(s[0]+s[1:-1][::-1]+s[-1],z) for s,z in queries]
 pp=predict(model,reversed_internal,info);y=np.array([float(r['observed_im']) for r in selected])
 result=dict(public_test=json.loads((OUT/'test-metrics.json').read_text()),rust_torch_max_abs_error=delta,
  order_diagnostic=dict(n=len(selected),original=metrics(y,oracle),internal_sequence_reversed=metrics(y,pp),median_prediction_change=float(np.median(np.abs(oracle-pp)))),
  A1_same_run_affine_not_independent_accuracy=diag,
  A1_note='Multiplicative calibration and error summaries on the same diagnostic A1 observations; optimistic calibration diagnostic, not independent accuracy. No A1 labels train or select the model; public sequence overlap is permitted across independent acquisition sources. A1 is excluded from searches.',
  architecture=info,training_report=meta['report'])
 (OUT/'diagnostics.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
 print(json.dumps({k:v for k,v in result.items() if k not in ['training_report']},indent=2))
if __name__=='__main__':main()
