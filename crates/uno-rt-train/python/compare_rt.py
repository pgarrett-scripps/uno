#!/usr/bin/env python3
"""Head-to-head RT comparison. Usage:
   compare_rt.py TRUE.tsv name=pred.tsv [name=pred.tsv ...]
Each pred file is seq\tpred. Reports raw + affine-aligned MAE/dt95/R2 plus
Pearson/Spearman (scale-invariant) against the true iRT, on the shared seq set.
"""
import sys
import numpy as np

def load(path):
    d = {}
    for line in open(path):
        a = line.rstrip("\n").split("\t")
        d[a[0]] = float(a[1])
    return d

true = load(sys.argv[1])
models = [(arg.split("=", 1)[0], load(arg.split("=", 1)[1])) for arg in sys.argv[2:]]

def rank(x):
    o = np.argsort(np.argsort(x))
    return o.astype(float)

def metrics(y, p):
    e = y - p
    r2 = 1 - np.sum(e * e) / np.sum((y - y.mean()) ** 2)
    return np.mean(np.abs(e)), np.percentile(np.abs(e), 95), r2

print(f"{'model':14s} {'n':>6s} | {'rawMAE':>7s} {'rawdt95':>8s} | "
      f"{'aMAE':>6s} {'adt95':>6s} {'aR2':>7s} | {'Pear':>6s} {'Spear':>6s}")
print("-" * 88)
for name, pred in models:
    keys = [k for k in true if k in pred]
    y = np.array([true[k] for k in keys])
    p = np.array([pred[k] for k in keys])
    rmae, rdt95, _ = metrics(y, p)
    # affine align p -> y
    A = np.vstack([p, np.ones_like(p)]).T
    a, b = np.linalg.lstsq(A, y, rcond=None)[0]
    pa = a * p + b
    amae, adt95, ar2 = metrics(y, pa)
    pear = np.corrcoef(y, p)[0, 1]
    spear = np.corrcoef(rank(y), rank(p))[0, 1]
    print(f"{name:14s} {len(keys):6d} | {rmae:7.3f} {rdt95:8.3f} | "
          f"{amae:6.3f} {adt95:6.3f} {ar2:7.4f} | {pear:6.4f} {spear:6.4f}")
