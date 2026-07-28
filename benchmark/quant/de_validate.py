#!/usr/bin/env python3
"""Downstream differential-expression test on uno's protein abundance matrix.

uno outputs a protein-group x sample abundance matrix (`--out-abundance`) and stops
there — it is an ID + quant tool, not a stats framework. This script is the kind of
downstream analysis a real user would run (MSstats/limma-style): it takes uno's
matrix plus an experimental design and calls what changed. It is a VALIDATION
harness (it also scores against the IonStar species ground truth), not part of uno.

Usage:
  de_validate.py abundance.tsv --dilution A:1,B:2,C:3,D:4,E:5   # monotonic trend
  de_validate.py abundance.tsv --contrast A,E                    # pairwise

Design is inferred from run-name substrings (Condition letter). No design file.
"""
import sys, argparse, re
import numpy as np
from scipy.stats import spearmanr, ttest_ind


def bh(pvals):
    """Benjamini-Hochberg q-values."""
    p = np.asarray(pvals)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    running = 1.0
    for rank in range(n - 1, -1, -1):
        i = order[rank]
        running = min(running, p[i] * n / (rank + 1))
        q[i] = running
    return q


def load(path):
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        runs = hdr[4:]
        accs, M = [], []
        for ln in f:
            r = ln.rstrip("\n").split("\t")
            accs.append(r[0])
            M.append([float(x) for x in r[4:]])
    return accs, runs, np.array(M)


def cond_of(run):
    m = re.search(r"[_/]([A-Z])[_/]", run) or re.search(r"ecoli_([A-E])_", run) or re.search(r"Condition_([A-Z])", run)
    return m.group(1) if m else "?"


def species(acc):
    a = acc.upper()
    return "ECOLI" if "ECOLI" in a else ("HUMAN" if "HUMAN" in a else "?")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("matrix")
    ap.add_argument("--dilution")           # cond:level,... (order only)
    ap.add_argument("--contrast")           # A,B
    ap.add_argument("--fdr", type=float, default=0.05)
    a = ap.parse_args()

    accs, runs, M = load(a.matrix)
    conds = np.array([cond_of(r) for r in runs])
    # median-normalize each run (log2), impute missing at per-run 5th pct
    L = np.where(M > 0, np.log2(M, where=M > 0), np.nan)
    floor = np.array([np.nanpercentile(L[:, j], 5) if np.isfinite(L[:, j]).any() else 0 for j in range(L.shape[1])])
    L = np.where(np.isnan(L), floor[None, :], L)
    L = L - np.nanmedian(L, 0) + np.nanmean(np.nanmedian(L, 0))

    pvals, stats = [], []
    if a.dilution:
        levels = {p.split(":")[0]: float(p.split(":")[1]) for p in a.dilution.split(",")}
        runlvl = np.array([levels.get(c, np.nan) for c in conds])
        keep = np.isfinite(runlvl)
        for i in range(len(accs)):
            rho, p = spearmanr(runlvl[keep], L[i, keep])
            pvals.append(1.0 if np.isnan(p) else p)
            stats.append(0.0 if np.isnan(rho) else rho)
        label = "rho"
    else:
        ca, cb = a.contrast.split(",")
        A = conds == ca; B = conds == cb
        for i in range(len(accs)):
            t, p = ttest_ind(L[i, B], L[i, A], equal_var=False)
            pvals.append(1.0 if np.isnan(p) else p)
            stats.append(float(np.mean(L[i, B]) - np.mean(L[i, A])))
        label = "log2fc"
    pvals = np.array(pvals); stats = np.array(stats)
    q = bh(pvals)

    # ground-truth scoring (ecoli = true DE, human = flat)
    sp = np.array([species(x) for x in accs])
    ne = (sp == "ECOLI").sum(); nh = (sp == "HUMAN").sum()
    print(f"{len(accs)} protein groups: {ne} ecoli, {nh} human")
    for thr in [a.fdr, 0.10, 0.20]:
        sig = q <= thr
        # spike-tracking direction (positive) for dilution; both signs for contrast
        d = sig & (stats > 0) if a.dilution else sig
        e = (d & (sp == "ECOLI")).sum(); h = (d & (sp == "HUMAN")).sum()
        print(f"  FDR<={thr}: ecoli {e}/{ne} ({100*e/max(ne,1):.0f}%)  human FP {h}  precision {100*e/max(e+h,1):.0f}%  ({label} dir)")


if __name__ == "__main__":
    main()
