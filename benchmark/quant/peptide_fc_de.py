#!/usr/bin/env python3
"""DirectMS1Quant-style peptide-fold-change differential test on uno's LFQ matrix.

uno emits a per-precursor LFQ matrix (--out-lfq); this is the kind of downstream
test a user runs on it (it is NOT part of uno). For a hi-vs-lo dilution contrast,
per protein it collects the per-peptide, per-replicate-pair log2 fold changes
(so n = peptides x replicate-pairs, as in DirectMS1Quant's n_arr), then a
one-sample t-test vs 0, a fold-change gate, and Benjamini-Hochberg at 1%.

Scores by species for the DirectMS1Quant UPS1-into-E.coli series (PXD033350):
UPS = spike (true positives), E. coli = flat (false positives), human = absent.

  python peptide_fc_de.py lfq.tsv 50 0_1        # hi=50 fmol vs lo=0.1 fmol

A peptide present in >=2 hi replicates but missing at lo is a genuine off event;
its lo cell is left-censored to a detection floor (1st-percentile intensity) so
the pair contributes a large fold change instead of being dropped.
"""
import sys, re, math
from statistics import mean, pstdev

UPS_OF = {"50": 1, "25": 2, "10": 3, "5": 4, "2_5": 5, "1": 6, "0_25": 7, "0_1": 8}


def betacf(a, b, x):
    qab, qap, qam = a + b, a + 1, a - 1
    c = 1.0; d = 1 - qab * x / qap
    d = 1 / (d if abs(d) > 1e-30 else 1e-30); h = d
    for m in range(1, 400):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1 / ((1 + aa * d) if abs(1 + aa * d) > 1e-30 else 1e-30); c = 1 + aa / c
        if abs(c) < 1e-30: c = 1e-30
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1 / ((1 + aa * d) if abs(1 + aa * d) > 1e-30 else 1e-30); c = 1 + aa / c
        if abs(c) < 1e-30: c = 1e-30
        de = d * c; h *= de
        if abs(de - 1) < 3e-12: break
    return h


def betai(a, b, x):
    if x <= 0: return 0.0
    if x >= 1: return 1.0
    bt = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log(1 - x))
    return bt * betacf(a, b, x) / a if x < (a + 1) / (a + b + 2) else 1 - bt * betacf(b, a, 1 - x) / b


def species(a):
    a = a.upper()
    return "UPS" if "UPS" in a else ("ECOLI" if "ECOLI" in a else "HUMAN")


def main():
    path, hi, lo = sys.argv[1], sys.argv[2], sys.argv[3]
    gate = 1.0  # |mean log2 FC| >= 1 (2-fold) to enter the BH set, DirectMS1Quant-style FC_pass
    hdr = open(path).readline().rstrip("\n").split("\t")

    def repmap(n):
        return {m.group(1): i for i, c in enumerate(hdr) for m in [re.search(rf"UPS_{n}_(\d\d)$", c)] if m}

    hm, lm = repmap(UPS_OF[hi]), repmap(UPS_OF[lo])
    pairs = [(hm[r], lm[r]) for r in hm if r in lm]
    lines = open(path).read().splitlines()[1:]

    allcols = [i for i, c in enumerate(hdr) if "UPS_" in c]
    allv = sorted(float(c.split("\t")[i]) for c in lines for i in allcols if float(c.split("\t")[i]) > 0)
    floor = allv[len(allv) // 100]  # 1st percentile detection floor

    prot = {}
    for ln in lines:
        c = ln.split("\t")
        if c[2] != "0":  # skip decoys
            continue
        acc = c[3].split(";")[0]
        his = [float(c[hi_i]) for hi_i, _ in pairs]
        if sum(1 for x in his if x > 0) < 2:  # need a confident up-state
            continue
        for hi_i, lo_i in pairs:
            a, b = float(c[hi_i]), float(c[lo_i])
            if a <= 0:
                continue
            b = b if b > 0 else floor  # left-censor the off event
            prot.setdefault(acc, []).append(math.log2(a / b))

    rows = []
    for acc, fcs in prot.items():
        n = len(fcs)
        if n < 3:
            continue
        m = mean(fcs); s = pstdev(fcs) * math.sqrt(n / (n - 1))
        if abs(m) < gate:
            continue
        se = s / math.sqrt(n) if s > 0 else 1e-9; t = m / se
        p = betai((n - 1) / 2, 0.5, (n - 1) / ((n - 1) + t * t))
        rows.append((p, species(acc)))

    rows.sort(); M = len(rows); thr = 0
    for k, (p, _) in enumerate(rows, 1):
        if p <= 0.01 * k / M:
            thr = k
    u = e = o = tu = te = 0
    for k, (p, s) in enumerate(rows):
        tu += s == "UPS"; te += s == "ECOLI"
        if k < thr:
            u += s == "UPS"; e += s == "ECOLI"; o += s == "HUMAN"
    tot = u + e + o
    print(f"  {hi}_vs_{lo}: UPS(TP)={u}  ECOLI_FP={e}  human={o}  "
          f"precision={100 * u / tot if tot else 0:.0f}%  (gated {M} of {len(prot)} prot)")


if __name__ == "__main__":
    main()
