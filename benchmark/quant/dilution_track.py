#!/usr/bin/env python3
"""Dilution-tracking validation for uno's protein abundance matrix.

Given uno's --out-abundance matrix over a spike-in dilution series, score how
well each protein's abundance tracks the spike amount, by species. On a clean
result the SPIKE species correlates positively with the dilution while the
constant background and any absent (entrapment) species stay flat.

Used for the DirectMS1Quant UPS1-into-E.coli series (PXD033350): UPS1 = spike,
E. coli = flat background, human = known-absent entrapment.

  python dilution_track.py abundance.tsv

Run columns are identified by the "UPS_<level>_" substring; the level->fmol map
below is DESCENDING (UPS_1 = 50 fmol, the highest; UPS_8 = 0.1 fmol, the lowest),
established empirically from spike-peptide signal.
"""
import sys, re, math
from statistics import mean, median

FMOL = {1: 50, 2: 25, 3: 10, 4: 5, 5: 2.5, 6: 1, 7: 0.25, 8: 0.1}


def rankdata(a):
    order = sorted(range(len(a)), key=lambda i: a[i]); r = [0.0] * len(a); i = 0
    while i < len(a):
        j = i
        while j + 1 < len(a) and a[order[j + 1]] == a[order[i]]:
            j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r


def spearman(x, y):
    rx, ry = rankdata(x), rankdata(y); n = len(x)
    mx, my = mean(rx), mean(ry)
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    dx = math.sqrt(sum((v - mx) ** 2 for v in rx)); dy = math.sqrt(sum((v - my) ** 2 for v in ry))
    return num / (dx * dy) if dx * dy else 0.0


def species(acc):
    """UPS1 / E. coli / human, by the deposit's own accession tags.

    The UPS1 spike-ins are tagged `<accession>ups|<NAME>_HUMAN_UPS`. Testing for a
    bare "UPS" substring after case-folding, which this did, also matches ordinary
    human accessions that merely contain those three letters (`sp|Q9UPS8|ANR26_HUMAN`).
    No such accession is present in this dataset, so no published number changes,
    but the test has to be on the tag rather than on the letters.
    """
    if "ups|" in acc or acc.upper().endswith("_UPS"):
        return "UPS"
    return "ECOLI" if "ECOLI" in acc.upper() else "HUMAN"


def main():
    path = sys.argv[1]
    hdr = open(path).readline().rstrip("\n").split("\t")
    run_cols = [(i, c) for i, c in enumerate(hdr) if "UPS_" in c]
    logf = [math.log10(FMOL[int(re.search(r"UPS_(\d)_", c).group(1))]) for _, c in run_cols]

    stats = {"UPS": [], "ECOLI": [], "HUMAN": []}
    for ln in open(path).read().splitlines()[1:]:
        c = ln.split("\t")
        vals = [float(c[i]) for i, _ in run_cols]
        det = [(v, f) for v, f in zip(vals, logf) if v > 0]
        if len(det) < 8:
            continue
        rho = spearman([f for _, f in det], [v for v, _ in det])
        stats[species(c[0])].append(rho)

    print("protein dilution response (Spearman rho of abundance vs log10 fmol, >=8 runs):")
    for sp in ("UPS", "ECOLI", "HUMAN"):
        r = stats[sp]
        if not r:
            print(f"  {sp}: 0 proteins"); continue
        # true median (mean of the two middle values when n is even). The old
        # `sorted(r)[len(r) // 2]` returned the upper-middle order statistic,
        # which reported +0.74 for the n=16 UPS set whose median is +0.71.
        med = median(r)
        frac = sum(1 for x in r if x > 0.5) / len(r)
        print(f"  {sp:6s} n={len(r):3d}  median_rho={med:+.2f}  mean={mean(r):+.2f}  frac(rho>0.5)={frac:.0%}")


if __name__ == "__main__":
    main()
