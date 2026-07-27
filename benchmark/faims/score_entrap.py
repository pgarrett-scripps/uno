#!/usr/bin/env python3
"""Honest species-entrapment protein FDP for uno OR ms1searchpy output.

Entrapment methodology: the sample is one species (here human, HeLa); the FASTA
adds an ABSENT entrapment species (here E. coli). Any accepted entrapment-species
protein is a false positive. Scale-corrected FDP at rank k (targets ranked by
score desc):

    FDP(k) = entrap_in_top_k * (n_target_db / n_entrap_db) / target_in_top_k

Reads both formats:
  - uno   --out-proteins: cols accession, decoy, ... , binomial_score(8), q_value
  - ms1searchpy *_proteins_full.tsv: cols dbname, score, matched, theoretical, decoy

Usage:  score_entrap.py uno   uno_proteins.tsv
        score_entrap.py ms1   LUM2_x.features_proteins_full.tsv
Env:    N_TARGET_DB / N_ENTRAP_DB (default 20193 / 4431 = human / E. coli here),
        TARGET_TAG / ENTRAP_TAG   (default HUMAN / ECOLI, matched in the accession)
"""
import sys, os

N_TARGET = int(os.environ.get("N_TARGET_DB", 20193))
N_ENTRAP = int(os.environ.get("N_ENTRAP_DB", 4431))
SCALE = N_TARGET / N_ENTRAP
TARGET_TAG = os.environ.get("TARGET_TAG", "HUMAN")
ENTRAP_TAG = os.environ.get("ENTRAP_TAG", "ECOLI")


def rows_uno(path):
    for ln in open(path).read().splitlines()[1:]:
        c = ln.split("\t")
        if c[1] != "0":  # skip decoys
            continue
        yield float(c[7]), c[0]           # binomial_score, accession


def rows_ms1(path):
    for ln in open(path).read().splitlines()[1:]:
        c = ln.split("\t")
        if c[4].strip().lower() in ("true", "1"):  # skip decoys
            continue
        yield float(c[1]), c[0]           # score, dbname


def curve(rows, thresholds=(0.01, 0.05)):
    rows = sorted(rows, reverse=True)
    tgt = ent = 0
    best = {t: 0 for t in thresholds}
    for _, acc in rows:
        a = acc.upper()
        if ENTRAP_TAG in a:
            ent += 1
        elif TARGET_TAG in a:
            tgt += 1
        else:
            continue
        fdp = ent * SCALE / max(tgt, 1)
        for t in thresholds:
            if fdp <= t:
                best[t] = tgt
    return best


def main():
    fmt, path = sys.argv[1], sys.argv[2]
    rows = rows_uno(path) if fmt == "uno" else rows_ms1(path)
    b = curve(list(rows))
    print(f"{path}: {TARGET_TAG} @1%={b[0.01]}  @5%={b[0.05]}  (scale {SCALE:.2f})")


if __name__ == "__main__":
    main()
