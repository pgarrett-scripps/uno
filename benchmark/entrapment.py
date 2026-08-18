#!/usr/bin/env python3
"""Unified yeast-entrapment protein FDR — reads uno OR ms1searchpy protein output.

Entrapment methodology: the sample is human (HeLa); the FASTA is human + yeast.
Yeast is absent from the sample, so any accepted _YEAST protein is a false
positive. Scale-corrected entrapment FDP at rank k (targets ranked by score desc):

    FDP(k) = yeast_in_top_k * (n_human_db / n_yeast_db) / human_in_top_k

q-value = running-min of FDP scanning from the deepest cutoff upward. Reported
depth at a threshold = the most human proteins whose q <= threshold.

Schemas auto-detected:
  uno : columns `accession`, `decoy` (0/1), `binomial_score`
  ms1 : columns `dbname`,    `decoy` (True/False), `score`

Usage: entrapment.py LABEL1 path1 [LABEL2 path2 ...]
Env:   N_HUMAN_DB / N_YEAST_DB override the DB sizes (defaults 20647 / 6066).
"""
import sys, csv, os

N_HUMAN_DB = int(os.environ.get("N_HUMAN_DB", 20647))
N_YEAST_DB = int(os.environ.get("N_YEAST_DB", 6066))
SCALE = N_HUMAN_DB / N_YEAST_DB


def load(path):
    """Return list of (score, is_yeast) for target human/yeast proteins."""
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        cols = r.fieldnames or []
        # detect schema
        if "binomial_score" in cols:
            acc_c, score_c = "accession", "binomial_score"
        elif "score" in cols and "dbname" in cols:
            acc_c, score_c = "dbname", "score"
        else:
            raise SystemExit(f"unrecognized schema in {path}: {cols}")
        rows = []
        for row in r:
            if str(row.get("decoy", "0")) not in ("0", "False", "false"):
                continue
            name = row[acc_c]
            is_y, is_h = "_YEAST" in name, "_HUMAN" in name
            if not (is_y or is_h):
                continue
            try:
                s = float(row[score_c])
            except (ValueError, KeyError):
                continue
            rows.append((s, is_y))
        return rows


def curve(rows, thresholds=(0.01, 0.05)):
    # Ties are broken ENTRAPMENT-FIRST, deterministically.
    #
    # Two things were wrong with the previous `key=lambda r: -r[0]`. Python's sort
    # is stable, so equal scores kept their order in the input FILE -- and the two
    # engines write their protein tables in different orders, so the same evidence
    # gave different curves. uno's binomial protein score is quantized, so ties are
    # common rather than exotic.
    #
    # Sorting entrapment hits ahead of true ones within a tie makes the estimate
    # CONSERVATIVE at every tie boundary: a tied block is charged its false
    # positives before it is credited its true ones, so the reported depth at a
    # given FDP can only be understated. The alternative (`r[1]`, target-first) is
    # equally deterministic but optimistic at every tie, which is the wrong
    # direction for a paper whose claim is that its FDR is honest. Applied
    # identically to both engines either way.
    rows.sort(key=lambda r: (-r[0], not r[1]))
    human = yeast = 0
    steps = []
    for _, is_y in rows:
        if is_y:
            yeast += 1
        else:
            human += 1
        if human:
            steps.append((human, (yeast * SCALE) / human))
    q = [0.0] * len(steps)
    rmin = float("inf")
    for i in range(len(steps) - 1, -1, -1):
        rmin = min(rmin, steps[i][1])
        q[i] = rmin
    out = {}
    for t in thresholds:
        best = 0
        for (h, _), qi in zip(steps, q):
            if qi <= t and h > best:
                best = h
        out[t] = best
    return out, human, yeast


if __name__ == "__main__":
    args = sys.argv[1:]
    for label, path in zip(args[0::2], args[1::2]):
        rows = load(path)
        out, th, ty = curve(rows)
        print(f"=== {label} ({os.path.basename(path)}) ===")
        print(f"  targets scored: human={th} yeast={ty}")
        for t, n in out.items():
            print(f"  depth @{int(t*100)}% entrapment FDR: {n} human proteins")
