#!/usr/bin/env python3
"""Run DeepLC (pretrained, no calibration) on a seq\tirt TSV and emit seq\tpred.

DeepLC predicts on its own abstract scale; compare via affine alignment / rank
correlation downstream (rt_eval does this in Rust, or compare_rt.py in numpy).
AA-only, C = carbamidomethyl (DeepLC modification 'Carbamidomethyl' on every C).
"""
import sys, time
import pandas as pd
from deeplc import DeepLC

inp, out = sys.argv[1], sys.argv[2]
seqs, ys = [], []
with open(inp) as f:
    for line in f:
        s, y = line.rstrip("\n").split("\t")
        seqs.append(s); ys.append(float(y))

# fixed carbamidomethyl on every cysteine, no other mods
def mods(seq):
    locs = [str(i + 1) + "|Carbamidomethyl" for i, c in enumerate(seq) if c == "C"]
    return "|".join(locs)

df = pd.DataFrame({"seq": seqs, "modifications": [mods(s) for s in seqs]})
t = time.time()
dlc = DeepLC()
preds = dlc.make_preds(seq_df=df, calibrate=False)
dt = time.time() - t
print(f"DeepLC predicted {len(preds)} in {dt:.1f}s ({len(preds)/dt:.0f}/s)", file=sys.stderr)

with open(out, "w") as o:
    for s, p in zip(seqs, preds):
        o.write(f"{s}\t{p}\n")
