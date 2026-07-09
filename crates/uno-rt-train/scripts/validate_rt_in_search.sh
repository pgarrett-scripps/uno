#!/usr/bin/env bash
# validate_rt_in_search.sh — judge an RT model by what actually matters:
# HONEST in-search protein IDs on run10 + PWL-calibrated RT MAE on run10 peptides.
#
# WHY THIS EXISTS (read before trusting any RT model):
#   ProteomeTools holdout MAE is NOT a proxy for search performance. The URTC CNN
#   beat DeepLC ~2x on holdout yet lost ~2x in search, because it overfits synthetic
#   iRT and does not generalize to a real LC gradient. See memory
#   uno-rt-search-validation. ALWAYS validate a new RT model with this harness.
#
#   Also: an entrapment RT library MUST cover yeast peptides. A human+ecoli-only
#   library silently drops yeast from RT gating -> fake 0% FDR. This harness predicts
#   the FULL entrapment candidate set (incl. yeast) so the FDR is honest.
#
# Usage:
#   validate_rt_in_search.sh <model.urtc|model.gbdt>
# Env overrides (defaults are the run10 benchmark in /tmp/cmp):
#   MATRIX, ENTRAP_FASTA, CANDS (all entrap candidate peptides, 1 per line, incl yeast),
#   RTCALIB (peptide<TAB>observed_rt confident set, from UNO_DUMP_RT_CALIB), CONFIG, UNO, TRAIN
set -euo pipefail

MODEL=${1:?usage: validate_rt_in_search.sh <model>}
REPO=$(cd "$(dirname "$0")/../../.." && pwd)
UNO=${UNO:-$REPO/target/release/uno}
TRAIN=${TRAIN:-$REPO/target/release/uno-rt-train}
MATRIX=${MATRIX:-/tmp/cmp/run10_uno_shape.tsv}
ENTRAP_FASTA=${ENTRAP_FASTA:-/tmp/entrap_db.fasta}
CANDS=${CANDS:-/tmp/cmp/entrap_cands_all.txt}      # 1.26M peps incl. yeast
RTCALIB=${RTCALIB:-/tmp/cmp/rtcalib.tsv}           # observed RT of confident peptides
CONFIG=${CONFIG:-/tmp/cmp/rtexp.toml}
# entrapment scaling: real DB = human+ecoli proteins, decoy reporter = yeast
NR=${NR:-24819}; NY=${NY:-6068}
WORK=${WORK:-/tmp/cmp}
tag=$(basename "$MODEL" | sed 's/\.[^.]*$//')
LIB=$WORK/val_${tag}_lib.tsv
OUT=$WORK/val_${tag}_prot.tsv
LOG=$WORK/val_${tag}.log

echo "== predicting $CANDS with $MODEL (batched) =="
"$TRAIN" --labels "$RTCALIB" --out "$MODEL" --predict-seqs "$CANDS" > "$LIB" 2>/dev/null
echo "   library rows: $(wc -l < "$LIB")"

echo "== uno entrapment search =="
RUST_LOG=info "$UNO" --intensity-matrix "$MATRIX" --fasta "$ENTRAP_FASTA" \
  --rt-library "$LIB" --config "$CONFIG" --out-proteins "$OUT" > "$LOG" 2>&1 || { tail -5 "$LOG"; exit 1; }
grep -iE "two-stage RT|proteins:" "$LOG" | sed 's/^/   /'

awk -F'\t' -v Nr="$NR" -v Ny="$NY" '
  NR>1 && $2==0 && $8<=0.01 { if ($1~/YEAST/) y++; else r++ }
  END { printf "   HONEST: real=%d  yeast=%d  emp_FDR=%.1f%%\n", r, y, (r>0 ? y*(Nr/Ny)/r*100 : 0) }' "$OUT"

echo "== PWL-calibrated RT MAE on run10 confident peptides =="
python3 - "$LIB" "$RTCALIB" <<'PY'
import sys, numpy as np
lib, calib = sys.argv[1], sys.argv[2]
obs={s:float(r) for s,r in (l.rstrip().split('\t') for l in open(calib) if '\t' in l)}
pred={}
for l in open(lib):
    p=l.rstrip().split('\t')
    if len(p)>=2:
        try: pred[p[0]]=float(p[1])
        except: pass
peps=[s for s in obs if s in pred]
if len(peps) < 100: print(f"   too few overlap ({len(peps)})"); sys.exit()
y=np.array([obs[s] for s in peps]); x=np.array([pred[s] for s in peps])
o=np.argsort(x); xs=x[o]; ys=y[o]; nb=60; knx=[]; kny=[]
for b in range(nb):
    lo=b*len(xs)//nb; hi=max((b+1)*len(xs)//nb, lo+1)
    knx.append(np.median(xs[lo:hi])); kny.append(np.median(ys[lo:hi]))
vals=[]; wts=[]
for v in kny:
    vals.append(v); wts.append(1.0)
    while len(vals)>1 and vals[-2]>vals[-1]:
        v2,w2=vals.pop(),wts.pop(); v1,w1=vals.pop(),wts.pop()
        vals.append((v1*w1+v2*w2)/(w1+w2)); wts.append(w1+w2)
ky=[]
for v,w in zip(vals,wts): ky += [v]*int(w)
kny=np.array(ky[:len(knx)]); knx=np.array(knx)
res=y-np.interp(x,knx,kny)
fr=[np.mean(np.abs(res)<=g)*100 for g in (0.5,1,2)]
print(f"   n={len(peps)}  PWL-MAE={np.mean(np.abs(res)):.3f} min  within[0.5,1,2]min={fr[0]:.0f}%,{fr[1]:.0f}%,{fr[2]:.0f}%")
print("   reference: DeepLC generic = 0.76 MAE / 89% within 1min / 827 honest IDs")
PY
