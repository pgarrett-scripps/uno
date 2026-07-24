#!/usr/bin/env python3
"""DirectMS1Quant-style protein-level DE test, validated on human/ecoli ground truth.

Inputs: uno --out-peptides table (peptide, decoy, proteins, charge, mass, rt, <runs>),
and the ID-confident protein set (uno --out-proteins, q<=0.01 targets).

Logic (co-direction + decoy-calibrated protein FDR):
  - median-normalize each run's log2 intensities (detected cells only)
  - per peptide: log2FC = mean(log2 E) - mean(log2 A), require >=MINREP detected per condition
  - a peptide is "significant" if |log2FC| >= Z * SD(null), null = decoy-peptide FC spread
  - per protein: n = quantifiable peptides, k = max co-directed significant peptides
  - protein DE score = k; FDR from decoy proteins run through the identical pipeline
Validation: at 5% FDR, DE calls should be ~all ecoli (true DE), few human (true flat).
"""
import sys, numpy as np, re
from collections import defaultdict

PEP = sys.argv[1]; PROT = sys.argv[2]
MINREP = 2; Z = 2.0

# --- load confident target proteins (q<=0.01) ---
conf = set()
with open(PROT) as f:
    h = f.readline().rstrip("\n").split("\t"); qi = h.index("q_value"); di = h.index("decoy")
    for ln in f:
        r = ln.rstrip("\n").split("\t")
        if r[di] == "0" and float(r[qi]) <= 0.01:
            conf.add(r[0])

# --- load peptide table ---
with open(PEP) as f:
    hdr = f.readline().rstrip("\n").split("\t")
    runs = hdr[6:]
    condA = [i for i,rn in enumerate(runs) if "_A_" in rn]
    condE = [i for i,rn in enumerate(runs) if "_E_" in rn]
    print(f"runs: {len(runs)}  A={len(condA)} E={len(condE)}")
    rows = []
    for ln in f:
        r = ln.rstrip("\n").split("\t")
        inten = np.array([float(x) for x in r[6:]])
        rows.append((r[0], r[1]=="1", r[2], inten))

# --- per-run median normalization on log2 (detected cells) ---
M = np.vstack([r[3] for r in rows])
logM = np.where(M > 0, np.log2(M), np.nan)
run_med = np.nanmedian(logM, axis=0)
logM = logM - run_med + np.nanmean(run_med)   # center runs to common median

def fc(row_log):
    a = row_log[condA]; e = row_log[condE]
    na = np.sum(~np.isnan(a)); ne = np.sum(~np.isnan(e))
    if na < MINREP or ne < MINREP: return None
    return np.nanmean(e) - np.nanmean(a)

# --- per-peptide FC ---
pep_fc = []
for i,(seq,dec,prots,_) in enumerate(rows):
    v = fc(logM[i])
    if v is not None: pep_fc.append((i, dec, prots, v))

# null SD from decoy-peptide FCs
dec_fcs = np.array([v for _,dec,_,v in pep_fc if dec])
sd_null = np.std(dec_fcs)
thr = Z * sd_null
print(f"quantifiable peptides: {len(pep_fc)} (decoy {dec_fcs.size}); null SD={sd_null:.3f}, sig |log2FC|>={thr:.3f}")

# --- roll up to proteins (protein = each accession; peptide may map to several) ---
# targets: only confident proteins; decoys: group by decoy accession as null
prot_peps = defaultdict(list)   # accession -> list of (fc, is_decoy_protein)
for i, dec, prots, v in pep_fc:
    for acc in prots.split(";"):
        is_dec = acc.startswith("DECOY") or dec
        if is_dec:
            prot_peps[("D",acc)].append(v)
        elif acc in conf:
            prot_peps[("T",acc)].append(v)

def stat(fcs):
    fcs = np.array(fcs)
    n = len(fcs)
    sig = fcs[np.abs(fcs) >= thr]
    kpos = np.sum(sig > 0); kneg = np.sum(sig < 0)
    k = max(kpos, kneg)
    med = np.median(sig) if len(sig) else 0.0
    return n, k, med

records = []  # (is_target, acc, n, k, med)
for (tag,acc), fcs in prot_peps.items():
    n,k,med = stat(fcs)
    records.append((tag=="T", acc, n, k, med))

# --- decoy-calibrated FDR by descending k (tie: |med|) ---
records.sort(key=lambda r: (r[3], abs(r[4])), reverse=True)
nt=nd=0; out=[]
for is_t, acc, n, k, med in records:
    if k < 1: continue
    if is_t: nt+=1
    else: nd+=1
    fdr = nd/max(nt,1)
    out.append((is_t, acc, n, k, med, fdr))

# proteins at 5% FDR (running-max q)
q=0.0; passing=[]
for rec in reversed(out):
    q=max(q, rec[5])
    if q<=0.05 and rec[0]: passing.append(rec)
def org(acc): return "ECOLI" if "ECOLI" in acc.upper() else ("HUMAN" if "HUMAN" in acc.upper() else "?")
from collections import Counter
c=Counter(org(a) for _,a,_,_,_,_ in passing)
up=[r for r in passing if r[4]>0]; dn=[r for r in passing if r[4]<0]
print(f"\nDE proteins @5% FDR: {len(passing)}  ->  {dict(c)}")
print(f"  direction: up(E>A)={len(up)} down(E<A)={len(dn)}")
print(f"  ecoli median log2FC: {np.median([r[4] for _,a,_,_,r4,_ in [(x[0],x[1],x[2],x[3],x[4],x[5]) for x in passing] for r4 in [x[4]] if org(a)=='ECOLI']) if c.get('ECOLI') else 'NA'}")
# background: how many confident human vs ecoli proteins total (quantifiable)
tot=Counter(org(acc) for is_t,acc,_,_,_ in records if is_t)
print(f"\nquantifiable confident proteins total: {dict(tot)}")
print(f"DE sensitivity (ecoli called / ecoli total): {c.get('ECOLI',0)}/{tot.get('ECOLI',0)}")
print(f"false-DE rate (human called / human total): {c.get('HUMAN',0)}/{tot.get('HUMAN',0)}")
