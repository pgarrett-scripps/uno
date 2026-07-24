#!/usr/bin/env python3
"""Build a koth-schema `intensity_matrix.tsv` from per-run koth features.parquet.

koth_align needs hills (deleted with the raw data), so this is a lightweight
stand-in: cross-run consensus by (charge, mass ppm, RT window), detected-only
(no MBR fill). Output columns match koth's writer exactly so uno's ingest is
identical to the real koth_align output:

  massCalib  mz  charge  rtApex  im  combined_score  seed_run  n_contributing_runs  <run1> <run2> ...

Also writes `runs.tsv`: column-order run names + condition labels for the
downstream quant stats.
"""
import sys, glob, os, re
import numpy as np
import pyarrow.parquet as pq

BATCH_DATA = "/home/patrick-garrett/Repos/koth_rust/benchmark/data/koth_sulfur_on.pre_sigma4_20260709_121922"
PPM = 8.0
RT_TOL = 0.5   # minutes; 2 h gradient, replicate runs -> stable RT

def cond_of(name):
    m = re.search(r"human_ecoli_([A-E])_", name)
    return m.group(1) if m else "?"

def load_run(pdir):
    t = pq.read_table(pdir, columns=["massCalib","mz","charge","rtApex","intensityApex","combined_score","im"])
    d = t.to_pydict()
    n = len(d["massCalib"])
    return {k: np.asarray(d[k], dtype=float) for k in ["massCalib","mz","charge","rtApex","intensityApex","combined_score","im"]}, n

def main():
    conds = sys.argv[1].split(",") if len(sys.argv) > 1 else ["A","E"]
    reps  = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    out   = sys.argv[3] if len(sys.argv) > 3 else "consensus"
    # pick runs
    dirs = sorted(glob.glob(f"{BATCH_DATA}/*.mzML"))
    chosen = []
    for c in conds:
        cd = [d for d in dirs if cond_of(os.path.basename(d)) == c][:reps]
        chosen += [(os.path.basename(d), c, os.path.join(d, "features.parquet")) for d in cd]
    print(f"runs ({len(chosen)}):")
    for nm,c,_ in chosen: print(f"  {c}  {nm[:40]}")

    # gather all features tagged by run index
    mass=[]; mz=[]; ch=[]; rt=[]; inten=[]; cs=[]; im=[]; run=[]
    for ri,(nm,c,pf) in enumerate(chosen):
        d,n = load_run(pf)
        mass.append(d["massCalib"]); mz.append(d["mz"]); ch.append(d["charge"])
        rt.append(d["rtApex"]); inten.append(d["intensityApex"]); cs.append(d["combined_score"]); im.append(d["im"])
        run.append(np.full(n, ri))
        print(f"  loaded {nm[:30]}: {n} features")
    mass=np.concatenate(mass); mz=np.concatenate(mz); ch=np.concatenate(ch).astype(int)
    rt=np.concatenate(rt); inten=np.concatenate(inten); cs=np.concatenate(cs); im=np.concatenate(im); run=np.concatenate(run)
    R=len(chosen)

    # cluster within charge -> mass consecutive-gap -> rt consecutive-gap
    clusters=[]  # list of index arrays
    for z in np.unique(ch):
        idx=np.where(ch==z)[0]
        o=idx[np.argsort(mass[idx])]
        m=mass[o]
        # new mass-group when consecutive relative gap > PPM
        gap=np.empty(len(m)); gap[0]=1e9; gap[1:]=(m[1:]-m[:-1])/m[1:]*1e6
        mgrp=np.cumsum(gap>PPM)
        for g in np.unique(mgrp):
            sel=o[mgrp==g]
            # split by RT gaps
            so=sel[np.argsort(rt[sel])]
            r=rt[so]
            if len(r)==1:
                clusters.append(so); continue
            rgap=np.empty(len(r)); rgap[0]=1e9; rgap[1:]=r[1:]-r[:-1]
            rgrp=np.cumsum(rgap>RT_TOL)
            for rg in np.unique(rgrp):
                clusters.append(so[rgrp==rg])

    # emit
    os.makedirs(out, exist_ok=True)
    mpath=os.path.join(out,"intensity_matrix.tsv")
    names=[nm for nm,_,_ in chosen]
    with open(mpath,"w") as fh:
        fh.write("massCalib\tmz\tcharge\trtApex\tim\tcombined_score\tseed_run\tn_contributing_runs\t"+"\t".join(names)+"\n")
        n_written=0
        for cl in clusters:
            iw=inten[cl]
            # per-run intensity (max within run if multiple)
            row=np.zeros(R)
            for j in cl:
                r=run[j]
                if inten[j]>row[r]: row[r]=inten[j]
            present=int((row>0).sum())
            seed=int(np.argmax(row))
            wsum=iw.sum() if iw.sum()>0 else 1.0
            cmass=float((mass[cl]*iw).sum()/wsum)
            cmz=float((mz[cl]*iw).sum()/wsum)
            crt=float((rt[cl]*iw).sum()/wsum)
            cz=int(np.bincount(ch[cl]).argmax())
            ccs=float(cs[cl].max())
            cim=im[cl][np.argmax(iw)]
            imstr=f"{cim:.6f}" if cim and cim>0 else ""
            cells="\t".join(f"{v:.5e}" if v>0 else "0" for v in row)
            fh.write(f"{cmass:.6f}\t{cmz:.6f}\t{cz}\t{crt:.6f}\t{imstr}\t{ccs:.6f}\t{seed}\t{present}\t{cells}\n")
            n_written+=1
    with open(os.path.join(out,"runs.tsv"),"w") as fh:
        fh.write("col\trun\tcondition\n")
        for i,(nm,c,_) in enumerate(chosen): fh.write(f"{i}\t{nm}\t{c}\n")
    print(f"wrote {n_written} consensus features x {R} runs -> {mpath}")

if __name__=="__main__":
    main()
