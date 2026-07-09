#!/usr/bin/env python3
"""chronologer_to_labels.py — convert the Chronologer harmonized RT corpus into
uno's RT training format (SEQUENCE<TAB>RT), in Hydrophobic-Index (HI, %ACN) space.

WHY: uno's URTC CNN was trained only on ProteomeTools/Prosit synthetic iRT and does
NOT generalize to real LC gradients (DeepLC 0.76 vs CNN 2.65 min PWL-MAE on run10 —
see memory uno-rt-search-validation). Chronologer's corpus is 2.25M peptides from 11
real-LC community datasets harmonized into a single physicochemical scale (HI = %ACN
at elution), which is exactly the diversity + common-scale needed to generalize.

DATA: searlelab/chronologer (Apache-2.0), file data/Chronologer_DB_220308.gz.
  git clone https://github.com/searlelab/chronologer
  python chronologer_to_labels.py chronologer/data/Chronologer_DB_220308.gz out_dir/

uno's CNN is AA-only (C is implicitly carbamidomethyl, uno's fixed mod). So this
keeps unmodified peptides + carbamidomethyl-C, and DROPS peptides carrying any other
modification. (To exploit Chronologer's 17 mod types you must first extend the URTC
featurizer to accept mods — a separate, larger change.)

Peptide encoding in the DB is EncyclopeDIA-style, e.g. ALSVLGC[+57.021464]GHTSSTK or
[42.010565]KGSPTPGFSTR. We detect columns by name (PeptideModSeq + an HI/RT column);
verify on first run with --peek.

Output: out_dir/{train,val,holdout}.tsv  (peptide-disjoint, deterministic hash split).
The HI label is a real common scale — uno's per-run PWL/RANSAC calibration then maps
HI -> observed minutes per run, same as it does for DeepLC/iRT today.
"""
import sys, os, gzip, re, hashlib, argparse

CAM_C = re.compile(r'C\[\+?57\.0214\d*\]')          # carbamidomethyl-C (keep the C)
ANY_MOD = re.compile(r'\[[^\]]*\]')                  # any [..] mass/name annotation
NTERM_MOD = re.compile(r'^\[[^\]]*\]')               # leading n-term mod

def clean_seq(modseq):
    """Return bare AA sequence if peptide has only (optional) carbamidomethyl-C,
    else None (carries a mod the AA-only CNN can't represent)."""
    s = CAM_C.sub('C', modseq)        # strip the fixed Cam mass, keep C
    if ANY_MOD.search(s):
        return None                   # some OTHER modification remains -> drop
    s = s.strip()
    if not s.isalpha() or not s.isupper():
        return None
    if not (6 <= len(s) <= 50):
        return None
    return s

def open_any(p):
    return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)

def split_of(seq, val=0.08, holdout=0.04):
    h = int(hashlib.md5(seq.encode()).hexdigest()[:8], 16) / 0xffffffff
    if h < holdout: return 'holdout'
    if h < holdout + val: return 'val'
    return 'train'

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('db'); ap.add_argument('outdir')
    ap.add_argument('--seq-col', default=None, help='override peptide column name')
    ap.add_argument('--rt-col', default=None, help='override HI/RT column name')
    ap.add_argument('--peek', action='store_true', help='print header + 3 rows and exit')
    a = ap.parse_args()

    f = open_any(a.db)
    header = f.readline().rstrip('\n')
    # Chronologer DB is tab- or comma-separated; sniff.
    sep = '\t' if header.count('\t') >= header.count(',') else ','
    cols = header.split(sep)
    if a.peek:
        print('sep=%r' % sep); print('columns:', cols)
        for _ in range(3): print(f.readline().rstrip())
        return
    def find(cands, override):
        if override: return cols.index(override)
        for i, c in enumerate(cols):
            cl = c.lower()
            if any(k in cl for k in cands): return i
        return None
    si = find(['peptidemodseq', 'modseq', 'peptide', 'sequence'], a.seq_col)
    ri = find(['hi', 'hydrophob', 'rt', 'tr', 'retention'], a.rt_col)
    if si is None or ri is None:
        sys.exit(f'could not find seq/rt columns in {cols}; use --seq-col/--rt-col')
    print(f'using seq col {cols[si]!r}, rt col {cols[ri]!r}, sep {sep!r}', file=sys.stderr)

    os.makedirs(a.outdir, exist_ok=True)
    outs = {k: open(os.path.join(a.outdir, f'{k}.tsv'), 'w') for k in ('train','val','holdout')}
    n_in = n_kept = n_dropmod = 0
    seen = {}   # seq -> list of rt (dedup later by mean)
    for line in f:
        n_in += 1
        parts = line.rstrip('\n').split(sep)
        if len(parts) <= max(si, ri): continue
        seq = clean_seq(parts[si])
        if seq is None:
            n_dropmod += 1; continue
        try: rt = float(parts[ri])
        except ValueError: continue
        seen.setdefault(seq, []).append(rt)
    for seq, rts in seen.items():
        rt = sum(rts) / len(rts)        # average replicates
        outs[split_of(seq)].write(f'{seq}\t{rt}\n')
        n_kept += 1
    for o in outs.values(): o.close()
    print(f'rows in={n_in} dropped(mod)={n_dropmod} unique-peptides-kept={n_kept}', file=sys.stderr)
    print(f'wrote {a.outdir}/train.tsv val.tsv holdout.tsv', file=sys.stderr)

if __name__ == '__main__':
    main()
