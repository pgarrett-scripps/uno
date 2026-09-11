#!/usr/bin/env python3
"""Build ms1searchpy `-pl` RT calibration lists from the PXD066701 DIA-NN report.

Diagnostic only. The deposit's own DIA-NN 2.0.2 analysis (CW0710_DN_all_.zip,
report.parquet) supplies MS2-verified peptide retention times for each raw
file. ms1searchpy's documented `-pl` option calibrates DeepLC on such a list
instead of self-calibrating. uno never sees these lists.

Filters: precursor Q.Value <= 0.01, no variable modifications (UniMod:4 on C
is the fixed mod both engines use), length 7-30, no missed cleavage, median RT
per stripped sequence.

usage: build_diann_calib_lists.py REPORT.parquet OUT_DIR STEM [STEM ...]
"""
import sys, pandas as pd
rep, out, stems = sys.argv[1], sys.argv[2], sys.argv[3:]
d = pd.read_parquet(rep, columns=["Run", "Stripped.Sequence", "Modified.Sequence", "RT", "Q.Value"])
for stem in stems:
    r = d[d.Run.str.contains(stem) & (d["Q.Value"] <= 0.01)]
    r = r[~r["Modified.Sequence"].str.replace("(UniMod:4)", "", regex=False).str.contains(r"\(")]
    r = r[r["Stripped.Sequence"].str.len().between(7, 30)]
    r = r[~r["Stripped.Sequence"].str[:-1].str.contains(r"[KR](?!P)")]
    g = (r.groupby("Stripped.Sequence").RT.median().reset_index()
           .rename(columns={"Stripped.Sequence": "peptide", "RT": "RT exp"}))
    g.to_csv(f"{out}/pl_diann_{stem}.tsv", sep="\t", index=False)
    print(stem, len(g), "peptides")
