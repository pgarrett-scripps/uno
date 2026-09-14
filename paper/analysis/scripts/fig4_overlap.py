#!/usr/bin/env python3
"""Protein-list stability within and between engines over all primary draws.

Panel a shows every pairwise Jaccard index and its median/range. Pairs reuse
searches, so they are descriptive comparisons, not independent replicates or
confidence intervals. Panel b shows each uno seed's recovery of the proteins
accepted by every available ms1searchpy run for that file.
"""
from __future__ import annotations

import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import matplotlib.pyplot as plt
import numpy as np
import entrapment_data as ed
import seed_data
import style
from _assets import record

THR = 0.01
OUT_DATA = style.FIGDATA / "overlap_analysis.tsv"
COLS = ("analysis", "dataset", "a", "b", "n_a", "n_b", "shared", "union", "fraction")


def main():
    data = seed_data.load()
    rows = ["\t".join(COLS)]
    values = {}

    def row(analysis, ds, a, b, aa, bb, denominator=None):
        union, shared = len(aa | bb), len(aa & bb)
        den = union if denominator is None else denominator
        fraction = shared / den if den else 0.0
        rows.append("\t".join(map(str, (analysis, ds, a, b, len(aa), len(bb), shared,
                                         union, f"{fraction:.8f}"))))
        values.setdefault((ds, analysis), []).append(100 * fraction)

    for ds, runs in data.items():
        uno = [(f"seed{r.rep}", r.accepted_human(THR)) for r in runs["uno"]]
        ms1 = [(f"rep{r.rep}", r.accepted_human(THR)) for r in runs["ms1"]]
        for analysis, pairs in (
            ("jaccard_within_uno", itertools.combinations(uno, 2)),
            ("jaccard_within_ms1", itertools.combinations(ms1, 2)),
            ("jaccard_uno_vs_ms1", itertools.product(uno, ms1)),
        ):
            for (a, aa), (b, bb) in pairs:
                row(analysis, ds, a, b, aa, bb)
        core = set.intersection(*(v for _, v in ms1))
        union = set.union(*(v for _, v in ms1))
        if not core:
            raise ValueError(f"empty comparator core on {ds}")
        for a, aa in uno:
            row("uno_recovers_ms1_core", ds, a, "ms1_core", aa, core, len(core))
            row("uno_confirmed_by_any_ms1", ds, a, "ms1_union", aa, union, len(aa))
            row("uno_in_no_ms1", ds, a, "uno_only", aa, aa - union, len(aa))
        for a, aa in ms1:
            others = set.union(*(v for b, v in ms1 if a != b))
            row("ms1_confirmed_by_other_ms1", ds, a, "ms1_others", aa, others, len(aa))

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(style.W2, 3.2),
                                     gridspec_kw={"width_ratios": [1.6, 1]})
    series = [("jaccard_within_uno", -0.26, style.C["uno"], "uno vs uno"),
              ("jaccard_uno_vs_ms1", 0, style.C["neutral"], "uno vs ms1searchpy"),
              ("jaccard_within_ms1", 0.26, style.C["ms1searchpy"], "ms1searchpy vs itself")]
    for i, ds in enumerate(ed.DATASETS):
        for analysis, offset, color, label in series:
            vv = values[ds, analysis]
            style.spread_marker(ax_a, i + offset, vv, color, jitter=0.16 / max(1, len(vv)-1),
                                half_width=0.105, ms=1.8, label=label if i == 0 else None)
        vv = values[ds, "uno_recovers_ms1_core"]
        style.spread_marker(ax_b, i, vv, style.C["uno"], jitter=0.038,
                            half_width=0.23, ms=3)
    for ax in (ax_a, ax_b):
        ax.set_xticks(range(len(ed.DATASETS)), [ed.SHORT[d] for d in ed.DATASETS],
                      fontsize=6.5, rotation=40, ha="right")
        ax.set_xlim(-0.6, len(ed.DATASETS) - 0.4)
        ax.set_ylim(0, 105)
        ax.set_xlabel("dataset")
    ax_a.set_ylabel("Jaccard index of 1% lists (%)")
    ax_b.set_ylabel("ms1searchpy core recovered by uno (%)")
    ax_a.legend(loc="lower left", fontsize=6.1)
    style.panel_label(ax_a, "a", x=-0.16)
    style.panel_label(ax_b, "b", x=-0.26)
    fig.tight_layout()
    OUT_DATA.write_text("\n".join(rows) + "\n")
    style.save(fig, "fig4_overlap", tight=False)
    for suffix, extension in (("", "png"), ("-pdf", "pdf")):
        record(f"fig.fig4-overlap{suffix}", f"figures/fig4_overlap.{extension}", kind="figure",
               inputs=seed_data.inputs(), desc="Within- and between-engine protein-list agreement over all primary draws")
    print(f"wrote {OUT_DATA}")


if __name__ == "__main__":
    main()
