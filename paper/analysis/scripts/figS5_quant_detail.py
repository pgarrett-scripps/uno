#!/usr/bin/env python3
"""figS5_quant_detail.png -- dilution-series detail behind fig5.

Same run as fig5_quant.py (see that file's header for the exact uno command and
the saved outputs under ``benchmark/quant/results/``). Nothing here is
re-searched; both panels read files on disk.

Panels
    a  per-species abundance trajectories across all eight dilution levels.
       Each protein is normalised to its own 50 fmol abundance; the line is the
       species median and the band is the interquartile range. The dashed guide
       is the ideal 1:1 dilution response. Source: ``benchmark/quant/results/ab.tsv``.

    b  detection-loss trace. UPS1 precursors confident at 50 fmol, and the
       fraction of them still detected at each lower level. Two definitions of
       "detected", because they answer different questions:

         * matrix cell filled -- non-zero intensity in >=2 of 3 replicates in
           uno's LFQ matrix. koth's alignment transfers (match-between-runs)
           into runs where the peptide was not independently found, so this is
           an upper bound.
         * independently detected -- koth per-cell q <= 0.01 in >=2 of 3
           replicates, from the aligner's own
           ``uno_proc/align/qvalue_matrix.tsv``. This is the honest detection
           rate and it is what limits low-abundance differential sensitivity.

       Both traces are computed here from
       ``benchmark/quant/results/lfq.tsv`` joined to
       ``/mnt/labdata/pgarrett/PXD033350_DirectMS1Quant/uno_proc_newkoth/align_avg/\
qvalue_matrix.tsv``
       on (massCalib, charge, rtApex) -- an exact join, all 11520 LFQ rows match.

       The two traces have DIFFERENT denominators and must be reported as such:
       the filled trace is over precursors whose cell is filled at 50 fmol
       (n=184), the q-confident trace over precursors confidently detected at
       50 fmol (n=82). The figure legend prints both.

    c  volcano for the 50 vs 0.1 fmol contrast (peptide_fc_de.analyze), moved
       here from the main-text quant figure together with the
       differential-testing discussion it illustrates (SI Section S7).
"""
from __future__ import annotations

import math
import pathlib
import re
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import style  # noqa: E402

sys.path.insert(0, str(style.BENCH / "quant"))
import dilution_track as dt  # noqa: E402
import peptide_fc_de as de  # noqa: E402

RESULTS = style.BENCH / "quant" / "results"
AB = RESULTS / "ab.tsv"
LFQ = RESULTS / "lfq.tsv"
ALIGN = pathlib.Path("/mnt/labdata/pgarrett/PXD033350_DirectMS1Quant/uno_proc_newkoth/align_avg")
QMAT = ALIGN / "qvalue_matrix.tsv"

FMOL = dt.FMOL
LEVELS = sorted(FMOL)
QTHR = 0.01          # koth per-cell q-value for "independently detected"
MIN_REPS = 2         # of 3 replicates

SPECIES_ORDER = ["UPS", "ECOLI", "HUMAN"]
SPECIES_LABEL = {"UPS": "UPS1 (spike)", "ECOLI": "E. coli (background)",
                 "HUMAN": "human (absent)"}
SPECIES_COLOR = {"UPS": style.C["uno"], "ECOLI": style.C["target"],
                 "HUMAN": style.C["entrapment"]}


def run_columns(header):
    cols = {}
    for i, c in enumerate(header):
        m = re.search(r"UPS_(\d)_(\d\d)$", c)
        if m:
            cols.setdefault(int(m.group(1)), []).append(i)
    if sorted(cols) != LEVELS:
        raise SystemExit(f"expected levels {LEVELS}, found {sorted(cols)}")
    return cols


# --------------------------------------------------------------------------
# panel a -- per-species trajectories
# --------------------------------------------------------------------------
def panel_a(ax):
    lines = AB.read_text().splitlines()
    cols = run_columns(lines[0].split("\t"))
    fmol = np.array([FMOL[lv] for lv in LEVELS], dtype=float)

    tr = {s: [] for s in SPECIES_ORDER}
    for ln in lines[1:]:
        c = ln.split("\t")
        means = []
        for lv in LEVELS:
            v = [float(c[i]) for i in cols[lv] if float(c[i]) > 0]
            means.append(sum(v) / len(v) if v else math.nan)
        means = np.array(means)
        if not np.isfinite(means[0]) or means[0] <= 0:
            continue
        tr[dt.species(c[0])].append(np.log2(means / means[0]))

    summary = {}
    for sp in SPECIES_ORDER:
        a = np.array(tr[sp])
        if not len(a):
            continue
        med = np.nanmedian(a, axis=0)
        summary[sp] = (len(a), med)
        if len(a) >= 4:
            lo = np.nanpercentile(a, 25, axis=0)
            hi = np.nanpercentile(a, 75, axis=0)
            ax.fill_between(fmol, lo, hi, color=SPECIES_COLOR[sp], alpha=0.18,
                            lw=0, zorder=2)
        ax.plot(fmol, med, color=SPECIES_COLOR[sp], lw=1.8, marker="o", ms=3.2,
                zorder=4, label=f"{SPECIES_LABEL[sp]}, n={len(a)}")

    ax.plot(fmol, np.log2(fmol / fmol[0]), color=style.C["neutral"], lw=1.0,
            ls="--", zorder=3, label="ideal 1:1 response")
    ax.axhline(0, color=style.C["light"], lw=0.8, zorder=1)
    ax.set_xscale("log")
    ax.set_xticks(fmol)
    ax.set_xticklabels([("%g" % f) for f in fmol], rotation=45, fontsize=6.6)
    ax.tick_params(axis="x", which="minor", length=0)
    ax.set_xlabel("UPS1 spike (fmol)")
    ax.set_ylabel("abundance, log$_2$ vs 50 fmol")
    ax.set_title("abundance trajectory by species", pad=4)
    ax.legend(loc="lower right", fontsize=6.6)
    return summary


# --------------------------------------------------------------------------
# panel b -- detection loss
# --------------------------------------------------------------------------
def panel_b(ax):
    qlines = QMAT.read_text().splitlines()
    qcols = run_columns(qlines[0].split("\t"))
    qrow = {}
    for ln in qlines[1:]:
        c = ln.split("\t")
        qrow[(round(float(c[0]), 5), c[2], round(float(c[3]), 4))] = c

    lines = LFQ.read_text().splitlines()
    lcols = run_columns(lines[0].split("\t"))
    ups = []
    for ln in lines[1:]:
        c = ln.split("\t")
        if c[2] == "0" and "ups|" in c[3]:
            key = (round(float(c[4]), 5), c[1], round(float(c[5]), 4))
            if key not in qrow:
                raise SystemExit(f"LFQ row {key} not in {QMAT}")
            ups.append((c, qrow[key]))

    def n_filled(row, lv):
        return sum(1 for i in lcols[lv] if float(row[i]) > 0)

    def n_conf(row, lv):
        return sum(1 for i in qcols[lv] if float(row[i]) <= QTHR)

    top = LEVELS[0]  # level 1 == 50 fmol
    base_fill = [p for p in ups if n_filled(p[0], top) >= MIN_REPS]
    base_conf = [p for p in ups if n_conf(p[1], top) >= MIN_REPS]

    fmol = np.array([FMOL[lv] for lv in LEVELS], dtype=float)
    fill = np.array([sum(1 for p in base_fill if n_filled(p[0], lv) >= MIN_REPS)
                     for lv in LEVELS], dtype=float)
    conf = np.array([sum(1 for p in base_conf if n_conf(p[1], lv) >= MIN_REPS)
                     for lv in LEVELS], dtype=float)

    series = [
        ("matrix cell filled", fill / len(base_fill), len(base_fill),
         style.C["accent"], "o", "-"),
        (f"independently detected (q$\\leq${QTHR:g})", conf / len(base_conf),
         len(base_conf), style.C["uno"], "s", "-"),
    ]
    for label, frac, n, col, mk, ls in series:
        ax.plot(fmol, 100 * frac, color=col, lw=1.6, ls=ls, marker=mk, ms=3.4,
                label=f"{label}, n={n}", zorder=4)

    ax.set_xscale("log")
    ax.set_xticks(fmol)
    ax.set_xticklabels([("%g" % f) for f in fmol], rotation=45, fontsize=6.6)
    ax.tick_params(axis="x", which="minor", length=0)
    ax.set_xlabel("UPS1 spike (fmol)")
    ax.set_ylabel("UPS1 precursors retained (%)")
    # headroom above 100% so the legend never sits on top of a curve
    ax.set_ylim(0, 134)
    ax.set_yticks([0, 20, 40, 60, 80, 100])
    ax.set_title("detection loss down the dilution", pad=4)
    ax.legend(loc="upper left", fontsize=6.6)
    return fmol, fill, len(base_fill), conf, len(base_conf)


# --------------------------------------------------------------------------
# panel c -- volcano, 50 vs 0.1 fmol
# --------------------------------------------------------------------------
def panel_c(ax):
    import math as _math

    # gate=0 -> every testable protein enters the cloud; the significance call
    # below is the gated (|log2 FC| >= 1) BH-1% set that peptide_fc_de prints.
    cloud, n_prot, _ = de.analyze(str(LFQ), "50", "0_1", gate=0.0)
    gated, _, thr = de.analyze(str(LFQ), "50", "0_1", gate=1.0)
    sig = {acc for _, _, acc, _, _ in gated[:thr]}

    ymax = 0.0
    counts = {s: 0 for s in SPECIES_ORDER}
    for p, sp, acc, m, n in cloud:
        y = -_math.log10(max(p, 1e-300))
        ymax = max(ymax, y)
        hit = acc in sig
        if hit:
            counts[sp] += 1
        ax.plot(m, y, ls="none", marker="o",
                ms=4.2 if hit else 2.4,
                mfc=SPECIES_COLOR[sp] if hit else "none",
                mec=SPECIES_COLOR[sp],
                mew=0.7 if hit else 0.5,
                alpha=1.0 if hit else 0.35,
                zorder=5 if hit else 2)

    pthr = gated[thr - 1][0] if thr else None
    if pthr:
        ax.axhline(-_math.log10(pthr), color=style.C["neutral"], lw=0.9, ls=":",
                   zorder=3)
        ax.annotate("BH 1%", xy=(0.015, -_math.log10(pthr)),
                    xycoords=("axes fraction", "data"),
                    va="bottom", ha="left", fontsize=6.6,
                    color=style.C["neutral"])
    for x in (-1, 1):
        ax.axvline(x, color=style.C["light"], lw=0.8, ls="--", zorder=1)

    handles = [plt.Line2D([], [], ls="none", marker="o", ms=4.2,
                          mfc=SPECIES_COLOR[s], mec=SPECIES_COLOR[s],
                          label=f"{SPECIES_LABEL[s]} — {counts[s]}")
               for s in SPECIES_ORDER]
    ax.legend(handles=handles, loc="upper left", title="called changed",
              title_fontsize=6.8, fontsize=6.6, alignment="left")
    ax.set_xlabel("mean log$_2$ FC (50 / 0.1 fmol)")
    ax.set_ylabel(r"$-\log_{10}\ p$")
    ax.set_ylim(-0.05 * ymax, ymax * 1.62)
    ax.set_title("50 vs 0.1 fmol contrast", pad=4)
    return cloud, counts, thr, n_prot


def main():
    fig, axes = plt.subplots(1, 3, figsize=(style.W2, 2.7))
    summary = panel_a(axes[0])
    fmol, fill, nf, conf, nc = panel_b(axes[1])
    cloud, counts, thr, n_prot = panel_c(axes[2])
    for ax, letter in zip(axes, "abc"):
        style.panel_label(ax, letter, x=-0.19)
    fig.tight_layout(w_pad=2.2)
    style.save(fig, "figS5_quant_detail", tight=False)

    for sp, (n, med) in summary.items():
        print(f"panel a: {sp:6s} n={n:3d}  median log2 vs 50 fmol = "
              + ", ".join(f"{FMOL[lv]}:{v:+.2f}" for lv, v in zip(LEVELS, med)))
    print(f"panel b: filled base n={nf}, q-confident base n={nc}")
    for lv, f, c in zip(LEVELS, fill, conf):
        print(f"  {FMOL[lv]:>5} fmol  filled {int(f):3d}/{nf} = {f / nf:.0%}   "
              f"q-confident {int(c):3d}/{nc} = {c / nc:.0%}")
    print(f"panel c: {len(cloud)} testable proteins of {n_prot}; BH-1% + 2-fold "
          f"accepts {thr} -> " + ", ".join(f"{s}={counts[s]}" for s in SPECIES_ORDER))


if __name__ == "__main__":
    main()
