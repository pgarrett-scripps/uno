#!/usr/bin/env python3
"""fig4_quant.png -- spike-in dilution response on DirectMS1Quant's own series.

Data: PXD033350, UPS1 spiked into an E. coli background, 8 levels
(50/25/10/5/2.5/1/0.25/0.1 fmol) x 3 replicates = 24 runs. Human is in the
search DB but absent from the sample, so it is a known-absent probe.

PROVENANCE -- every number below traces to files produced by these commands
(run 2026-07-28 from the repo root, uno built with `cargo build --release`):

    # target-only DB: 48 UPS1 + 4431 E. coli + 20193 human
    awk '/^>/{p=(!/^>DECOY_/)} p' \
        /mnt/labdata/pgarrett/PXD033350_DirectMS1Quant/part1/\
sprot_ecoli_human_with_ups_replaced_shuffled.fasta \
        > benchmark/quant/results/ups_ecoli_human.fasta

    ./target/release/uno \
        --matrix /mnt/labdata/pgarrett/PXD033350_DirectMS1Quant/uno_proc_newkoth/\
align_avg/intensity_matrix.tsv \
        --matrix-min-runs 6 \
        --fasta benchmark/quant/results/ups_ecoli_human.fasta \
        --rt-model models/rt_cnn_chrono_global_v3_ch32.urtc \
        --profile-reassign \
        --out-abundance benchmark/quant/results/ab.tsv \
        --out-lfq benchmark/quant/results/lfq.tsv \
        --out-proteins benchmark/quant/results/proteins.tsv

    python benchmark/quant/dilution_track.py benchmark/quant/results/ab.tsv
    python benchmark/quant/peptide_fc_de.py benchmark/quant/results/lfq.tsv 50 0_1

Saved outputs: benchmark/quant/results/{ab,lfq,proteins}.tsv and
analysis_stdout.txt. NOTE these paths are gitignored (benchmark/.gitignore), so
they are NOT in the repository; they ship in the Zenodo archive instead.

This figure plots the AVERAGINE-PROJECTION matrix (koth_align with
`[lfq] averagine_projection = true`, under uno_proc_newkoth/align_avg), which is
uno's standard quant path and the source of the 326-protein result. The raw
box-sum variant is retained for comparison only, as {ab,lfq,proteins}_boxsum.tsv
(243 proteins under the current koth); the older *_avg.tsv files are the
superseded old-koth run.

Panels
    a  UPS1 protein abundance vs spike amount, each protein normalised to its
       own 50 fmol value; the dashed guide is the ideal 1:1 dilution response.
    b  Spearman rho of abundance vs log10(fmol), by species.

The volcano for the 50 vs 0.1 fmol contrast moved to the SI
(figS5_quant_detail.py panel c) together with the differential-testing
discussion it illustrates.
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

RESULTS = style.BENCH / "quant" / "results"
AB = RESULTS / "ab.tsv"

# level -> fmol, exactly as benchmark/quant/dilution_track.py defines it
FMOL = dt.FMOL
LEVELS = sorted(FMOL)  # 1..8, descending in fmol

SPECIES_ORDER = ["UPS", "ECOLI", "HUMAN"]
SPECIES_LABEL = {
    "UPS": "UPS1 (spike)",
    "ECOLI": "E. coli (flat background)",
    "HUMAN": "human (absent probe)",
}
SPECIES_COLOR = {
    "UPS": style.C["uno"],
    "ECOLI": style.C["target"],
    "HUMAN": style.C["entrapment"],
}


# --------------------------------------------------------------------------
# loading
# --------------------------------------------------------------------------
def run_columns(header):
    """{level: [column index per replicate]} from the ``UPS_<level>_<rep>`` names."""
    cols = {}
    for i, c in enumerate(header):
        m = re.search(r"UPS_(\d)_(\d\d)$", c)
        if m:
            cols.setdefault(int(m.group(1)), []).append(i)
    if sorted(cols) != LEVELS:
        raise SystemExit(f"expected levels {LEVELS}, found {sorted(cols)}")
    return cols


def load_abundance(path):
    """(rows, level_columns); each row is (accession, species, values-by-index)."""
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    cols = run_columns(header)
    rows = []
    for ln in lines[1:]:
        c = ln.split("\t")
        rows.append((c[0], dt.species(c[0]), c))
    return rows, cols


def level_means(cells, cols):
    """Mean of the detected (>0) replicates at each level; NaN if none detected."""
    out = []
    for lv in LEVELS:
        v = [float(cells[i]) for i in cols[lv] if float(cells[i]) > 0]
        out.append(sum(v) / len(v) if v else math.nan)
    return np.array(out)


# --------------------------------------------------------------------------
# panel a -- UPS1 abundance vs spike amount
# --------------------------------------------------------------------------
def panel_a(ax, rows, cols):
    fmol = np.array([FMOL[lv] for lv in LEVELS], dtype=float)
    traces = []
    for acc, sp, cells in rows:
        if sp != "UPS":
            continue
        m = level_means(cells, cols)
        if not np.isfinite(m[0]) or m[0] <= 0:
            continue
        traces.append(np.log2(m / m[0]))  # relative to the 50 fmol level
    traces = np.array(traces)

    for t in traces:
        ax.plot(fmol, t, color=SPECIES_COLOR["UPS"], lw=0.7, alpha=0.35,
                solid_capstyle="round", zorder=2)
    med = np.nanmedian(traces, axis=0)
    ax.plot(fmol, med, color=SPECIES_COLOR["UPS"], lw=2.0, zorder=4,
            marker="o", ms=3.4, label=f"UPS1 median (n={len(traces)})")

    ideal = np.log2(fmol / fmol[0])
    ax.plot(fmol, ideal, color=style.C["neutral"], lw=1.0, ls="--", zorder=3,
            label="ideal 1:1 response")

    ax.set_xscale("log")
    ax.set_xticks(fmol)
    ax.set_xticklabels([("%g" % f) for f in fmol], rotation=45, fontsize=6.6)
    ax.tick_params(axis="x", which="minor", length=0)
    ax.set_xlabel("UPS1 spike (fmol)")
    ax.set_ylabel("abundance, log$_2$ vs 50 fmol")
    ax.axhline(0, color=style.C["light"], lw=0.8, zorder=1)
    ax.legend(loc="lower right")
    ax.set_title("UPS1 tracks the dilution", pad=4)
    return traces, med


# --------------------------------------------------------------------------
# panel b -- Spearman rho by species
# --------------------------------------------------------------------------
def panel_b(ax, rows, cols):
    """Reproduces benchmark/quant/dilution_track.py's per-protein rho exactly."""
    logf_all = [math.log10(FMOL[lv]) for lv in LEVELS for _ in cols[lv]]
    idx_all = [i for lv in LEVELS for i in cols[lv]]

    rho = {s: [] for s in SPECIES_ORDER}
    for acc, sp, cells in rows:
        vals = [float(cells[i]) for i in idx_all]
        det = [(v, f) for v, f in zip(vals, logf_all) if v > 0]
        if len(det) < 8:  # dilution_track.py's minimum-run gate
            continue
        rho[sp].append(dt.spearman([f for _, f in det], [v for v, _ in det]))

    rng = np.random.default_rng(0)  # jitter only; deterministic seed
    for k, sp in enumerate(SPECIES_ORDER):
        v = np.array(rho[sp])
        y = k + (rng.random(len(v)) - 0.5) * 0.42 if len(v) > 1 else np.full(len(v), k)
        ax.plot(v, y, ls="none", marker="o", ms=2.6 if len(v) > 1 else 4.0,
                mfc=SPECIES_COLOR[sp], mec="none",
                alpha=0.55 if len(v) > 1 else 1.0, zorder=3)
        if len(v):
            m = float(np.median(v))
            ax.plot([m, m], [k - 0.34, k + 0.34], color=SPECIES_COLOR[sp],
                    lw=2.0, zorder=5, solid_capstyle="butt")
            ax.annotate(f"n={len(v)}, med {m:+.2f}", xy=(-1.0, k - 0.46),
                        va="bottom", ha="left", fontsize=6.6,
                        color=style.C["neutral"])

    ax.axvline(0, color=style.C["light"], lw=0.8, zorder=1)
    ax.set_yticks(range(len(SPECIES_ORDER)))
    ax.set_yticklabels(["UPS1", "E. coli", "human"])
    ax.set_ylim(-0.6, len(SPECIES_ORDER) - 0.4)
    ax.invert_yaxis()
    ax.set_xlim(-1.05, 1.05)
    ax.set_xlabel(r"Spearman $\rho$ (abundance vs log$_{10}$ fmol)")
    ax.set_title("spike responds, background flat", pad=4)
    return rho


# --------------------------------------------------------------------------
def main():
    rows, cols = load_abundance(AB)
    fig, axes = plt.subplots(1, 2, figsize=(style.W2, 2.75))

    traces, med = panel_a(axes[0], rows, cols)
    rho = panel_b(axes[1], rows, cols)

    for ax, letter in zip(axes, "ab"):
        style.panel_label(ax, letter)

    fig.tight_layout(w_pad=2.4)
    style.save(fig, "fig4_quant", tight=False)

    # echo every plotted summary number so the log is auditable
    print(f"panel a: {len(traces)} UPS1 proteins; median log2 vs 50 fmol = "
          + ", ".join(f"{FMOL[lv]}fmol {v:+.2f}" for lv, v in zip(LEVELS, med)))
    for sp in SPECIES_ORDER:
        v = rho[sp]
        if v:
            print(f"panel b: {sp:6s} n={len(v):3d} median_rho={float(np.median(v)):+.2f} "
                  f"frac(rho>0.5)={sum(1 for x in v if x > 0.5) / len(v):.0%}")


if __name__ == "__main__":
    main()
