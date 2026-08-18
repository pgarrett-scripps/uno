#!/usr/bin/env python3
"""Fig 1 -- the uno pipeline and scoring cascade. Schematic; no data.

Drawn with matplotlib patches so it restyles with the rest of the figure set.
Stage names and parameters are taken from the code, not from memory:

  crates/uno/src/main.rs             CLI surface (--features / --raw / --matrix,
                                     --fasta, --rt-model, --out-proteins,
                                     --out-lfq, --out-abundance)
  crates/uno-core/src/database.rs    "Phase 2: in-silico peptide database"
                                     sage-core digest; defaults cleave_at "KR",
                                     restrict "P", Cys +57.02146,
                                     decoy_mass_shift 11.0 Da, decoy_shuffle
  crates/uno-core/src/decoy.rs       decoy generators (default = peptide-interior
                                     global-pool shuffle, handled in database.rs)
  crates/uno-core/src/matching.rs    "Phase 3: match single-run features to
                                     theoretical peptides by precursor mass";
                                     MatchConfig defaults 10 ppm,
                                     isotope_offsets [-1, 0, 1]
  crates/uno-core/src/rt_cnn.rs      distilled pure-Rust RT CNN
  crates/uno-core/src/calibrate.rs   monotone (PAVA) piecewise-linear calibration
  crates/uno-core/src/score.rs       "Phase 5: discriminant scoring, bootstrap RT
                                     calibration, and target-decoy FDR";
                                     top_k default 50, 3-fold CV split BY FEATURE
  crates/uno-core/src/protein.rs     "Phase 6: protein-level aggregation + FDR"
                                     binomial surprise + greedy parsimony
  benchmark/uno.toml                 the digest used in this paper: 7-30 aa, 0 MC

Honesty: uno does NOT read raw spectra or find features on the path benchmarked
in this paper -- it consumes a biosaur2 feature list, exactly as ms1searchpy
does. Feature detection is drawn as an external (dashed) stage.

Design notes:
  * Everything the single static binary does sits inside one outlined
    container, which IS the figure's message: the external world contributes a
    FASTA and one feature table, and one process does the rest.
  * Each cascade stage carries a small wordless glyph of its mechanism (mass
    window, RT calibration fit, target/decoy score separation, peptide-to-
    protein aggregation). Target/decoy in the glyphs is never color-alone:
    decoys are open markers and both humps are direct-labeled, because the
    Wong green/grey pair is weak under deuteranopia.
  * Boxes are laid out top-down with an explicit height budget per text line
    (TITLE_H / LINE_H below), so adding a line of text grows the box instead
    of overflowing it. UNITS_PER_INCH fixes the vertical data-unit -> inch
    mapping; LINE_H / UNITS_PER_INCH must stay >= the rendered line height of
    DETAIL_SIZE text (6.0 pt x 1.42 = 0.118 in; 2.0 / 16.8 = 0.119 in).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402

import style  # noqa: E402

# fills / edges by box role
ROLE = {
    "uno": dict(fc="#E4EFF7", ec=style.C["uno"], ls="solid", lw=1.0),
    "ext": dict(fc="#FFFFFF", ec=style.C["decoy"], ls=(0, (3, 2)), lw=0.9),
    "out": dict(fc="#E6F4EF", ec=style.C["target"], ls="solid", lw=1.0),
}
CONTAINER_FC = "#F5F9FC"
GLYPH_NOTE = "#5A7A8C"

TITLE_SIZE, DETAIL_SIZE, NOTE_SIZE = 7.2, 6.0, 5.2
PAD_TOP, TITLE_H, LINE_H, PAD_BOT = 1.6, 2.5, 2.0, 0.7
GAP = 2.6            # vertical space between stacked boxes (the arrow lives here)
UNITS_PER_INCH = 16.8
CAPTION_IN = 0.42    # strip reserved at the bottom for the figure's footnote

# lane geometry (x, data units 0..100 over the 7 in width)
CONT_X0, CONT_X1 = 3.0, 69.5          # the uno container
LX, LW = 22.0, 34.0                   # database lane, inside the container
MX, MW = 36.25, 62.0                  # cascade spine, inside the container
RX, RW = 85.5, 28.0                   # external feature lane, outside
GLYPH_W = 12.0                        # width reserved at a box's right edge


def box_h(detail: str | None, min_h: float = 0.0) -> float:
    n = len(detail.split("\n")) if detail else 0
    return max(PAD_TOP + TITLE_H + n * LINE_H + PAD_BOT, min_h)


def box(ax, cx, top, w, title, detail=None, role="uno", glyph=None,
        min_h: float = 0.0) -> float:
    """Draw a box whose TOP edge is at `top`; return its bottom edge.

    With `glyph`, the right GLYPH_W units are reserved and the callable is
    handed that region's (x0, x1, y0, y1); text centers in the remainder.
    """
    h = box_h(detail, min_h)
    st = ROLE[role]
    ax.add_patch(FancyBboxPatch(
        (cx - w / 2, top - h), w, h,
        boxstyle="round,pad=0,rounding_size=1.0",
        facecolor=st["fc"], edgecolor=st["ec"], linewidth=st["lw"],
        linestyle=st["ls"], zorder=2, mutation_aspect=0.85,
    ))
    tcx = cx - GLYPH_W / 2 if glyph else cx
    n = len(detail.split("\n")) if detail else 0
    text_h = TITLE_H + n * LINE_H
    y = top - h / 2 + text_h / 2  # vertically centre the text block
    ax.text(tcx, y, title, ha="center", va="top", fontsize=TITLE_SIZE,
            fontweight="bold", color=style.C["neutral"], zorder=3)
    if detail:
        ax.text(tcx, y - TITLE_H, detail, ha="center", va="top",
                fontsize=DETAIL_SIZE, color=style.C["neutral"], zorder=3,
                linespacing=1.42)
    if glyph:
        pad = 1.0
        glyph(ax, cx + w / 2 - GLYPH_W, cx + w / 2 - pad,
              top - h + pad, top - pad)
    return top - h


def arrow(ax, x0, y0, x1, y1, lw=1.0):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="-|>", color=style.C["neutral"],
                                linewidth=lw, shrinkA=0, shrinkB=0,
                                mutation_scale=8), zorder=4)


# ---------------------------------------------------------------------------
# stage glyphs -- wordless mechanism sketches, drawn into (x0, x1, y0, y1).
# Markers are scatter points (size in pt^2) so the x/y unit anisotropy of the
# canvas cannot squash them.
# ---------------------------------------------------------------------------
def glyph_matching(ax, x0, x1, y0, y1):
    """A feature's mass window: candidates inside the +-10 ppm band."""
    axis_y = y0 + 0.20 * (y1 - y0)
    bx0, bx1 = x0 + 0.34 * (x1 - x0), x0 + 0.62 * (x1 - x0)
    ax.add_patch(Rectangle((bx0, axis_y), bx1 - bx0, y1 - axis_y - 0.8,
                           facecolor=style.C["uno"], alpha=0.13, zorder=2.5,
                           linewidth=0))
    ax.plot([x0, x1], [axis_y, axis_y], color=style.C["decoy"], lw=0.8,
            zorder=3, solid_capstyle="butt")
    fx = (bx0 + bx1) / 2
    ax.scatter([fx], [y1 - 2.0], marker="D", s=11, color=style.C["uno"],
               zorder=4)
    ax.plot([fx, fx], [axis_y, y1 - 2.5], color=style.C["uno"], lw=0.7,
            zorder=3.5)
    # candidates on the mass axis: filled = target, open = decoy, outside = cut
    ax.scatter([bx0 + 0.7, fx + 0.9], [axis_y, axis_y], s=9,
               color=style.C["target"], zorder=4)
    ax.scatter([bx1 - 0.4], [axis_y], s=9, facecolor="white",
               edgecolor=style.C["decoy"], linewidth=0.8, zorder=4)
    ax.scatter([x0 + 0.9, x1 - 0.9], [axis_y, axis_y], marker="x", s=8,
               color=style.C["decoy"], alpha=0.65, linewidth=0.9, zorder=4)
    ax.text((bx0 + bx1) / 2, y0, "±10 ppm", ha="center", va="bottom",
            fontsize=NOTE_SIZE, style="italic", color=GLYPH_NOTE, zorder=4)


def glyph_rt(ax, x0, x1, y0, y1):
    """Monotone calibration: predicted RT scatter, PAVA piecewise fit."""
    px0, py0 = x0 + 0.8, y0 + 1.0
    ax.plot([px0, px0, x1 - 0.4], [y1 - 0.6, py0, py0],
            color=style.C["decoy"], lw=0.8, zorder=3, solid_capstyle="butt")
    w, h = (x1 - 0.9) - (px0 + 0.5), (y1 - 1.0) - (py0 + 0.5)
    pts = [(0.06, 0.10), (0.18, 0.30), (0.30, 0.18), (0.42, 0.44),
           (0.55, 0.52), (0.66, 0.72), (0.80, 0.66), (0.92, 0.90)]
    ax.scatter([px0 + 0.5 + u * w for u, v in pts],
               [py0 + 0.5 + v * h for u, v in pts],
               s=3.5, color=style.C["decoy"], alpha=0.75, zorder=3.5,
               linewidth=0)
    seg = [(0.02, 0.08), (0.34, 0.26), (0.62, 0.64), (0.98, 0.92)]
    ax.plot([px0 + 0.5 + u * w for u, v in seg],
            [py0 + 0.5 + v * h for u, v in seg],
            color=style.C["uno"], lw=1.3, zorder=4, solid_capstyle="round")


def glyph_scoring(ax, x0, x1, y0, y1):
    """Two score distributions separating; dashes mark the FDR threshold."""
    base = y0 + 0.24 * (y1 - y0)
    top = y1 - 2.3
    w = x1 - x0

    def hump(cu, su, color, alpha):
        c, s = x0 + cu * w, su * w
        xs = [c + f * s for f in (-2.2, -1.5, -0.8, 0.0, 0.8, 1.5, 2.2)]
        ys = [base, base + 0.25 * (top - base), base + 0.75 * (top - base),
              top, base + 0.75 * (top - base), base + 0.25 * (top - base),
              base]
        ax.fill(xs, ys, color=color, alpha=alpha, linewidth=0, zorder=3)

    hump(0.32, 0.085, style.C["decoy"], 0.45)
    hump(0.64, 0.10, style.C["target"], 0.55)
    ax.plot([x0, x1], [base, base], color=style.C["decoy"], lw=0.8, zorder=3.5,
            solid_capstyle="butt")
    thr = x0 + 0.80 * w
    ax.plot([thr, thr], [base, y1 - 0.6], color=style.C["neutral"], lw=0.8,
            ls=(0, (2, 1.4)), zorder=4)
    ax.text(x0 + 0.30 * w, y1 - 0.4, "decoy", ha="center", va="top",
            fontsize=NOTE_SIZE, style="italic", color=GLYPH_NOTE, zorder=4)
    ax.text(x0 + 0.63 * w, y0, "target", ha="center", va="bottom",
            fontsize=NOTE_SIZE, style="italic", color=GLYPH_NOTE, zorder=4)


def glyph_protein(ax, x0, x1, y0, y1):
    """Peptide hits aggregating onto one protein; an unmatched decoy stays out."""
    pep_x = x0 + 0.22 * (x1 - x0)
    prot = (x0 + 0.74 * (x1 - x0), (y0 + y1) / 2)
    ys = [y0 + f * (y1 - y0) for f in (0.78, 0.55, 0.32)]
    for py in ys:
        ax.plot([pep_x, prot[0]], [py, prot[1]], color=style.C["target"],
                lw=0.7, alpha=0.8, zorder=3)
    ax.scatter([pep_x] * 3, ys, s=8, color=style.C["target"], zorder=4)
    ax.scatter([pep_x], [y0 + 0.06 * (y1 - y0)], s=8, facecolor="white",
               edgecolor=style.C["decoy"], linewidth=0.8, zorder=4)
    ax.scatter([prot[0]], [prot[1]], s=30, color=style.C["uno"], zorder=4)
    ax.text(prot[0], y0, "protein", ha="center", va="bottom",
            fontsize=NOTE_SIZE, style="italic", color=GLYPH_NOTE, zorder=4)


def main() -> int:
    fig, ax = plt.subplots(figsize=(style.W2, 8.0))
    ax.set_xlim(0, 100)
    ax.set_axis_off()

    TOP = 98.0

    # ---------------- external inputs, top row ------------------------------
    fasta_bot = box(ax, LX, TOP, LW, "protein FASTA", role="ext")

    y = TOP
    lane_r = [
        ("one LC–MS run", None),
        ("biosaur2 feature detection", None),
        ("single-run feature table",
         "mass, retention time, charge,\nintensity and feature quality"),
    ]
    for i, (title, detail) in enumerate(lane_r):
        if i:
            arrow(ax, RX, y, RX, y - GAP)
            y -= GAP
        y = box(ax, RX, y, RW, title, detail, "ext")
    feat_bottom = y

    # ---------------- the uno container -------------------------------------
    # Drawn before its children would be tidier, but its bottom depends on the
    # cascade height; record top now, draw at the end under everything (zorder 1).
    cont_top = fasta_bot - GAP / 2

    # database lane, inside the container
    arrow(ax, LX, fasta_bot, LX, fasta_bot - GAP)
    y = fasta_bot - GAP
    y = box(ax, LX, y, LW, "digest and decoy generation",
            "tryptic target peptides and\nshuffled, mass-shifted decoys")
    arrow(ax, LX, y, LX, y - GAP)
    y = box(ax, LX, y - GAP, LW, "mass-sorted peptide database")
    db_bottom = y

    # container label, in the free area right of the database lane
    ax.text(CONT_X1 - 2.0, cont_top - 1.6, "uno", ha="right", va="top",
            fontsize=9.5, fontweight="bold", color=style.C["uno"], zorder=3)
    ax.text(CONT_X1 - 2.0, cont_top - 4.4, "one static Rust binary",
            ha="right", va="top", fontsize=6.0, style="italic",
            color=GLYPH_NOTE, zorder=3)

    # ---------------- scoring cascade ----------------------------------------
    stages = [
        ("precursor-mass candidate matching",
         "features are matched to mass-compatible target and decoy peptides",
         glyph_matching),
        ("retention-time prediction and calibration",
         "an in-process model maps candidate sequences onto this run's gradient",
         glyph_rt),
        ("peptide–feature match scoring",
         "mass, retention time and available feature evidence are combined",
         glyph_scoring),
        ("protein inference and FDR control",
         "binomial protein score → parsimony → target/decoy competition",
         glyph_protein),
    ]
    arrow(ax, LX, db_bottom, LX, db_bottom - GAP)
    y = db_bottom - GAP
    for i, (title, detail, glyph) in enumerate(stages):
        if i:
            arrow(ax, MX, y, MX, y - GAP)
            y -= GAP
        y = box(ax, MX, y, MW, title, detail, "uno", glyph=glyph, min_h=9.4)
        if i == 0:
            match_top, match_bot = y + box_h(detail, 9.4), y
    spine_bottom = y
    cont_bot = spine_bottom - 1.8

    # feature table enters the binary at the matching stage
    entry_y = (match_top + match_bot) / 2
    ax.plot([RX, RX], [feat_bottom, entry_y], color=style.C["neutral"], lw=1.0,
            zorder=1.5, solid_capstyle="butt")
    arrow(ax, RX, entry_y, MX + MW / 2, entry_y)

    # ---------------- outputs ------------------------------------------------
    oy = cont_bot - GAP
    arrow(ax, LX, spine_bottom, LX, oy)
    box(ax, LX, oy, LW, "proteins.tsv",
        "accession, group members, peptide\ncounts, binomial score, q-value",
        "out")
    # the quantitation product exits sideways, under the data lane
    side_y = spine_bottom + box_h("x", 9.4) / 2
    ax.plot([MX + MW / 2, RX, RX], [side_y, side_y, oy + 1.0],
            color=style.C["neutral"], lw=1.0, zorder=1.5,
            solid_capstyle="butt")
    arrow(ax, RX, oy + 1.0, RX, oy)
    bottom = box(ax, RX, oy, RW + 3.0, "optional abundance matrices",
                 "identified peptides mapped\nto aligned intensities", "out")

    # the container itself, under everything inside it
    ax.add_patch(FancyBboxPatch(
        (CONT_X0, cont_bot), CONT_X1 - CONT_X0, cont_top - cont_bot,
        boxstyle="round,pad=0,rounding_size=1.4",
        facecolor=CONTAINER_FC, edgecolor=style.C["uno"], linewidth=1.2,
        zorder=1, mutation_aspect=0.85,
    ))

    fig.text(0.5, 0.030,
             "Solid blue boxes inside the outlined panel are stages of the single uno "
             "binary; dashed boxes are inputs and external tools; green boxes are outputs.",
             ha="center", va="bottom", fontsize=6.1, color=style.C["neutral"])
    fig.text(0.5, 0.004,
             "For the benchmark, biosaur2 performs feature detection and both engines "
             "receive the identical single-run feature table.",
             ha="center", va="bottom", fontsize=6.1, color=style.C["neutral"])

    # Size the canvas to the drawing: keep UNITS_PER_INCH fixed (so the height
    # budget above stays valid) and let the figure grow/shrink instead.
    ax.set_ylim(bottom - 2.0, 100.0)
    yrange = 100.0 - (bottom - 2.0)
    fig_h = yrange / UNITS_PER_INCH + CAPTION_IN
    fig.set_size_inches(style.W2, fig_h)
    fig.subplots_adjust(left=0.004, right=0.996, top=1.0, bottom=CAPTION_IN / fig_h)
    style.save(fig, "fig1_pipeline", tight=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
