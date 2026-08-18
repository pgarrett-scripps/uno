#!/usr/bin/env python3
"""Shared figure style for the uno JPR manuscript.

EVERY figure script in ``paper/scripts/figures/`` imports this module so the
whole set reads as one system. Import it *before* creating any figure -- the
rcParams are applied at import time.

Quick start
-----------
::

    import sys, pathlib
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

    import matplotlib.pyplot as plt
    import style

    fig, ax = plt.subplots(1, 2, figsize=(style.W2, 2.8))
    ax[0].bar(x, y, color=style.C["uno"], label="uno")
    style.panel_label(ax[0], "a")
    style.panel_label(ax[1], "b")
    style.save(fig, "fig2_headtohead")      # -> paper/figures/fig2_headtohead.{png,pdf}

API
---
Sizes (inches; ACS two-column page)
    ``W1`` = 3.4   single column
    ``W2`` = 7.0   full width
    ``figsize(ncol=2, h=3.0)`` -> ``(W2, 3.0)`` convenience

Colors -- ``C`` is a dict of *named roles*, not raw hexes. Use the role, never
the hex, so a palette change propagates. All entries are from the Wong
(Nat. Methods 2011) colorblind-safe set; no red/green pairing anywhere.

    C["uno"]          blue        the engine this paper describes
    C["ms1searchpy"]  orange      the reference engine
    C["target"]       green       target hits / true positives
    C["decoy"]        grey        decoy hits
    C["entrapment"]   purple      entrapment (yeast) false positives
    C["accent"]       sky         a third series when two roles are not enough
    C["accent2"]      vermillion  a fourth series (safe against blue/orange)
    C["neutral"]      dark grey   guide lines, annotations, box edges
    C["light"]        pale grey   fills, shading, minor gridlines

    ``ENGINE`` maps the two engines to (color, marker) for consistent scatter
    styling: ``ENGINE["uno"] == (C["uno"], "o")``.

Paths (all ``pathlib.Path``, absolute)
    ``REPO``    repository root
    ``PAPER``   ``<repo>/paper``
    ``FIGDIR``  ``<repo>/paper/figures``     (created on import)
    ``FIGDATA`` ``<repo>/paper/figdata``     (created on import)
    ``BENCH``   ``<repo>/benchmark``

Helpers
    ``panel_label(ax, "a")``       bold lowercase letter, top-left, outside the
                                   axes; pass ``x=``/``y=`` in axes coords to
                                   nudge it.
    ``save(fig, "fig2_headtohead")``  ``tight_layout`` (unless
                                   ``tight=False``) then write BOTH
                                   ``figures/<name>.png`` at 300 dpi and
                                   ``figures/<name>.pdf``. Returns the png path.
                                   Do NOT pass an extension in ``name``.
    ``read_tsv(path)``             ``pandas.read_csv`` with ``sep='\\t'`` and
                                   ``#`` comment lines skipped -- the format
                                   every ``paper/figdata/*.tsv`` uses.
    ``median_range(values)``       ``(median, min, max)`` for replicate points.
    ``spread_marker(ax, x, values, color)``  draw replicate points + a median
                                   tick at ``x``; the house style for "this tool
                                   is non-deterministic, here are its reruns".

House rules
    * Never invent a number. Every plotted value traces to a file on disk
      (see ``paper/FIGURE_MANIFEST.md``).
    * No top/right spines, no gridlines by default, no chartjunk.
    * Panel letters are lowercase, bold, at the top-left of each panel.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# --------------------------------------------------------------------------
# paths
# --------------------------------------------------------------------------
PAPER = Path(__file__).resolve().parents[2]
REPO = PAPER.parent
FIGDIR = PAPER / "figures"
FIGDATA = PAPER / "figdata"
BENCH = REPO / "benchmark"

FIGDIR.mkdir(parents=True, exist_ok=True)
FIGDATA.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------------------------------
# canvas sizes -- ACS / JPR two-column page
# --------------------------------------------------------------------------
W1 = 3.4  # single column
W2 = 7.0  # full width
DPI = 300


def figsize(ncol: int = 2, h: float = 3.0) -> tuple[float, float]:
    """``(width, h)`` for a 1- or 2-column figure."""
    return (W1 if ncol == 1 else W2, h)


# --------------------------------------------------------------------------
# palette -- Wong (Nat. Methods 2011) colorblind-safe; no red/green pair
# --------------------------------------------------------------------------
C = {
    "uno": "#0072B2",          # blue
    "ms1searchpy": "#E69F00",  # orange
    "target": "#009E73",       # bluish green
    "decoy": "#7F7F7F",        # grey
    "entrapment": "#CC79A7",   # reddish purple
    "accent": "#56B4E9",       # sky blue
    "accent2": "#D55E00",      # vermillion
    "neutral": "#333333",
    "light": "#D9D9D9",
}

ENGINE = {
    "uno": (C["uno"], "o"),
    "ms1searchpy": (C["ms1searchpy"], "s"),
}

# --------------------------------------------------------------------------
# rcParams -- applied at import
# --------------------------------------------------------------------------
plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": DPI,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "pdf.fonttype": 42,      # embed TrueType, editable text in the PDF
    "ps.fonttype": 42,
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 8,
    "axes.titlesize": 8.5,
    "axes.labelsize": 8,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": False,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size": 3.0,
    "ytick.major.size": 3.0,
    "legend.fontsize": 7.5,
    "legend.frameon": False,
    "legend.handlelength": 1.4,
    "legend.handletextpad": 0.5,
    "legend.labelspacing": 0.3,
    "legend.borderaxespad": 0.3,
    "lines.linewidth": 1.2,
    "lines.markersize": 4.0,
    "patch.linewidth": 0.6,
    "figure.titlesize": 9,
})


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def panel_label(ax, letter: str, x: float = -0.16, y: float = 1.06, **kw):
    """Bold lowercase panel letter at the top-left, in axes coordinates.

    Nudge with ``x``/``y`` when a panel has a wide y-axis label. Returns the
    ``Text`` artist.
    """
    opts = dict(fontsize=10, fontweight="bold", va="bottom", ha="left",
                transform=ax.transAxes, color="black")
    opts.update(kw)
    return ax.text(x, y, letter, **opts)


def save(fig, name: str, tight: bool = True, pdf: bool = True) -> Path:
    """Write ``paper/figures/<name>.png`` (300 dpi) and ``<name>.pdf``.

    ``name`` must NOT carry an extension. Returns the png path.
    """
    if Path(name).suffix:
        raise ValueError(f"pass a bare figure name without an extension, got {name!r}")
    if tight:
        fig.tight_layout()
    png = FIGDIR / f"{name}.png"
    fig.savefig(png, dpi=DPI)
    if pdf:
        fig.savefig(FIGDIR / f"{name}.pdf")
    plt.close(fig)
    print(f"wrote {png}" + (f" and {FIGDIR / (name + '.pdf')}" if pdf else ""))
    return png


def read_tsv(path):
    """Read a ``paper/figdata`` TSV: tab-separated, ``#`` provenance lines skipped."""
    import pandas as pd

    return pd.read_csv(path, sep="\t", comment="#")


def median_range(values) -> tuple[float, float, float]:
    """``(median, min, max)`` of an iterable of replicate measurements."""
    import numpy as np

    v = np.asarray(list(values), dtype=float)
    return float(np.median(v)), float(v.min()), float(v.max())


def spread_marker(ax, x: float, values, color: str, *, jitter: float = 0.055,
                  half_width: float = 0.16, ms: float = 3.4, zorder: int = 4,
                  label: str | None = None, median_label: str | None = None):
    """House style for a non-deterministic tool's replicate measurements.

    Draws every replicate as an open circle (spread symmetrically around ``x``
    by ``jitter``), a thin vertical connector spanning min..max, and a solid
    horizontal tick of half-width ``half_width`` at the median. The reader sees
    both the individual reruns and the summary statistic.
    """
    import numpy as np

    v = np.asarray(list(values), dtype=float)
    n = len(v)
    offs = (np.arange(n) - (n - 1) / 2) * jitter if n > 1 else np.zeros(1)
    ax.plot([x, x], [v.min(), v.max()], color=color, lw=0.8, zorder=zorder - 1,
            solid_capstyle="butt")
    ax.plot(x + offs, np.sort(v), ls="none", marker="o", mfc="white", mec=color,
            mew=0.9, ms=ms, zorder=zorder, label=label)
    med = float(np.median(v))
    ax.plot([x - half_width, x + half_width], [med, med], color=color, lw=1.6,
            solid_capstyle="butt", zorder=zorder + 1, label=median_label)
    return med
