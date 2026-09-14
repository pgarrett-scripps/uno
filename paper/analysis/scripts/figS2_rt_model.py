#!/usr/bin/env python3
"""Fig S4 -- what the RT model is and how it was trained.

  a  the network: residue tokens -> learned embedding -> conv stem -> three
     residual dilated blocks -> masked mean+max pool -> concatenate the global
     physicochemical features -> two-layer head -> one abstract RT. Shapes and
     hyperparameters are READ FROM THE SHIPPED MODEL FILE's header, not typed.
  b  the training corpus: peptides kept per source dataset after the
     amino-acid-only filter, and the fraction of each source that survives it.
  c  holdout accuracy by peptide length, the axis along which the model is
     weakest.
  d  holdout accuracy per source dataset, showing the model is not carried by
     one contributor.

DATA: models/rt_cnn_chrono_global_v3_ch32.urtc (header),
      data/rt/chronologer/{train,val,holdout}.tsv (the prepared corpus),
      data/rt/chronologer/Chronologer_DB_220308.gz (source labels), and
      predictions from target/release/rt_eval on the holdout partition.
      Written by paper/analysis/scripts/figS2_rt_model.py.
"""

from __future__ import annotations

import collections
import csv
import gzip
import os
import pathlib
import re
import struct
import subprocess
import sys
import tempfile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

import style  # noqa: E402
from _assets import record  # noqa: E402

MODEL = style.REPO / "models" / "rt_cnn_chrono_global_v3_ch32.urtc"
CORPUS = style.REPO / "data" / "rt" / "chronologer"
DB = CORPUS / "Chronologer_DB_220308.gz"
RT_EVAL = style.REPO / "target" / "release" / "rt_eval"
LABELS_PY = style.REPO / "crates" / "uno-rt-train" / "python" / "chronologer_to_labels.py"
LEN_BINS = [(6, 8), (9, 11), (12, 15), (16, 20), (21, 30), (31, 50)]


def model_header() -> dict:
    b = MODEL.read_bytes()
    if b[:4] != b"URTC":
        raise RuntimeError(f"{MODEL} is not a URTC model")
    p = 4

    def u32():
        nonlocal p
        v = struct.unpack_from("<I", b, p)[0]
        p += 4
        return v

    def f32():
        nonlocal p
        v = struct.unpack_from("<f", b, p)[0]
        p += 4
        return v

    version = u32()
    y_mean, y_std = f32(), f32()
    emb, ch, k, n_blocks = u32(), u32(), u32(), u32()
    dils = [u32() for _ in range(n_blocks)]
    head = u32()
    n_global = u32() if version >= 2 else 0
    return {"version": version, "y_mean": y_mean, "y_std": y_std, "emb": emb, "ch": ch,
            "k": k, "dilations": dils, "head": head, "n_global": n_global,
            "n_params": (len(b) - p) // 4}


def clean_seq_fn():
    """The corpus filter, imported from the preparation script itself."""
    src = LABELS_PY.read_text()
    start = src.index("CAM_C = ")
    end = src.index("def open_any")
    ns: dict = {"re": re}
    exec(compile(src[start:end], str(LABELS_PY), "exec"), ns)  # noqa: S102
    return ns["clean_seq"]


def read_labels(path: pathlib.Path) -> dict[str, float]:
    out = {}
    with path.open() as handle:
        for line in handle:
            if line.startswith("#") or line.startswith("sequence"):
                continue
            seq, value = line.rstrip("\n").split("\t")
            out[seq] = float(value)
    return out


def holdout_predictions(holdout: pathlib.Path) -> dict[str, float]:
    if not RT_EVAL.is_file():
        raise RuntimeError(f"missing {RT_EVAL}; build rt_eval in release mode first")
    with tempfile.TemporaryDirectory(prefix="uno-rt-fig-") as tmp:
        dump = pathlib.Path(tmp) / "preds.tsv"
        subprocess.run([str(RT_EVAL), str(MODEL), str(holdout), str(dump)],
                       cwd=style.REPO, check=True, capture_output=True, text=True)
        out = {}
        with dump.open() as handle:
            for line in handle:
                seq, value = line.rstrip("\n").split("\t")
                out[seq] = float(value)
    return out


def corpus_sources(clean_seq, holdout: set[str]):
    """Rows per source, rows surviving the filter, and holdout peptides per source."""
    total: collections.Counter = collections.Counter()
    kept: collections.Counter = collections.Counter()
    hold: dict[str, list[str]] = collections.defaultdict(list)
    with gzip.open(DB, "rt") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader)
        for source, modseq, _prosit_rt, _hi in reader:
            total[source] += 1
            seq = clean_seq(modseq)
            if seq is None:
                continue
            kept[source] += 1
            if seq in holdout:
                hold[source].append(seq)
    return total, kept, hold


def box(ax, x, y, w, h, text, color, fontsize=6.2, text_color="white"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.02",
                                facecolor=color, edgecolor="none", zorder=2))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fontsize,
            color=text_color, zorder=3, linespacing=1.35)


def arrow(ax, x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=7,
                                 lw=0.9, color=style.C["neutral"], zorder=1,
                                 shrinkA=0, shrinkB=0))


def panel_architecture(ax, hdr):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    blue, grey, dark = style.C["uno"], style.C["light"], style.C["neutral"]
    dil = ", ".join(str(d) for d in hdr["dilations"])

    xs, w = 0.02, 0.152
    ys, h = 0.66, 0.30
    box(ax, xs, ys, w, h, "peptide\nD F E A …  K", grey, text_color=dark)
    arrow(ax, xs + w, ys + h / 2, xs + w + 0.036, ys + h / 2)
    x = xs + w + 0.036
    box(ax, x, ys, w, h, f"embedding\n21 × {hdr['emb']}", blue)
    arrow(ax, x + w, ys + h / 2, x + w + 0.036, ys + h / 2)
    x += w + 0.036
    box(ax, x, ys, w, h, f"convolution\n{hdr['ch']} channels, k = {hdr['k']}", blue)
    arrow(ax, x + w, ys + h / 2, x + w + 0.036, ys + h / 2)
    x += w + 0.036
    box(ax, x, ys, w * 1.34, h,
        f"{len(hdr['dilations'])} residual blocks\ndilations {dil}", blue)
    # down into the pooling row
    ax.add_patch(FancyArrowPatch((x + w * 0.67, ys), (x + w * 0.67, ys - 0.08),
                                 arrowstyle="-|>", mutation_scale=7, lw=0.9,
                                 color=dark, zorder=1, shrinkA=0, shrinkB=0))

    y2 = ys - 0.08 - h
    xa = 0.10
    box(ax, xa, y2, w * 1.5, h, f"masked mean + max\n2 × {hdr['ch']} = {2 * hdr['ch']}", blue)
    box(ax, xa + w * 1.5 + 0.03, y2, w * 1.35, h,
        f"global features\n{hdr['n_global']} values", style.C["accent"])
    ax.text(xa + w * 1.5 + 0.015, y2 + h / 2, "+", ha="center", va="center", fontsize=9,
            color=dark, zorder=3)
    xh = xa + w * 1.5 + 0.03 + w * 1.35 + 0.036
    arrow(ax, xh - 0.036, y2 + h / 2, xh, y2 + h / 2)
    box(ax, xh, y2, w * 1.15, h, f"output layers\n{2 * hdr['ch'] + hdr['n_global']} → {hdr['head']} → 1",
        blue)
    arrow(ax, xh + w * 1.15, y2 + h / 2, xh + w * 1.15 + 0.036, y2 + h / 2)
    ax.text(xh + w * 1.15 + 0.045, y2 + h / 2, "hydrophobic\nindex", ha="left", va="center",
            fontsize=6.2, color=dark, linespacing=1.35)

    ax.text(0.0, 0.055,
            f"{hdr['n_params']:,} parameters   ·   ReLU   ·   trained on hydrophobic index "
            f"(mean {hdr['y_mean']:.1f}, SD {hdr['y_std']:.1f})\n"
            "amino acids only, cysteine implicitly carbamidomethyl",
            ha="left", va="center", fontsize=6.2, color=dark, linespacing=1.5)


def main() -> int:
    hdr = model_header()
    clean_seq = clean_seq_fn()
    holdout_lab = read_labels(CORPUS / "holdout.tsv")
    n_train = len(read_labels(CORPUS / "train.tsv"))
    n_val = len(read_labels(CORPUS / "val.tsv"))
    pred = holdout_predictions(CORPUS / "holdout.tsv")
    total, kept, hold = corpus_sources(clean_seq, set(holdout_lab))

    seqs = sorted(pred)
    y = np.array([holdout_lab[s] for s in seqs])
    p = np.array([pred[s] for s in seqs])
    lengths = np.array([len(s) for s in seqs])

    fig = plt.figure(figsize=(style.W2, 5.6))
    gs = fig.add_gridspec(3, 2, height_ratios=[0.62, 1.0, 1.0], hspace=0.62, wspace=0.34,
                          left=0.09, right=0.98, top=0.96, bottom=0.08)
    ax_arch = fig.add_subplot(gs[0, :])
    ax_corp = fig.add_subplot(gs[1, :])
    ax_len = fig.add_subplot(gs[2, 0])
    ax_src = fig.add_subplot(gs[2, 1])

    panel_architecture(ax_arch, hdr)
    style.panel_label(ax_arch, "a", x=-0.055, y=0.99)

    # ---- b: corpus composition ----
    order = [s for s, _ in kept.most_common()]
    x = np.arange(len(order), dtype=float)
    kept_v = np.array([kept[s] for s in order], dtype=float)
    drop_v = np.array([total[s] - kept[s] for s in order], dtype=float)
    ax_corp.bar(x, kept_v / 1000, 0.66, color=style.C["uno"], zorder=2,
                label="kept: unmodified or carbamidomethyl cysteine")
    ax_corp.bar(x, drop_v / 1000, 0.66, bottom=kept_v / 1000, color=style.C["light"],
                zorder=2, label="dropped: carries another modification")
    ax_corp.set_xticks(x)
    ax_corp.set_xticklabels([s.replace("_", "\n") for s in order], fontsize=5.6)
    ax_corp.set_xlim(-0.6, len(order) - 0.4)
    ax_corp.set_ylabel("peptide rows (thousands)")
    ax_corp.legend(loc="upper right", fontsize=6.3, handlelength=1.0)
    ax_corp.set_title(
        f"training corpus: {n_train:,} train / {n_val:,} validation / {len(holdout_lab):,} "
        "holdout peptides, split by sequence hash", fontsize=7.5, pad=3)
    style.panel_label(ax_corp, "b", x=-0.062, y=1.05)

    # ---- c: accuracy by peptide length ----
    xs = np.arange(len(LEN_BINS), dtype=float)
    mae, counts = [], []
    for lo, hi in LEN_BINS:
        m = (lengths >= lo) & (lengths <= hi)
        mae.append(float(np.abs(y[m] - p[m]).mean()))
        counts.append(int(m.sum()))
    ax_len.bar(xs, mae, 0.62, color=style.C["uno"], zorder=2)
    for xi, v, c in zip(xs, mae, counts):
        ax_len.text(xi, v + 0.012, f"{c:,}", ha="center", va="bottom", fontsize=5.6,
                    color=style.C["neutral"])
    ax_len.set_xticks(xs)
    ax_len.set_xticklabels([f"{lo} to {hi}" for lo, hi in LEN_BINS], fontsize=6.3)
    ax_len.set_xlabel("peptide length (residues)")
    ax_len.set_ylabel("holdout mean absolute error\n(hydrophobic index)")
    ax_len.set_ylim(0, max(mae) * 1.25)
    ax_len.set_title("accuracy by length, bar labels are peptide counts",
                     fontsize=7.5, pad=3)
    style.panel_label(ax_len, "c", x=-0.27, y=1.05)

    # ---- d: accuracy by source dataset ----
    rows = []
    for source, members in hold.items():
        uniq = sorted(set(members))
        if len(uniq) < 200:
            continue
        e = np.array([pred[s] - holdout_lab[s] for s in uniq])
        rows.append((source, len(uniq), float(np.abs(e).mean())))
    rows.sort(key=lambda r: r[2])
    xs = np.arange(len(rows), dtype=float)
    ax_src.bar(xs, [r[2] for r in rows], 0.62, color=style.C["uno"], zorder=2)
    overall = float(np.abs(y - p).mean())
    ax_src.axhline(overall, color=style.C["neutral"], lw=0.9, ls=(0, (4, 2)), zorder=3)
    ax_src.text(len(rows) - 0.45, overall, f" all {overall:.2f}", ha="right", va="bottom",
                fontsize=6, color=style.C["neutral"])
    ax_src.set_xticks(xs)
    ax_src.set_xticklabels([r[0].split("_")[0] for r in rows], fontsize=6.3, rotation=40,
                           ha="right", rotation_mode="anchor")
    ax_src.set_xlim(-0.6, len(rows) - 0.4)
    ax_src.set_ylim(0, max(r[2] for r in rows) * 1.25)
    ax_src.set_ylabel("holdout mean absolute error\n(hydrophobic index)")
    ax_src.set_title("accuracy by source dataset, at least 200 holdout peptides",
                     fontsize=7.5, pad=3)
    style.panel_label(ax_src, "d", x=-0.27, y=1.05)

    style.save(fig, "figS2_rt_model", tight=False)
    provenance = [MODEL, CORPUS / "train.tsv", CORPUS / "val.tsv", CORPUS / "holdout.tsv",
                  DB, LABELS_PY,
                  style.REPO / "crates" / "uno-rt-train" / "python" / "train_rt.py",
                  style.REPO / "crates" / "uno-core" / "src" / "rt_cnn.rs"]
    inputs = [pathlib.Path(os.path.relpath(q, style.PAPER)).as_posix() for q in provenance]
    record("fig.figs2-rt-model", "figures/figS2_rt_model.png", kind="figure",
           inputs=inputs,
           desc="Architecture, training corpus and holdout accuracy of uno's RT model")
    record("fig.figs2-rt-model-pdf", "figures/figS2_rt_model.pdf", kind="figure",
           inputs=inputs, desc="Vector twin of the RT model figure")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
