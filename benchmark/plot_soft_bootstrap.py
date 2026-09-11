#!/usr/bin/env python3
"""Plot the frozen soft-bootstrap Stage-A reruns from the verified, complete arm reports."""
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/uno-soft-matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "benchmark/results/soft_bootstrap"
ARMS = ("soft",)
LABELS = ("V3", "Soft confidence\nbootstrap")
SEEDS = (0, 1, 2, 5, 8)
FILES = ("astral_1ng_rep1", "astral_1ng_rep2", "astral_1ng_rep4",
         "astral_200pg_rep1", "astral_200pg_rep2", "astral_200pg_rep3")


def main():
    panels = {}
    hashes = {}
    for arm in ARMS:
        path = RESULTS / f"{arm}-stage-a-report.json"
        report = json.loads(path.read_text())
        assert report["arm"] == arm and report["stage"] == "lowinput-screen"
        panels[arm] = report["panels"]["lowinput"]
        rows = panels[arm]["individual_results"]
        assert {(r["dataset"], r["seed"], r["variant"]) for r in rows} == {
            (f, s, v) for f in FILES for s in SEEDS for v in ("baseline", "candidate")}
        assert len(rows) == 60
        hashes[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    baseline = {(r["dataset"], r["seed"]): r["internal_q1_targets"]
                for r in panels[ARMS[0]]["individual_results"] if r["variant"] == "baseline"}
    for panel in panels.values():
        assert baseline == {(r["dataset"], r["seed"]): r["internal_q1_targets"]
                            for r in panel["individual_results"] if r["variant"] == "baseline"}
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5), constrained_layout=True)
    colors = plt.get_cmap("tab10")(np.arange(len(SEEDS)))
    offsets = np.linspace(-.17, .17, len(SEEDS))
    for ax, dataset in zip(axes.flat, FILES):
        values = [[baseline[dataset, s] for s in SEEDS]]
        for arm in ARMS:
            lookup = {r["seed"]: r["internal_q1_targets"]
                      for r in panels[arm]["individual_results"]
                      if r["variant"] == "candidate" and r["dataset"] == dataset}
            values.append([lookup[s] for s in SEEDS])
        for x, counts in enumerate(values):
            ax.vlines(x, min(counts), max(counts), color="0.7", linewidth=1.3, zorder=1)
            ax.hlines(np.median(counts), x-.23, x+.23, color="0.2", linewidth=2, zorder=2)
            for i, count in enumerate(counts):
                ax.scatter(x+offsets[i], count, color=colors[i], s=28, zorder=3,
                           label=f"Seed {SEEDS[i]}" if x == 0 else None)
        ax.set_title(dataset.replace("astral_", "").replace("_", " "))
        ax.set_xticks(range(len(LABELS)), LABELS, fontsize=9)
        ax.set_ylabel("Reported target proteins (q ≤ 1%)")
        ax.set_ylim(bottom=-35)
        ax.grid(axis="y", alpha=.2)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("UNO soft calibration: all five reruns per file", fontsize=16)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center",
               ncol=5, bbox_to_anchor=(.5, -.045))
    for suffix in ("png", "pdf"):
        fig.savefig(RESULTS / f"stage-a-reported-counts.{suffix}", dpi=180, bbox_inches="tight")
    plt.close(fig)
    (RESULTS / "stage-a-figure-provenance.json").write_text(json.dumps({
        "input_sha256": hashes,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "caption": "Each dot is one parent seed; bars mark the median and lines span the range. Zero means an empty reported q≤1% target list. These are correlated algorithmic reruns on six acquisition files, not independent biological replicates. Equal seeds across arms do not imply identical internal random draws."
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
