#!/usr/bin/env python3
"""Plot every paired low-input seed, without modifying manuscript figures.

Run with the paper/analysis environment after stability_study.py completes.
"""
from pathlib import Path
import argparse
import hashlib
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from stability_report import load_study


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    study, rows, fingerprints = load_study(args.study)
    datasets = sorted({r["dataset"] for r in rows}, key=lambda d: ("200pg" not in d, d))
    if len(datasets) != 6:
        raise ValueError("expected the six-file low-input panel")
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42})
    colors = ("#56616D", "#007F86")
    args.output.mkdir(parents=True, exist_ok=True)
    generated = []
    def seed_label(values):
        values = sorted(values)
        return (f"{values[0]}–{values[-1]}" if len(values) > 1 and values == list(range(values[0], values[-1]+1))
                else ", ".join(map(str, values)))
    unseen = sorted(set(study["seeds"]) & set(study["held_out_seeds"]))
    development = sorted(set(study["seeds"]) - set(unseen))
    metrics = (("depth_1pct", "Human proteins at 1% entrapment threshold"),
               ("depth_5pct", "Human proteins at 5% entrapment threshold"),
               ("internal_q1_human", "Human proteins reported at q ≤ 1%"))
    for metric, ylabel in metrics:
        fig, axes = plt.subplots(2, 3, figsize=(8.4, 6.1), sharey=True)
        ymax = max(r[metric] for r in rows)
        for ax, ds in zip(axes.flat, datasets):
            data = {(r["variant"], r["seed"]): r for r in rows if r["dataset"] == ds}
            seeds = sorted({seed for _, seed in data})
            for j, seed in enumerate(seeds):
                pair = [data[v, seed] for v in ("baseline", "candidate")]
                offset = (j - (len(seeds) - 1) / 2) * 0.055
                xx, yy = [offset, 1 + offset], [r[metric] for r in pair]
                ax.plot(xx, yy, color="#CDD3D6", lw=0.7, zorder=1)
                for i, row in enumerate(pair):
                    new_seed = row["seed_role"] == "held_out"
                    ax.scatter(xx[i], yy[i], marker="D" if new_seed else "o", s=30,
                               facecolor="white" if new_seed else colors[i],
                               edgecolor=colors[i], linewidth=1.2, zorder=3)
            for i, variant in enumerate(("baseline", "candidate")):
                values = sorted(data[variant, seed][metric] for seed in seeds)
                n = len(values)
                med = (values[(n-1)//2] + values[n//2]) / 2
                ax.plot([i-0.17, i+0.17], [med, med], color=colors[i], lw=2.2, zorder=4)
            label = ds.replace("astral_", "").replace("200pg", "200 pg").replace("1ng", "1 ng").replace("_rep", " · replicate ")
            ax.set_title(label, loc="left", fontsize=10)
            ax.set_xticks([0, 1], ["Baseline", "Candidate"])
            ax.set_xlim(-0.35, 1.35)
            ax.set_ylim(-ymax*0.04, ymax*1.10)
            ax.grid(axis="y", color="#ECEFF1", lw=0.7)
            ax.set_axisbelow(True)
        fig.supylabel(ylabel, x=0.015, fontsize=10)
        fig.suptitle("UNO low-input seed stability", fontsize=14, x=0.08, ha="left")
        legend = [Line2D([], [], marker="o", ls="", color="#56616D", label=f"Development seeds {seed_label(development)}"),
                  Line2D([], [], marker="D", ls="", markerfacecolor="white", color="#56616D", label=f"Unseen seeds {seed_label(unseen)}"),
                  Line2D([], [], color="#56616D", lw=2.2, label="Median")]
        fig.legend(handles=legend, loc="upper center", bbox_to_anchor=(0.54, 0.943), ncol=3, frameon=False)
        fig.text(0.08, 0.03, f"Lines connect the same seed. All {len(study['seeds'])} draws are shown, including zero-depth searches.\n"
                 "Unseen seeds test new randomness on the same files. Reported-list FDP is evaluated separately.", fontsize=8)
        fig.subplots_adjust(top=0.84, bottom=0.14, left=0.09, right=0.98, hspace=0.34, wspace=0.17)
        for ext in ("png", "pdf"):
            path = args.output / f"lowinput_{metric}.{ext}"
            metadata = {"CreationDate": None, "ModDate": None} if ext == "pdf" else None
            fig.savefig(path, dpi=180, metadata=metadata)
            generated.append(path)
        plt.close(fig)
    sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
    (args.output / "plot_provenance.json").write_text(json.dumps({
        "script_sha256": sha(__file__), "input_hashes": fingerprints,
        "outputs": {str(p.resolve()): sha(p) for p in generated},
        "matplotlib": matplotlib.__version__, "seeds": study["seeds"],
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
