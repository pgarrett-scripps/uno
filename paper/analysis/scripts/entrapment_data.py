#!/usr/bin/env python3
"""Entrapment FDR calibration data for fig4 / figS1.

Two DIFFERENT quantities are computed here and must never be conflated:

  reported q   -- the engine's OWN target/decoy FDR estimate.
                  uno : the ``q_value`` column it writes (running-min of
                        (d+1)/t over its shuffled-decoy protein list).
                  ms1searchpy : reproduced EXACTLY as ms1searchpy computes it --
                        target/decoy pair competition on the protein list, then
                        ``pyteomics.auxiliary.filter(..., formula=1,
                        correction=1)`` i.e. q = (d_cum + 1) / t_cum with the
                        monotonicity/tie fixups of
                        ``pyteomics.auxiliary.target_decoy._calculate_qvalues``.
                        Validated: the protein count at reported q <= 0.01
                        reproduces the row count of ms1searchpy's own
                        ``*_proteins.tsv`` exactly for all 18 runs.

  entrapment FDP -- the ground-truth-ish empirical false discovery proportion
                  from the yeast entrapment database. Computed with the
                  functions COPIED VERBATIM from benchmark/entrapment.py (see
                  ``load`` / ``curve`` below): the sample is human HeLa, the
                  FASTA is human + yeast, so every accepted _YEAST target is a
                  false positive, DB-size corrected by n_human_db/n_yeast_db,
                  and turned into a q-value by a running minimum scanned from
                  the deepest cutoff upward.

Data sources (all on disk, nothing transcribed):
  uno : benchmark/results/uno_<dataset>.tsv
  ms1 : benchmark/results/ms1_<dataset>_rep{1,2,3}/*_proteins_full.tsv
"""

from __future__ import annotations

import csv
import glob
import os
from bisect import bisect_right
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
RESULTS = REPO / "benchmark" / "results"

DATASETS = [
    "hela5min_09746",
    "hela5min_09766",
    "hela5min_09734",
    "hela5min_09778",
    "hela5min_09950",
    "faims_LUM2_05433_VG",
]

# Pretty labels for panel titles / annotations.
NICE = {
    "hela5min_09746": "HeLa 5 min #09746",
    "hela5min_09766": "HeLa 5 min #09766",
    "hela5min_09734": "HeLa 5 min #09734",
    "hela5min_09778": "HeLa 5 min #09778",
    "hela5min_09950": "HeLa 5 min #09950",
    "faims_LUM2_05433_VG": "FAIMS LUM2_05433",
}
SHORT = {
    "hela5min_09746": "09746",
    "hela5min_09766": "09766",
    "hela5min_09734": "09734",
    "hela5min_09778": "09778",
    "hela5min_09950": "09950",
    "faims_LUM2_05433_VG": "FAIMS",
}

MS1_DECOY_PREFIX = "DECOY_"

# --------------------------------------------------------------------------
# COPIED VERBATIM from benchmark/entrapment.py (the canonical scorer).
# Only change: `load` also returns the accession, needed to join the engine's
# own reported q onto the same rows. The FDP arithmetic is byte-identical.
# --------------------------------------------------------------------------
N_HUMAN_DB = int(os.environ.get("N_HUMAN_DB", 20647))
N_YEAST_DB = int(os.environ.get("N_YEAST_DB", 6066))
SCALE = N_HUMAN_DB / N_YEAST_DB


def load(path):
    """Return list of (score, is_yeast, accession) for target human/yeast proteins.

    Verbatim from benchmark/entrapment.py::load, plus the accession in the tuple.
    """
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        cols = r.fieldnames or []
        # detect schema
        if "binomial_score" in cols:
            acc_c, score_c = "accession", "binomial_score"
        elif "score" in cols and "dbname" in cols:
            acc_c, score_c = "dbname", "score"
        else:
            raise SystemExit(f"unrecognized schema in {path}: {cols}")
        rows = []
        for row in r:
            if str(row.get("decoy", "0")) not in ("0", "False", "false"):
                continue
            name = row[acc_c]
            is_y, is_h = "_YEAST" in name, "_HUMAN" in name
            if not (is_y or is_h):
                continue
            try:
                s = float(row[score_c])
            except (ValueError, KeyError):
                continue
            rows.append((s, is_y, name))
        return rows


def curve(rows):
    """Running-min entrapment q per accepted-target rank.

    Verbatim arithmetic from benchmark/entrapment.py::curve; instead of only
    returning the depth at 1%/5% it returns the full step list so the figures
    can plot it. ``rows`` are (score, is_yeast, ...) tuples.

    Returns (scores, human, yeast, q) as parallel lists, one entry per target
    protein in descending score order.
    """
    rows = sorted(rows, key=lambda r: -r[0])
    human = yeast = 0
    scores, hs, ys, fdp = [], [], [], []
    for row in rows:
        is_y = row[1]
        if is_y:
            yeast += 1
        else:
            human += 1
        if human:
            scores.append(row[0])
            hs.append(human)
            ys.append(yeast)
            fdp.append((yeast * SCALE) / human)
    q = [0.0] * len(fdp)
    rmin = float("inf")
    for i in range(len(fdp) - 1, -1, -1):
        rmin = min(rmin, fdp[i])
        q[i] = rmin
    return scores, hs, ys, q


# --------------------------------------------------------------------------
# engine-reported q
# --------------------------------------------------------------------------
class StepQ:
    """Reported q as a monotone step function of the score threshold."""

    def __init__(self, scores_desc, qs):
        self.s = list(scores_desc)          # descending
        self.q = list(qs)
        self._neg = [-x for x in self.s]    # ascending, for bisect

    def __call__(self, score):
        i = bisect_right(self._neg, -score) - 1
        if i < 0:
            return self.q[0] if self.q else 1.0
        return self.q[i]


def uno_reported_q(path) -> StepQ:
    """uno's own q_value column, as a step function of binomial_score."""
    rows = []
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            rows.append((float(row["binomial_score"]), float(row["q_value"])))
    rows.sort(key=lambda r: -r[0])
    return StepQ([r[0] for r in rows], [r[1] for r in rows])


def ms1_reported_q(path, prefix: str = MS1_DECOY_PREFIX) -> StepQ:
    """ms1searchpy's own protein FDR, reproduced exactly.

    ms1searchpy/main.py::final_iteration does target/decoy *pair competition*
    (the lower-scoring member of a target/DECOY_target pair is dropped) and then
    ``aux.filter(..., formula=1, correction=1)``. pyteomics' formula 1 with
    correction 1 and ratio 1 is q = (decoys_so_far + 1) / targets_so_far,
    followed by the tie/monotonicity sweep in
    ``pyteomics.auxiliary.target_decoy._calculate_qvalues``.
    """
    raw = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            raw[row["dbname"]] = float(row["score"])
    keep = dict(raw)
    for k, v in raw.items():                      # pair competition
        if k.startswith(prefix):
            if raw.get(k[len(prefix):], -1e6) > v:
                keep.pop(k, None)
        else:
            if raw.get(prefix + k, -1e6) > v:
                keep.pop(k, None)
    items = sorted(keep.items(), key=lambda x: -x[1])
    d = t = 0
    q = []
    for k, _ in items:
        if k.startswith(prefix):
            d += 1
        else:
            t += 1
        q.append((d + 1) / t if t else float("inf"))
    sc = [v for _, v in items]
    for i in range(len(q) - 1, 0, -1):             # pyteomics tie/monotone sweep
        if sc[i] == sc[i - 1] or q[i - 1] > q[i]:
            q[i - 1] = q[i]
    return StepQ(sc, q)


# --------------------------------------------------------------------------
# assembly
# --------------------------------------------------------------------------
class Run:
    """One engine run on one dataset: paired reported-q and entrapment-FDP steps."""

    def __init__(self, dataset, engine, rep, path):
        self.dataset, self.engine, self.rep, self.path = dataset, engine, rep, path
        rows = load(path)
        self.scores, self.human, self.yeast, self.entq = curve(rows)
        step = uno_reported_q(path) if engine == "uno" else ms1_reported_q(path)
        self.repq = [step(s) for s in self.scores]
        # reported q is monotone in rank by construction, but enforce it so
        # `depth_at_reported` can use a simple scan.
        run = 0.0
        for i, v in enumerate(self.repq):
            run = max(run, v)
            self.repq[i] = run

    # -- depth lookups -----------------------------------------------------
    def depth_at_entrapment(self, thr):
        """Most human proteins accepted with entrapment q <= thr (entrapment.py rule)."""
        best = 0
        for h, q in zip(self.human, self.entq):
            if q <= thr and h > best:
                best = h
        return best

    def _last_index_reported(self, thr):
        idx = -1
        for i, q in enumerate(self.repq):
            if q <= thr:
                idx = i
        return idx

    def at_reported(self, thr):
        """(n_human, n_yeast, entrapment_q) at the engine's own q <= thr cutoff."""
        i = self._last_index_reported(thr)
        if i < 0:
            return 0, 0, 0.0
        return self.human[i], self.yeast[i], self.entq[i]

    def calibration(self):
        """(reported_q, entrapment_q) step pairs, deduplicated on reported q."""
        xs, ys = [], []
        last = None
        for rq, eq in zip(self.repq, self.entq):
            if rq != last:
                xs.append(rq)
                ys.append(eq)
                last = rq
            else:
                ys[-1] = eq          # keep the deepest cutoff at this reported q
        return xs, ys

    def depth_curve(self):
        """(entrapment_q, n_human) monotone accepted-depth curve."""
        xs, ys = [], []
        best = 0
        for h, q in zip(self.human, self.entq):
            if h > best:
                best = h
                xs.append(q)
                ys.append(h)
        return xs, ys


def ms1_path(dataset, rep):
    hits = sorted(glob.glob(str(RESULTS / f"ms1_{dataset}_rep{rep}" / "*_proteins_full.tsv")))
    hits = [h for h in hits if "noexclusion" not in os.path.basename(h)]
    if not hits:
        raise SystemExit(f"no ms1 protein output for {dataset} rep{rep}")
    return hits[0]


def load_all():
    """{dataset: {'uno': Run, 'ms1': [Run, Run, Run]}} for all six datasets."""
    out = {}
    for ds in DATASETS:
        uno_p = RESULTS / f"uno_{ds}.tsv"
        if not uno_p.exists():
            raise SystemExit(f"missing {uno_p}")
        out[ds] = {
            "uno": Run(ds, "uno", 1, str(uno_p)),
            "ms1": [Run(ds, "ms1", r, ms1_path(ds, r)) for r in (1, 2, 3)],
        }
    return out


TABLE_GRID = (0.001, 0.002, 0.005, 0.01, 0.02, 0.03, 0.05, 0.10)


def write_table(path, data=None):
    """Emit paper/figdata/entrapment_curves.tsv -- the numbers fig4/figS1 plot.

    Long format, one row per (dataset, engine, replicate, cutoff_type, cutoff):

      cutoff_type = reported_q    the cutoff is applied to the ENGINE'S OWN
                                  q-value; entrapment_fdp is then the empirical
                                  truth at that cutoff (fig4a / figS1 y-axis).
      cutoff_type = entrapment_q  the cutoff is applied to the entrapment FDP
                                  itself; n_human is then the honest depth
                                  quoted in benchmark/RESULTS_canonical.md
                                  (fig4b y-axis).
    """
    data = data or load_all()
    lines = [
        "# entrapment_curves.tsv -- generated by paper/scripts/figures/entrapment_data.py",
        "# Provenance: benchmark/results/uno_<dataset>.tsv and",
        "#             benchmark/results/ms1_<dataset>_rep{1,2,3}/*_proteins_full.tsv",
        "# FDP definition copied verbatim from benchmark/entrapment.py:",
        f"#   FDP(k) = yeast_in_top_k * ({N_HUMAN_DB}/{N_YEAST_DB}) / human_in_top_k,",
        "#   q = running minimum of FDP scanned from the deepest cutoff upward.",
        "# reported_q: uno = its own `q_value` column; ms1searchpy = its own protein FDR",
        "#   (target/decoy pair competition + pyteomics aux.filter formula=1 correction=1),",
        "#   reproduced exactly -- protein count at reported_q<=0.01 matches ms1searchpy's",
        "#   own *_proteins.tsv row count for all 18 runs.",
        "# n_human / n_yeast are TARGET proteins accepted at the cutoff.",
        "\t".join(["dataset", "engine", "rep", "cutoff_type", "cutoff",
                   "n_human", "n_yeast", "entrapment_fdp"]),
    ]
    for ds in DATASETS:
        runs = [("uno", data[ds]["uno"])] + [("ms1searchpy", m) for m in data[ds]["ms1"]]
        for eng, run in runs:
            for t in TABLE_GRID:
                h, y, e = run.at_reported(t)
                lines.append(f"{ds}\t{eng}\t{run.rep}\treported_q\t{t:g}\t{h}\t{y}\t{e:.6f}")
            for t in TABLE_GRID:
                h = run.depth_at_entrapment(t)
                lines.append(f"{ds}\t{eng}\t{run.rep}\tentrapment_q\t{t:g}\t{h}\t\t{t:.6f}")
    Path(path).write_text("\n".join(lines) + "\n")
    print(f"wrote {path} ({len(lines) - 12} data rows)")


if __name__ == "__main__":
    data = load_all()
    for ds in DATASETS:
        u = data[ds]["uno"]
        h, y, e = u.at_reported(0.01)
        print(f"{ds:22s} uno  rep-q<=1% -> {h:5d} human {y:4d} yeast  entrapFDP={e:.4f}"
              f"  depth@1%entrap={u.depth_at_entrapment(0.01)}")
        for m in data[ds]["ms1"]:
            h, y, e = m.at_reported(0.01)
            print(f"{'':22s} ms1{m.rep} rep-q<=1% -> {h:5d} human {y:4d} yeast  "
                  f"entrapFDP={e:.4f}  depth@1%entrap={m.depth_at_entrapment(0.01)}")
