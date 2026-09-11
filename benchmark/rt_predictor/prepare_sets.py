#!/usr/bin/env python3
"""Build the sequence-disjoint, unmodified-peptide RT benchmark sets.

The large public source files live under ``data/rt`` and are intentionally
gitignored.  This script distills them into small, versioned benchmark sets.
Peptides seen in the Chronologer train or validation partitions are excluded
from every set, repeated measurements are collapsed to their median, and a
stable SHA-256 ordering selects at most 10,000 sequences per source.

"Unmodified" here permits the fixed carbamidomethyl cysteine assumed by both
uno and the DeepLC invocation in ``run.py``; variable modifications are absent.
"""

from __future__ import annotations

import csv
import hashlib
import statistics
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RT_DATA = ROOT / "data" / "rt"
OUT = Path(__file__).resolve().parent / "sets"
AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
MAX_ROWS = 10_000


def load_seen() -> set[str]:
    seen: set[str] = set()
    for name in ("train.tsv", "val.tsv"):
        with (RT_DATA / "chronologer" / name).open() as handle:
            for line in handle:
                seq = line.split("\t", 1)[0]
                if seq:
                    seen.add(seq)
    return seen


def valid(seq: str, seen: set[str]) -> bool:
    return 7 <= len(seq) <= 30 and set(seq) <= AA and seq not in seen


def stable_sample(name: str, values: dict[str, float]) -> list[tuple[str, float]]:
    def key(item: tuple[str, float]) -> bytes:
        seq = item[0]
        return hashlib.sha256(f"uno-rt-benchmark-v1\0{name}\0{seq}".encode()).digest()

    return sorted(values.items(), key=key)[:MAX_ROWS]


def write_set(name: str, rows: list[tuple[str, float]], source: str, label: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"{name}.tsv"
    with path.open("w") as out:
        out.write("# uno RT predictor benchmark v1\n")
        out.write(f"# source: {source}\n")
        out.write("# filter: canonical residues; length 7-30; no variable modifications; "
                  "sequence absent from Chronologer train+validation\n")
        out.write("# repeated sequence labels: median; selection: lowest deterministic "
                  "SHA-256 ranks, maximum 10000\n")
        out.write(f"# label: {label}\n")
        out.write("sequence\tretention_time\n")
        for seq, rt in rows:
            out.write(f"{seq}\t{rt:.10g}\n")
    print(f"wrote {path} ({len(rows)} peptides)")


def chronologer(seen: set[str]) -> None:
    values: dict[str, float] = {}
    source = RT_DATA / "chronologer" / "holdout.tsv"
    with source.open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) >= 2 and valid(fields[0], seen):
                values[fields[0]] = float(fields[1])
    write_set(
        "chronologer_holdout",
        stable_sample("chronologer_holdout", values),
        "searlelab/chronologer Chronologer_DB_220308.gz, held-out partition",
        "harmonized hydrophobic index",
    )


def prosit(seen: set[str]) -> None:
    grouped: dict[str, list[float]] = defaultdict(list)
    source = RT_DATA / "prosit" / "prosit_irt.csv"
    with source.open(newline="") as handle:
        for row in csv.DictReader(handle):
            seq = row["sequence"]
            # mod_sequence == sequence excludes oxidation and every variable mod.
            if row["state"] == "holdout" and row["mod_sequence"] == seq and valid(seq, seen):
                grouped[seq].append(float(row["irt"]))
    values = {seq: statistics.median(rts) for seq, rts in grouped.items()}
    write_set(
        "prosit_unmodified_holdout",
        stable_sample("prosit_unmodified_holdout", values),
        "ProteomeTools--Prosit iRT data (figshare article 10043132), source holdout rows",
        "median indexed retention time",
    )


def deeplc_unmodified(seen: set[str]) -> None:
    grouped: dict[str, list[float]] = defaultdict(list)
    source = RT_DATA / "deeplc" / "unmod.csv"
    with source.open(newline="") as handle:
        for row in csv.DictReader(handle):
            seq = row["seq"]
            if not row["modifications"].strip() and valid(seq, seen):
                grouped[seq].append(float(row["tr"]))
    values = {seq: statistics.median(rts) for seq, rts in grouped.items()}
    write_set(
        "deeplc_unmodified",
        stable_sample("deeplc_unmodified", values),
        "CompOmics/DeepLC examples/datasets/unmod.csv",
        "median experimental retention time",
    )


def main() -> None:
    seen = load_seen()
    chronologer(seen)
    prosit(seen)
    deeplc_unmodified(seen)


if __name__ == "__main__":
    main()
