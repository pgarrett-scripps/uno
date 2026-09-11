"""Strict, headered RT data and deterministic peptide-disjoint splits."""
import csv
import hashlib
import math
import random
import statistics
from dataclasses import dataclass
from pathlib import Path

AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
UNITS = ("minutes", "seconds", "irt", "hi")


@dataclass(frozen=True)
class Row:
    sequence: str
    rt: float
    run_id: str = ""


def peptide_hash(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_data(path, units):
    if units not in UNITS:
        raise ValueError(f"units must be one of {UNITS}")
    grouped = {}
    with open(path, newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fields = reader.fieldnames or []
        if not {"sequence", "rt"}.issubset(fields) or len(fields) != len(set(fields)):
            raise ValueError(f"{path}: require unique column headers including sequence and rt")
        for line, record in enumerate(reader, 2):
            seq = (record.get("sequence") or "").strip()
            if not 1 <= len(seq) <= 50 or not set(seq) <= AA:
                raise ValueError(f"{path}:{line}: sequence must contain 1–50 canonical amino acids; variable modifications are unsupported")
            try:
                rt = float(record["rt"])
            except (ValueError, TypeError) as error:
                raise ValueError(f"{path}:{line}: invalid RT") from error
            if not math.isfinite(rt) or (units in ("minutes", "seconds") and rt < 0):
                raise ValueError(f"{path}:{line}: RT must be finite and time values nonnegative")
            run = (record.get("run_id") or "").strip()
            if "run_id" in fields and not run:
                raise ValueError(f"{path}:{line}: run_id must be populated for every row")
            if None in record or any(value is None for value in record.values()):
                raise ValueError(f"{path}:{line}: malformed TSV row")
            grouped.setdefault((seq, run), []).append(rt)
    if not grouped:
        raise ValueError(f"{path}: no RT observations")
    if units in ("minutes", "seconds") and len({run for _, run in grouped}) > 1:
        raise ValueError("multiple runs with raw time units require prior normalization to a shared irt/hi scale; raw gradients must not be pooled")
    return [Row(seq, statistics.median(values), run) for (seq, run), values in sorted(grouped.items())]


def split_data(rows, validation=None, holdout_run=None, seed=0, val_fraction=0.2):
    if validation is not None and holdout_run is not None:
        raise ValueError("validation and holdout_run are mutually exclusive")
    if validation is not None:
        train, val = list(rows), list(validation)
        if {r.sequence for r in train} & {r.sequence for r in val}:
            raise ValueError("training and validation contain overlapping peptides")
    elif holdout_run is not None:
        val = [r for r in rows if r.run_id == holdout_run]
        held = {r.sequence for r in val}
        train = [r for r in rows if r.run_id != holdout_run and r.sequence not in held]
    else:
        if not 0 < val_fraction < 1:
            raise ValueError("val_fraction must be between zero and one")
        seqs = sorted({r.sequence for r in rows})
        random.Random(seed).shuffle(seqs)
        count = max(2, round(len(seqs) * val_fraction))
        held = set(seqs[:count])
        train = [r for r in rows if r.sequence not in held]
        val = [r for r in rows if r.sequence in held]
    if min(len({r.sequence for r in train}), len({r.sequence for r in val})) < 2:
        raise ValueError("need at least two distinct peptides in each split; add data or revise the holdout")
    return train, val
