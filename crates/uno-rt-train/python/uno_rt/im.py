"""Charge-aware IM training using the RT CNN and the same adaptation workflow."""
import argparse
import csv
import json
import math
import statistics
from dataclasses import dataclass

from .dataset import AA

UNITS = ("inverse_mobility", "proportional_inverse_mobility")

@dataclass(frozen=True)
class Row:
    sequence: str
    charge: int
    rt: float  # shared regression trainer's target field
    run_id: str = ""


def read_data(path, units):
    if units not in UNITS:
        raise ValueError(f"IM units must be one of {UNITS}; CCS must be converted explicitly")
    grouped = {}
    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fields = reader.fieldnames or []
        if not {"sequence", "charge", "im"} <= set(fields) or len(fields) != len(set(fields)):
            raise ValueError("require unique sequence, charge, im headers")
        for line, r in enumerate(reader, 2):
            s = (r.get("sequence") or "").strip()
            if not 1 <= len(s) <= 50 or not set(s) <= AA:
                raise ValueError(f"{path}:{line}: require 1–50 canonical residues with fixed CAM C")
            try:
                z, value = int(r["charge"]), float(r["im"])
            except (ValueError, TypeError) as e:
                raise ValueError(f"{path}:{line}: invalid charge or IM") from e
            if not 1 <= z <= 8 or not math.isfinite(value) or value <= 0:
                raise ValueError(f"{path}:{line}: charge must be 1–8 and IM finite and positive")
            run = (r.get("run_id") or "").strip()
            if ("run_id" in fields and not run) or None in r or any(v is None for v in r.values()):
                raise ValueError(f"{path}:{line}: malformed row or empty run_id")
            grouped.setdefault((s, z, run), []).append(value)
    if not grouped:
        raise ValueError("no IM observations")
    return [Row(s, z, statistics.median(v), run) for (s, z, run), v in sorted(grouped.items())]


def train(data, out, **kwargs):
    from .training import train as train_shared
    return train_shared(data, out, _target="im", **kwargs)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data", required=True); p.add_argument("--out", required=True)
    p.add_argument("--units", required=True, choices=UNITS); p.add_argument("--method", required=True)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--init"); g.add_argument("--from-scratch", action="store_true")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--validation"); g.add_argument("--holdout-run")
    p.add_argument("--seed", type=int, default=0); p.add_argument("--epochs", type=int, default=40)
    p.add_argument("--config")
    print(json.dumps(train(**vars(p.parse_args())), indent=2, allow_nan=False))

if __name__ == "__main__":
    main()
