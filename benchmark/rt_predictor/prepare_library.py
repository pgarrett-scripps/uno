#!/usr/bin/env python3
"""Build the synthetic tryptic peptide library used for the throughput benchmark.

Throughput depends on the peptides, not on which experimental dataset they came
from, so the timing library is an in-silico tryptic digest of the human entries
in the benchmark FASTA: cleavage after K/R (no proline rule, no missed
cleavages), canonical residues only, length 7--30, unique sequences, ordered by
SHA-256 so the file is deterministic. The first N rows are the N-peptide batch.
"""

from __future__ import annotations

import hashlib
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
FASTA = HERE.parent / "data" / "entrap_human_yeast.fasta"
OUT = HERE / "sets" / "synthetic_tryptic.tsv"
MAX_PEPTIDES = 100_000
CANONICAL = set("ACDEFGHIKLMNPQRSTVWY")


def proteins(path: Path):
    header, chunks = None, []
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:].strip(), []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def main() -> None:
    peptides: set[str] = set()
    for header, seq in proteins(FASTA):
        if "_HUMAN" not in header.split()[0]:
            continue
        for pep in re.split(r"(?<=[KR])", seq):
            if 7 <= len(pep) <= 30 and set(pep) <= CANONICAL:
                peptides.add(pep)
    ranked = sorted(peptides, key=lambda p: hashlib.sha256(p.encode()).digest())[:MAX_PEPTIDES]
    with OUT.open("w") as handle:
        handle.write("# uno RT predictor throughput library v1\n")
        handle.write(f"# source: in-silico tryptic digest of human entries in {FASTA.name}\n")
        handle.write("# filter: cleave after K/R, no missed cleavages; canonical residues; length 7-30; unique\n")
        handle.write(f"# selection: lowest deterministic SHA-256 ranks, maximum {MAX_PEPTIDES}\n")
        handle.write("sequence\n")
        handle.write("\n".join(ranked) + "\n")
    print(f"wrote {OUT} ({len(ranked)} of {len(peptides)} unique tryptic peptides)")


if __name__ == "__main__":
    main()
