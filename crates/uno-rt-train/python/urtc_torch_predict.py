#!/usr/bin/env python3
"""PyTorch/Rust parity oracle, supports URTC v1/v2/v3.
Usage: urtc_torch_predict.py MODEL.urtc SEQS.tsv OUT.tsv
"""
import sys
from uno_rt.artifact import load_model
from uno_rt.training import predict

path, sequences_file, out = sys.argv[1:]
model, info, _ = load_model(path)
with open(sequences_file) as stream:
    sequences = [line.strip().split("\t")[0] for line in stream if line.strip()]
values = predict(model, sequences, info)
with open(out, "w") as stream:
    for sequence, value in zip(sequences, values):
        stream.write(f"{sequence}\t{value}\n")
