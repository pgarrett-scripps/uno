#!/usr/bin/env python3
"""Reconstruct the trained RTCnn from a URTC binary and predict with PyTorch.

This is the parity oracle: the Rust URTC inference must reproduce these numbers.
Reads weights straight out of the URTC file into the same nn.Module that
produced the training metrics, so it exercises torch's exact conv/pool ops.
Usage: urtc_torch_predict.py MODEL.urtc SEQS.tsv OUT.tsv
"""
import struct, sys
import numpy as np
import torch
from train_rt import RTCnn, feat_tokens, NAA

path, seqs_file, out = sys.argv[1], sys.argv[2], sys.argv[3]
b = open(path, "rb").read()
p = 0
def u32():
    global p; v = struct.unpack_from("<I", b, p)[0]; p += 4; return v
def f32():
    global p; v = struct.unpack_from("<f", b, p)[0]; p += 4; return v
def arr(n):
    global p; v = np.frombuffer(b, "<f4", n, p).copy(); p += 4 * n; return v

assert b[:4] == b"URTC"; p = 4
ver = u32(); y_mean = f32(); y_std = f32()
emb = u32(); ch = u32(); k = u32()
nb = u32(); dils = [u32() for _ in range(nb)]; head = u32()
n_global = u32() if ver >= 2 else 0
head_in = 2 * ch + n_global
print(f"config: ver={ver} emb={emb} ch={ch} k={k} blocks={dils} head={head} n_global={n_global}", file=sys.stderr)
from train_rt import feat_global, N_GLOBAL
assert n_global in (0, N_GLOBAL), f"n_global {n_global} != {N_GLOBAL}"

model = RTCnn(emb=emb, ch=ch, k=k, blocks=dils, head=head)
with torch.no_grad():
    model.emb.weight.copy_(torch.from_numpy(arr((NAA + 1) * emb).reshape(NAA + 1, emb)))
    model.stem.weight.copy_(torch.from_numpy(arr(ch * emb * k).reshape(ch, emb, k)))
    model.stem.bias.copy_(torch.from_numpy(arr(ch)))
    for blk in model.blocks:
        blk.c1.weight.copy_(torch.from_numpy(arr(ch * ch * k).reshape(ch, ch, k)))
        blk.c1.bias.copy_(torch.from_numpy(arr(ch)))
        blk.c2.weight.copy_(torch.from_numpy(arr(ch * ch * k).reshape(ch, ch, k)))
        blk.c2.bias.copy_(torch.from_numpy(arr(ch)))
    h = list(model.head)
    h[0].weight.copy_(torch.from_numpy(arr(head * head_in).reshape(head, head_in)))
    h[0].bias.copy_(torch.from_numpy(arr(head)))
    h[2].weight.copy_(torch.from_numpy(arr(1 * head).reshape(1, head)))
    h[2].bias.copy_(torch.from_numpy(arr(1)))
assert p == len(b), f"trailing bytes: {p} != {len(b)}"

seqs = [l.split("\t")[0] for l in open(seqs_file)]
T, L = feat_tokens(seqs)
G = feat_global(seqs)
model.eval()
with torch.no_grad():
    preds = []
    for i in range(0, len(seqs), 8192):
        out_b = model(torch.from_numpy(T[i:i+8192]), torch.from_numpy(L[i:i+8192]), torch.from_numpy(G[i:i+8192]))
        preds.append(out_b.numpy())
preds = np.concatenate(preds) * y_std + y_mean
with open(out, "w") as o:
    for s, pr in zip(seqs, preds):
        o.write(f"{s}\t{pr}\n")
print(f"wrote {len(preds)} torch preds -> {out}", file=sys.stderr)
