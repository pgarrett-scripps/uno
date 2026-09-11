"""URTC v1/v2 compatibility and v3 self-describing, inference-only artifacts.

v3 extends the v2 header after n_global with u32 JSON length + UTF-8 JSON.
Weights remain little-endian float32, in PyTorch tensor order. No pickle.
"""
import json
import math
import os
import struct
import tempfile
from pathlib import Path

import numpy as np
import torch

from .network import RTCnn

FEATURE_SCHEMA = "uno-aa-global-v1"
ENCODING = "canonical-aa-fixed-cam"


def rebase_output(model, info, mean, std, slope=1.0, intercept=0.0):
    """Preserve physical predictions (or an explicit affine transform) on a new label scale."""
    with torch.no_grad():
        final = model.head[-1]
        factor = slope * info["y_std"] / std
        final.weight.mul_(factor)
        final.bias.mul_(factor).add_((slope * info["y_mean"] + intercept - mean) / std)
    return dict(info, y_mean=mean, y_std=std)


def load_model(path):
    blob = Path(path).read_bytes()
    if blob[:4] != b"URTC":
        raise ValueError("only URTC CNN models can be trained")
    pos = 4

    def unpack(fmt):
        nonlocal pos
        size = struct.calcsize(fmt)
        if pos + size > len(blob):
            raise ValueError("truncated URTC model")
        value = struct.unpack_from(fmt, blob, pos)
        pos += size
        return value

    version, mean, std, emb, ch, k, nb = unpack("<IffIIII")
    if version not in (1, 2, 3):
        raise ValueError(f"unsupported URTC version {version}")
    if not (math.isfinite(mean) and math.isfinite(std) and std > 0):
        raise ValueError("invalid URTC normalization")
    if not (1 <= emb <= 1024 and 1 <= ch <= 1024 and 1 <= k <= 31 and k % 2 and nb <= 32):
        raise ValueError("invalid or unsupported CNN dimensions")
    blocks = list(unpack(f"<{nb}I"))
    head, = unpack("<I")
    ng = unpack("<I")[0] if version >= 2 else 0
    if not 1 <= head <= 4096 or any(not 1 <= d <= 1024 for d in blocks) or ng not in (0, 5, 11):
        raise ValueError("invalid or unsupported CNN architecture/feature schema")
    metadata = {}
    if version == 3:
        size, = unpack("<I")
        if size > 64 * 1024 * 1024 or pos + size > len(blob):
            raise ValueError("invalid URTC metadata length")
        metadata = json.loads(blob[pos:pos + size])
        pos += size
        if not isinstance(metadata, dict) or metadata.get("feature_schema") != FEATURE_SCHEMA or metadata.get("sequence_encoding") != ENCODING:
            raise ValueError("unsupported URTC feature schema or sequence encoding")
    cfg = dict(emb=emb, ch=ch, k=k, blocks=blocks, head=head, n_global=ng)
    # Validate weight length before allocating the network.
    count = 21 * emb + ch * emb * k + ch + nb * 2 * (ch * ch * k + ch)
    count += head * (2 * ch + ng) + head + head + 1
    if len(blob) - pos != count * 4:
        raise ValueError("URTC weights are truncated or have trailing bytes")
    model = RTCnn(**cfg)
    with torch.no_grad():
        for parameter in model.parameters():
            size = parameter.numel()
            values = np.frombuffer(blob, dtype="<f4", count=size, offset=pos).copy()
            if not np.isfinite(values).all():
                raise ValueError("nonfinite CNN weights")
            parameter.copy_(torch.from_numpy(values.reshape(tuple(parameter.shape))))
            pos += 4 * size
    return model.eval(), dict(cfg, version=version, y_mean=mean, y_std=std), metadata


def save_model(path, model, info, metadata):
    path = Path(path)
    meta = dict(metadata, feature_schema=FEATURE_SCHEMA, sequence_encoding=ENCODING)
    encoded = json.dumps(meta, allow_nan=False, separators=(",", ":")).encode()
    if len(encoded) > 64 * 1024 * 1024:
        raise ValueError("URTC metadata exceeds 64 MiB")
    # Write completely before publishing; never replace an existing model.
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as stream:
            temporary = stream.name
            stream.write(b"URTC")
            stream.write(struct.pack("<IffIIII", 3, info["y_mean"], info["y_std"], info["emb"], info["ch"], info["k"], len(info["blocks"])))
            stream.write(struct.pack(f"<{len(info['blocks'])}I", *info["blocks"]))
            stream.write(struct.pack("<III", info["head"], info["n_global"], len(encoded)))
            stream.write(encoded)
            for parameter in model.parameters():
                stream.write(parameter.detach().cpu().numpy().astype("<f4").tobytes())
        # The reader verifies the exact artifact before it becomes user-visible.
        load_model(temporary)
        os.link(temporary, path)
    finally:
        if temporary is not None:
            os.unlink(temporary)
