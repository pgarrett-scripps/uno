"""CNN adaptation with peptide-disjoint validation and training-only alignment."""
import copy
import json
import math
import os
import warnings
from datetime import datetime, timezone
from importlib.resources import files
from pathlib import Path

import numpy as np
# Deterministic CUDA linear layers require this before their first invocation.
os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
import torch

from . import __version__
from .artifact import load_model, save_model, rebase_output
from .dataset import file_hash, peptide_hash, read_data, split_data
from .network import RTCnn, feat_global, feat_tokens


def predict(model, seqs, info, batch_size=512, device="cpu"):
    model.eval()
    result = []
    with torch.no_grad():
        for start in range(0, len(seqs), batch_size):
            inputs = features(seqs[start:start + batch_size], device)
            result.append(model(*inputs).cpu().numpy())
    return np.concatenate(result).astype(np.float64) * info["y_std"] + info["y_mean"]


def features(seqs, device):
    is_im = bool(seqs) and isinstance(seqs[0], tuple)
    sequences = [s for s, _ in seqs] if is_im else seqs
    tokens, lengths = feat_tokens(sequences)
    glob = feat_global(sequences)
    if is_im:
        glob = np.column_stack((glob, np.array([z / 4 for _, z in seqs], dtype=np.float32)))
    return tuple(torch.from_numpy(v).to(device) for v in (tokens, lengths, glob))


def metrics(y, p):
    if not np.isfinite(p).all():
        raise ValueError("training produced nonfinite predictions; check labels and learning rate")
    error = y - p
    total = float(np.sum((y - y.mean()) ** 2))
    return dict(mae=float(np.mean(np.abs(error))), rmse=float(np.sqrt(np.mean(error ** 2))),
                p95_absolute_error=float(np.percentile(np.abs(error), 95)),
                r2=1 - float(np.sum(error ** 2)) / total if total > 0 else None)


def align_on_training(y, predictions):
    variance = float(np.sum((predictions - predictions.mean()) ** 2))
    slope = float(np.sum((predictions - predictions.mean()) * (y - y.mean()))) / variance if variance > 1e-12 else 1.0
    return slope, float(y.mean() - slope * predictions.mean())


def train(data, out, *, units, method, init=None, from_scratch=False, validation=None,
          holdout_run=None, seed=0, epochs=40, config=None, _target="rt"):
    """Train a method-specific CNN and atomically export a URTC v3 model.

    ``config`` is an optional dict (or JSON file) for batch_size, lr, patience,
    val_fraction, device, threads and architecture. Architecture overrides are
    only accepted with from_scratch. Validation is for model selection; use an
    independent dataset with `uno rt evaluate` to estimate generalization.
    """
    if init is not None and from_scratch:
        raise ValueError("init and from_scratch are mutually exclusive")
    if not method or not method.strip():
        raise ValueError("provide a method name describing this column/separation")
    if Path(out).exists():
        raise ValueError(f"output already exists: {out}; choose a new model path")
    if epochs < 1 or seed < 0 or seed >= 2 ** 32:
        raise ValueError("epochs must be positive and seed must be in 0..2^32")
    options = json.loads(Path(config).read_text()) if isinstance(config, (str, Path)) else dict(config or {})
    allowed = {"batch_size", "lr", "patience", "val_fraction", "device", "threads", "architecture"}
    if set(options) - allowed:
        raise ValueError(f"unknown advanced options: {sorted(set(options) - allowed)}")
    bs, patience = options.get("batch_size", 512), options.get("patience", 6)
    lr = options.get("lr", 1e-3 if from_scratch else 1e-4)
    threads = options.get("threads", 4)
    if any(not isinstance(v, int) or v <= 0 for v in (bs, patience, threads)) or not math.isfinite(lr) or lr <= 0:
        raise ValueError("batch_size, patience, threads and lr must be positive")
    if "architecture" in options and not from_scratch:
        raise ValueError("architecture is read from --init; overrides require --from-scratch")
    reader = read_data
    if _target == "im":
        from .im import read_data as reader
        if init is None and not from_scratch:
            raise ValueError("IM training requires --init or --from-scratch; no bundled IM model is assumed")
    rows = reader(data, units)
    val_rows = reader(validation, units) if validation else None
    tr, va = split_data(rows, val_rows, holdout_run, seed, options.get("val_fraction", .2))
    if len(tr) < 100:
        warnings.warn("fewer than 100 training observations: use an independent test set before adopting this model")
    tr_seqs, va_seqs = [r.sequence for r in tr], [r.sequence for r in va]
    tr_inputs = [(r.sequence, r.charge) for r in tr] if _target == "im" else tr_seqs
    va_inputs = [(r.sequence, r.charge) for r in va] if _target == "im" else va_seqs
    y = np.array([r.rt for r in tr], dtype=np.float64)
    vy = np.array([r.rt for r in va], dtype=np.float64)
    mean, std = float(y.mean()), float(y.std())
    if not math.isfinite(std) or std <= 1e-8 or max(abs(mean), std) > 1e30:
        raise ValueError("training RTs need finite, representable nonzero variance")
    device = options.get("device", "cuda" if torch.cuda.is_available() else "cpu")
    if device not in ("cpu", "cuda") or (device == "cuda" and not torch.cuda.is_available()):
        raise ValueError("device must be cpu or an available cuda device")
    torch.manual_seed(seed)
    np.random.seed(seed)
    torch.set_num_threads(threads)
    torch.use_deterministic_algorithms(True)
    baseline = None
    ancestry = set()
    parent_hash = None
    if from_scratch:
        info = dict(emb=32, ch=32, k=5, blocks=[1, 2, 4], head=64, n_global=12 if _target == "im" else 11)
        architecture = options.get("architecture", {})
        if set(architecture) - set(info):
            raise ValueError("unknown CNN architecture fields")
        info.update(architecture)
        if not (1 <= info["emb"] <= 1024 and 1 <= info["ch"] <= 1024 and 1 <= info["head"] <= 4096
                and 1 <= info["k"] <= 31 and info["k"] % 2 and len(info["blocks"]) <= 32
                and all(1 <= d <= 1024 for d in info["blocks"]) and info["n_global"] in ((12,) if _target == "im" else (0, 5, 11))):
            raise ValueError("invalid CNN architecture")
        model = RTCnn(**info)
    else:
        source = Path(init) if init else Path(str(files("uno_rt").joinpath("data/default.urtc")))
        model, info, parent = load_model(source, target=_target)
        parent_hash = file_hash(source)
        ancestry = set(parent.get("seen_peptide_hashes", []))
        if ancestry & {peptide_hash(s) for s in va_seqs}:
            raise ValueError("validation peptides were used to train/select the starting model; use unseen validation data")
        training_preds = predict(model, tr_inputs, info, bs)
        slope, intercept = align_on_training(y, training_preds)
        # Re-express the old physical predictions on the new normalization.
        # Calibration is fitted ONLY on training observations, never validation.
        info = rebase_output(model, info, mean, std, slope, intercept)
    info.update(y_mean=mean, y_std=std)
    model.to(device)
    if not from_scratch:
        baseline = metrics(vy, predict(model, va_inputs, info, bs, device))
    # Include epoch zero: unsuccessful adaptation retains the calibrated parent.
    best_mae = baseline["mae"] if baseline else float("inf")
    best_state = copy.deepcopy(model.state_dict()) if baseline else None
    best_epoch, stale = 0, 0
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    target = (y - mean) / std
    history = []
    for epoch in range(1, epochs + 1):
        model.train()
        order = np.random.permutation(len(tr))
        for start in range(0, len(tr), bs):
            indices = order[start:start + bs]
            inputs = features([tr_inputs[i] for i in indices], device)
            labels = torch.tensor(target[indices], dtype=torch.float32, device=device)
            optimizer.zero_grad()
            loss = torch.nn.functional.smooth_l1_loss(model(*inputs), labels)
            if not torch.isfinite(loss):
                raise ValueError("nonfinite training loss")
            loss.backward()
            optimizer.step()
        result = metrics(vy, predict(model, va_inputs, info, bs, device))
        history.append(dict(epoch=epoch, **result))
        if result["mae"] < best_mae:
            best_mae, best_epoch, stale = result["mae"], epoch, 0
            best_state = copy.deepcopy(model.state_dict())
        else:
            stale += 1
        print(f"epoch {epoch}: validation MAE={result['mae']:.5g} {units}", flush=True)
        if stale >= patience:
            break
    model.load_state_dict(best_state)
    report = dict(units=units, train_observations=len(tr), validation_observations=len(va),
                  split="explicit" if validation else "heldout_run" if holdout_run else "peptide",
                  holdout_run=holdout_run, best_epoch=best_epoch, baseline=baseline,
                  validation=metrics(vy, predict(model, va_inputs, info, bs, device)), history=history,
                  interpretation="Validation selected the checkpoint; evaluate independent data and search performance before adoption.")
    metadata = dict(rt_units=units, method=method.strip(), created_at=datetime.now(timezone.utc).isoformat(),
                    trainer_version=__version__, torch_version=torch.__version__, numpy_version=np.__version__,
                    seed=seed, epochs=epochs, options=dict(options, device=device, lr=lr, batch_size=bs, patience=patience, threads=threads),
                    source_sha256=file_hash(data), validation_sha256=file_hash(validation) if validation else None,
                    parent_sha256=parent_hash, report=report,
                    ancestry_coverage="complete" if from_scratch else "recorded lineage only; legacy/pretrained corpus may be unknown",
                    seen_peptide_hashes=sorted(ancestry | {peptide_hash(s) for s in tr_seqs + va_seqs}))
    if _target == "im":
        metadata["im_units"] = metadata.pop("rt_units")
        metadata["charge_range"] = [min(r.charge for r in tr), max(r.charge for r in tr)]
    save_model(out, model.cpu(), info, metadata, target=_target)
    return report
