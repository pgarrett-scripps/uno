#!/usr/bin/env python3
"""Train an RT model for uno on real experimental iRT (ProteomeTools/Prosit).

Two architectures:
  posmlp : positional one-hot (MAXLEN*20) -> ReLU MLP. Exports to the existing
           `UMLP` binary (feature_kind=1) -> drop-in, no Rust changes.
  cnn    : learned AA embedding -> residual dilated Conv1d stack -> masked
           mean+max pool -> MLP head. Exports to a new `URTC` binary (needs a
           Rust inference path). Chronologer-style; the architecture the
           literature says beats DeepLC.

Labels are iRT (run-independent); uno applies per-run calibration downstream, so
the absolute scale only needs to be monotone+linear with observed RT. AA-only:
C is implicitly carbamidomethyl (uno's fixed mod), no variable mods.
"""
import argparse, struct, sys, time
import numpy as np
import torch
import torch.nn as nn

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_IDX = {c: i for i, c in enumerate(AA)}      # 0..19
PAD = 20                                        # embedding pad index (cnn)
MAXLEN = 50
NAA = 20


def load_tsv(path):
    seqs, ys = [], []
    with open(path) as f:
        for line in f:
            s, y = line.rstrip().split("\t")
            seqs.append(s)
            ys.append(float(y))
    return seqs, np.asarray(ys, dtype=np.float32)


# ---------- featurizers ----------
def feat_posmlp(seqs):
    """Dense positional one-hot, (N, MAXLEN*NAA) float32 (matches Rust)."""
    X = np.zeros((len(seqs), MAXLEN * NAA), dtype=np.float32)
    for i, s in enumerate(seqs):
        for pos, c in enumerate(s[:MAXLEN]):
            a = AA_IDX.get(c)
            if a is not None:
                X[i, pos * NAA + a] = 1.0
    return X


def feat_tokens(seqs):
    """Integer token ids (N, MAXLEN) with PAD; plus lengths."""
    T = np.full((len(seqs), MAXLEN), PAD, dtype=np.int64)
    L = np.zeros(len(seqs), dtype=np.int64)
    for i, s in enumerate(seqs):
        n = min(len(s), MAXLEN)
        L[i] = n
        for pos, c in enumerate(s[:MAXLEN]):
            T[i, pos] = AA_IDX.get(c, PAD)
    return T, L


# Explicit global physicochemical features appended to the pooled conv vector
# (DeepLC-style). Order + formulas MUST match `global_feats` / `N_GLOBAL` in
# crates/uno-core/src/rt_cnn.rs EXACTLY. Computed from the FULL sequence.
N_GLOBAL = 11
_ACIDIC, _BASIC, _PHOB, _AROM = set("DE"), set("KR"), set("FWYLIVM"), set("FWY")
# Kyte-Doolittle hydropathy; keep identical to kd() in rt_cnn.rs.
_KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
       "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
       "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}

def feat_global(seqs):
    """(N, N_GLOBAL) float32, order matching global_feats() in rt_cnn.rs (v3):
    [len/50, (D+E)/len, (K+R)/len, H/len, hydrophobic/len, KDsum/50, GRAVY,
     aromatic/len, Pro/len, Nterm_KD/4.5, Cterm_KD/4.5]."""
    G = np.zeros((len(seqs), N_GLOBAL), dtype=np.float32)
    for i, s in enumerate(seqs):
        n = max(len(s), 1)
        de = sum(c in _ACIDIC for c in s)
        kr = sum(c in _BASIC for c in s)
        h = s.count("H")
        phob = sum(c in _PHOB for c in s)
        arom = sum(c in _AROM for c in s)
        pro = s.count("P")
        kd_sum = sum(_KD.get(c, 0.0) for c in s)
        nterm = _KD.get(s[0], 0.0) if s else 0.0
        cterm = _KD.get(s[-1], 0.0) if s else 0.0
        G[i] = (len(s) / 50.0, de / n, kr / n, h / n, phob / n,
                kd_sum / 50.0, kd_sum / n, arom / n, pro / n,
                nterm / 4.5, cterm / 4.5)
    return G


# ---------- models ----------
class PosMLP(nn.Module):
    def __init__(self, hidden=(512, 128)):
        super().__init__()
        dims = [MAXLEN * NAA, *hidden]
        layers = []
        for a, b in zip(dims[:-1], dims[1:]):
            layers += [nn.Linear(a, b), nn.ReLU()]
        layers += [nn.Linear(dims[-1], 1)]
        self.net = nn.Sequential(*layers)

    def forward(self, x):
        return self.net(x).squeeze(-1)


class ResBlock(nn.Module):
    def __init__(self, ch, k, dil):
        super().__init__()
        pad = (k - 1) // 2 * dil
        self.c1 = nn.Conv1d(ch, ch, k, padding=pad, dilation=dil)
        self.c2 = nn.Conv1d(ch, ch, k, padding=pad, dilation=dil)
        self.act = nn.ReLU()

    def forward(self, x, m):
        # mask after every conv so PAD positions stay exactly zero -> inference
        # can run over the real peptide length (no fixed-MAXLEN padding work).
        h = self.act(self.c1(x)) * m
        h = self.c2(h) * m
        return self.act(x + h) * m


class RTCnn(nn.Module):
    """Embedding -> conv stem -> residual dilated blocks -> masked mean+max -> head."""
    def __init__(self, emb=32, ch=64, blocks=(1, 2, 4), k=5, head=64):
        super().__init__()
        self.emb = nn.Embedding(NAA + 1, emb, padding_idx=PAD)
        self.stem = nn.Conv1d(emb, ch, k, padding=(k - 1) // 2)
        self.blocks = nn.ModuleList([ResBlock(ch, k, d) for d in blocks])
        self.act = nn.ReLU()
        # head input = pooled (2*ch) ++ explicit global physicochemical features
        self.head = nn.Sequential(
            nn.Linear(2 * ch + N_GLOBAL, head), nn.ReLU(), nn.Linear(head, 1)
        )

    def forward(self, tok, length, glob):
        m = (tok != PAD).float().unsqueeze(1)          # (N,1,L) mask
        x = self.emb(tok).transpose(1, 2)              # (N,emb,L)
        x = self.act(self.stem(x)) * m
        for b in self.blocks:
            x = b(x, m)
        # masked mean + max pool over length
        s = x.sum(-1)
        mean = s / m.sum(-1).clamp(min=1)
        x_masked = x.masked_fill(m == 0, float("-inf"))
        mx = x_masked.max(-1).values
        mx = torch.nan_to_num(mx, neginf=0.0)
        return self.head(torch.cat([mean, mx, glob], -1)).squeeze(-1)


# ---------- metrics ----------
def metrics(y, p):
    e = y - p
    ss_res = float(np.sum(e * e))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot
    mae = float(np.mean(np.abs(e)))
    rmse = (ss_res / len(y)) ** 0.5
    dt95 = float(np.percentile(np.abs(e), 95))
    return r2, mae, rmse, dt95


# ---------- export ----------
def export_umlp(path, model, y_mean, y_std):
    """Write PosMLP to the existing UMLP v2 binary (feature_kind=1, unstandardized
    input so the Rust sparse positional path activates)."""
    lin = [m for m in model.net if isinstance(m, nn.Linear)]
    in_dim = MAXLEN * NAA
    with open(path, "wb") as f:
        f.write(b"UMLP")
        f.write(struct.pack("<I", 2))            # version
        f.write(struct.pack("<I", 1))            # feature_kind = positional
        f.write(struct.pack("<f", y_mean))
        f.write(struct.pack("<f", y_std))
        f.write(struct.pack("<I", in_dim))
        f.write(struct.pack(f"<{in_dim}f", *([0.0] * in_dim)))   # x_mean=0
        f.write(struct.pack(f"<{in_dim}f", *([1.0] * in_dim)))   # x_std=1
        f.write(struct.pack("<I", len(lin)))
        for i, L in enumerate(lin):
            relu = 1 if i < len(lin) - 1 else 0
            w = L.weight.detach().numpy().astype("<f4")   # (out,in) row-major
            b = L.bias.detach().numpy().astype("<f4")
            f.write(struct.pack("<I", w.shape[0]))
            f.write(struct.pack("<I", w.shape[1]))
            f.write(struct.pack("<I", relu))
            f.write(w.tobytes())
            f.write(b.tobytes())
    print(f"exported UMLP -> {path}")


def export_urtc(path, model, y_mean, y_std, cfg):
    """Write RTCnn to a new URTC binary. Layout documented in rt_cnn.rs."""
    def wf(f, arr):
        f.write(np.asarray(arr, dtype="<f4").tobytes())
    with open(path, "wb") as f:
        f.write(b"URTC")
        f.write(struct.pack("<I", 2))                    # version (2 = global feats)
        f.write(struct.pack("<f", y_mean))
        f.write(struct.pack("<f", y_std))
        f.write(struct.pack("<I", cfg["emb"]))
        f.write(struct.pack("<I", cfg["ch"]))
        f.write(struct.pack("<I", cfg["k"]))
        f.write(struct.pack("<I", len(cfg["blocks"])))
        f.write(struct.pack(f"<{len(cfg['blocks'])}I", *cfg["blocks"]))
        f.write(struct.pack("<I", cfg["head"]))
        f.write(struct.pack("<I", N_GLOBAL))             # v2: count of appended global feats
        # embedding (NAA+1, emb)
        wf(f, model.emb.weight.detach().numpy())
        # stem conv (ch, emb, k) + bias
        wf(f, model.stem.weight.detach().numpy()); wf(f, model.stem.bias.detach().numpy())
        # blocks: each c1,c2 (ch,ch,k)+bias
        for b in model.blocks:
            wf(f, b.c1.weight.detach().numpy()); wf(f, b.c1.bias.detach().numpy())
            wf(f, b.c2.weight.detach().numpy()); wf(f, b.c2.bias.detach().numpy())
        # head: Linear(2ch,head),Linear(head,1)
        h = list(model.head)
        wf(f, h[0].weight.detach().numpy()); wf(f, h[0].bias.detach().numpy())
        wf(f, h[2].weight.detach().numpy()); wf(f, h[2].bias.detach().numpy())
    print(f"exported URTC -> {path}")


def load_urtc_weights(model, path):
    """Legacy adapter; use uno_rt.train to infer architecture and RT normalization."""
    from uno_rt.artifact import load_model
    loaded, info, _ = load_model(path)
    model.load_state_dict(loaded.state_dict())
    return info


# ---------- train ----------
def main():
    print("DEPRECATED research trainer: use `uno rt train` for validated CNN training and fine-tuning.", file=sys.stderr)
    ap = argparse.ArgumentParser()
    ap.add_argument("--arch", choices=["posmlp", "cnn"], default="cnn")
    ap.add_argument("--train", required=True)
    ap.add_argument("--val", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--epochs", type=int, default=40)
    ap.add_argument("--bs", type=int, default=2048)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--max-train", type=int, default=0)
    ap.add_argument("--patience", type=int, default=6)
    ap.add_argument("--ch", type=int, default=64)
    ap.add_argument("--emb", type=int, default=32)
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--blocks", type=str, default="1,2,4")
    ap.add_argument("--head", type=int, default=64)
    ap.add_argument("--init", type=str, default="", help="warm-start from a URTC file (full-network fine-tune)")
    args = ap.parse_args()
    blocks = [int(x) for x in args.blocks.split(",")]

    torch.manual_seed(0)
    t0 = time.time()
    tr_seq, tr_y = load_tsv(args.train)
    va_seq, va_y = load_tsv(args.val)
    if args.max_train and len(tr_seq) > args.max_train:
        tr_seq, tr_y = tr_seq[: args.max_train], tr_y[: args.max_train]
    print(f"loaded train={len(tr_seq)} val={len(va_seq)} in {time.time()-t0:.1f}s")

    y_mean, y_std = float(tr_y.mean()), float(tr_y.std())
    tr_ys = (tr_y - y_mean) / y_std
    va_ys = (va_y - y_mean) / y_std

    dev = "cuda" if torch.cuda.is_available() else "cpu"
    if args.arch == "posmlp":
        tr_inputs = (torch.from_numpy(feat_posmlp(tr_seq)),)
        va_inputs = (torch.from_numpy(feat_posmlp(va_seq)),)
        model = PosMLP()
    else:
        Ttr, Ltr = feat_tokens(tr_seq)
        Tva, Lva = feat_tokens(va_seq)
        tr_inputs = (torch.from_numpy(Ttr), torch.from_numpy(Ltr), torch.from_numpy(feat_global(tr_seq)))
        va_inputs = (torch.from_numpy(Tva), torch.from_numpy(Lva), torch.from_numpy(feat_global(va_seq)))
        model = RTCnn(emb=args.emb, ch=args.ch, blocks=blocks, k=args.k, head=args.head)
        if args.init:
            info = load_urtc_weights(model, args.init)
            # Preserve physical predictions when changing target normalization.
            with torch.no_grad():
                model.head[-1].weight.mul_(info["y_std"] / y_std)
                model.head[-1].bias.mul_(info["y_std"] / y_std).add_((info["y_mean"] - y_mean) / y_std)
            print(f"warm-started from {args.init}")
    model = model.to(dev)
    tr_inputs = tuple(t.to(dev) for t in tr_inputs)
    va_inputs = tuple(t.to(dev) for t in va_inputs)
    print(f"device={dev}")
    n_params = sum(p.numel() for p in model.parameters())
    print(f"arch={args.arch} params={n_params:,}  featurized in {time.time()-t0:.1f}s")

    ytr = torch.from_numpy(tr_ys).to(dev)
    opt = torch.optim.Adam(model.parameters(), lr=args.lr)
    sched = torch.optim.lr_scheduler.ReduceLROnPlateau(opt, factor=0.5, patience=2)
    lossf = nn.SmoothL1Loss()
    n = len(tr_seq)
    best = (1e9, None)
    bad = 0
    for ep in range(args.epochs):
        model.train()
        perm = torch.randperm(n)
        te = time.time()
        tot = 0.0
        for i in range(0, n, args.bs):
            idx = perm[i : i + args.bs]
            batch = [t[idx] for t in tr_inputs]
            opt.zero_grad()
            out = model(*batch)
            loss = lossf(out, ytr[idx])
            loss.backward()
            opt.step()
            tot += loss.item() * len(idx)
        model.eval()
        with torch.no_grad():
            preds = []
            for i in range(0, len(va_seq), 8192):
                batch = [t[i : i + 8192] for t in va_inputs]
                preds.append(model(*batch).cpu().numpy())
        pv = np.concatenate(preds) * y_std + y_mean
        r2, mae, rmse, dt95 = metrics(va_y, pv)
        sched.step(mae)
        print(f"ep{ep:02d} loss={tot/n:.4f} | val R2={r2:.4f} MAE={mae:.3f} "
              f"RMSE={rmse:.3f} dt95={dt95:.3f} | {time.time()-te:.1f}s")
        if mae < best[0]:
            best = (mae, {k: v.clone() for k, v in model.state_dict().items()})
            bad = 0
        else:
            bad += 1
            if bad >= args.patience:
                print(f"early stop (no val MAE improvement in {args.patience})")
                break
    model.load_state_dict(best[1])
    model = model.to("cpu")
    print(f"best val MAE={best[0]:.3f} (iRT units)")

    if args.arch == "posmlp":
        export_umlp(args.out, model, y_mean, y_std)
    else:
        cfg = dict(emb=args.emb, ch=args.ch, k=args.k, blocks=blocks, head=args.head)
        export_urtc(args.out, model, y_mean, y_std, cfg)


if __name__ == "__main__":
    main()
