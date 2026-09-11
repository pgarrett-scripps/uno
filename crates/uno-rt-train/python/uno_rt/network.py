"""CNN and features shared by training and the PyTorch parity oracle."""
import numpy as np
import torch
from torch import nn

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_IDX = {c: i for i, c in enumerate(AA)}      # 0..19
PAD = 20                                        # embedding pad index (cnn)
MAXLEN = 50
NAA = 20


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
    def __init__(self, emb=32, ch=64, blocks=(1, 2, 4), k=5, head=64, n_global=N_GLOBAL):
        super().__init__()
        self.n_global = n_global
        self.emb = nn.Embedding(NAA + 1, emb, padding_idx=PAD)
        self.stem = nn.Conv1d(emb, ch, k, padding=(k - 1) // 2)
        self.blocks = nn.ModuleList([ResBlock(ch, k, d) for d in blocks])
        self.act = nn.ReLU()
        # head input = pooled (2*ch) ++ explicit global physicochemical features
        self.head = nn.Sequential(
            nn.Linear(2 * ch + n_global, head), nn.ReLU(), nn.Linear(head, 1)
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
        return self.head(torch.cat([mean, mx, glob[:, :self.n_global]], -1)).squeeze(-1)


