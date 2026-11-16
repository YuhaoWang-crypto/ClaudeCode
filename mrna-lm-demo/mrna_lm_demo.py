"""
mRNA-LM DEMO  (faithful, self-contained, CPU-friendly reimplementation)
=======================================================================

Reproduces the core idea of Sanofi's mRNA-LM (https://github.com/Sanofi-Public/mRNA-LM):

  An mRNA molecule has three functional regions:
        5'UTR ---- CDS (coding sequence) ---- 3'UTR
  mRNA-LM encodes each region with its OWN BERT tower, then fuses the three
  region embeddings and predicts a molecular property (here: mRNA half-life).

What is faithful to the original repo:
  * Tokenization is IDENTICAL to dataload.py:
        - 5'UTR / 3'UTR : single-nucleotide tokens  (k=1, stride 1) over {A,U,G,C,N}
        - CDS          : codon tokens               (k=3, stride 3) -> 125-token vocab
        - all T -> U before tokenizing
  * Architecture mirrors FullModel.py (non-CLIP path):
        3 independent Transformer encoders -> masked mean-pool per region
        -> concat(5utr, cds, 3utr) -> dense+GELU+LayerNorm+dropout -> linear head
        -> MSE regression   (num_labels == 1 branch)
  * Trained & evaluated on the REAL data shipped in the repo:
        data/mrna_half-life.csv   (12,919 human mRNAs, z-scored half-life)
  * Reported with Spearman / Pearson, the same metrics the paper uses.

What is DIFFERENT (so it runs in minutes on a CPU, no 1GB downloads, no GPU):
  * The three towers are SMALL and RANDOMLY INITIALISED and trained from scratch,
    instead of loading the pretrained CodonBERT / 5UTRBERT / 3UTRBERT checkpoints
    and LoRA-finetuning them. This demonstrates the *method*, not the paper's SOTA
    numbers. (See codonbert_real_demo.py for loading the real pretrained weights.)

Run:
    python mrna_lm_demo.py
"""

import os, math, random, argparse
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from scipy.stats import pearsonr, spearmanr

SEED = 0
random.seed(SEED); np.random.seed(SEED); torch.manual_seed(SEED)

# ----------------------------------------------------------------------------
# 1. Tokenization  (copied in spirit from the repo's dataload.py / *Model.py)
# ----------------------------------------------------------------------------
NUC = list("AUGCN")
SPECIAL = ["[PAD]", "[UNK]", "[CLS]", "[SEP]", "[MASK]"]

# 5'/3'UTR vocab: single nucleotides
UTR_VOCAB = {t: i for i, t in enumerate(SPECIAL + NUC)}
# CDS vocab: all 3-mers over {A,U,G,C,N}  -> 5^3 = 125 codons
CODONS = [a + b + c for a in NUC for b in NUC for c in NUC]
CDS_VOCAB = {t: i for i, t in enumerate(SPECIAL + CODONS)}

def mytok(seq, k, s):
    """Exactly the repo's tokenizer: T->U, then sliding k-mers of stride s."""
    seq = str(seq).upper().replace("T", "U")
    return [seq[j:j + k] for j in range(0, len(seq) - k + 1, s)]

def encode(seq, k, s, vocab, max_len):
    toks = mytok(seq, k, s)[:max_len]
    ids = [vocab.get(t, vocab["[UNK]"]) for t in toks]
    mask = [1] * len(ids)
    pad = max_len - len(ids)
    ids += [vocab["[PAD]"]] * pad
    mask += [0] * pad
    return ids, mask

# region-specific max token counts (truncated for CPU speed; attention is O(L^2))
L5, LC, L3 = 100, 200, 200

# ----------------------------------------------------------------------------
# 2. One BERT-style tower  (small Transformer encoder)
# ----------------------------------------------------------------------------
class Tower(nn.Module):
    def __init__(self, vocab_size, max_len, d=96, layers=2, heads=4, ff=192, drop=0.1):
        super().__init__()
        self.emb = nn.Embedding(vocab_size, d, padding_idx=0)
        self.pos = nn.Embedding(max_len, d)
        enc = nn.TransformerEncoderLayer(d, heads, ff, drop, batch_first=True,
                                         activation="gelu", norm_first=True)
        self.encoder = nn.TransformerEncoder(enc, layers)
        self.d = d

    def forward(self, ids, mask):
        pos = torch.arange(ids.size(1), device=ids.device).unsqueeze(0)
        x = self.emb(ids) + self.pos(pos)
        x = self.encoder(x, src_key_padding_mask=(mask == 0))
        # masked mean pool
        m = mask.unsqueeze(-1).float()
        return (x * m).sum(1) / m.sum(1).clamp(min=1e-9)

# ----------------------------------------------------------------------------
# 3. mRNA-LM: three towers -> fuse -> regression head   (cf. FullModel.py)
# ----------------------------------------------------------------------------
class MRNALM(nn.Module):
    def __init__(self, d=96, head_dim=192, drop=0.3):
        super().__init__()
        self.t5 = Tower(len(UTR_VOCAB), L5, d)
        self.tc = Tower(len(CDS_VOCAB), LC, d)
        self.t3 = Tower(len(UTR_VOCAB), L3, d)
        self.fuse = nn.Linear(d * 3, head_dim)
        self.act = nn.GELU()
        self.norm = nn.LayerNorm(head_dim)
        self.drop = nn.Dropout(drop)
        self.head = nn.Linear(head_dim, 1)

    def forward(self, b):
        e5 = self.t5(b["ids5"], b["m5"])
        ec = self.tc(b["idsc"], b["mc"])
        e3 = self.t3(b["ids3"], b["m3"])
        h = torch.cat([e5, ec, e3], dim=1)
        h = self.drop(self.norm(self.act(self.fuse(h))))
        return self.head(h).squeeze(-1)

# ----------------------------------------------------------------------------
# 4. Data: real half-life CSV, fold-based split
# ----------------------------------------------------------------------------
def build(df, max_n=None):
    if max_n and len(df) > max_n:
        df = df.sample(max_n, random_state=SEED)
    ids5, m5, idsc, mc, ids3, m3, ys = [], [], [], [], [], [], []
    for _, r in df.iterrows():
        a, am = encode(r.UTR5, 1, 1, UTR_VOCAB, L5)
        c, cm = encode(r.CDS,  3, 3, CDS_VOCAB, LC)
        t, tm = encode(r.UTR3, 1, 1, UTR_VOCAB, L3)
        ids5.append(a); m5.append(am); idsc.append(c); mc.append(cm)
        ids3.append(t); m3.append(tm); ys.append(float(r.y))
    T = lambda x: torch.tensor(x, dtype=torch.long)
    return {"ids5": T(ids5), "m5": T(m5), "idsc": T(idsc), "mc": T(mc),
            "ids3": T(ids3), "m3": T(m3),
            "y": torch.tensor(ys, dtype=torch.float)}

def batches(data, bs, shuffle):
    n = len(data["y"]); idx = list(range(n))
    if shuffle: random.shuffle(idx)
    for i in range(0, n, bs):
        j = idx[i:i + bs]
        yield {k: v[j] for k, v in data.items()}

def evaluate(model, data, bs=64):
    model.eval(); preds = []
    with torch.no_grad():
        for b in batches(data, bs, False):
            preds.append(model(b).numpy())
    p = np.concatenate(preds); y = data["y"].numpy()
    return spearmanr(y, p)[0], pearsonr(y, p)[0], p

# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="data/mrna_half-life.csv")
    ap.add_argument("--epochs", type=int, default=4)
    ap.add_argument("--train_n", type=int, default=4000, help="subsample train for speed")
    ap.add_argument("--bs", type=int, default=32)
    ap.add_argument("--lr", type=float, default=3e-4)
    args = ap.parse_args()

    torch.set_num_threads(os.cpu_count() or 4)
    df = pd.read_csv(args.csv).dropna(subset=["y"])
    # folds 0-7 train, 8 valid, 9 test  (the CSV ships 10 CV folds in `split`)
    tr = df[df.split <= 7]; va = df[df.split == 8]; te = df[df.split == 9]
    print(f"data: {len(df)} mRNAs | train {len(tr)} valid {len(va)} test {len(te)}")
    print("encoding (truncate 5utr<=%d nt, cds<=%d codons, 3utr<=%d nt) ..." % (L5, LC, L3))
    dtr = build(tr, args.train_n); dva = build(va, 1000); dte = build(te, 1000)

    model = MRNALM()
    nparam = sum(p.numel() for p in model.parameters())
    print(f"model: tri-tower mRNA-LM, {nparam/1e6:.2f}M params\n")

    opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=1e-2)
    lossf = nn.MSELoss()

    # baseline before training
    s, p, _ = evaluate(model, dte)
    print(f"[init ] test Spearman {s:+.3f}  Pearson {p:+.3f}  (random model)")

    for ep in range(1, args.epochs + 1):
        model.train(); tot = 0; nb = 0
        for b in batches(dtr, args.bs, True):
            opt.zero_grad()
            loss = lossf(model(b), b["y"])
            loss.backward(); opt.step()
            tot += loss.item(); nb += 1
        sv, pv, _ = evaluate(model, dva)
        print(f"[ep {ep}] train MSE {tot/nb:.3f} | valid Spearman {sv:+.3f} Pearson {pv:+.3f}")

    s, p, preds = evaluate(model, dte)
    print(f"\n==== TEST (held-out fold 9) ====")
    print(f"Spearman {s:+.3f}   Pearson {p:+.3f}   (n={len(preds)})")
    print("Random baseline is ~0.0; a positive correlation means the tri-tower")
    print("model learned to read half-life signal from 5'UTR + CDS + 3'UTR.\n")

    # show a few predictions
    y = dte["y"].numpy()
    order = np.argsort(-preds)
    print("example predictions (z-scored half-life):")
    print("  rank   predicted    actual")
    for r in list(order[:3]) + list(order[-3:]):
        print(f"   {'' :>2}      {preds[r]:+.2f}      {y[r]:+.2f}")

if __name__ == "__main__":
    main()
