"""CodonBERT demo: extract mRNA embeddings and benchmark on mRFP expression.

Uses the pretrained Sanofi CodonBERT (BertForPreTraining) to turn codon
sequences into embeddings, then trains a ridge regressor to predict protein
expression, comparing CLS-pooled vs mean-pooled embeddings against a
codon-frequency baseline. This is a lightweight *frozen-feature probe* (no
fine-tuning) so it runs on CPU in a few minutes; the paper's fine-tuned numbers
are higher.

Paths are taken from environment variables so the script is portable:
  CODONBERT_SRC    = <clone>/benchmarks        (holds utils/ and data/)
  CODONBERT_MODEL  = <dir with config.json + pytorch_model.bin>

Example:
  export CODONBERT_SRC=$PWD/CodonBERT/benchmarks
  export CODONBERT_MODEL=$PWD/cb_model/codonbert
  python -I codonbert_mrfp_demo.py
"""
import os, sys, time
import numpy as np, pandas as pd
import torch
from transformers import BertForPreTraining
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score

SRC = os.environ["CODONBERT_SRC"]            # .../CodonBERT/benchmarks
MODEL_DIR = os.environ["CODONBERT_MODEL"]    # .../codonbert
DATA = os.path.join(SRC, "CodonBERT/data/fine-tune/mRFP_Expression.csv")
sys.path.append(SRC)
from utils.tokenizer import get_tokenizer, mytok

tok = get_tokenizer()
model = BertForPreTraining.from_pretrained(MODEL_DIR)
model.eval()
print(f"Loaded CodonBERT: {sum(p.numel() for p in model.parameters())/1e6:.1f}M params")

df = pd.read_csv(DATA)
print(f"mRFP dataset: {df.shape[0]} seqs, columns={list(df.columns)}")
train_df = df[df.Split == "train"].reset_index(drop=True)
test_df  = df[df.Split == "test"].reset_index(drop=True)
if len(test_df) == 0:  # fallback if no official split
    train_df, test_df = df.iloc[:int(.8*len(df))], df.iloc[int(.8*len(df)):]
print(f"train={len(train_df)}  test={len(test_df)}")

@torch.no_grad()
def embed(seqs, max_len=1024, bs=16):
    cls_all, mean_all = [], []
    for i in range(0, len(seqs), bs):
        batch = seqs[i:i+bs]
        toks = [tok.encode(" ".join(mytok(s, 3, 3))).ids[:max_len] for s in batch]
        L = max(len(t) for t in toks)
        ids = torch.zeros(len(toks), L, dtype=torch.long)
        mask = torch.zeros(len(toks), L, dtype=torch.long)
        for j, t in enumerate(toks):
            ids[j, :len(t)] = torch.tensor(t); mask[j, :len(t)] = 1
        hs = model(ids, attention_mask=mask, output_hidden_states=True).hidden_states[-1]
        cls_all.append(hs[:, 0].numpy())
        m = (hs * mask.unsqueeze(-1)).sum(1) / mask.sum(1, keepdim=True)
        mean_all.append(m.numpy())
        print(f"  embedded {min(i+bs,len(seqs))}/{len(seqs)}", end="\r")
    print()
    return np.vstack(cls_all), np.vstack(mean_all)

def codon_freq(seqs):  # 64-dim baseline
    vocab = [a+b+c for a in "AUGC" for b in "AUGC" for c in "AUGC"]
    idx = {k: i for i, k in enumerate(vocab)}
    X = np.zeros((len(seqs), 64))
    for r, s in enumerate(seqs):
        km = mytok(s, 3, 3)
        for k in km:
            if k in idx: X[r, idx[k]] += 1
        if km: X[r] /= len(km)
    return X

yr, yt = train_df.Value.values, test_df.Value.values
t0 = time.time()
Xr_cls, Xr_mean = embed(train_df.Sequence.tolist())
Xt_cls, Xt_mean = embed(test_df.Sequence.tolist())
print(f"embedding time: {time.time()-t0:.0f}s")

def evaluate(name, Xr, Xt):
    m = Ridge(alpha=1.0).fit(Xr, yr)
    p = m.predict(Xt)
    print(f"{name:28s}  Pearson={pearsonr(p,yt)[0]:.3f}  Spearman={spearmanr(p,yt)[0]:.3f}  R2={r2_score(yt,p):.3f}")

print("\n=== mRFP expression prediction (test set) ===")
evaluate("codon-freq baseline (64d)", codon_freq(train_df.Sequence.tolist()), codon_freq(test_df.Sequence.tolist()))
evaluate("CodonBERT CLS (768d)", Xr_cls, Xt_cls)
evaluate("CodonBERT mean-pool (768d)", Xr_mean, Xt_mean)
