#!/usr/bin/env python3
"""
Baseline aptamer-protein interaction scorer (Stage 0 pre-filter before Boltz).

Two modes, auto-selected by what's installed:
  * If scikit-learn is available: train LogisticRegression on known binders (label=1) +
    target-shuffled synthetic negatives (label=0), report P(bind) for each query.
  * Else (pure-Python fallback): k-NN cosine similarity in pair-feature space — score a
    query by its similarity to the most similar known binder, up-weighted when the known
    binder's target resembles the query target.

This is a transfer-friendly BASELINE. For the real model, replace featurize.* with
RNA-FM (aptamer) + ESM-2 (protein) embeddings and swap LogisticRegression for the
two-tower cross-attention head in README.md — the I/O here stays the same.

Usage:
  python predict_baseline.py --train seed_dataset.csv --query queries.csv
"""
import argparse, csv, math, sys
from featurize import pair_features, aptamer_features, protein_features, _norm_nuc

def read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))

def cosine(a, b):
    dot = sum(x*y for x, y in zip(a, b))
    na = math.sqrt(sum(x*x for x in a)); nb = math.sqrt(sum(y*y for y in b))
    return dot / (na*nb) if na and nb else 0.0

def target_sim(t1, t2):
    """Cheap protein similarity via AA-composition cosine (proxy for family relatedness)."""
    return cosine(protein_features(t1), protein_features(t2))

# ---------- pure-Python kNN fallback ----------
def knn_score(train, apt, tgt, k=3):
    q = pair_features(apt, tgt)
    sims = []
    for r in train:
        if r.get("label", "1") not in ("1", 1, "", None):
            continue
        s = cosine(q, pair_features(r["sequence"], r.get("target_seq", "")))
        ts = target_sim(tgt, r.get("target_seq", "")) if r.get("target_seq") else 0.5
        sims.append(0.5*s + 0.5*ts)   # blend pair-similarity and target-relatedness
    sims.sort(reverse=True)
    return sum(sims[:k]) / max(1, min(k, len(sims)))

# ---------- optional logistic-regression path ----------
def try_sklearn(train, queries):
    try:
        import numpy as np
        from sklearn.linear_model import LogisticRegression
    except Exception:
        return None
    import random
    random.seed(7)
    pos = [r for r in train if str(r.get("label", "1")) == "1"]
    if len(pos) < 4:
        return None  # too few to train; use fallback
    tgts = [r.get("target_seq", "") for r in pos]
    X, y = [], []
    for r in pos:
        X.append(pair_features(r["sequence"], r.get("target_seq", ""))); y.append(1)
        # synthetic negative: same aptamer, a different (shuffled) target
        alt = random.choice([t for t in tgts if t != r.get("target_seq", "")] or tgts)
        X.append(pair_features(r["sequence"], alt)); y.append(0)
    clf = LogisticRegression(max_iter=1000).fit(np.array(X), np.array(y))
    out = []
    for q in queries:
        p = clf.predict_proba([pair_features(q["sequence"], q.get("target_seq", ""))])[0][1]
        out.append((q, float(p)))
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train", required=True)
    ap.add_argument("--query", required=True)
    args = ap.parse_args()
    train = read_csv(args.train); queries = read_csv(args.query)

    sk = try_sklearn(train, queries)
    mode = "logreg (sklearn)" if sk else "kNN cosine (pure-Python fallback)"
    if sk:
        scored = sorted(sk, key=lambda t: t[1], reverse=True)
        rows = [(q["apt_id"], q.get("target_name", ""), q["chem"], p) for q, p in scored]
        col = "P(bind)"
    else:
        tmp = [(q, knn_score(train, q["sequence"], q.get("target_seq", "")))
               for q in queries]
        tmp.sort(key=lambda t: t[1], reverse=True)
        rows = [(q["apt_id"], q.get("target_name", ""), q["chem"], s) for q, s in tmp]
        col = "sim-score"

    print(f"\nStage-0 aptamer prioritization  [mode: {mode}]")
    print("=" * 62)
    print(f"{'#':>2}  {'apt_id':<16} {'target':<10} {'chem':<5} {col:>9}")
    print("-" * 62)
    for i, (aid, tgt, chem, sc) in enumerate(rows, 1):
        print(f"{i:>2}  {aid:<16} {tgt:<10} {chem:<5} {sc:>9.3f}")
    print("=" * 62)
    print("Baseline pre-filter only. Pass survivors to Boltz-2.1 (Stage 1) for structural"
          " scoring. Scores are relative, NOT affinities.")

if __name__ == "__main__":
    main()
