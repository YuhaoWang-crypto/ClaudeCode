#!/usr/bin/env python3
"""
Nucleotide "language model" — a demonstrator generator for de novo aptamer candidates.

An order-k Markov chain over nucleotides IS the simplest sequence language model:
P(x_i | x_{i-k..i-1}). Trained on known binders (optionally only those against
targets similar to your target), it samples NOVEL candidate sequences that share the
compositional/motif statistics of real aptamers, and scores each by log-likelihood.

This is the honest baseline for the "aptamer language model" idea. Scale-up path
(README.md): fine-tune RNA-FM / a small GPT on the full UTexas corpus for a real,
context-aware generator conditioned on the target embedding (cf. AptaBLE).

Usage:
  python generate_lm.py --train seed_dataset.csv --k 3 --n 10 --lmin 26 --lmax 34 --chem DNA
"""
import argparse, csv, random
from collections import defaultdict

def norm(seq, chem):
    s = seq.upper().replace(" ", "")
    return s.replace("U", "T") if chem == "DNA" else s.replace("T", "U")

def train_markov(seqs, k):
    trans = defaultdict(lambda: defaultdict(int))
    starts = defaultdict(int)
    for s in seqs:
        if len(s) <= k:
            continue
        starts[s[:k]] += 1
        for i in range(len(s) - k):
            trans[s[i:i+k]][s[i+k]] += 1
    return trans, starts

def sample(trans, starts, k, lmin, lmax, alphabet, rng):
    if starts:
        seeds = list(starts.keys()); w = list(starts.values())
        seq = rng.choices(seeds, weights=w)[0]
    else:
        seq = "".join(rng.choice(alphabet) for _ in range(k))
    target_len = rng.randint(lmin, lmax)
    while len(seq) < target_len:
        ctx = seq[-k:]
        nxt = trans.get(ctx)
        if not nxt:
            seq += rng.choice(alphabet)          # backoff to uniform
        else:
            ch = list(nxt.keys()); w = list(nxt.values())
            seq += rng.choices(ch, weights=w)[0]
    return seq[:target_len]

def loglik(seq, trans, k):
    import math
    lp = 0.0; n = 0
    for i in range(len(seq) - k):
        ctx, ch = seq[i:i+k], seq[i+k]
        nxt = trans.get(ctx)
        tot = sum(nxt.values()) if nxt else 0
        p = (nxt.get(ch, 0) + 0.5) / (tot + 0.5*4)   # add-0.5 smoothing
        lp += math.log(p); n += 1
    return lp / n if n else -99.0

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train", required=True)
    ap.add_argument("--k", type=int, default=3)
    ap.add_argument("--n", type=int, default=10)
    ap.add_argument("--lmin", type=int, default=26)
    ap.add_argument("--lmax", type=int, default=34)
    ap.add_argument("--chem", choices=["DNA", "RNA"], default="DNA")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    rng = random.Random(args.seed)

    with open(args.train, newline="") as fh:
        rows = list(csv.DictReader(fh))
    seqs = [norm(r["sequence"], args.chem) for r in rows
            if r.get("sequence") and str(r.get("label", "1")) == "1"]
    if len(seqs) < 2:
        print("Not enough training sequences — load the UTexas export for real use.")
    trans, starts = train_markov(seqs, args.k)
    alphabet = "ACGT" if args.chem == "DNA" else "ACGU"

    cands = set()
    tries = 0
    while len(cands) < args.n and tries < args.n * 50:
        s = sample(trans, starts, args.k, args.lmin, args.lmax, alphabet, rng)
        if s not in seqs:            # keep only novel
            cands.add(s)
        tries += 1
    scored = sorted(cands, key=lambda s: loglik(s, trans, args.k), reverse=True)

    print(f"\nGenerated {len(scored)} novel {args.chem} candidates "
          f"(order-{args.k} Markov LM, trained on {len(seqs)} binders)")
    print("=" * 60)
    print(f"{'#':>2}  {'avg logP':>8}  sequence (5'->3')")
    print("-" * 60)
    for i, s in enumerate(scored, 1):
        print(f"{i:>2}  {loglik(s, trans, args.k):>8.3f}  {s}")
    print("=" * 60)
    print("Novel candidates -> fold-filter (ViennaRNA) -> Boltz-2.1 co-fold (Stage 1).")
    print("NOTE: a k-mer Markov LM has no target conditioning; for target-aware generation")
    print("use the RNA-FM / GPT fine-tune described in README.md.")

if __name__ == "__main__":
    main()
