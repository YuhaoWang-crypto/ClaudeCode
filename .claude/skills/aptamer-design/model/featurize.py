#!/usr/bin/env python3
"""
Featurize aptamer (nucleic-acid) and protein sequences with pure-Python, dependency-free
descriptors. These are the "poor-man's embeddings" for the baseline; swap in RNA-FM /
ESM-2 embeddings later for the full model (see model/README.md).

Aptamer features: k-mer composition (k=1..3, RNA U folded to T), length (log), GC fraction.
Protein features: amino-acid composition (20) + a few physicochemical fractions.
"""
from itertools import product
import math

NUC = "ACGT"
AA = "ACDEFGHIKLMNPQRSTVWY"

def _norm_nuc(seq):
    return seq.upper().replace("U", "T").replace(" ", "")

def kmer_vector(seq, k):
    seq = _norm_nuc(seq)
    keys = ["".join(p) for p in product(NUC, repeat=k)]
    counts = {kk: 0 for kk in keys}
    n = 0
    for i in range(len(seq) - k + 1):
        km = seq[i:i+k]
        if km in counts:
            counts[km] += 1; n += 1
    return [counts[kk] / n if n else 0.0 for kk in keys]

def aptamer_features(seq):
    s = _norm_nuc(seq)
    gc = (s.count("G") + s.count("C")) / len(s) if s else 0.0
    feats = []
    for k in (1, 2, 3):
        feats += kmer_vector(s, k)
    feats += [math.log1p(len(s)), gc]
    return feats

def protein_features(seq):
    s = (seq or "").upper()
    s = "".join(c for c in s if c in AA)
    L = len(s) or 1
    comp = [s.count(a) / L for a in AA]
    hydro = sum(s.count(a) for a in "AILMFWVY") / L      # hydrophobic
    pos = sum(s.count(a) for a in "KRH") / L             # basic (matters for aptamers!)
    neg = sum(s.count(a) for a in "DE") / L
    return comp + [hydro, pos, neg, math.log1p(L)]

def pair_features(apt_seq, prot_seq):
    """Concatenated features for an (aptamer, protein) pair."""
    return aptamer_features(apt_seq) + protein_features(prot_seq)

if __name__ == "__main__":
    v = aptamer_features("GGGTAGGGTAGGGTAGGGACTGCAAGT")
    print(f"aptamer feature dim = {len(v)} (4+16+64 kmer + len + gc)")
    print(f"protein feature dim = {len(protein_features('VCNGIGIGEF'))}")
