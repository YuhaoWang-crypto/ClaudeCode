"""Gene-pair (ternary) state code. Python standard library only.

Each panel entry is a gene pair (i, j). For a target log2FC profile b the pair
gives one ternary symbol:
    +1 if b_i - b_j >  tau
    -1 if b_i - b_j < -tau
     0 otherwise (dead zone: difference too small to trust)
A k-pair panel turns a profile into a k-symbol code, stored as two bitmasks
(positive, negative) so similarity is a popcount.

Only within-sample comparisons are used, so the code does not change when every
gene is shifted by the same amount (e.g. a global normalization offset) and only
the comparisons near tau react to monotone distortions. It still needs log2FC
against a matched control: the reference data are contrasts, not raw expression.

Classification: under a state with reference pair difference d, measured
difference ~ N(s * d, 2 sigma^2), where s is the unknown response strength.
"none" has s = 0. Each stimulus class takes the best s on a fixed grid.
Missing genes drop their pairs; they are never imputed.
"""
import math

# Response strength is unknown, so each stimulus class may take any strength from
# 0.1x to 3x its reference (log-spaced). A narrow grid lets "weak IFNa" and
# "IFNg" swap at the grid edge when the two references are nearly collinear.
STRENGTH_GRID = tuple(round(0.1 * 30 ** (i / 11), 4) for i in range(12))
_EPS = 1e-9


def _phi(x):
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))


def symbol(diff, tau):
    return 1 if diff > tau else (-1 if diff < -tau else 0)


def symbol_probs(d, sigma, tau):
    """P(symbol = +1, 0, -1) when the true pair difference is d."""
    s2 = sigma * math.sqrt(2)
    p_pos = 1 - _phi((tau - d) / s2)
    p_neg = _phi((-tau - d) / s2)
    p_zero = max(1 - p_pos - p_neg, 0.0)
    return {1: p_pos, 0: p_zero, -1: p_neg}


def encode(panel, target, tau):
    """Return (symbols, pos_mask, neg_mask, missing_mask). Symbol None = pair not measured."""
    symbols, pos, neg, miss = [], 0, 0, 0
    for r, p in enumerate(panel):
        if p["gene_i"] in target and p["gene_j"] in target:
            s = symbol(target[p["gene_i"]] - target[p["gene_j"]], tau)
            pos |= (s == 1) << r
            neg |= (s == -1) << r
        else:
            s = None
            miss |= 1 << r
        symbols.append(s)
    return symbols, pos, neg, miss


def code_string(symbols):
    return "".join({1: "+", 0: "0", -1: "-", None: "."}[s] for s in symbols)


def hamming(code_a, code_b):
    """Distance between two (pos_mask, neg_mask) codes: count of differing symbols."""
    (pa, na), (pb, nb) = code_a, code_b
    return bin((pa ^ pb) | (na ^ nb)).count("1")


def class_loglik(panel, symbols, d_key, sigma, tau, strengths):
    best = -math.inf
    for s in strengths:
        ll = 0.0
        for p, sym in zip(panel, symbols):
            if sym is None:
                continue
            ll += math.log(symbol_probs(s * p[d_key], sigma, tau)[sym] + _EPS)
        best = max(best, ll)
    return best


def score_pairs(panel, target, sigma, tau, classes, strengths=STRENGTH_GRID):
    """classes: subset of ["IFNA", "IFNG"] that are identifiable for this cell type."""
    symbols, pos, neg, miss = encode(panel, target, tau)
    measured = sum(s is not None for s in symbols)
    if measured == 0:
        raise ValueError("None of the panel gene pairs were measured.")
    ll = {"none": class_loglik(panel, symbols, "d_ifna", sigma, tau, (0.0,))}
    for c in classes:
        ll[c] = class_loglik(panel, symbols, {"IFNA": "d_ifna", "IFNG": "d_ifng"}[c], sigma, tau, strengths)
    best_signal = max((ll[c] for c in classes), default=-math.inf)
    return {
        "code": code_string(symbols),
        "pos_mask": pos, "neg_mask": neg, "missing_mask": miss,
        "pairs_measured": measured, "pairs_total": len(panel),
        "loglik": ll,
        # log-likelihood ratio of the best stimulus class over "no response"
        "llr_signal_vs_none": best_signal - ll["none"],
        "best_class": max(ll, key=ll.get),
    }
