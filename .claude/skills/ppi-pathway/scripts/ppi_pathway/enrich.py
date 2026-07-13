"""
enrich.py — offline over-representation (hypergeometric) pathway enrichment.

The STRING API (``string_api.enrichment``) is the easiest route and needs no
gene-set files. Use *this* module when you want full control / reproducibility:
supply your own gene sets as a GMT file (KEGG, Reactome, MSigDB Hallmark, GO —
all distributed in ``.gmt`` format) and run the classic hypergeometric test with
Benjamini-Hochberg FDR — no network calls, fully deterministic.

GMT format: one gene set per line, tab-separated:
    <set_name>\t<description>\t<gene1>\t<gene2>\t...
"""
from __future__ import annotations

from math import lgamma, log, exp
from pathlib import Path
from typing import Iterable


def read_gmt(path: str | Path) -> dict[str, set[str]]:
    sets: dict[str, set[str]] = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                sets[parts[0]] = {g.strip().upper() for g in parts[2:] if g.strip()}
    return sets


def _log_binom(n: int, k: int) -> float:
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)


def _hypergeom_sf(k: int, M: int, n: int, N: int) -> float:
    """P(X >= k) for the hypergeometric distribution.

    M = population size, n = successes in population (gene-set size within
    background), N = sample size (query size within background), k = observed
    overlap. Pure-python; uses scipy if available for speed/accuracy.
    """
    try:
        from scipy.stats import hypergeom
        return float(hypergeom.sf(k - 1, M, n, N))
    except Exception:
        pass
    # log-space fallback
    lo, hi = k, min(n, N)
    denom = _log_binom(M, N)
    total = 0.0
    for i in range(lo, hi + 1):
        lp = _log_binom(n, i) + _log_binom(M - n, N - i) - denom
        total += exp(lp)
    return min(1.0, total)


def _bh_fdr(pvals: list[float]) -> list[float]:
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    fdr = [0.0] * m
    prev = 1.0
    for rank, idx in enumerate(reversed(order), start=1):
        i = m - rank  # position from the top
        q = pvals[idx] * m / (m - i)
        prev = min(prev, q)
        fdr[idx] = min(1.0, prev)
    return fdr


def enrich(query: Iterable[str], gene_sets: dict[str, set[str]],
           background: Iterable[str] | None = None,
           min_set: int = 3, max_set: int = 500,
           fdr_max: float = 0.05) -> list[dict]:
    """Hypergeometric over-representation of ``query`` against ``gene_sets``.

    ``background``: the universe of testable genes (e.g. all genes measured in
    your assay). Defaults to the union of all gene-set members — pass your real
    background for a statistically honest test.
    """
    q = {g.strip().upper() for g in query if str(g).strip()}
    if background is None:
        universe = set().union(*gene_sets.values()) if gene_sets else set()
    else:
        universe = {g.strip().upper() for g in background}
    q &= universe
    M = len(universe)
    N = len(q)
    if M == 0 or N == 0:
        return []

    rows = []
    for name, members in gene_sets.items():
        s = members & universe
        if not (min_set <= len(s) <= max_set):
            continue
        overlap = q & s
        k = len(overlap)
        if k == 0:
            continue
        p = _hypergeom_sf(k, M, len(s), N)
        rows.append({
            "term": name,
            "overlap": k,
            "set_size": len(s),
            "query_size": N,
            "p_value": p,
            "genes": sorted(overlap),
        })
    if not rows:
        return []
    fdr = _bh_fdr([r["p_value"] for r in rows])
    for r, q_ in zip(rows, fdr):
        r["fdr"] = q_
    rows = [r for r in rows if r["fdr"] <= fdr_max]
    rows.sort(key=lambda r: r["p_value"])
    return rows
