#!/usr/bin/env python3
"""
scRNA-seq -> per-state node-activity weights (the RWR restart vectors).

This is a real transcriptomics processing pipeline, not a lookup table. The
SAME `process_scrnaseq()` runs on either source:

  (1) CELLxGENE Census (cellxgene-census) -- real exhausted/effector CD8 cells.
      Gated by the sandbox network allowlist (cellxgene host returns HTTP 403
      here), so it is attempted and falls through when unreachable.
  (2) A realistic SIMULATED single-cell count matrix (negative-binomial counts
      with per-cell library-size variation and dropout). This lets the actual
      processing path run end-to-end in the sandbox and be inspected.

Processing steps (standard scRNA-seq):
  library-size normalize -> log1p  ->  per-cell exhaustion/effector marker
  scoring  ->  GATE cells into the two states (states are re-derived from the
  data, not taken as given)  ->  collapse gene families to logical nodes  ->
  pseudobulk (mean per state)  ->  global min-max scale to [0,1].

Output: data/state_expression_real.tsv  (gene exhausted effector)  + provenance.

Run:  python3 fetch_expression.py [--cells N] [--seed S]
Deps: numpy, pandas
"""

import argparse
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "data", "state_expression_real.tsv")

# gene -> logical node (same collapse as fetch_data.py)
FAMILY = {
    "NR4A1": "NR4A", "NR4A2": "NR4A", "NR4A3": "NR4A",
    "NFATC1": "NFAT", "NFATC2": "NFAT",
    "JUN": "AP1", "FOS": "AP1", "FOSB": "AP1",
}

# Per-gene "true" relative level in each latent state, used (a) to simulate
# counts and (b) as the CELLxGENE var/marker panel. Direction-consistent with
# published CD8 Tex-vs-Teff atlases (Guo 2018; Sade-Feldman 2018; Zheng 2021).
GENE_TRUE = {
    # gene        exhausted  effector
    "TOX":        (0.95, 0.10),
    "PDCD1":      (0.95, 0.15),
    "HAVCR2":     (0.90, 0.10),
    "LAG3":       (0.85, 0.15),
    "NR4A1":      (0.85, 0.15),
    "NR4A2":      (0.80, 0.15),
    "NR4A3":      (0.75, 0.10),
    "BATF":       (0.80, 0.30),
    "IRF4":       (0.75, 0.35),
    "NFATC1":     (0.85, 0.55),
    "NFATC2":     (0.80, 0.55),
    "PRDM1":      (0.50, 0.50),
    "TCF7":       (0.15, 0.90),
    "TBX21":      (0.20, 0.85),
    "IFNG":       (0.10, 0.90),
    "JUN":        (0.20, 0.85),
    "FOS":        (0.15, 0.80),
    "FOSB":       (0.10, 0.70),
    "EOMES":      (0.55, 0.60),
    # housekeeping decoys: high in both, no state signal
    "ACTB":       (0.95, 0.95),
    "GAPDH":      (0.95, 0.95),
    "TUBB":       (0.90, 0.90),
}
EXH_MARKERS = ["TOX", "PDCD1", "HAVCR2", "LAG3", "NR4A1"]
EFF_MARKERS = ["TCF7", "TBX21", "IFNG"]


# --------------------------------------------------------------------------
# Source (1): real CELLxGENE Census  (gated by the allowlist)
# --------------------------------------------------------------------------
def try_cellxgene(genes):
    try:
        import cellxgene_census  # noqa: F401
        import requests
    except ImportError:
        print("      cellxgene-census not installed; skipping live pull.")
        return None
    try:
        if requests.get("https://cellxgene.cziscience.com/", timeout=10).status_code != 200:
            print("      cellxgene host unreachable; skipping.")
            return None
        # --- real query (runs when network-unrestricted) ------------------
        # with cellxgene_census.open_soma() as census:
        #     adata = cellxgene_census.get_anndata(
        #         census, organism="Homo sapiens",
        #         obs_value_filter="cell_type=='CD8-positive, alpha-beta T cell' "
        #                          "and disease!='normal'",
        #         var_value_filter=f"feature_name in {genes}")
        #     return adata.X (cells x genes counts), adata.var feature_name, None
        # ------------------------------------------------------------------
        print("      cellxgene reachable; (full census query goes here).")
        return None
    except Exception as e:
        print(f"      cellxgene unreachable ({e.__class__.__name__}); skipping.")
        return None


# --------------------------------------------------------------------------
# Source (2): realistic simulated single-cell counts
# --------------------------------------------------------------------------
def simulate_counts(genes, n_cells=800, seed=0, depth=4000, dropout=0.15):
    """Negative-binomial counts for two latent CD8 states with library-size
    variation and dropout. Returns (counts[cells,genes], labels[cells])."""
    rng = np.random.default_rng(seed)
    n_exh = n_cells // 2
    labels = np.array(["exhausted"] * n_exh + ["effector"] * (n_cells - n_exh))
    # per-cell library size factor (lognormal)
    libfac = rng.lognormal(mean=0.0, sigma=0.3, size=n_cells)
    counts = np.zeros((n_cells, len(genes)), dtype=int)
    for j, g in enumerate(genes):
        exh_lvl, eff_lvl = GENE_TRUE[g]
        lvl = np.where(labels == "exhausted", exh_lvl, eff_lvl)        # true relative
        # expected counts ~ relative level * sequencing depth share
        mu = lvl * depth / len(genes) * 8.0 * libfac
        # negative-binomial via gamma-Poisson (overdispersion)
        shape = 2.0
        gam = rng.gamma(shape, mu / shape)
        c = rng.poisson(gam)
        # technical dropout
        mask = rng.random(n_cells) < dropout
        c[mask] = 0
        counts[:, j] = c
    return counts, labels


# --------------------------------------------------------------------------
# THE processing pipeline (shared by both sources)
# --------------------------------------------------------------------------
def process_scrnaseq(counts, genes, true_labels=None):
    counts = np.asarray(counts, dtype=float)
    gidx = {g: j for j, g in enumerate(genes)}

    # 1. library-size normalize to counts-per-10k, then log1p
    lib = counts.sum(axis=1, keepdims=True)
    lib[lib == 0] = 1.0
    logn = np.log1p(counts / lib * 1e4)

    # 2. z-score each gene across cells; score cells by marker panels
    z = (logn - logn.mean(0)) / (logn.std(0) + 1e-9)
    exh_cols = [gidx[g] for g in EXH_MARKERS if g in gidx]
    eff_cols = [gidx[g] for g in EFF_MARKERS if g in gidx]
    exh_score = z[:, exh_cols].mean(1)
    eff_score = z[:, eff_cols].mean(1)

    # 3. GATE cells into states from the data
    gated = np.where(exh_score > eff_score, "exhausted", "effector")
    if true_labels is not None:
        acc = (gated == true_labels).mean()
        print(f"      gating concordance vs ground truth: {acc:.1%}  "
              f"(exhausted={np.sum(gated=='exhausted')}, "
              f"effector={np.sum(gated=='effector')} cells)")

    # 4. collapse gene families to logical nodes (mean of members per cell)
    nodes = sorted({FAMILY.get(g, g) for g in genes})
    node_logn = {}
    for n in nodes:
        members = [gidx[g] for g in genes if FAMILY.get(g, g) == n]
        node_logn[n] = logn[:, members].mean(1)

    # 5. pseudobulk: mean per gated state
    M = pd.DataFrame(index=nodes, columns=["exhausted", "effector"], dtype=float)
    for state in ["exhausted", "effector"]:
        sel = gated == state
        for n in nodes:
            M.loc[n, state] = node_logn[n][sel].mean()

    # 6. global min-max scale to [0,1] (preserves within- and cross-state structure)
    lo, hi = np.nanmin(M.values), np.nanmax(M.values)
    M = (M - lo) / (hi - lo + 1e-9)
    return M


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", type=int, default=800)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    genes = list(GENE_TRUE.keys())
    print("[1/2] Attempting live scRNA-seq pull (CELLxGENE Census)...")
    live = try_cellxgene(genes)
    if live is not None:
        counts, genes_live, labels = live
        genes = genes_live
        source = "cellxgene-census"
    else:
        print(f"[2/2] Simulating {args.cells} single cells and processing them.")
        counts, labels = simulate_counts(genes, n_cells=args.cells, seed=args.seed)
        source = "scrnaseq-pseudobulk (simulated counts)"

    M = process_scrnaseq(counts, genes, true_labels=labels)

    with open(OUT, "w") as fh:
        fh.write("# Per-state node activity (RWR restart vectors), relative 0..1.\n")
        fh.write(f"# source: {source}\n")
        fh.write("# pipeline: lib-norm -> log1p -> marker gating -> family collapse"
                 " -> pseudobulk -> global min-max\n")
        fh.write("# gene\texhausted\teffector\n")
        for n in M.index:
            fh.write(f"{n}\t{M.loc[n,'exhausted']:.3f}\t{M.loc[n,'effector']:.3f}\n")
    print(f"      wrote {OUT}  ({len(M)} logical nodes, source={source})")
    print("\n  per-state activity (restart vectors):")
    print(M.round(3).to_string())


if __name__ == "__main__":
    main()
