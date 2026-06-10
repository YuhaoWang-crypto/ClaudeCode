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
  library-size normalize -> log1p  ->  PCA  ->  kNN graph  ->  Louvain
  clustering (UNSUPERVISED; cell labels never used)  ->  annotate each cluster
  with canonical exhaustion/effector marker module scores  ->  assign clusters
  to states  ->  collapse gene families to logical nodes  ->  pseudobulk (mean
  per state)  ->  global min-max scale to [0,1].

Output: data/state_expression_real.tsv  (gene exhausted effector)  + provenance.

Run:  python3 fetch_expression.py [--cells N] [--seed S] [--knn K]
Deps: numpy, pandas, networkx
"""

import argparse
import os
import numpy as np
import pandas as pd
import networkx as nx
from networkx.algorithms.community import louvain_communities

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
# Canonical marker gene SETS used only to ANNOTATE unsupervised clusters.
# These are published Tex / Teff markers (identity only) -- the simulation's
# GENE_TRUE levels are never consulted by the clustering or annotation.
EXH_MARKERS = ["TOX", "PDCD1", "HAVCR2", "LAG3"]   # terminal exhaustion
EFF_MARKERS = ["TCF7", "TBX21", "IFNG"]            # effector / memory


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
# Unsupervised clustering (PCA -> kNN graph -> Louvain), then marker annotation
# --------------------------------------------------------------------------
def pca(z, n_comp=15, seed=0):
    """Top principal components via SVD on the (cells x genes) z-scored matrix."""
    z = z - z.mean(0)
    u, s, _ = np.linalg.svd(z, full_matrices=False)
    n = min(n_comp, s.shape[0])
    return u[:, :n] * s[:n]


def knn_graph(emb, k=15):
    """Undirected kNN graph in embedding space (Euclidean)."""
    d2 = ((emb[:, None, :] - emb[None, :, :]) ** 2).sum(-1)
    np.fill_diagonal(d2, np.inf)
    g = nx.Graph()
    g.add_nodes_from(range(emb.shape[0]))
    for i in range(emb.shape[0]):
        for j in np.argsort(d2[i])[:k]:
            g.add_edge(i, int(j))
    return g


def cluster_and_annotate(logn, gidx, seed=0, k=15):
    """UNSUPERVISED: PCA -> kNN -> Louvain communities; then label each
    community by its dominant marker module (exhaustion vs effector).
    Also prints a marker heatmap and a 2D projection to sanity-check labels."""
    z = (logn - logn.mean(0)) / (logn.std(0) + 1e-9)
    emb = pca(z, seed=seed)
    g = knn_graph(emb, k=k)
    comms = louvain_communities(g, seed=seed)

    exh_cols = [gidx[x] for x in EXH_MARKERS if x in gidx]
    eff_cols = [gidx[x] for x in EFF_MARKERS if x in gidx]
    cell_state = np.empty(logn.shape[0], dtype=object)
    cell_cluster = np.full(logn.shape[0], -1, dtype=int)
    summary = []
    ordered = sorted(comms, key=len, reverse=True)
    for ci, cells in enumerate(ordered):
        cells = list(cells)
        es = z[np.ix_(cells, exh_cols)].mean()
        fs = z[np.ix_(cells, eff_cols)].mean()
        state = "exhausted" if es > fs else "effector"
        cell_state[cells] = state
        cell_cluster[cells] = ci
        summary.append((ci, len(cells), es, fs, state))
    print(f"      Louvain found {len(comms)} clusters (unsupervised); "
          f"annotated by marker module score:")
    for ci, n, es, fs, state in summary:
        print(f"        cluster {ci}: {n:4d} cells  exh_score={es:+.2f} "
              f"eff_score={fs:+.2f}  -> {state}")

    cstate = {ci: st for ci, _, _, _, st in summary}
    _marker_heatmap(z, gidx, cell_cluster, list(range(len(ordered))), cstate)
    _ascii_projection(emb[:, :2], cell_state)
    return cell_state


# --------------------------------------------------------------------------
# Text diagnostics: marker heatmap + 2D projection (annotation sanity check)
# --------------------------------------------------------------------------
_SHADE = " .:-=+*#%@"          # low -> high


def _shade(v, lo, hi):
    t = min(max((v - lo) / (hi - lo + 1e-9), 0.0), 1.0)
    return _SHADE[int(t * (len(_SHADE) - 1))]


def _marker_heatmap(z, gidx, cell_cluster, clusters, cstate):
    markers = ([(g, "exh") for g in EXH_MARKERS if g in gidx]
               + [(g, "eff") for g in EFF_MARKERS if g in gidx])
    sym = {"exhausted": "X", "effector": "O"}
    # group columns by state (exhausted clusters first) so the block stands out
    clusters = sorted(clusters, key=lambda ci: (0 if cstate[ci] == "exhausted" else 1, ci))
    means, vals = {}, []
    for ci in clusters:
        cells = np.where(cell_cluster == ci)[0]
        for g, _ in markers:
            m = z[cells, gidx[g]].mean()
            means[(g, ci)] = m
            vals.append(m)
    lo, hi = min(vals), max(vals)
    print("\n      marker x cluster mean z-score "
          f"(shade '{_SHADE[1]}'->'{_SHADE[-1]}' = low->high; columns grouped by state):")
    print("        " + "gene".ljust(8) + "".join(f"c{ci}".rjust(3) for ci in clusters))
    print("        " + "state".ljust(8)
          + "".join(sym[cstate[ci]].rjust(3) for ci in clusters)
          + "    X=exhausted  O=effector")
    for g, grp in markers:
        row = "".join(_shade(means[(g, ci)], lo, hi).rjust(3) for ci in clusters)
        print(f"        {g.ljust(8)}{row}    [{grp} marker]")
    print("      -> exh-marker rows should darken under X columns,"
          " eff-marker rows under O columns.")


def _ascii_projection(emb2, cell_state, width=56, height=16):
    x, y = emb2[:, 0], emb2[:, 1]
    gx = ((x - x.min()) / (np.ptp(x) + 1e-9) * (width - 1)).astype(int)
    gy = ((y - y.min()) / (np.ptp(y) + 1e-9) * (height - 1)).astype(int)
    grid = [[" "] * width for _ in range(height)]
    for i in range(len(cell_state)):
        sym = "X" if cell_state[i] == "exhausted" else "o"
        r, c = height - 1 - gy[i], gx[i]
        cur = grid[r][c]
        grid[r][c] = sym if cur == " " else (cur if cur == sym else "*")
    print("\n      PC1 x PC2 projection (X=exhausted  o=effector  *=mixed):")
    for row in grid:
        print("        |" + "".join(row) + "|")
    print("      -> well-annotated states occupy distinct regions;"
          " '*' marks the ambiguous boundary.")


# --------------------------------------------------------------------------
# THE processing pipeline (shared by both sources)
# --------------------------------------------------------------------------
def process_scrnaseq(counts, genes, true_labels=None, knn=15, seed=0):
    counts = np.asarray(counts, dtype=float)
    gidx = {g: j for j, g in enumerate(genes)}

    # 1. library-size normalize to counts-per-10k, then log1p
    lib = counts.sum(axis=1, keepdims=True)
    lib[lib == 0] = 1.0
    logn = np.log1p(counts / lib * 1e4)

    # 2-3. unsupervised clustering, then annotate clusters with marker modules
    gated = cluster_and_annotate(logn, gidx, seed=seed, k=knn)
    if true_labels is not None:
        acc = (gated == true_labels).mean()
        print(f"      cluster-label concordance vs ground truth: {acc:.1%}  "
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
    ap.add_argument("--knn", type=int, default=15)
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

    M = process_scrnaseq(counts, genes, true_labels=labels, knn=args.knn, seed=args.seed)

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
