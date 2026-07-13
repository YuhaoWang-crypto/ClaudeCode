"""
subnetwork.py — turn an omics hit list into an interpretable module.

Given a loaded PPI graph (from ``network.py``) and a gene list (e.g. your
differentially expressed / mutated genes), extract the induced subnetwork,
optionally grow it by one neighbour shell, and rank the genes by network
centrality so you can spot the hubs / bottlenecks that organise the module.

Outputs are plain dicts / files (edge list, node table) so you can load them
straight into Cytoscape, gephi, or a notebook.
"""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable


def induced(graph, genes: Iterable[str], expand: int = 0):
    """Subgraph on ``genes``. If ``expand`` > 0, add the first-neighbour shell
    ``expand`` times (connectors that link your hits together)."""
    seeds = {g for g in genes if g in graph}
    nodes = set(seeds)
    frontier = set(seeds)
    for _ in range(expand):
        shell = set()
        for n in frontier:
            shell.update(graph.neighbors(n))
        frontier = shell - nodes
        nodes.update(shell)
    sub = graph.subgraph(nodes).copy()
    for n in sub.nodes:
        sub.nodes[n]["is_seed"] = n in seeds
    sub.graph["n_seeds_found"] = len(seeds)
    sub.graph["n_seeds_missing"] = len(set(genes)) - len(seeds)
    return sub


def rank_genes(subgraph) -> list[dict]:
    """Centrality table: degree, betweenness, clustering. Hubs (high degree) and
    bottlenecks (high betweenness) are the module's control points."""
    import networkx as nx

    if subgraph.number_of_nodes() == 0:
        return []
    deg = dict(subgraph.degree())
    btw = nx.betweenness_centrality(subgraph)
    clu = nx.clustering(subgraph)
    rows = [{
        "gene": n,
        "is_seed": subgraph.nodes[n].get("is_seed", False),
        "degree": deg[n],
        "betweenness": round(btw[n], 5),
        "clustering": round(clu[n], 4),
    } for n in subgraph.nodes]
    rows.sort(key=lambda r: (r["degree"], r["betweenness"]), reverse=True)
    return rows


def largest_component(subgraph):
    import networkx as nx
    if subgraph.number_of_nodes() == 0:
        return subgraph
    cc = max(nx.connected_components(subgraph), key=len)
    return subgraph.subgraph(cc).copy()


def write_cytoscape(subgraph, out_prefix: str | Path) -> dict[str, Path]:
    """Write ``<prefix>.edges.tsv`` and ``<prefix>.nodes.tsv`` for Cytoscape."""
    out_prefix = Path(out_prefix)
    edges = out_prefix.with_suffix(".edges.tsv")
    nodes = out_prefix.with_suffix(".nodes.tsv")

    with open(edges, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["source", "target", "score", "in_string", "in_biogrid",
                    "evidence"])
        for u, v, d in subgraph.edges(data=True):
            w.writerow([u, v, d.get("score", ""), d.get("in_string", ""),
                        d.get("in_biogrid", ""), d.get("evidence", "")])

    ranked = {r["gene"]: r for r in rank_genes(subgraph)}
    with open(nodes, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene", "is_seed", "degree", "betweenness", "clustering"])
        for n in subgraph.nodes:
            r = ranked.get(n, {})
            w.writerow([n, subgraph.nodes[n].get("is_seed", False),
                        r.get("degree", ""), r.get("betweenness", ""),
                        r.get("clustering", "")])
    return {"edges": edges, "nodes": nodes}


def summary(subgraph) -> dict:
    import networkx as nx
    n, m = subgraph.number_of_nodes(), subgraph.number_of_edges()
    lcc = largest_component(subgraph)
    return {
        "nodes": n,
        "edges": m,
        "seeds_found": subgraph.graph.get("n_seeds_found"),
        "seeds_missing": subgraph.graph.get("n_seeds_missing"),
        "density": round(nx.density(subgraph), 5) if n > 1 else 0.0,
        "largest_component_nodes": lcc.number_of_nodes(),
        "components": nx.number_connected_components(subgraph) if n else 0,
    }
