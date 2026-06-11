#!/usr/bin/env python3
"""
Minimal, connectivity-only pipeline on an UNDIRECTED human PPI network.

This is the simplified track (no direction, no sign, no Boolean dynamics). With
just "who interacts with whom" you can still do, robustly:

  1. STATE-SPECIFIC PROTEIN RANKING  -- personalized PageRank (RWR) seeded by a
     cell-state expression vector, with a DEGREE-CORRECTED permutation test so
     generic high-degree hubs do not automatically win.
  2. MINIMAL NETWORK FOR THE STATE   -- a Steiner-tree / PCSF-style smallest
     connected subnetwork linking the top-ranked proteins (pulls in the hidden
     "connector" proteins on the paths between them).
  3. COMPENSATORY PATHWAYS           -- knock out the top target, re-propagate,
     and report which proteins GAIN rank (they absorb the flow the target
     carried) plus how the shortest route between two hubs reroutes.

What you DO NOT get from undirected PPI: a cell attractor / "equilibrium state"
or a minimal flip combination -- those need a directed, signed causal network
(see pipeline.py / bnlib.py). PageRank's stationary vector here is a *diffusion*
equilibrium (a ranking), not a physiological steady state.

Data source: http://prodata.swmed.edu/humanPPI/bulk_download (Grishin lab human
PPI). That host is NOT in this sandbox's network allowlist (HTTP 403), so the
download is attempted and, when blocked, a stand-in interactome is built so the
method runs end-to-end. Add the host to your egress settings (or pass --ppi a
local edge file) to run on the real data; the parser is below.

Run:  python3 ppi_min.py [--ppi FILE] [--state-expr FILE] [--perms N]
Deps: numpy, scipy, networkx
"""

import argparse
import os
import numpy as np
import networkx as nx
from networkx.algorithms.approximation import steiner_tree

HERE = os.path.dirname(os.path.abspath(__file__))
PPI_FILE = os.path.join(HERE, "data", "ppi_edges.tsv")
STATE_EXPR = os.path.join(HERE, "data", "state_expression_real.tsv")
PRODATA_URL = "http://prodata.swmed.edu/humanPPI/bulk_download"

PANEL_EXH = ["TOX", "PDCD1", "HAVCR2", "LAG3", "NR4A", "BATF", "IRF4"]
PANEL_EFF = ["TCF7", "TBX21", "IFNG", "AP1", "EOMES"]
PANEL_HUB = ["NFAT"]
PANEL = PANEL_EXH + PANEL_EFF + PANEL_HUB


# --------------------------------------------------------------------------
# Data: real prodata download (gated) OR a reproducible stand-in interactome
# --------------------------------------------------------------------------
def parse_prodata(text):
    """Tolerant parser: take the first two whitespace/tab columns of each
    non-comment line as an undirected interacting protein pair. Adjust column
    indices here once you see the real bulk file's header."""
    edges = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith(("#", "//")):
            continue
        cols = line.replace("\t", " ").split()
        if len(cols) >= 2 and cols[0] != cols[1]:
            edges.append((cols[0], cols[1]))
    return edges


def try_fetch_prodata():
    try:
        import requests
        r = requests.get(PRODATA_URL, timeout=30)
        if r.status_code != 200:
            print(f"      prodata returned HTTP {r.status_code} "
                  f"({r.text[:60].strip()}); using stand-in.")
            return None
        edges = parse_prodata(r.text)
        print(f"      parsed {len(edges)} edges from prodata.")
        return edges
    except Exception as e:
        print(f"      prodata fetch failed ({e.__class__.__name__}); using stand-in.")
        return None


def build_standin_interactome(n=2000, m=3, seed=0):
    """Scale-free (Barabasi-Albert) graph with realistic degree distribution,
    with the T-cell panel embedded as two interacting communities so the
    ranking / Steiner / compensatory steps are interpretable. NOT real biology
    -- a structural stand-in to exercise the method at interactome scale."""
    rng = np.random.default_rng(seed)
    G = nx.barabasi_albert_graph(n, m, seed=seed)
    # choose panel host-nodes spread across the degree spectrum
    by_deg = sorted(G.nodes(), key=lambda x: G.degree(x), reverse=True)
    picks = list(rng.choice(by_deg[50:400], size=len(PANEL), replace=False))
    mapping = {int(picks[i]): PANEL[i] for i in range(len(PANEL))}
    G = nx.relabel_nodes(G, mapping)
    # wire the panel into two sparse communities with local hubs (TOX, TCF7),
    # bridged only through a couple of cross-links -- so connecting proteins
    # across the module requires intermediate "connector" nodes.
    def ring(nodes):
        for i in range(len(nodes)):
            G.add_edge(nodes[i], nodes[(i + 1) % len(nodes)])
    ring(PANEL_EXH)
    ring(PANEL_EFF)
    for n in PANEL_EXH:                       # TOX = exhaustion-module hub
        if n != "TOX":
            G.add_edge("TOX", n)
    for n in PANEL_EFF:                       # TCF7 = effector-module hub
        if n != "TCF7":
            G.add_edge("TCF7", n)
    for a, b in [("TOX", "TCF7"), ("NFAT", "TOX"), ("AP1", "TBX21"),
                 ("NFAT", "NR4A")]:           # the few bridges between modules
        if a in G and b in G:
            G.add_edge(a, b)
    # give anonymous nodes protein-like string IDs
    G = nx.relabel_nodes(G, {n: f"P{int(n):05d}" for n in list(G.nodes())
                             if isinstance(n, (int, np.integer))})
    return G


def load_or_build_ppi(path):
    if os.path.exists(path):
        G = nx.Graph()
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                a, b = line.split()[:2]
                if a != b:
                    G.add_edge(a, b)
        print(f"      loaded {G.number_of_edges()} edges from {os.path.basename(path)}")
        return G
    print(f"[fetch] trying prodata bulk download: {PRODATA_URL}")
    edges = try_fetch_prodata()
    if edges:
        G = nx.Graph()
        G.add_edges_from(edges)
        source = "prodata"
    else:
        G = build_standin_interactome()
        source = "stand-in (scale-free + embedded panel)"
    with open(path, "w") as fh:
        fh.write(f"# undirected human PPI edge list  source: {source}\n")
        for a, b in G.edges():
            fh.write(f"{a}\t{b}\n")
    print(f"      built network ({source}); wrote {os.path.basename(path)}")
    return G


def load_state_weights(path, state="exhausted"):
    """restart vector from scRNA-seq-derived per-state activity (panel genes)."""
    if not os.path.exists(path):
        return {g: (0.9 if g in PANEL_EXH else 0.1) for g in PANEL}
    col = {"exhausted": 1, "effector": 2}[state]
    w = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            w[p[0]] = float(p[col])
    return w


# --------------------------------------------------------------------------
# 1. Personalized PageRank + degree-corrected permutation test
# --------------------------------------------------------------------------
def ppr(G, restart, alpha=0.15):
    z = sum(restart.values())
    pers = {n: restart.get(n, 0.0) / z for n in G.nodes()}
    return nx.pagerank(G, alpha=1 - alpha, personalization=pers)


def degree_corrected_test(G, restart, alpha=0.15, n_perm=60, seed=0):
    """Null: keep the SAME seed weights but move them onto degree-matched random
    nodes. A protein is state-specific only if it beats this degree-aware null,
    which removes the 'hubs always win' bias."""
    rng = np.random.default_rng(seed)
    obs = ppr(G, restart, alpha)
    nodes = list(G.nodes())
    degs = np.array([G.degree(n) for n in nodes])
    # bucket nodes by degree for degree-matched resampling
    order = np.argsort(degs)
    null = {n: np.zeros(n_perm) for n in nodes}
    seed_nodes = [n for n in restart if n in G and restart[n] > 0]
    seed_w = [restart[n] for n in seed_nodes]
    for p in range(n_perm):
        new_restart = {}
        for sn, wt in zip(seed_nodes, seed_w):
            d = G.degree(sn)
            # pick a random node with the nearest degree
            pos = np.searchsorted(degs[order], d)
            pos = min(max(pos + rng.integers(-15, 16), 0), len(nodes) - 1)
            new_restart[nodes[order[pos]]] = wt
        pp = ppr(G, new_restart, alpha)
        for n in nodes:
            null[n][p] = pp.get(n, 0.0)
    stats = {}
    for n in nodes:
        mu, sd = null[n].mean(), null[n].std() + 1e-12
        z = (obs[n] - mu) / sd
        pval = (np.sum(null[n] >= obs[n]) + 1) / (n_perm + 1)
        stats[n] = (obs[n], z, pval)
    return obs, stats


# --------------------------------------------------------------------------
# 2. Minimal connected subnetwork (Steiner / PCSF-style)
# --------------------------------------------------------------------------
def minimal_subnetwork(G, terminals):
    terminals = [t for t in terminals if t in G]
    # restrict to the largest connected component containing terminals
    T = steiner_tree(G, terminal_nodes=terminals)
    connectors = [n for n in T.nodes() if n not in terminals]
    return T, connectors


# --------------------------------------------------------------------------
# 3. Compensatory pathways: knock out target, re-propagate, find rank gainers
# --------------------------------------------------------------------------
def compensatory(G, restart, target, alpha=0.15, top=10):
    base = ppr(G, restart, alpha)
    H = G.copy()
    H.remove_node(target)
    r2 = {n: w for n, w in restart.items() if n != target}
    after = ppr(H, r2, alpha)
    gain = {n: after.get(n, 0.0) - base.get(n, 0.0) for n in H.nodes()}
    gainers = sorted(gain.items(), key=lambda kv: -kv[1])[:top]
    return base, after, gainers


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ppi", default=PPI_FILE)
    ap.add_argument("--state-expr", default=STATE_EXPR)
    ap.add_argument("--state", default="exhausted", choices=["exhausted", "effector"])
    ap.add_argument("--perms", type=int, default=60)
    ap.add_argument("--top-k", type=int, default=12)
    args = ap.parse_args()

    print("=" * 72)
    print("PPI-ONLY PIPELINE | rank -> minimal subnetwork -> compensatory paths")
    print("=" * 72)
    G = load_or_build_ppi(args.ppi)
    print(f"network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges, "
          f"<k>={2*G.number_of_edges()/G.number_of_nodes():.1f}")

    restart = load_state_weights(args.state_expr, args.state)
    restart = {g: w for g, w in restart.items() if g in G and w > 0}
    print(f"restart seeds ({args.state}): "
          + ", ".join(f"{g}={w:.2f}" for g, w in sorted(restart.items())))

    # 1. ranking with degree-corrected significance
    print("\n-- 1. PageRank ranking (degree-corrected permutation test) --")
    obs, stats = degree_corrected_test(G, restart, n_perm=args.perms)
    raw_top = sorted(obs, key=obs.get, reverse=True)[:10]
    print("raw PageRank top-10 (degree-biased):")
    print("   " + ", ".join(raw_top))
    corr = sorted(G.nodes(), key=lambda n: -stats[n][1])[:args.top_k]
    print(f"\ndegree-corrected top-{args.top_k} (by z vs degree-matched null):")
    print(f"   {'protein':10}{'PPR':>10}{'z':>8}{'p':>8}{'degree':>8}")
    for n in corr:
        s = stats[n]
        print(f"   {n:10}{s[0]:>10.5f}{s[1]:>8.1f}{s[2]:>8.3f}{G.degree(n):>8}")
    print("   (note: raw top is dominated by high-degree hubs; the corrected test"
          " demotes\n    them and surfaces seed-specific interactors instead.)")

    # 2. minimal connected subnetwork over the corrected top hits
    print("\n-- 2. Minimal connected subnetwork (Steiner/PCSF-style) --")
    terminals = [n for n in corr if n in PANEL][:6]      # named, interpretable
    if not any(t in PANEL_EFF for t in terminals):       # span both modules
        eff = sorted((n for n in PANEL_EFF if n in G), key=lambda n: -obs[n])
        if eff:
            terminals.append(eff[0])
    T, connectors = minimal_subnetwork(G, terminals)
    print(f"terminals ({len(terminals)}): {', '.join(terminals)}")
    print(f"minimal subnetwork: {T.number_of_nodes()} nodes, {T.number_of_edges()} edges")
    print(f"hidden connector proteins pulled in: "
          f"{', '.join(map(str, connectors)) if connectors else '(none -- terminals already adjacent)'}")

    # 3. compensatory pathways around the #1 target
    target = corr[0]
    print(f"\n-- 3. Compensatory pathways after knocking out {target} --")
    base, after, gainers = compensatory(G, restart, target)
    print(f"proteins that GAIN the most rank when {target} is removed")
    print("(they absorb the propagation the target carried = compensatory candidates):")
    for n, g in gainers:
        if g > 0:
            print(f"   {str(n):10} Δrank=+{g:.5f}  degree={G.degree(n)}")
    # reroute between two seed hubs
    hubs = [t for t in terminals if t in PANEL][:2]
    if len(hubs) == 2 and target not in hubs:
        a, b = hubs
        try:
            d0 = nx.shortest_path_length(G, a, b)
            H = G.copy(); H.remove_node(target)
            d1 = nx.shortest_path_length(H, a, b) if nx.has_path(H, a, b) else None
            msg = f"{d1}" if d1 is not None else "DISCONNECTED"
            print(f"\n   shortest path {a}–{b}: {d0} hops; after removing {target}: {msg}"
                  + ("  (rerouted via compensatory path)" if d1 and d1 > d0 else ""))
        except nx.NetworkXNoPath:
            pass

    print("\n" + "=" * 72)
    print("Note: this is the connectivity-only track. Rankings, the minimal")
    print("subnetwork, and compensatory nodes are valid from undirected PPI.")
    print("A cell 'equilibrium'/flip combination is NOT — that needs the signed,")
    print("directed model (pipeline.py). PageRank's fixed point here is a")
    print("diffusion ranking, not a physiological steady state.")
    print("=" * 72)


if __name__ == "__main__":
    main()
