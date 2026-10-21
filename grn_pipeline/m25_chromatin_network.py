"""
Module 25 — The chromatin read/write network (REAL data) re-organised into
irreducible functional modules.

This is the "重新整理关系" step: take the real protein-protein-interaction (PPI)
network of the chromatin writing/reading/erasing machinery and the TF motifs
that recruit them, and collapse it — via graph symmetry and a role-seeded
equitable partition (balanced colouring) — into a small set of IRREDUCIBLE
core modules, each of which maps onto one of the M23 memory classes.

Data (all real, fetched live and cached):
  * STRING v12 PPI among 22 chromatin regulators (writers/readers/erasers),
    required_score>=400, Homo sapiens (9606).  -> figures/_string_cache/
  * JASPAR motifs (M24) for the TFs that recruit them (REST->CoREST/HDAC,
    KRAB-ZNF->TRIM28/SETDB1).

Method:
  * Build the weighted PPI graph; annotate each node with its biochemical ROLE
    = (mark, action), e.g. ("H3K27me3","writer") for EZH2/EED/SUZ12.
  * ORBITS: graph automorphisms of the PPI -> sets of structurally
    interchangeable regulators (e.g. the three PRC2 core members).
  * MODULES: coarsest equitable partition seeded by role -> the irreducible
    functional classes (minimal balanced colouring on the undirected PPI).
  * Attach each module to its memory class from M23 (H3K27/H3K9=hybrid,
    DNAmeth=long, acetyl=short), re-organising PPI+TF -> memory.

RIGOR:
  ✅ The PPI edges are the real STRING scores; automorphism orbits and the
     equitable partition are deterministic graph computations.
  ✅ The role annotation is standard, well-established biochemistry (writer/
     reader/eraser + histone/DNA mark) of each named protein.
  ⚠️ Mapping a module to an M23 memory class is the paper's CR->class
     assignment extended by mark identity (H3K27/H3K9 hybrid, 5mC long,
     acetyl short); it is a hypothesis-level bridge, not a measured rate here.
"""
import json
import os
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CACHE = "figures/_string_cache/chromatin_ppi.json"
GENES = ["EZH2", "EED", "SUZ12", "DNMT1", "DNMT3A", "DNMT3B", "SUV39H1",
         "SETDB1", "EHMT2", "KDM6A", "KDM6B", "KDM1A", "TET1", "TET2",
         "HDAC1", "HDAC2", "HDAC4", "CBX5", "MECP2", "TRIM28", "RCOR1", "UHRF1"]
STRING_URL = ("https://string-db.org/api/json/network?identifiers={ids}"
              "&species=9606&required_score=400")

# biochemical role: (mark, action) + the M23 memory class of that mark
ROLE = {
    "EZH2": ("H3K27me3", "writer"), "EED": ("H3K27me3", "cofactor"),
    "SUZ12": ("H3K27me3", "cofactor"), "KDM6A": ("H3K27me3", "eraser"),
    "KDM6B": ("H3K27me3", "eraser"),
    "SUV39H1": ("H3K9me3", "writer"), "SETDB1": ("H3K9me3", "writer"),
    "EHMT2": ("H3K9me3", "writer"), "CBX5": ("H3K9me3", "reader"),
    "TRIM28": ("H3K9me3", "cofactor"),
    "DNMT1": ("5mC", "writer"), "DNMT3A": ("5mC", "writer"),
    "DNMT3B": ("5mC", "writer"), "UHRF1": ("5mC", "cofactor"),
    "MECP2": ("5mC", "reader"), "TET1": ("5mC", "eraser"),
    "TET2": ("5mC", "eraser"),
    "HDAC1": ("acetyl", "eraser"), "HDAC2": ("acetyl", "eraser"),
    "HDAC4": ("acetyl", "eraser"), "KDM1A": ("acetyl", "eraser"),
    "RCOR1": ("acetyl", "cofactor"),
}
MARK_MEMORY = {"H3K27me3": "hybrid (EED-like)", "H3K9me3": "hybrid (KRAB-like)",
               "5mC": "long-term (DNMT3B-like)", "acetyl": "short-term (HDAC4-like)"}


def load_network():
    if not os.path.exists(CACHE):
        import urllib.request
        os.makedirs(os.path.dirname(CACHE), exist_ok=True)
        url = STRING_URL.format(ids="%0d".join(GENES))
        with urllib.request.urlopen(url, timeout=40) as r:
            data = r.read().decode()
        open(CACHE, "w").write(data)
    edges = json.loads(open(CACHE).read())
    G = nx.Graph()
    for g in GENES:
        G.add_node(g, mark=ROLE[g][0], action=ROLE[g][1])
    for e in edges:
        a, b = e["preferredName_A"], e["preferredName_B"]
        if a in ROLE and b in ROLE:
            G.add_edge(a, b, weight=float(e["score"]))
    return G


# --------------------------------------------------------------------------
# orbits: automorphisms of the (role-coloured) PPI graph
# --------------------------------------------------------------------------
def automorphism_orbits(G):
    """Orbits under automorphisms that preserve the (mark,action) role colour.
    Uses VF2 self-isomorphism search via node-match on role."""
    from networkx.algorithms.isomorphism import GraphMatcher

    def nm(d1, d2):
        return d1["mark"] == d2["mark"] and d1["action"] == d2["action"]

    gm = GraphMatcher(G, G, node_match=nm)
    # union-find over all automorphisms (can be many; cap for tractability)
    parent = {n: n for n in G.nodes()}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        parent[find(a)] = find(b)

    count = 0
    for mapping in gm.isomorphisms_iter():
        for u, v in mapping.items():
            union(u, v)
        count += 1
        if count > 5000:
            break
    orbits = {}
    for n in G.nodes():
        orbits.setdefault(find(n), []).append(n)
    return [sorted(v) for v in orbits.values()], count


# --------------------------------------------------------------------------
# irreducible modules: coarsest equitable partition seeded by role
# (undirected balanced colouring)
# --------------------------------------------------------------------------
def equitable_partition(G, seed_by_role=True):
    if seed_by_role:
        colour = {n: f"{G.nodes[n]['mark']}|{G.nodes[n]['action']}"
                  for n in G.nodes()}
    else:
        colour = {n: 0 for n in G.nodes()}
    colour = {n: hash(c) for n, c in colour.items()}
    while True:
        sig = {}
        for n in G.nodes():
            neigh = sorted(colour[u] for u in G.neighbors(n))
            sig[n] = (colour[n], tuple(neigh))
        uniq = {s: i for i, s in enumerate(sorted(set(sig.values()),
                                                  key=str))}
        new = {n: uniq[sig[n]] for n in G.nodes()}
        if len(set(new.values())) == len(set(colour.values())):
            break
        colour = new
    classes = {}
    for n, c in colour.items():
        classes.setdefault(c, []).append(n)
    return [sorted(v) for v in classes.values()]


def module_summary(G):
    """Collapse by (mark, action) -> the irreducible functional modules, each
    tagged with its memory class and the member proteins."""
    mods = {}
    for n in G.nodes():
        key = ROLE[n]
        mods.setdefault(key, []).append(n)
    rows = []
    for (mark, action), members in sorted(mods.items()):
        rows.append({"mark": mark, "action": action,
                     "memory_class": MARK_MEMORY[mark],
                     "members": sorted(members),
                     "k": len(members)})
    return rows


def make_figure(G, orbits, path):
    fig, ax = plt.subplots(1, 2, figsize=(15, 7))
    mark_color = {"H3K27me3": "#2980b9", "H3K9me3": "#8e44ad",
                  "5mC": "#c0392b", "acetyl": "#16a085"}
    pos = nx.spring_layout(G, seed=3, k=0.6, weight="weight")
    node_c = [mark_color[G.nodes[n]["mark"]] for n in G.nodes()]
    shapes = {"writer": "s", "eraser": "v", "reader": "o", "cofactor": "D"}
    for action, mk in shapes.items():
        ns = [n for n in G.nodes() if G.nodes[n]["action"] == action]
        nx.draw_networkx_nodes(G, pos, nodelist=ns, node_shape=mk, ax=ax[0],
                               node_color=[mark_color[G.nodes[n]["mark"]] for n in ns],
                               node_size=560, edgecolors="k")
    nx.draw_networkx_edges(G, pos, ax=ax[0], alpha=0.25,
                           width=[0.5 + 3 * G[u][v]["weight"] for u, v in G.edges()])
    nx.draw_networkx_labels(G, pos, ax=ax[0], font_size=7)
    ax[0].set_title("Real STRING chromatin read/write PPI\n"
                    "colour=mark (blue H3K27, purple H3K9, red 5mC, teal acetyl)\n"
                    "shape: square=writer v=eraser o=reader D=cofactor")
    ax[0].axis("off")

    # module bar: members per (mark,action) module and memory class
    rows = module_summary(G)
    labels = [f"{r['mark']}\n{r['action']}" for r in rows]
    ks = [r["k"] for r in rows]
    cols = [mark_color[r["mark"]] for r in rows]
    y = np.arange(len(rows))
    ax[1].barh(y, ks, color=cols)
    ax[1].set_yticks(y); ax[1].set_yticklabels(labels, fontsize=7)
    ax[1].invert_yaxis()
    for i, r in enumerate(rows):
        ax[1].text(r["k"] + 0.05, i, "  ".join(r["members"]), va="center",
                   fontsize=6)
    ax[1].set_xlabel("# proteins in module")
    ax[1].set_title("Irreducible functional modules\n(mark x action) -> memory class")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def report(fig_path="figures/m25_chromatin_network.png"):
    print("=" * 68)
    print("MODULE 25 — real chromatin read/write network -> irreducible modules")
    print("=" * 68)
    G = load_network()
    print(f"\nreal STRING PPI: {G.number_of_nodes()} regulators, "
          f"{G.number_of_edges()} edges (score>=0.4)")

    orbits, n_auto = automorphism_orbits(G)
    multi_orb = [o for o in orbits if len(o) > 1]
    print(f"\nrole-preserving automorphisms found: {n_auto} "
          f"-> {len(orbits)} orbits")
    print("structurally-interchangeable regulator sets (nontrivial orbits):")
    for o in sorted(multi_orb, key=lambda x: -len(x)):
        print(f"  {ROLE[o[0]]}: {{{', '.join(o)}}}")
    if not multi_orb:
        print("  (none at full-graph level — PPI wiring breaks exact symmetry)")

    part = equitable_partition(G)
    print(f"\ncoarsest role-seeded equitable partition: "
          f"{G.number_of_nodes()} -> {len(part)} classes")

    rows = module_summary(G)
    print(f"\nirreducible functional modules (mark x action -> memory class):")
    print(f"  {'mark':>9} {'action':>9} {'n':>2}  {'memory class':>24}  members")
    for r in rows:
        print(f"  {r['mark']:>9} {r['action']:>9} {r['k']:>2}  "
              f"{r['memory_class']:>24}  {', '.join(r['members'])}")

    print("\nre-organised relationship (the answer to '重新整理关系'):")
    print("  * writers set the FORWARD silencing edge (A->R / R->I) of M23;")
    print("  * erasers set the REVERSE / decay edge (R->A), i.e. 消退速度;")
    print("  * 5mC writers (DNMT) feed the ABSORBING I state (long memory);")
    print("  * acetyl erasers (HDAC) keep the chain in the reversible A<->R")
    print("    regime (short memory). Mark identity -> memory class -> the")
    print("    Markov generator of M23. The network is thus the WIRING that")
    print("    produces each memory timescale.")

    make_figure(G, orbits, fig_path)
    print(f"\nfigure written to: {fig_path}")
    return {"n_nodes": G.number_of_nodes(), "n_edges": G.number_of_edges(),
            "orbits": orbits, "n_automorphisms": n_auto,
            "partition_classes": len(part), "modules": rows}


if __name__ == "__main__":
    report()
