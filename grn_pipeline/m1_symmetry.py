"""
Module 1 — Symmetry reduction of a regulatory network.

Idea (rigorous):
  A network's structural symmetries are its graph automorphism group Aut(G).
  Nodes that map to each other under Aut(G) form ORBITS (an equitable
  partition). Collapsing each orbit to a single node gives the QUOTIENT
  network, which is smaller but dynamically compatible (Golubitsky-Stewart
  balanced colourings; MacArthur & Sanchez-Garcia 2008). The nodes that
  survive quotienting are the IRREDUCIBLE CORE.

Concrete system:
  A canonical RTK -> RAS -> RAF -> MEK -> ERK cascade in which the three
  human RAS paralogues (HRAS/KRAS/NRAS) share identical wiring. They are a
  genuine automorphism orbit -> the group S_3 (order 6) -> three upstream
  nodes collapse to one effective node. This is the mathematical realisation
  of "different upstream copies collapse onto a single downstream node".
"""
import networkx as nx
from networkx.algorithms.isomorphism import DiGraphMatcher
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def build_mapk_network():
    """RTK/RAS/RAF/MEK/ERK cascade with 3 redundant RAS paralogues.
    All edges are activating (sign=+1); edge signs are respected by the
    automorphism search so an activating edge can never map onto a
    repressing one."""
    G = nx.DiGraph()
    edges = [
        ("Ligand", "RTK", +1),
        ("RTK", "SOS", +1),
        ("SOS", "HRAS", +1),
        ("SOS", "KRAS", +1),
        ("SOS", "NRAS", +1),
        ("HRAS", "RAF", +1),
        ("KRAS", "RAF", +1),
        ("NRAS", "RAF", +1),
        ("RAF", "MEK", +1),
        ("MEK", "ERK", +1),
    ]
    for u, v, s in edges:
        G.add_edge(u, v, sign=s)
    return G


def automorphisms(G):
    """Return every automorphism of G as a node->node dict.
    Edge signs must be preserved (edge_match)."""
    def em(e1, e2):
        return e1.get("sign") == e2.get("sign")
    gm = DiGraphMatcher(G, G, edge_match=em)
    return list(gm.isomorphisms_iter())


def orbits(G, autos):
    """Orbit partition of the node set under the automorphism group."""
    orb = {}
    for v in G.nodes():
        img = frozenset(phi[v] for phi in autos)
        orb.setdefault(img, set()).add(v)
    # canonicalise: map each node to a representative orbit label
    partition = []
    seen = set()
    for v in G.nodes():
        if v in seen:
            continue
        img = frozenset(phi[v] for phi in autos)
        partition.append(sorted(img))
        seen |= img
    return partition


def quotient(G, partition):
    """Collapse each orbit block to one node -> irreducible-core network."""
    blocks = [frozenset(b) for b in partition]
    Q = nx.quotient_graph(G, blocks, relabel=False, create_using=nx.DiGraph)
    # give the blocks readable labels
    label = {}
    for b in blocks:
        b = frozenset(b)
        if len(b) == 1:
            label[b] = next(iter(b))
        else:
            # e.g. {HRAS,KRAS,NRAS} -> "RAS*"
            common = _common_stem(b)
            label[b] = (common + "*") if common else "|".join(sorted(b))
    return nx.relabel_nodes(Q, label)


def _common_stem(names):
    names = list(names)
    # crude common-suffix detector for paralogues like *RAS
    for k in range(min(len(n) for n in names), 0, -1):
        suf = names[0][-k:]
        if all(n.endswith(suf) for n in names):
            return suf
    return ""


def analyse():
    G = build_mapk_network()
    autos = automorphisms(G)
    part = orbits(G, autos)
    Q = quotient(G, part)
    nontrivial = [b for b in part if len(b) > 1]
    return {
        "G": G,
        "autos": autos,
        "group_order": len(autos),
        "partition": part,
        "nontrivial_orbits": nontrivial,
        "n_nodes": G.number_of_nodes(),
        "n_core_nodes": Q.number_of_nodes(),
        "quotient": Q,
    }


def make_figure(r, path="figures/m1_symmetry.png"):
    G, Q = r["G"], r["quotient"]
    orbit = set(r["nontrivial_orbits"][0]) if r["nontrivial_orbits"] else set()
    pos_full = {
        "Ligand": (0, 4), "RTK": (0, 3), "SOS": (0, 2),
        "HRAS": (-1, 1), "KRAS": (0, 1), "NRAS": (1, 1),
        "RAF": (0, 0), "MEK": (0, -1), "ERK": (0, -2)}
    pos_q = {"Ligand": (0, 4), "RTK": (0, 3), "SOS": (0, 2),
             "RAS*": (0, 1), "RAF": (0, 0), "MEK": (0, -1), "ERK": (0, -2)}
    fig, ax = plt.subplots(1, 2, figsize=(9, 6))
    cols = ["#e74c3c" if n in orbit else "#5dade2" for n in G.nodes()]
    nx.draw(G, pos_full, ax=ax[0], with_labels=True, node_color=cols,
            node_size=1400, font_size=8, arrows=True,
            edgecolors="k", font_weight="bold")
    ax[0].set_title(f"Full network (9 nodes)\n"
                    f"red = automorphism orbit  |Aut|={r['group_order']}=S$_3$")
    qcols = ["#e67e22" if n == "RAS*" else "#5dade2" for n in Q.nodes()]
    nx.draw(Q, pos_q, ax=ax[1], with_labels=True, node_color=qcols,
            node_size=1400, font_size=8, arrows=True,
            edgecolors="k", font_weight="bold")
    ax[1].set_title("Quotient / irreducible core (7 nodes)\n"
                    "3 RAS paralogues collapse to one node")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def report():
    r = analyse()
    print("=" * 68)
    print("MODULE 1 — Automorphism symmetry reduction (RTK/RAS/RAF/MEK/ERK)")
    print("=" * 68)
    print(f"nodes in full network        : {r['n_nodes']}")
    print(f"|Aut(G)| (automorphism group): {r['group_order']}  "
          f"(= |S_3| = 6 -> permutations of the 3 RAS paralogues)")
    print("non-trivial orbits (symmetric / redundant nodes):")
    for b in r["nontrivial_orbits"]:
        print(f"    {{{', '.join(b)}}}  -> collapses to ONE core node")
    print(f"irreducible-core network size: {r['n_core_nodes']} nodes")
    print("core reaction chain:")
    for u, v in r["quotient"].edges():
        print(f"    {u:>6}  ->  {v}")
    p = make_figure(r)
    print(f"figure written to: {p}")
    return r


if __name__ == "__main__":
    report()
