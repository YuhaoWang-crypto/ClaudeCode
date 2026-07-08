"""
Module 11 — Input-tree fibration & fiber representatives (Morone-Leifer-Makse).

Generalises M1/M5. Graph automorphism (M1) requires a GLOBAL symmetry; the
fibration symmetry of Morone, Leifer & Makse (PNAS 2020) only requires nodes
to have isomorphic INPUT TREES. Two nodes fall in the same 'fiber' if they
receive structurally equivalent upstream signals, so they synchronise and can
be read out by a single 'fiber representative' (e.g. pERK1/2 for MAPK1/MAPK3).

Algorithm: the minimal balanced colouring = coarsest equitable partition under
incoming edges. We seed colours by node ROLE (role-constrained fibration, as
in the uploaded report) and refine by the multiset of (edge-sign, in-neighbour
colour) until the partition stops splitting. Fibers = final colour classes.

We run it on the expanded MAPK/ERK paralogue graph and recover fiber
representatives (RAS, RAF, MEK, ERK, RSK, DUSP ...), reproducing the report's
input-fibration compression.
"""
import networkx as nx


def build_expanded_mapk():
    """Paralogue-level EGFR/RAS/RAF/MEK/ERK graph with ERK->DUSP feedback."""
    G = nx.DiGraph()
    roles = {}

    def add(nodes, role):
        for n in nodes:
            roles[n] = role

    ligands = ["EGF", "TGFA"]
    ras = ["HRAS", "KRAS", "NRAS"]
    raf = ["ARAF", "BRAF", "RAF1"]
    mek = ["MAP2K1", "MAP2K2"]
    erk = ["MAPK1", "MAPK3"]
    rsk = ["RPS6KA1", "RPS6KA2", "RPS6KA3"]
    sos = ["SOS1", "SOS2"]
    dusp = ["DUSP5", "DUSP6"]
    tf = ["ELK1", "MYC", "CREB1", "FOS", "EGR1", "CCND1"]
    add(ligands, "ligand"); add(["EGFR"], "RTK"); add(["GRB2"], "adaptor")
    add(sos, "GEF"); add(ras, "GTPase"); add(raf, "MAP3K"); add(mek, "MAP2K")
    add(erk, "MAPK"); add(rsk, "RSK"); add(dusp, "DUSP"); add(tf, "TF")

    def link(us, vs, sign=+1):
        for u in us:
            for v in vs:
                G.add_edge(u, v, sign=sign)

    link(ligands, ["EGFR"]); link(["EGFR"], ["GRB2"]); link(["GRB2"], sos)
    link(sos, ras); link(ras, raf); link(raf, mek); link(mek, erk)
    link(erk, rsk); link(erk, dusp); link(erk, tf)
    link(dusp, erk, sign=-1)                       # negative feedback
    for n in G.nodes():
        G.nodes[n]["role"] = roles.get(n, "other")
    return G


def fibration_partition(G, role_seed=True):
    """Coarsest equitable partition under incoming edges (minimal balanced
    colouring). Returns dict node->fiber_id."""
    if role_seed:
        colour = {n: G.nodes[n].get("role", "x") for n in G.nodes()}
    else:
        colour = {n: 0 for n in G.nodes()}

    def signature(n):
        ins = sorted((G[u][n]["sign"], colour[u]) for u in G.predecessors(n))
        return (colour[n], tuple(ins))

    while True:
        sig = {n: signature(n) for n in G.nodes()}
        # relabel signatures to compact colours
        uniq = {s: i for i, s in enumerate(sorted(set(sig.values()),
                                                  key=lambda x: str(x)))}
        new = {n: uniq[sig[n]] for n in G.nodes()}
        if all(_same_partition(colour, new, a, b)
               for a in G.nodes() for b in G.nodes()):
            break
        if _npart(colour) == _npart(new):
            break
        colour = new
    return colour


def _npart(colour):
    return len(set(colour.values()))


def _same_partition(c1, c2, a, b):
    return (c1[a] == c1[b]) == (c2[a] == c2[b])


def fibers(G, colour):
    f = {}
    for n, c in colour.items():
        f.setdefault(c, []).append(n)
    return [sorted(v) for v in f.values()]


def report():
    print("=" * 68)
    print("MODULE 11 — Input-tree fibration & fiber representatives")
    print("=" * 68)
    G = build_expanded_mapk()
    colour = fibration_partition(G, role_seed=True)
    fbs = fibers(G, colour)
    multi = [f for f in fbs if len(f) > 1]
    n0, n1 = G.number_of_nodes(), len(fbs)
    print(f"expanded MAPK graph: {n0} nodes, {G.number_of_edges()} edges")
    print(f"fibration -> {n1} fibers  (compression ratio "
          f"{n0/n1:.2f}x)")
    print(f"\nmulti-node fibers (paralogues collapsing to one representative):")
    for f in sorted(multi, key=lambda x: -len(x)):
        role = G.nodes[f[0]]["role"]
        rep = _representative(f)
        print(f"  [{role:7s}] {{{', '.join(f)}}}  ->  read out via {rep}")
    print(f"\ncontrast with M1 automorphism (global symmetry, found only the")
    print(f"strict RAS S_3 orbit). Fibration finds ALL input-equivalent")
    print(f"paralogue layers at once -> {len(multi)} multi-node fibers.")
    print("fiber representative = the single pan-assay readout for that layer,")
    print("e.g. MAPK1/MAPK3 -> pERK1/2, SOS1/SOS2 -> pan-SOS, etc.")
    return {"n_nodes": n0, "n_fibers": n1, "multi_fibers": multi}


def _representative(members):
    reps = {"MAPK1": "pERK1/2", "SOS1": "pan-SOS", "HRAS": "pan-RAS",
            "ARAF": "pan-RAF", "MAP2K1": "pMEK1/2", "RPS6KA1": "pan-RSK",
            "DUSP5": "pan-DUSP", "EGF": "pan-ligand", "ELK1": "ERK-TF panel"}
    for m in members:
        if m in reps:
            return reps[m]
    return members[0]


if __name__ == "__main__":
    report()
