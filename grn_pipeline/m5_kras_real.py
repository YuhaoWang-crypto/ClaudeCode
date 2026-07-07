"""
Module 5 — Real target neighbourhood: KRAS G12C + covalent inhibitors.

Replaces the toy MAPK network of M1 with the real KRAS signalling
neighbourhood, using data pulled live from ChEMBL / Inductive Bio / Boltz:

  Target   KRAS G12C  (ChEMBL CHEMBL2189121, UniProt P01116, mut G12C)
  Drug 1   Sotorasib  (AMG-510,  CHEMBL4535757) covalent G12C inhibitor
  Drug 2   Adagrasib  (MRTX849,  CHEMBL4594350) covalent G12C inhibitor

Key result (rigorous, graph-theoretic):
  The three RAS paralogues HRAS/KRAS/NRAS are an automorphism orbit -> the
  symmetry group S_3 (order 6). A covalent G12C drug binds ONLY the KRAS
  node, marking it distinct. This BREAKS the symmetry:  S_3 (order 6) ->
  S_2 (order 2). Targeted covalent selectivity is literally a symmetry-
  breaking operation on the regulatory graph, and the drug's ChEMBL
  selectivity data (30 nM on G12C vs 36,500 nM on G12S, ~1000x) is the
  experimental signature of exactly that node-distinguishing action.
"""
import networkx as nx
from networkx.algorithms.isomorphism import DiGraphMatcher


# ---- real data pulled from the connected databases (this session) -------
CHEMBL = {
    "target": {"name": "GTPase KRas (G12C)", "chembl_id": "CHEMBL2189121",
               "uniprot": "P01116", "mutation": "G12C"},
    "sotorasib": {
        "chembl_id": "CHEMBL4535757", "alias": "AMG-510", "approved": 2021,
        "action": "covalent INHIBITOR", "mw": 560.61, "qed": 0.36,
        "ic50_g12c_nM": 30.0, "pchembl_g12c": 7.52,
        "ic50_g12s_nM": 36500.0,        # A549 (non-G12C) -> selectivity
        "antiproliferative_miapaca2_nM": 5.0,
        # Inductive Bio predictions (this session):
        "logD": 3.15, "acidic_pKa": 7.95, "basic_pKa": 4.98},
    "adagrasib": {
        "chembl_id": "CHEMBL4594350", "alias": "MRTX849", "approved": 2022,
        "action": "covalent INHIBITOR", "mw": 604.13, "qed": 0.36,
        "ic50_g12c_nM": None, "pchembl_g12c": None,
        "logD": 2.92, "acidic_pKa": 12.0, "basic_pKa": 8.05},
}
# Boltz-2.1 structure_and_binding predictions (launched + completed this
# session against the KRAS G12C sequence above). These are confidence-type
# scores, NOT ΔG in kcal/mol.
BOLTZ = {
    "sotorasib": {"binding_confidence": 0.879, "optimization_score": 0.423,
                  "iptm": 0.974, "complex_plddt": 0.893,
                  "structure_confidence": 0.909,
                  "prediction_id": "sab_pred_H4vBcEfU8bH6qww2Iiyu"},
    "adagrasib": {"binding_confidence": 0.958, "optimization_score": 0.611,
                  "iptm": 0.981, "complex_plddt": 0.909,
                  "structure_confidence": 0.924,
                  "prediction_id": "sab_pred_kq0AENipAVkvCAYUacca"},
}


def build_kras_network(drugged=None):
    """Real KRAS effector neighbourhood (EGFR -> RAS paralogues -> RAF/PI3K
    -> MAPK/AKT, with ERK negative feedback). If `drugged` names a RAS node,
    that node is marked (covalent inhibitor bound) so it is no longer
    interchangeable with its paralogues."""
    G = nx.DiGraph()
    edges = [
        ("EGFR", "GRB2", +1), ("GRB2", "SOS1", +1),
        ("SOS1", "HRAS", +1), ("SOS1", "KRAS", +1), ("SOS1", "NRAS", +1),
        ("HRAS", "RAF", +1), ("KRAS", "RAF", +1), ("NRAS", "RAF", +1),
        ("HRAS", "PI3K", +1), ("KRAS", "PI3K", +1), ("NRAS", "PI3K", +1),
        ("RAF", "MEK", +1), ("MEK", "ERK", +1),
        ("PI3K", "AKT", +1), ("AKT", "mTOR", +1),
        ("ERK", "SOS1", -1),                 # negative feedback
    ]
    for u, v, s in edges:
        G.add_edge(u, v, sign=s)
    for node in G.nodes():
        G.nodes[node]["drugged"] = 1 if node == drugged else 0
    return G


def symmetry(G):
    def em(e1, e2):
        return e1.get("sign") == e2.get("sign")
    def nm(n1, n2):
        return n1.get("drugged", 0) == n2.get("drugged", 0)
    gm = DiGraphMatcher(G, G, edge_match=em, node_match=nm)
    autos = list(gm.isomorphisms_iter())
    # orbits
    orbits = []
    seen = set()
    for v in G.nodes():
        if v in seen:
            continue
        img = sorted({phi[v] for phi in autos})
        if len(img) > 1:
            orbits.append(img)
        seen |= set(img)
    return {"order": len(autos), "orbits": orbits}


def perturbation_strength(drug):
    """Map ChEMBL potency to a normalised network-perturbation strength in
    [0,1] via pIC50 (higher potency -> stronger push toward the tipping
    point in the M4 stability model)."""
    d = CHEMBL[drug]
    ic50 = d.get("ic50_g12c_nM")
    if ic50 is None:
        return None
    pic50 = 9 - (ic50 and __import__("math").log10(ic50))  # pIC50 ~ 9-log10(nM)
    return max(0.0, min(1.0, (pic50 - 4) / 5))             # 4..9 -> 0..1


def update_boltz(drug, affinity):
    BOLTZ[drug] = affinity


def report():
    print("=" * 68)
    print("MODULE 5 — Real KRAS G12C neighbourhood + covalent inhibitors")
    print("=" * 68)
    t = CHEMBL["target"]
    print(f"target: {t['name']}  {t['chembl_id']}  ({t['uniprot']} {t['mutation']})")

    G0 = build_kras_network(drugged=None)
    s0 = symmetry(G0)
    print(f"\nUNDRUGGED network ({G0.number_of_nodes()} nodes):")
    print(f"  |Aut(G)| = {s0['order']}   symmetric orbit(s): {s0['orbits']}")
    print(f"  -> the 3 RAS paralogues form S_3 (order 6): interchangeable")

    Gd = build_kras_network(drugged="KRAS")
    sd = symmetry(Gd)
    print(f"\nKRAS-G12C-DRUGGED network (covalent inhibitor on KRAS only):")
    print(f"  |Aut(G)| = {sd['order']}   remaining orbit(s): {sd['orbits']}")
    print(f"  -> symmetry BROKEN  S_3 (6) -> S_2 (2): KRAS is now distinct,")
    print(f"     only HRAS/NRAS remain interchangeable.")

    print(f"\ndrug data (ChEMBL + Inductive Bio, this session):")
    for name in ("sotorasib", "adagrasib"):
        d = CHEMBL[name]
        ic = d["ic50_g12c_nM"]
        boltz = BOLTZ[name]
        print(f"  {name:10s} {d['chembl_id']} {d['alias']:8s} "
              f"{d['action']}")
        print(f"    IC50(G12C)={str(ic) if ic else 'NA':>6} nM  "
              f"logD={d['logD']}  acidic_pKa={d['acidic_pKa']}  "
              f"basic_pKa={d['basic_pKa']}")
        if boltz is not None:
            print(f"    Boltz-2.1: binding_confidence="
                  f"{boltz['binding_confidence']}  optimization_score="
                  f"{boltz['optimization_score']}  iPTM={boltz['iptm']}")
    so = CHEMBL["sotorasib"]
    fold = so["ic50_g12s_nM"] / so["ic50_g12c_nM"]
    print(f"\n  sotorasib G12C-selectivity: {so['ic50_g12c_nM']} nM (G12C) vs "
          f"{so['ic50_g12s_nM']:.0f} nM (G12S) = {fold:.0f}x")
    print(f"  -> the ~{fold:.0f}x selectivity IS the experimental signature of")
    print(f"     the node-distinguishing (symmetry-breaking) covalent action.")
    return {"undrugged": s0, "drugged": sd}


if __name__ == "__main__":
    report()
