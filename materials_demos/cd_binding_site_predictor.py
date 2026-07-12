"""
Demo #5 — Geometric prediction of Cd(2+) binding sites in a protein.

Motivation (from the chat context)
-----------------------------------
Client (bizlikery): has a protein, wants to know WHICH sites bind divalent
cadmium (Cd2+), in order to design point mutants and measure binding by MST.
They mentioned AutoDock / Discovery Studio. Before running docking, a fast and
physically-motivated step is to predict candidate metal sites directly from the
structure and rank residues to mutate.

WHAT THIS IS
------------
A runnable predictor that:
    1. downloads a protein structure (PDB),
    2. collects candidate metal-donor atoms,
    3. clusters donors that could jointly coordinate one Cd2+ ion,
    4. scores each cluster with a Cd2+-specific (soft-acid, HSAB) weighting,
    5. ranks candidate sites and lists the residues to mutate for MST.

Cd2+ is a soft/borderline Lewis acid, so it prefers soft donors:
    thiolate/thioether S (Cys, Met) > imidazole N (His) > carboxylate O
    (Asp, Glu) > amide/hydroxyl O (backbone, Ser/Thr/Asn/Gln).

VALIDATION
----------
Run on carbonic anhydrase II (1CA2): the predictor should surface the known
catalytic metal site (His94, His96, His119) as a top-ranked cluster, which is a
sanity check that the geometry/scoring is sound.

WHAT THIS IS *NOT*
------------------
This is geometric/heuristic prediction, NOT docking or QM. A production job
would follow up with AutoDock (with a Cd2+ parameterisation), MetalPDB/CHED
comparison, and ideally QM/MM refinement of the coordination sphere. Use this
to shortlist sites and mutation targets cheaply, then dock/refine the top hits.
"""

import os
import sys
import urllib.request
import numpy as np
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))

# Cd2+ donor-atom dictionary: (resname, atomname) -> (donor_type, weight)
# weights follow HSAB softness preference for Cd2+ (soft S > N > O)
DONORS = {
    ("CYS", "SG"): ("S(thiolate)", 3.0),
    ("MET", "SD"): ("S(thioether)", 2.5),
    ("HIS", "ND1"): ("N(imidazole)", 2.0),
    ("HIS", "NE2"): ("N(imidazole)", 2.0),
    ("ASP", "OD1"): ("O(carboxylate)", 1.5),
    ("ASP", "OD2"): ("O(carboxylate)", 1.5),
    ("GLU", "OE1"): ("O(carboxylate)", 1.5),
    ("GLU", "OE2"): ("O(carboxylate)", 1.5),
    ("ASN", "OD1"): ("O(amide)", 1.0),
    ("GLN", "OE1"): ("O(amide)", 1.0),
    ("SER", "OG"): ("O(hydroxyl)", 0.8),
    ("THR", "OG1"): ("O(hydroxyl)", 0.8),
    ("TYR", "OH"): ("O(phenol)", 0.8),
}

# two ligands coordinating the SAME metal (M-L ~2.4 A, L-M-L ~90-180 deg)
# sit roughly 3.0-4.8 A apart. Use a tolerant window for donor-donor linking.
DD_MIN, DD_MAX = 2.6, 5.0


def fetch_pdb(pdb_id):
    cache = os.path.join(HERE, f"{pdb_id}.pdb")
    if not os.path.exists(cache):
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        urllib.request.urlretrieve(url, cache)
    return cache


def parse_donors(pdb_path):
    """Return list of donor atoms: dict(resname,resid,chain,atom,type,weight,xyz)."""
    donors = []
    seen = set()
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            altloc = line[16]
            if altloc not in (" ", "A"):
                continue
            resname = line[17:20].strip()
            atom = line[12:16].strip()
            key = (resname, atom)
            if key not in DONORS:
                continue
            chain = line[21]
            resid = line[22:26].strip()
            uid = (chain, resid, atom)
            if uid in seen:
                continue
            seen.add(uid)
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            dtype, w = DONORS[key]
            donors.append(dict(resname=resname, resid=resid, chain=chain,
                               atom=atom, type=dtype, weight=w,
                               xyz=np.array([x, y, z])))
    return donors


def _bron_kerbosch(R, P, X, adj, cliques):
    """Maximal-clique enumeration; a metal site needs donors MUTUALLY close
    (all coordinating one central ion), which a clique captures and simple
    connected-components (chaining) does not."""
    if not P and not X:
        if len(R) >= 2:
            cliques.append(set(R))
        return
    pivot = max(P | X, key=lambda u: len(adj[u] & P)) if (P | X) else None
    for v in list(P - (adj[pivot] if pivot is not None else set())):
        _bron_kerbosch(R | {v}, P & adj[v], X & adj[v], adj, cliques)
        P = P - {v}
        X = X | {v}


def cluster_sites(donors):
    """Find sets of donor atoms that could jointly coordinate one Cd2+ (cliques
    in the donor-donor coordination graph), then score them."""
    n = len(donors)
    xyz = np.array([d["xyz"] for d in donors])
    adj = {i: set() for i in range(n)}
    for i in range(n):
        r = np.sqrt(((xyz - xyz[i]) ** 2).sum(1))
        for j in np.where((r >= DD_MIN) & (r <= DD_MAX))[0]:
            if j != i:
                adj[i].add(int(j))

    cliques = []
    _bron_kerbosch(set(), set(range(n)), set(), adj, cliques)

    sites, seen = [], set()
    for members in cliques:
        members = tuple(sorted(members))
        residues = {(donors[m]["chain"], donors[m]["resid"]) for m in members}
        if len(residues) < 2:                    # need >=2 distinct residues
            continue
        pts = xyz[list(members)]
        centroid = pts.mean(0)
        spread = np.sqrt(((pts - centroid) ** 2).sum(1).max())  # max radius
        if spread > 3.2:                         # too diffuse to be one metal site
            continue
        key = tuple(sorted(residues))
        if key in seen:
            continue
        seen.add(key)
        weight = sum(donors[m]["weight"] for m in members)
        # score: soft-donor weight, reward >=3 residues, penalise diffuse geometry
        score = weight * (1.0 + 0.3 * (len(residues) - 2)) / (1.0 + 0.4 * spread)
        sites.append(dict(members=list(members), residues=residues,
                          centroid=centroid, spread=spread,
                          weight=weight, score=score))
    return sorted(sites, key=lambda s: -s["score"])


def report(donors, sites, pdb_id, topn=6):
    print(f"\n=== Predicted Cd2+ binding sites for {pdb_id} "
          f"({len(donors)} candidate donor atoms) ===\n")
    lines = []
    for rank, s in enumerate(sites[:topn], 1):
        res_list = sorted({(donors[m]["resname"], donors[m]["chain"],
                            donors[m]["resid"]) for m in s["members"]},
                          key=lambda t: (t[1], int(t[2])))
        res_str = ", ".join(f"{r[0]}{r[2]}/{r[1]}" for r in res_list)
        donor_types = sorted({donors[m]["type"] for m in s["members"]})
        print(f"  #{rank}  score={s['score']:5.2f}  spread={s['spread']:.2f} Å")
        print(f"       residues : {res_str}")
        print(f"       donors   : {', '.join(donor_types)}")
        print(f"       -> MST mutation targets: "
              f"{', '.join(f'{r[0]}{r[2]}' for r in res_list)}\n")
        lines.append((rank, s, res_list))
    return lines


def main(pdb_id="1CA2"):
    path = fetch_pdb(pdb_id)
    donors = parse_donors(path)
    sites = cluster_sites(donors)
    ranked = report(donors, sites, pdb_id)

    # figure: score bars + donor-composition
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        top = sites[:6]
        labels = [f"#{i+1}" for i in range(len(top))]
        scores = [s["score"] for s in top]
        weights = [s["weight"] for s in top]
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
        ax[0].barh(labels[::-1], scores[::-1], color="#2c6fbb")
        ax[0].set_xlabel("site score (Cd²⁺ affinity proxy)")
        ax[0].set_title("Ranked candidate Cd²⁺ sites")
        ax[0].grid(alpha=0.3, axis="x")
        # donor-type composition of the top site
        best = sites[0]
        comp = defaultdict(float)
        for m in best["members"]:
            comp[donors[m]["type"]] += donors[m]["weight"]
        ax[1].bar(list(comp.keys()), list(comp.values()), color="#e08e0b")
        ax[1].set_ylabel("summed donor weight")
        ax[1].set_title("Top site donor composition")
        ax[1].tick_params(axis="x", rotation=25, labelsize=8)
        ax[1].grid(alpha=0.3, axis="y")
        fig.suptitle(f"Demo #5 · Cd²⁺ binding-site prediction on {pdb_id}\n"
                     "geometric HSAB predictor (shortlist before AutoDock/QM)",
                     fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, 0.9))
        out = os.path.join(HERE, "cd_binding_sites.png")
        fig.savefig(out, dpi=150)
        print(f"saved figure -> {out}")
    except Exception as exc:
        print(f"(plot skipped: {exc})")


if __name__ == "__main__":
    pdb = sys.argv[1] if len(sys.argv) > 1 else "1CA2"
    main(pdb)
