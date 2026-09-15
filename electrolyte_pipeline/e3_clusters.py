"""
E3 — SSIP / CIP / AGG and the coordination network   [digest p6]

"配位状态的命名要落到明确的距离或连接判据."  Criteria used here, all reported
in the JSON next to the numbers:

  contact(cation, anion) :=  min over O(anion) of |r_cat - r_O|  <  r_cut
  r_cut                  :=  first minimum of g_cat-O(anion)(r) from E2
                             (falls back to 3.0 Å if E2 has not been run)

  Graph: nodes = cations + anions, edges = contacts.  Per frame, each cation
  is classified from the connected component it belongs to:
      no-contact  (free + SSIP; the two are not separated here — a second-shell
                   criterion would be needed, and is not defined in this module)
      CIP         component = exactly {1 cation, 1 anion}
      AGG         component has ≥ 3 ions
  Reported: state fractions with block errors, cluster-size distribution,
  anion participation (fraction of anions with ≥1 cation), cations per bound
  anion (how many Li+ one TFSI- bridges — "跨壳层连接").
"""
from __future__ import annotations

import json
import os

import numpy as np

from electrolyte_pipeline.e0_systems import WORK, FIG, SERIES
from electrolyte_pipeline.traj import Traj


def analyse(label: str, workdir: str = WORK, r_cut: float | None = None, stride: int = 1,
            nblocks: int = 5) -> dict:
    import networkx as nx
    from MDAnalysis.lib.distances import distance_array
    tr = Traj(label, workdir)
    if "TFSI" not in tr.counts():
        return {"label": label, "skipped": "no salt"}
    cat = tr.sel(species=tr.cation)
    o_an = tr.sel(species="TFSI", element="O")
    an_of_o = tr.molid[o_an]
    anions = np.unique(an_of_o)
    if r_cut is None:
        e2 = os.path.join(tr.dir, "e2_rdf.json")
        r_cut = json.load(open(e2))["pairs"]["O(TFSI)"]["r_min"] if os.path.exists(e2) else 3.0
        r_cut_source = "E2 first minimum of g_cat-O(TFSI)" if os.path.exists(e2) else "default 3.0 Å"
    else:
        r_cut_source = "user"

    states = []            # per frame: (n_free, n_cip, n_agg)
    sizes = {}             # cluster size -> count (ions, over all frames)
    anion_bound = []       # fraction of anions with >=1 cation
    cat_per_anion = []     # mean cations per bound anion
    nfr = 0
    for ts in tr.frames(stride):
        d = distance_array(tr.u.atoms.positions[cat], tr.u.atoms.positions[o_an], box=ts.dimensions)
        G = nx.Graph()
        G.add_nodes_from([("c", i) for i in range(len(cat))])
        G.add_nodes_from([("a", int(a)) for a in anions])
        ci, oi = np.where(d < r_cut)
        for c, o in zip(ci, oi):
            G.add_edge(("c", int(c)), ("a", int(an_of_o[o])))
        n_free = n_cip = n_agg = 0
        for comp in nx.connected_components(G):
            ncat = sum(1 for n in comp if n[0] == "c")
            nan = len(comp) - ncat
            if ncat == 0:
                continue
            if len(comp) == 1:
                n_free += 1
            elif len(comp) == 2:
                n_cip += 1
            else:
                n_agg += ncat
            if len(comp) >= 2:
                sizes[len(comp)] = sizes.get(len(comp), 0) + len(comp)
        states.append((n_free, n_cip, n_agg))
        deg = np.array([G.degree(("a", int(a))) for a in anions])
        anion_bound.append((deg > 0).mean())
        cat_per_anion.append(deg[deg > 0].mean() if (deg > 0).any() else 0.0)
        nfr += 1
    S = np.array(states, float) / len(cat)
    blocks = np.array_split(np.arange(nfr), nblocks)
    bm = np.array([S[b].mean(axis=0) for b in blocks])
    frac = S.mean(axis=0)
    sem = bm.std(axis=0, ddof=1) / np.sqrt(nblocks)
    total_ion_frames = sum(sizes.values()) + S[:, 0].sum() * len(cat)
    size_dist = {int(k): v / total_ion_frames for k, v in sorted(sizes.items())}
    res = {"label": label, "cation": tr.cation, "r_cut_A": float(r_cut), "r_cut_source": r_cut_source,
           "criteria": "contact = any O(TFSI) within r_cut of the cation; CIP = 2-ion component; "
                       "AGG = component with >=3 ions; free/SSIP not separated",
           "n_frames": nfr, "n_cations": int(len(cat)), "n_anions": int(len(anions)),
           "frac_no_contact": float(frac[0]), "frac_CIP": float(frac[1]), "frac_AGG": float(frac[2]),
           "sem": sem.round(4).tolist(), "block_means": bm.round(3).tolist(),
           "anion_participation": float(np.mean(anion_bound)),
           "cations_per_bound_anion": float(np.mean(cat_per_anion)),
           "cluster_size_dist_ionfrac": size_dist,
           "largest_cluster_ions": int(max(sizes) if sizes else 1)}
    with open(os.path.join(tr.dir, "e3_clusters.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    return res


def report(labels=None, workdir: str = WORK) -> dict:
    labels = labels or [c.label for c in SERIES if os.path.exists(os.path.join(workdir, c.label, "prod.dcd"))]
    out = {}
    print("E3  ion-pairing states (criteria in e3_clusters.json; free/SSIP not separated)")
    print(f"{'system':12s} {'r_cut':>5s} {'no-contact':>10s} {'CIP':>8s} {'AGG':>8s} {'anion-part':>10s} "
          f"{'cat/anion':>9s} {'largest':>7s}")
    for lab in labels:
        r = analyse(lab, workdir)
        if "skipped" in r:
            continue
        out[lab] = r
        print(f"{lab:12s} {r['r_cut_A']:5.2f} {r['frac_no_contact']:6.2f}±{r['sem'][0]:.2f} "
              f"{r['frac_CIP']:5.2f}±{r['sem'][1]:.2f} {r['frac_AGG']:5.2f}±{r['sem'][2]:.2f} "
              f"{r['anion_participation']:10.2f} {r['cations_per_bound_anion']:9.2f} {r['largest_cluster_ions']:7d}")
    _plot(out)
    return out


def _plot(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if not out:
        return
    from electrolyte_pipeline.e0_systems import SERIES_BY_LABEL
    labs = [l for l in out if l.startswith("Li")]
    x = [SERIES_BY_LABEL[l].counts["DME"] / SERIES_BY_LABEL[l].counts["Li"] if l in SERIES_BY_LABEL else 0 for l in labs]
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    for key, name, mk in (("frac_no_contact", "no contact (free+SSIP)", "s"), ("frac_CIP", "CIP", "o"), ("frac_AGG", "AGG", "^")):
        y = [out[l][key] for l in labs]
        e = [out[l]["sem"][["frac_no_contact", "frac_CIP", "frac_AGG"].index(key)] for l in labs]
        ax[0].errorbar(x, y, yerr=e, marker=mk, label=name, capsize=2)
    ax[0].set_xlabel("DME : Li ratio  (low ratio = high concentration)"); ax[0].set_ylabel("fraction of cations")
    ax[0].invert_xaxis(); ax[0].legend(fontsize=8); ax[0].set_ylim(0, 1)
    for l in out:
        sd = out[l]["cluster_size_dist_ionfrac"]
        ax[1].plot([int(k) for k in sd], list(sd.values()), "o-", ms=3, label=l)
    ax[1].set_xlabel("cluster size (ions)"); ax[1].set_ylabel("fraction of ions in clusters of this size")
    ax[1].set_yscale("log"); ax[1].legend(fontsize=7)
    fig.suptitle("E3  ion-pairing states vs concentration · cluster-size distribution", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "e3_clusters.png"), dpi=140); plt.close(fig)


if __name__ == "__main__":
    import sys
    report(sys.argv[1:] or None)
