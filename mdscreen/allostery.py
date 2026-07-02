"""Allostery analysis: dynamic cross-correlation + apo-vs-holo perturbation.

Two complementary signals for "does this ligand cause an allosteric effect?":

1. Dynamic cross-correlation matrix (DCCM) of Cα motions, turned into a
   residue interaction network whose betweenness centrality highlights
   residues that mediate long-range communication (allosteric relay points).

2. Apo-vs-holo comparison: per-residue RMSF and correlation changes induced by
   ligand binding. Residues far from the binding site whose dynamics shift
   substantially are candidate allosteric hotspots — the fingerprint of an
   allosteric (rather than purely orthosteric) modulator.
"""
from __future__ import annotations

import os
from typing import Optional

import numpy as np

from .utils import log, ensure_dir, write_json


def _ca_positions(topology_pdb: str, trajectory_dcd: str):
    import MDAnalysis as mda
    from MDAnalysis.analysis import align
    u = mda.Universe(topology_pdb, trajectory_dcd)
    ca = u.select_atoms("protein and name CA")
    if len(ca) == 0:
        return None, None
    # superpose to remove global translation/rotation before correlations
    ref = align.AverageStructure(u, u, select="protein and name CA",
                                 ref_frame=0).run().results.universe
    align.AlignTraj(u, ref, select="protein and name CA", in_memory=True).run()
    coords = np.array([ca.positions.copy() for _ in u.trajectory])  # (T, N, 3)
    return ca.resids, coords


def dccm(coords: np.ndarray) -> np.ndarray:
    """Normalised dynamic cross-correlation matrix from (T, N, 3) coords."""
    mean = coords.mean(axis=0)
    disp = coords - mean                          # (T, N, 3)
    T, N, _ = disp.shape
    flat = disp.reshape(T, N, 3)
    # covariance of displacement vectors between residues
    cov = np.einsum("tid,tjd->ij", flat, flat) / T          # <Δi·Δj>
    var = np.sqrt(np.einsum("tid,tid->i", flat, flat) / T)  # sqrt<Δi·Δi>
    denom = np.outer(var, var)
    denom[denom == 0] = 1e-9
    return cov / denom


def network_hotspots(resids, corr: np.ndarray, coords: np.ndarray,
                     contact_cutoff: float = 8.0, top: int = 15):
    """Build a Cα network (edges = spatially contacting, |corr|-weighted) and
    return residues ranked by betweenness centrality."""
    import networkx as nx
    mean = coords.mean(axis=0)
    N = len(resids)
    G = nx.Graph()
    G.add_nodes_from(range(N))
    dmat = np.linalg.norm(mean[:, None, :] - mean[None, :, :], axis=-1)
    for i in range(N):
        for j in range(i + 1, N):
            if dmat[i, j] <= contact_cutoff and abs(corr[i, j]) > 0.3:
                # edge length small when correlation strong (info flows easily)
                G.add_edge(i, j, weight=-np.log(abs(corr[i, j]) + 1e-6))
    if G.number_of_edges() == 0:
        return []
    bc = nx.betweenness_centrality(G, weight="weight")
    ranked = sorted(bc.items(), key=lambda kv: kv[1], reverse=True)[:top]
    return [{"resid": int(resids[i]), "betweenness": round(v, 4)} for i, v in ranked]


def analyze_single(topology_pdb: str, trajectory_dcd: str, outdir: str,
                   tag: str = "holo") -> dict:
    ensure_dir(outdir)
    resids, coords = _ca_positions(topology_pdb, trajectory_dcd)
    if coords is None:
        return {"error": "no Cα atoms found"}
    corr = dccm(coords)
    np.save(os.path.join(outdir, f"{tag}_dccm.npy"), corr)
    _plot_matrix(outdir, f"{tag}_dccm", corr)
    hotspots = network_hotspots(resids, corr, coords)
    out = {"tag": tag, "n_residues": len(resids),
           "communication_hotspots": hotspots}
    write_json(os.path.join(outdir, f"{tag}_allostery.json"), out)
    return out


def compare_apo_holo(apo_top: str, apo_dcd: str, holo_top: str, holo_dcd: str,
                     outdir: str, binding_site_resids: Optional[list] = None) -> dict:
    """Flag allosteric perturbation from apo->holo dynamics changes."""
    ensure_dir(outdir)
    apo_res, apo_xyz = _ca_positions(apo_top, apo_dcd)
    holo_res, holo_xyz = _ca_positions(holo_top, holo_dcd)
    if apo_xyz is None or holo_xyz is None:
        return {"error": "missing Cα coordinates"}

    # Match residues present in both (holo has extra ligand but same protein).
    common = sorted(set(apo_res.tolist()) & set(holo_res.tolist()))
    ai = {r: k for k, r in enumerate(apo_res.tolist())}
    hi = {r: k for k, r in enumerate(holo_res.tolist())}
    apo_idx = [ai[r] for r in common]
    holo_idx = [hi[r] for r in common]

    apo_rmsf = _rmsf(apo_xyz[:, apo_idx, :])
    holo_rmsf = _rmsf(holo_xyz[:, holo_idx, :])
    d_rmsf = holo_rmsf - apo_rmsf

    site = set(binding_site_resids or [])
    perturbed = []
    thr = np.mean(np.abs(d_rmsf)) + np.std(np.abs(d_rmsf))
    for r, dv in zip(common, d_rmsf):
        distal = r not in site
        if abs(dv) > thr and distal:
            perturbed.append({"resid": int(r), "d_rmsf_ang": round(float(dv), 3)})
    perturbed.sort(key=lambda x: -abs(x["d_rmsf_ang"]))

    # Global allosteric signal: rmsf change magnitude at distal residues.
    distal_mask = np.array([r not in site for r in common])
    allo_signal = float(np.mean(np.abs(d_rmsf[distal_mask]))) if distal_mask.any() else 0.0

    out = {
        "n_common_residues": len(common),
        "allosteric_signal_ang": round(allo_signal, 4),
        "allosteric_verdict": _verdict(allo_signal, len(perturbed)),
        "distal_perturbed_residues": perturbed[:25],
    }
    _save_delta(outdir, common, d_rmsf)
    write_json(os.path.join(outdir, "apo_holo_comparison.json"), out)
    log.info("Allosteric signal=%.3f Å, %d distal residues perturbed -> %s",
             allo_signal, len(perturbed), out["allosteric_verdict"])
    return out


# ---------------------------------------------------------------------------
def _rmsf(coords):
    mean = coords.mean(axis=0)
    return np.sqrt(((coords - mean) ** 2).sum(axis=-1).mean(axis=0))


def _verdict(signal, n_perturbed):
    if signal > 0.6 and n_perturbed >= 5:
        return "likely allosteric modulator (strong distal dynamic response)"
    if signal > 0.3 and n_perturbed >= 3:
        return "possible allosteric effect (moderate distal response)"
    return "no clear allosteric signal (effect appears local/orthosteric)"


def _plot_matrix(outdir, name, mat):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(5, 4.2))
        im = ax.imshow(mat, cmap="coolwarm", vmin=-1, vmax=1, origin="lower")
        fig.colorbar(im, ax=ax, label="cross-correlation")
        ax.set_xlabel("residue index")
        ax.set_ylabel("residue index")
        ax.set_title(name)
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, f"{name}.png"), dpi=120)
        plt.close(fig)
    except Exception as exc:  # pragma: no cover
        log.warning("DCCM plot skipped: %s", exc)


def _save_delta(outdir, resids, d_rmsf):
    import csv
    with open(os.path.join(outdir, "delta_rmsf.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["resid", "delta_rmsf_ang"])
        for r, d in zip(resids, d_rmsf):
            w.writerow([r, d])
