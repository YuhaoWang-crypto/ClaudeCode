"""
Structural anatomy of the two scaffolds named in the quote.

  7CBC  de novo designed switch caging a hemagglutinin binder (LucCage family,
        Quijano-Rubio 2021) -- the "LucCage (LOCKR) architecture" option.
  1OMP  maltose-binding protein, apo (open)      -- the "Type I clamshell"
  1ANF  maltose-binding protein, maltose-bound (closed)   option.

Everything here runs on downloaded PDB files with numpy only; no licensed
structural-biology package is required.
"""

from __future__ import annotations

import urllib.request
from pathlib import Path

import numpy as np

RCSB = "https://files.rcsb.org/download/{}.pdb"


def fetch(pdb_id: str, directory: str | Path = ".") -> Path:
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{pdb_id.lower()}.pdb"
    if not path.exists():
        urllib.request.urlretrieve(RCSB.format(pdb_id.upper()), path)
    return path


def load_ca(path, chain="A", hetatm=False):
    """Alpha-carbon coordinates keyed by residue number."""
    ca = {}
    starts = ("ATOM", "HETATM") if hetatm else ("ATOM",)
    for line in open(path):
        if not line.startswith(starts):
            continue
        if line[12:16].strip() != "CA" or line[21] != chain:
            continue
        if line[16] not in (" ", "A"):
            continue
        ca[int(line[22:26])] = np.array(
            [float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return ca


def helix_records(path, chain="A"):
    """(start, end) of every HELIX record for a chain, from the PDB header."""
    out = []
    for line in open(path):
        if line.startswith("HELIX") and line[19] == chain:
            out.append((int(line[21:25]), int(line[33:37])))
    return out


def kabsch(p, q):
    """Optimal superposition of p onto q. Returns rotation, centroids, RMSD."""
    pc, qc = p.mean(0), q.mean(0)
    h = (p - pc).T @ (q - qc)
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    r = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    rms = float(np.sqrt((((r @ (p - pc).T).T + qc - q) ** 2).sum(1).mean()))
    return r, pc, qc, rms


def rotation_angle_deg(r):
    return float(np.degrees(np.arccos(np.clip((np.trace(r) - 1.0) / 2.0, -1.0, 1.0))))


# --------------------------------------------------------------------------
# 1. LOCKR / LucCage cage-latch topology
# --------------------------------------------------------------------------

LUCCAGE_SEGMENTS = [
    ("H1", 2, 47), ("H2", 49, 94), ("H3", 98, 143), ("H4", 145, 191),
    ("H5", 194, 240), ("H6", 244, 269),
    ("b1", 273, 285), ("b2", 288, 299), ("b3", 302, 316),
]


def contact_matrix(ca, segments, cutoff=10.0):
    """Count of CA pairs closer than `cutoff` between every pair of segments."""
    n = len(segments)
    m = np.zeros((n, n), dtype=int)
    coords = [np.array([ca[r] for (_, a, b) in [seg] for r in range(a, b + 1)
                        if r in ca]) for seg in segments]
    for i in range(n):
        for j in range(i + 1, n):
            if len(coords[i]) == 0 or len(coords[j]) == 0:
                continue
            d = np.linalg.norm(coords[i][:, None, :] - coords[j][None, :, :], axis=-1)
            m[i, j] = m[j, i] = int((d < cutoff).sum())
    return m


def luccage_topology(pdb_path, chain="A"):
    # hetatm=True picks up selenomethionine (MSE), which is a real residue of the
    # chain recorded as HETATM; solvent and ethanol have no CA atom so are ignored
    ca = load_ca(pdb_path, chain, hetatm=True)
    m = contact_matrix(ca, LUCCAGE_SEGMENTS)
    names = [s[0] for s in LUCCAGE_SEGMENTS]
    lengths = {}
    for name, a, b in LUCCAGE_SEGMENTS:
        pts = np.array([ca[r] for r in range(a, b + 1) if r in ca])
        lengths[name] = float(np.linalg.norm(pts[0] - pts[-1])) if len(pts) > 1 else 0.0
    return {"names": names, "contacts": m, "end_to_end_ang": lengths,
            "n_residues": len(ca)}


# --------------------------------------------------------------------------
# 2. Type I periplasmic-binding-protein clamshell hinge motion
# --------------------------------------------------------------------------

# Sharff/Quiocho domain definition for E. coli MBP
MBP_N_DOMAIN = [(1, 109), (264, 309)]
MBP_C_DOMAIN = [(114, 258), (316, 370)]


def _expand(ranges, available):
    out = []
    for a, b in ranges:
        out.extend(r for r in range(a, b + 1) if r in available)
    return out


def clamshell_hinge(open_pdb, closed_pdb, chain="A",
                    n_domain=MBP_N_DOMAIN, c_domain=MBP_C_DOMAIN):
    """Rigid-body hinge angle between the two lobes of a Type I clamshell.

    Superposes each domain independently; the relative rotation between the two
    superpositions is the hinge closure angle.  Low internal RMSD for each
    domain is the evidence that the motion really is rigid-body hinge bending
    (which is what makes the scaffold usable as an allosteric input).
    """
    o = load_ca(open_pdb, chain)
    c = load_ca(closed_pdb, chain)
    common = sorted(set(o) & set(c))
    d1 = _expand(n_domain, set(common))
    d2 = _expand(c_domain, set(common))

    p_all = np.array([o[r] for r in common])
    q_all = np.array([c[r] for r in common])
    _, _, _, rms_global = kabsch(p_all, q_all)

    r1, pc1, qc1, rms1 = kabsch(np.array([o[r] for r in d1]),
                                np.array([c[r] for r in d1]))
    r2, _, _, rms2 = kabsch(np.array([o[r] for r in d2]),
                             np.array([c[r] for r in d2]))

    aligned = {r: (r1 @ (o[r] - pc1)) + qc1 for r in common}
    disp = np.array([np.linalg.norm(aligned[r] - c[r]) for r in d2])

    return {
        "n_common": len(common),
        "rmsd_global_ang": rms_global,
        "rmsd_n_domain_ang": rms1,
        "rmsd_c_domain_ang": rms2,
        "hinge_angle_deg": rotation_angle_deg(r2 @ r1.T),
        "c_domain_mean_displacement_ang": float(disp.mean()),
        "c_domain_max_displacement_ang": float(disp.max()),
    }
