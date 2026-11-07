"""
Shrake-Rupley SASA and the RASA / DeltaRASA glycan-burial metric.

This is a faithful re-implementation of the glycan-shielding score used in
Lunn-Halbert et al., "Programming N-linked glycan composition through de novo
protein design" (bioRxiv 2026.10.07.757460).

Definition, from the paper's Methods and Extended Data Fig. 7 legend
---------------------------------------------------------------------
RASA = relative accessible surface area, computed per atom with a 1.4 A probe
and normalised by the theoretical maximum accessibility of an *isolated* atom,

    SASA_max(i) = 4 * pi * (r_i + r_probe)^2

It is computed twice for every structure:

  pass 1 ("context")    glycan atoms WITH the protein present  -> total RASA
  pass 2 ("isolated")   the same glycan with PROTEIN ATOMS REMOVED, which
                        captures self-occlusion by neighbouring glycan atoms

    DeltaRASA(i) = | RASA_isolated(i) - RASA_context(i) |

Higher DeltaRASA = more occlusion caused by protein structure. The subtraction
is the whole point: it removes burial that the glycan inflicts on itself, so
what is left is protein-dependent shielding only.

Important ceiling (stated explicitly in the ED Fig. 7 legend): because bonded
glycan neighbours already occlude part of the isolated-atom reference,
DeltaRASA *cannot* approach 100% even for a completely buried atom. The paper's
per-atom colour bar tops out at 50%. Do not interpret DeltaRASA as "fraction
buried"; it is a relative shielding score.

Why a local implementation instead of Bio.PDB.SASA.ShrakeRupley
---------------------------------------------------------------
The two-pass metric needs independent control of *query* atoms (the ones whose
SASA is reported) and *occluder* atoms (the ones that block them). Biopython's
``ShrakeRupley.compute`` only takes a single entity, so expressing pass 1 vs
pass 2 would mean building throwaway Structure objects. We instead implement
the algorithm directly (same golden-spiral sphere, same probe, and Biopython's
own ``ATOMIC_RADII`` table) and verify numerical agreement with Biopython in
``validate_against_biopython``, which ``g1_delta_rasa`` runs as a self-check.

Dependencies: numpy + biopython only (no scipy).
"""

from __future__ import annotations

import math

import numpy as np
from Bio.PDB.SASA import ATOMIC_RADII

__all__ = [
    "PROBE_RADIUS",
    "N_POINTS",
    "golden_sphere",
    "radii_for_elements",
    "max_sasa",
    "sasa",
    "rasa",
    "delta_rasa",
    "validate_against_biopython",
]

#: Probe radius in Angstrom. The paper specifies 1.4 A (water).
PROBE_RADIUS = 1.40

#: Sphere points per atom. Biopython's default is 100; the paper does not say.
#: 960 is used here because the per-atom DeltaRASA is a small difference of two
#: similar numbers, so quadrature noise matters. See ``g1`` for the convergence
#: check that justifies this.
N_POINTS = 960

#: Fallback radius for elements missing from Biopython's table.
DEFAULT_RADIUS = 1.80

# Query atoms are processed in blocks so the point-vs-neighbour distance
# tensor stays small regardless of structure size.
_BLOCK = 128


def golden_sphere(n_points: int = N_POINTS) -> np.ndarray:
    """``n_points`` roughly-even points on the unit sphere (golden spiral).

    Same construction as Biopython's ``ShrakeRupley._compute_sphere`` so the
    two implementations quadrature-match rather than merely agree on average.
    """
    n = int(n_points)
    if n < 1:
        raise ValueError("n_points must be >= 1")
    dl = math.pi * (3.0 - math.sqrt(5.0))
    dz = 2.0 / n
    k = np.arange(n, dtype=np.float64)
    z = 1.0 - dz / 2.0 - k * dz
    r = np.sqrt(np.clip(1.0 - z * z, 0.0, None))
    lon = k * dl
    return np.column_stack((np.cos(lon) * r, np.sin(lon) * r, z))


def radii_for_elements(elements) -> np.ndarray:
    """Van-der-Waals radii for a sequence of element symbols.

    Uses Biopython's ``ATOMIC_RADII`` so radii match the tool the paper names.
    """
    out = np.empty(len(elements), dtype=np.float64)
    for i, el in enumerate(elements):
        key = (el or "").strip().capitalize()
        out[i] = ATOMIC_RADII.get(key, DEFAULT_RADIUS)
    return out


def max_sasa(radii: np.ndarray, probe: float = PROBE_RADIUS) -> np.ndarray:
    """Theoretical maximum accessibility of an isolated atom, 4*pi*(r+rp)^2."""
    return 4.0 * np.pi * (np.asarray(radii, dtype=np.float64) + probe) ** 2


def sasa(
    query_xyz: np.ndarray,
    query_radii: np.ndarray,
    occluder_xyz: np.ndarray | None = None,
    occluder_radii: np.ndarray | None = None,
    probe: float = PROBE_RADIUS,
    n_points: int = N_POINTS,
) -> np.ndarray:
    """Per-atom SASA (A^2) of the query atoms.

    Query atoms are occluded by *each other* and by any ``occluder_*`` atoms.
    Passing no occluders gives the "isolated" pass of the DeltaRASA metric.
    """
    q_xyz = np.asarray(query_xyz, dtype=np.float64).reshape(-1, 3)
    q_r = np.asarray(query_radii, dtype=np.float64).reshape(-1)
    if q_xyz.shape[0] != q_r.shape[0]:
        raise ValueError("query_xyz and query_radii length mismatch")
    if q_xyz.shape[0] == 0:
        return np.zeros(0)

    if occluder_xyz is None or len(occluder_xyz) == 0:
        all_xyz, all_r = q_xyz, q_r
    else:
        o_xyz = np.asarray(occluder_xyz, dtype=np.float64).reshape(-1, 3)
        o_r = np.asarray(occluder_radii, dtype=np.float64).reshape(-1)
        if o_xyz.shape[0] != o_r.shape[0]:
            raise ValueError("occluder_xyz and occluder_radii length mismatch")
        all_xyz = np.vstack((q_xyz, o_xyz))
        all_r = np.concatenate((q_r, o_r))

    sphere = golden_sphere(n_points)
    ext = all_r + probe          # extended radii of every potential occluder
    q_ext = q_r + probe
    n_q = q_xyz.shape[0]
    out = np.empty(n_q, dtype=np.float64)

    for start in range(0, n_q, _BLOCK):
        stop = min(start + _BLOCK, n_q)
        blk = np.arange(start, stop)
        # Distances from this block of query atoms to every atom in the system.
        d = np.linalg.norm(q_xyz[blk, None, :] - all_xyz[None, :, :], axis=2)
        # An atom can occlude atom i only if their extended spheres overlap.
        touching = d < (q_ext[blk, None] + ext[None, :])
        touching[np.arange(len(blk)), blk] = False       # never self-occlude

        for local, i in enumerate(blk):
            nb = np.flatnonzero(touching[local])
            if nb.size == 0:
                out[i] = 4.0 * np.pi * q_ext[i] ** 2
                continue
            pts = q_xyz[i] + sphere * q_ext[i]           # (n_points, 3)
            dd = np.linalg.norm(pts[:, None, :] - all_xyz[nb][None, :, :], axis=2)
            buried = (dd < ext[nb][None, :]).any(axis=1)
            frac = 1.0 - buried.mean()
            out[i] = 4.0 * np.pi * q_ext[i] ** 2 * frac

    return out


def rasa(
    query_xyz: np.ndarray,
    query_radii: np.ndarray,
    occluder_xyz: np.ndarray | None = None,
    occluder_radii: np.ndarray | None = None,
    probe: float = PROBE_RADIUS,
    n_points: int = N_POINTS,
) -> np.ndarray:
    """Per-atom RASA: SASA normalised by the isolated-atom maximum."""
    s = sasa(query_xyz, query_radii, occluder_xyz, occluder_radii, probe, n_points)
    return s / max_sasa(query_radii, probe)


def delta_rasa(
    glycan_xyz: np.ndarray,
    glycan_radii: np.ndarray,
    protein_xyz: np.ndarray,
    protein_radii: np.ndarray,
    probe: float = PROBE_RADIUS,
    n_points: int = N_POINTS,
    as_percent: bool = True,
) -> dict:
    """Two-pass DeltaRASA for one glycan against a protein context.

    ``protein_xyz`` is whatever should count as "context". Pass protein atoms
    only to reproduce the paper exactly; pass protein + *other* glycan trees to
    additionally score shielding by neighbouring glycans (see ``g2``).

    Returns a dict with per-atom arrays and the two summary statistics the
    paper ranks designs on: mean per-atom DeltaRASA, and the number of glycan
    atoms with non-zero DeltaRASA ("contact atoms").
    """
    ctx = rasa(glycan_xyz, glycan_radii, protein_xyz, protein_radii, probe, n_points)
    iso = rasa(glycan_xyz, glycan_radii, None, None, probe, n_points)
    d = np.abs(iso - ctx)
    scale = 100.0 if as_percent else 1.0
    # Quadrature noise floor: one sphere point out of n_points.
    tol = 1.0 / n_points
    return {
        "rasa_context": ctx * scale,
        "rasa_isolated": iso * scale,
        "delta_rasa": d * scale,
        "mean_delta_rasa": float(d.mean() * scale),
        "max_delta_rasa": float(d.max() * scale) if d.size else 0.0,
        "n_contact_atoms": int((d > tol).sum()),
        "n_atoms": int(d.size),
        "units": "percent" if as_percent else "fraction",
    }


def validate_against_biopython(n_points: int = 100, tol: float = 1e-6) -> dict:
    """Check this implementation against ``Bio.PDB.SASA.ShrakeRupley``.

    Builds a small synthetic structure, runs both codes with the same probe and
    sphere resolution, and returns the maximum absolute per-atom deviation.
    Called by ``g1_delta_rasa`` so the pipeline proves its own SASA engine
    rather than asserting it.
    """
    from Bio.PDB.Atom import Atom
    from Bio.PDB.Chain import Chain
    from Bio.PDB.Model import Model
    from Bio.PDB.Residue import Residue
    from Bio.PDB.SASA import ShrakeRupley
    from Bio.PDB.Structure import Structure

    rng = np.random.default_rng(0)
    elements = ["C", "N", "O", "C", "O", "C", "N", "C", "O", "C", "C", "O"]
    xyz = rng.normal(scale=2.6, size=(len(elements), 3))

    struct = Structure("t")
    model = Model(0)
    chain = Chain("A")
    res = Residue((" ", 1, " "), "TST", "")
    for i, (el, c) in enumerate(zip(elements, xyz)):
        res.add(Atom(f"{el}{i}", c, 0.0, 1.0, " ", f"{el}{i}", i + 1, element=el))
    chain.add(res)
    model.add(chain)
    struct.add(model)

    ShrakeRupley(probe_radius=PROBE_RADIUS, n_points=n_points).compute(struct, level="A")
    ref = np.array([a.sasa for a in struct.get_atoms()])

    mine = sasa(xyz, radii_for_elements(elements), n_points=n_points)

    max_abs = float(np.max(np.abs(ref - mine)))
    rel = float(max_abs / max(ref.max(), 1e-12))
    return {
        "n_atoms": len(elements),
        "n_points": n_points,
        "biopython_total_A2": float(ref.sum()),
        "local_total_A2": float(mine.sum()),
        "max_abs_atom_diff_A2": max_abs,
        "max_rel_atom_diff": rel,
        "agrees": bool(max_abs < tol),
    }
