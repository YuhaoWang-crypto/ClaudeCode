"""
Module 7 - turn the open state into something you can actually screen against.

Exports representative open and closed conformations and measures the pocket,
so the result is a docking target rather than a statistic.

HARD LIMITATION, STATED UP FRONT: the published dataset contains C-alpha
coordinates only (the authors exported C-alpha traces with JEDi).  You cannot
dock into a C-alpha trace.  What this module therefore delivers is
  (a) the pocket's residue list and C-alpha-level dimensions in the open state,
  (b) which experimental structure is the right docking template,
  (c) representative frames as PDB, for side-chain rebuilding or as restraint
      targets if you re-run the MD with all atoms (see run_gromacs.py).
"""
import os

import numpy as np

from .common import (FIG, amb_to_idx, apo_indices, banner, idx_to_amb,
                     crystal_ca, kabsch_rmsd, load)
from .m2_cryptic_pocket import (POCKET_RES, collective_variables,
                                opening_axis)

OUT = os.path.join(FIG, "m7_structures")
# the residues that line the 1PZO/1PZP non-competitive site
SITE_RES = [213, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226,
            227, 228, 229, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280]


def write_ca_pdb(path, coords, resnums, bfac=None):
    """Minimal C-alpha-only PDB (ALA placeholder residues)."""
    with open(path, "w") as fh:
        fh.write("REMARK  C-alpha trace from MSM pipeline; side chains absent\n")
        for n, (r, xyz) in enumerate(zip(resnums, coords), start=1):
            b = 0.0 if bfac is None else float(bfac[n - 1])
            fh.write(f"ATOM  {n:5d}  CA  ALA A{r:4d}    "
                     f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00{b:6.2f}           C\n")
        fh.write("END\n")


def main():
    banner("M7  the open state as a screening target")
    os.makedirs(OUT, exist_ok=True)
    X, names = load()
    proj, cv1, _ = collective_variables()
    _, _, L = opening_axis()
    resnums = [int(idx_to_amb(i)) for i in range(X.shape[2])]

    # Two different "best" frames, and they are NOT the same frame:
    #   (a) the furthest along the opening axis  -> shows how far it can go
    #   (b) the closest to the inhibitor-bound crystal -> the docking-relevant one
    apo = apo_indices(names)
    pr = [amb_to_idx(a) for a in POCKET_RES]
    ref_open = np.array([crystal_ca("1PZO")(a) for a in POCKET_RES])

    far_k, far_f, far_v = None, None, -1e9
    for k in apo:
        f = int(np.argmax(proj[k]))
        if proj[k][f] > far_v:
            far_k, far_f, far_v = k, f, float(proj[k][f])
    print(f"  (a) furthest along the axis: {names[far_k]} frame {far_f}"
          f"  axis {far_v:.2f} A (crystal axis is {L:.2f})")
    print(f"      pocket RMSD to 1PZO: "
          f"{kabsch_rmsd(np.asarray(X[far_k][far_f])[pr], ref_open):.2f} A"
          f"  -> it OVERSHOOTS; a wide-open loop is not a bound-like pocket")

    # The right question is not "which frame is closest to 1PZO" - a
    # well-ordered frame is close to EVERY crystal.  It is "which frame is
    # closer to the bound form THAN to the closed form", i.e. maximise
    # RMSD(->closed) - RMSD(->bound).  That is the discriminative measure.
    ref_closed = np.array([crystal_ca("1BTL")(a) for a in POCKET_RES])
    best = (-1e9, None, None, None, None)
    for k in apo:
        cand = np.where(proj[k] > 4.0)[0]
        for f in cand[::5]:
            fr = np.asarray(X[k][f])[pr]
            ro = kabsch_rmsd(fr, ref_open)
            rc = kabsch_rmsd(fr, ref_closed)
            if rc - ro > best[0]:
                best = (rc - ro, k, int(f), ro, rc)
    sel, bo_k, bo_f, bo_ro, bo_rc = best
    print(f"  (b) most BOUND-LIKE open frame (max RMSD_closed - RMSD_bound):")
    print(f"      {names[bo_k]} frame {bo_f}  axis {proj[bo_k][bo_f]:.2f} A")
    print(f"      RMSD to 1PZO {bo_ro:.2f} A   to 1BTL {bo_rc:.2f} A"
          f"   selectivity {sel:+.2f} A")

    # and the ensemble-level version of the same question
    n_sel, n_tot = 0, 0
    for k in apo:
        cand = np.where(proj[k] > 6.5)[0]
        for f in cand[::20]:
            fr = np.asarray(X[k][f])[pr]
            n_tot += 1
            n_sel += kabsch_rmsd(fr, ref_closed) > kabsch_rmsd(fr, ref_open)
    print(f"      across the open ensemble, {100*n_sel/max(n_tot,1):.0f} % of frames"
          f" (n={n_tot}) are closer to 1PZO than to 1BTL")

    closed_k = apo[0]
    closed_f = int(np.argmin(proj[closed_k]))
    for tag, k, f in [("axis_extreme", far_k, far_f),
                      ("bound_like", bo_k, bo_f),
                      ("closed", closed_k, closed_f)]:
        pth = os.path.join(OUT, f"tem_{tag}_{names[k]}_f{f}.pdb")
        write_ca_pdb(pth, np.asarray(X[k][f]), resnums)
        print(f"      wrote {os.path.relpath(pth, os.path.dirname(FIG))}")
    best_k, best_f, best_v = bo_k, bo_f, float(proj[bo_k][bo_f])

    # pocket geometry in open vs closed ensembles
    print(f"\n  pocket mouth geometry ({len(SITE_RES)} lining residues):")
    si = [amb_to_idx(a) for a in SITE_RES]
    h11 = [amb_to_idx(a) for a in range(216, 230)]
    h12 = [amb_to_idx(a) for a in range(271, 281)]
    sub = np.asarray(X[apo_indices(names, "TEM-10")])
    pv = proj[apo_indices(names, "TEM-10")]
    om, cm = pv > 6.5, pv < 1.5
    for tag, m in [("open  ", om), ("closed", cm)]:
        fr = sub[m]
        if len(fr) > 4000:
            fr = fr[:: len(fr) // 4000]
        sep = np.linalg.norm(fr[:, h11].mean(1) - fr[:, h12].mean(1), axis=-1)
        rg = np.sqrt(((fr[:, si] - fr[:, si].mean(1, keepdims=True)) ** 2).sum(-1).mean(-1))
        print(f"    {tag}  n={m.sum():6d}  H11-H12 centroid separation "
              f"{sep.mean():5.2f} +- {sep.std():.2f} A   site Rg {rg.mean():5.2f} A")

    # which experimental structure should you dock into?
    print("\n  docking template choice (pocket RMSD of the most 1PZO-like open frame):")
    q = np.asarray(X[best_k][best_f])[pr]
    for code in ["1BTL", "1ZG4", "1ERM", "1PZO", "1PZP"]:
        ref = np.array([crystal_ca(code)(a) for a in POCKET_RES])
        print(f"    {code}  {kabsch_rmsd(q, ref):5.2f} A")
    print("""
  THE HONEST READING - and it cuts against an easy story:

  1. The opening is real, it is variant-selective (M5), and it runs along the
     axis that crystallography says separates apo from inhibitor-bound.

  2. But the MD never builds a convincingly BOUND-LIKE pocket.  Selecting the
     frame closest to 1PZO just finds a well-ordered frame that is close to
     every crystal, including the closed one.  Using the discriminative
     measure instead, only a minority of open frames are closer to 1PZO than
     to 1BTL, and the furthest-open frames overshoot into a splayed loop that
     resembles neither.

  So: conformational selection gets the lid off; the bound geometry itself is
  substantially induced by the ligand.  That is the same split Bowman &
  Geissler reported for this protein
  (PNAS 2012, doi:10.1073/pnas.1209309109) - and it is a reason to dock into
  the crystal, not into these frames.

  PRACTICAL CONSEQUENCE
    * dock into 1PZO / 1PZP (the ligand-bound open form), not into an apo
      crystal and not into a C-alpha MD frame;
    * the MD's contribution is not the coordinates, it is the POPULATION
      argument: it says which ENZYME VARIANT presents that pocket most often.

  CHEMICAL STARTING POINT
    The 1PZO ligand (PDB chemical id CBT, N,N-bis(4-chlorobenzyl)-1H-tetrazol-
    5-amine) and the 1PZP ligand (FTA) are the known non-competitive chemotype
    for this site.  A ChEMBL similarity search (>=70 %, Morgan/FPSim2) against
    CHEMBL2065 (Beta-lactamase TEM) returns NO analogues of CBT, i.e. the
    series was never followed up in the deposited bioactivity literature.  For
    a pocket with a crystal structure and a 282x variant-selectivity handle,
    that is an unusually open field.""")
    return {"most_open": names[best_k], "axis": best_v}


if __name__ == "__main__":
    main()
