"""
E9 — reactions and SEI: what this pipeline deliberately does NOT claim   [digest p10]

The force field used in E1–E7 is a fixed-topology SMIRNOFF model: bonds are
defined once and can never break or form.  Therefore

  ❌ no EC/DME reduction, no radical chemistry, no Li2EDC/Li2BDC products,
     no SEI growth can appear in these trajectories — not because they do not
     happen, but because the model cannot represent them (digest p3: "不能据此
     直接研究新的断键、成键过程").

What *is* done here is a demonstration that the reaction-event bookkeeping of
the digest ("怎样定义一个反应事件": bond-length/lifetime criteria applied to a
trajectory) returns zero events on a fixed-topology trajectory, and the
reasoning for why a reactive study was not attempted within this project:

  * eReaxFF/ReaxFF parameter sets for Li/EC/DME/TFSI exist in the literature
    but none is bundled with the tools available here; an untested parameter
    set would violate the digest's own rule (p3 "反应力场需要单独的参数化与验证").
  * AIMD at the scale of Fig.26/28 (hundreds of atoms, 10^2–10^3 ps) is far
    outside the CPU/GPU budget of this session; a few ps of AIMD would show
    "no reaction" — which the digest warns is not evidence of absence.

A minimal, honest next step (not run): GFN2-xTB MD of one Li+(EC)_n cluster
with an added electron (doublet) to look for the C–O ring-opening of EC⁻ —
still a *model* observation, not an SEI prediction.
"""
from __future__ import annotations

import json
import os

import numpy as np

from electrolyte_pipeline.e0_systems import WORK

BOND_CRIT_A = {("C", "O"): 1.8, ("C", "C"): 1.9, ("S", "N"): 2.1, ("S", "O"): 1.9, ("C", "F"): 1.7, ("C", "S"): 2.2}


def bond_topology_check(label: str = "LiTFSI_r10", workdir: str = WORK, stride: int = 25) -> dict | None:
    """Apply a distance-based bond criterion frame by frame and count topology changes."""
    from electrolyte_pipeline.traj import Traj
    from MDAnalysis.lib.distances import self_distance_array
    if not os.path.exists(os.path.join(workdir, label, "prod.dcd")):
        return None
    tr = Traj(label, workdir)
    heavy = np.where(tr.element != "H")[0]
    ref = None
    changes, nfr, max_stretch = 0, 0, 0.0
    from scipy.spatial.distance import squareform
    for ts in tr.frames(stride):
        d = squareform(self_distance_array(tr.u.atoms.positions[heavy], box=ts.dimensions))
        el = tr.element[heavy]
        bonded = np.zeros_like(d, bool)
        for (a, b), cut in BOND_CRIT_A.items():
            m = (el[:, None] == a) & (el[None, :] == b) | (el[:, None] == b) & (el[None, :] == a)
            bonded |= m & (d < cut)
        np.fill_diagonal(bonded, False)
        if ref is None:
            ref = bonded.copy()
            ref_d = d[ref]
        else:
            changes += int((bonded != ref).sum() // 2)
            max_stretch = max(max_stretch, float((d[ref] / ref_d).max()))
        nfr += 1
    res = {"label": label, "frames_checked": nfr, "bond_criterion_A": {f"{a}-{b}": c for (a, b), c in BOND_CRIT_A.items()},
           "bonds_in_reference_frame": int(ref.sum() // 2), "topology_changes_detected": changes,
           "max_bond_stretch_ratio": max_stretch,
           "verdict": "❌ fixed-topology force field: zero bond-topology events by construction; "
                      "this is a property of the model, not a chemical result"}
    json.dump(res, open(os.path.join(tr.dir, "e9_reactive.json"), "w"), indent=1)
    return res


def report(workdir: str = WORK):
    print("E9  reactions / SEI — scope statement")
    r = bond_topology_check(workdir=workdir)
    if r:
        print(f"    {r['label']}: {r['bonds_in_reference_frame']} heavy-atom bonds tracked over {r['frames_checked']} frames; "
              f"topology changes = {r['topology_changes_detected']}; max stretch ratio {r['max_bond_stretch_ratio']:.2f}")
        print(f"    {r['verdict']}")
    print("    Reactive MD (ReaxFF/eReaxFF) and AIMD were NOT run: no validated Li/EC/DME/TFSI reactive parameter set "
          "is available here, and short AIMD showing 'no reaction' would not be evidence (digest p3/p10).")
    return r


if __name__ == "__main__":
    report()
