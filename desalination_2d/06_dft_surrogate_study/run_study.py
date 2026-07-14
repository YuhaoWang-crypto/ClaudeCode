#!/usr/bin/env python3
"""
End-to-end pore-crossing study for 2D RO membranes, run LOCALLY with the CHGNet
universal potential as a DFT surrogate.

  Materials : graphene, then MoS2
  Permeants : H2O, Na, Cl
  Observable: the energy profile E(z) as the permeant is dragged along z through
              the pore (a reaction-coordinate PMF). The crossing barrier
              E_a = max(E) - E(feed side) is the first-order selectivity number:
              a good desalination pore has LOW water barrier, HIGH ion barriers.

IMPORTANT -- what this is and isn't
-----------------------------------
CHGNet is a universal MLIP trained on Materials-Project PBE data. Using it here
makes the WHOLE pipeline run end-to-end on a laptop and gives real, reproducible
numbers of the right order of magnitude and the right qualitative trends
(pore-size scaling, water-vs-ion ordering, graphene-vs-MoS2). It is NOT a
substitute for the user's own vdW-DF2 DFT on charged, explicitly-hydrated ions:
absolute barriers and the ionic dehydration penalty require that. Every number
here is labelled "CHGNet surrogate" in the report for exactly this reason.

Run:
    python run_study.py --quick     # small smoke test
    python run_study.py             # full study -> results/*.json, *.png, table
"""
import argparse
import json
import os
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "01_build"))
import build_graphene_pore as bgp   # noqa: E402
import build_mos2_pore as bmp       # noqa: E402

from ase import Atom, Atoms          # noqa: E402
from ase.constraints import FixAtoms, FixedPlane  # noqa: E402
from ase.optimize import BFGS        # noqa: E402

_CALC = None


def get_calc():
    global _CALC
    if _CALC is None:
        from chgnet.model.dynamics import CHGNetCalculator
        _CALC = CHGNetCalculator()
    return _CALC


def make_permeant(kind):
    """Return (Atoms, anchor_index) for the permeant; anchor is the atom whose
    z we drive along the reaction coordinate."""
    if kind == "H2O":
        d, ang = 0.9572, np.deg2rad(104.52)
        mol = Atoms("OH2", positions=[
            [0, 0, 0],
            [d, 0, 0],
            [d * np.cos(ang), d * np.sin(ang), 0]])
        return mol, 0
    return Atoms(kind, positions=[[0, 0, 0]]), 0


def pore_centre(membrane):
    cx = membrane.cell[0, 0] / 2 + membrane.cell[1, 0] / 2
    cy = membrane.cell[1, 1] / 2
    z0 = membrane.positions[:, 2].mean()
    return cx, cy, z0


def scan_profile(membrane, permeant_kind, z_rel, fmax=0.2, max_steps=25):
    """Drive the permeant anchor to each z in z_rel (relative to the membrane
    plane); at each z relax the permeant in-plane (membrane fixed). Return the
    energy profile (eV, referenced to the first/feed point)."""
    calc = get_calc()
    cx, cy, z0 = pore_centre(membrane)
    n_mem = len(membrane)

    energies = []
    for zr in z_rel:
        system = membrane.copy()
        mol, anchor = make_permeant(permeant_kind)
        mol.translate([cx, cy, z0 + zr] - mol.positions[anchor])
        system += mol
        system.calc = calc

        # fix the membrane; constrain the anchor to its z-plane (free x,y),
        # let the rest of the permeant relax freely.
        anchor_global = n_mem + anchor
        cons = [FixAtoms(indices=list(range(n_mem))),
                FixedPlane(anchor_global, [0, 0, 1])]
        system.set_constraint(cons)
        opt = BFGS(system, logfile=None)
        try:
            opt.run(fmax=fmax, steps=max_steps)
            e = system.get_potential_energy()
        except Exception as exc:  # pragma: no cover
            print(f"    [warn] z={zr:+.1f} failed: {exc}")
            e = np.nan
        energies.append(e)
    energies = np.array(energies)
    energies = energies - energies[0]   # reference to feed side
    return energies


def barrier_stats(z_rel, energies):
    valid = ~np.isnan(energies)
    e = energies[valid]
    if len(e) < 2:
        return {"barrier_eV": None}
    imax = int(np.nanargmax(energies))
    imin = int(np.nanargmin(energies))
    return {
        "barrier_eV": round(float(np.nanmax(energies)), 3),
        "min_energy_eV": round(float(np.nanmin(energies)), 3),
        "z_at_barrier": round(float(z_rel[imax]), 2),
        "z_at_min": round(float(z_rel[imin]), 2),
    }


def run_material(name, membrane, permeants, z_rel, outdir, fmax, max_steps):
    print(f"\n=== {name}: {len(membrane)} atoms, cell "
          f"{[round(x,2) for x in membrane.cell.lengths()]} ===")
    results = {"material": name, "natoms": len(membrane), "permeants": {}}
    profiles = {}
    for p in permeants:
        print(f"  scanning {p} ...", flush=True)
        e = scan_profile(membrane, p, z_rel, fmax=fmax, max_steps=max_steps)
        stats = barrier_stats(z_rel, e)
        results["permeants"][p] = {**stats, "profile_eV": [
            round(float(x), 4) if not np.isnan(x) else None for x in e]}
        profiles[p] = e
        print(f"    E_a({p}) = {stats.get('barrier_eV')} eV "
              f"at z={stats.get('z_at_barrier')} A")
    _plot(name, z_rel, profiles, outdir)
    return results


def _plot(name, z_rel, profiles, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 4.2))
    colors = {"H2O": "steelblue", "Na": "darkorange", "Cl": "seagreen"}
    for p, e in profiles.items():
        ax.plot(z_rel, e, "o-", label=p, color=colors.get(p))
    ax.axvline(0, color="k", ls="--", lw=1, alpha=0.6)
    ax.set_xlabel("z relative to membrane plane (Å)  —  feed < 0 < permeate")
    ax.set_ylabel("energy relative to feed (eV)")
    ax.set_title(f"{name}: pore-crossing profiles (CHGNet surrogate)")
    ax.legend()
    fig.tight_layout()
    path = os.path.join(outdir, f"profile_{name}.png")
    fig.savefig(path, dpi=130)
    plt.close(fig)
    print(f"  wrote {path}")


def pore_size_scan(radii, z_rel, outdir, fmax, max_steps):
    """Part A: water crossing barrier vs graphene pore radius. The central
    engineering trade-off: smaller pore rejects more but throttles water."""
    print("\n########## PART A: water barrier vs graphene pore size ##########")
    rows = []
    for r in radii:
        membrane, info = bgp.build_pore(6, 4, r, passivate="H")
        e = scan_profile(membrane, "H2O", z_rel, fmax=fmax, max_steps=max_steps)
        stats = barrier_stats(z_rel, e)
        rows.append({"pore_radius_A": r, "n_removed_C": info["n_removed"],
                     **stats})
        print(f"  r={r:.1f} A ({info['n_removed']} C removed): "
              f"E_a(H2O) = {stats['barrier_eV']} eV")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 4.2))
    rr = [x["pore_radius_A"] for x in rows]
    ba = [x["barrier_eV"] for x in rows]
    ax.plot(rr, ba, "o-", color="steelblue")
    ax.set_xlabel("graphene pore radius (Å)")
    ax.set_ylabel("water crossing barrier E_a (eV)")
    ax.set_title("Part A: smaller pore → higher water barrier (CHGNet surrogate)")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    p = os.path.join(outdir, "poresize_water_barrier.png")
    fig.savefig(p, dpi=130)
    plt.close(fig)
    print(f"  wrote {p}")
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true",
                    help="Small/fast smoke test (fewer z points, Na only)")
    ap.add_argument("--pore-radius", type=float, default=2.5,
                    help="Pore radius for the Part-B selectivity comparison")
    ap.add_argument("--outdir", default=os.path.join(
        os.path.dirname(__file__), "results"))
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    if args.quick:
        z_rel = np.linspace(-4, 4, 5)
        permeants = ["Na"]
        radii = [2.5]
        fmax, max_steps = 0.3, 12
    else:
        z_rel = np.linspace(-5, 5, 11)
        permeants = ["H2O", "Na", "Cl"]
        radii = [2.0, 2.5, 3.0, 3.5]
        fmax, max_steps = 0.2, 25

    # Part A: pore-size scan (graphene, water)
    poresize_rows = pore_size_scan(radii, z_rel, args.outdir, fmax, max_steps)

    # Part B: 3-permeant selectivity at a moderate pore, both materials
    print("\n########## PART B: H2O/Na/Cl selectivity at r="
          f"{args.pore_radius} Å ##########")
    graphene, _ = bgp.build_pore(6, 4, args.pore_radius, passivate="H")
    mos2, _ = bmp.build_pore(6, 6, args.pore_radius, rim="auto")
    materials = [run_material("graphene", graphene, permeants, z_rel,
                              args.outdir, fmax, max_steps),
                 run_material("MoS2", mos2, permeants, z_rel,
                              args.outdir, fmax, max_steps)]

    out = {
        "method": "CHGNet universal potential (DFT surrogate)",
        "z_reaction_coord_A": [round(float(z), 2) for z in z_rel],
        "partA_poresize_scan_graphene_H2O": poresize_rows,
        "partB_selectivity_pore_radius_A": args.pore_radius,
        "partB_materials": materials,
    }
    jpath = os.path.join(args.outdir, "study_results.json")
    with open(jpath, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\nWrote {jpath}")

    print("\n============ PART B SELECTIVITY SUMMARY (CHGNet surrogate) ======")
    print(f"{'material':10s} {'permeant':8s} {'E_a (eV)':>9s} {'z* (A)':>7s}")
    for m in materials:
        for p, s in m["permeants"].items():
            print(f"{m['material']:10s} {p:8s} {str(s['barrier_eV']):>9s} "
                  f"{str(s['z_at_barrier']):>7s}")


if __name__ == "__main__":
    main()
