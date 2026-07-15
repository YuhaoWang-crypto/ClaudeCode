#!/usr/bin/env python3
"""
Step 2 -- HYDRATED-ion pore-crossing PMF (fixes the bare-ion defect).

The Part-B result in REPORT.md showed bare Na/Cl have NO crossing barrier -- the
opposite of desalination selectivity -- because real rejection comes from the
DEHYDRATION penalty, which only exists when the ion carries its water shell.

This script builds Na+.(H2O)6 and Cl-.(H2O)6 hydrated complexes, drives the ION
along z through the pore (reaction coordinate = ion z), and relaxes the shell at
each step (membrane fixed). As the ion enters the sub-nm pore its waters cannot
follow -> the shell is stripped -> a real barrier appears. We compare
bare-ion vs hydrated-ion vs lone-water on the same pore.

Method label: CHGNet universal-potential surrogate (offline, bundled weights).
Right physics and ordering; absolute eV need the client's QE + vdW-DF2 with
charged supercells (tot_charge = +/-1) and explicit water.

    python run_hydrated_pmf.py            # graphene, r=2.5 A
    python run_hydrated_pmf.py --quick    # fast smoke test
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

import run_study as rs               # reuse calc + pore_centre + plotting bits

# ion-O first-shell distances (A) and typical coordination numbers
ION_O = {"Na": 2.35, "Cl": 3.15}
ION_N = {"Na": 6, "Cl": 6}


def build_hydrated_ion(ion, n_water=6):
    """Ion at origin surrounded by n_water waters (octahedral for n=6), each
    water oxygen pointing at the ion at the first-shell distance."""
    d = ION_O[ion]
    # octahedral directions (+/-x, +/-y, +/-z); trim to n_water
    dirs = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0],
                     [0, -1, 0], [0, 0, 1], [0, 0, -1]], dtype=float)[:n_water]
    dOH, ang = 0.9572, np.deg2rad(104.52)

    atoms = Atoms(ion, positions=[[0, 0, 0]])
    for u in dirs:
        o = u * d
        # build a water with O at o, bisector pointing back toward the ion
        back = -u  # from O toward ion
        # pick an arbitrary perpendicular to place the two H
        perp = np.cross(back, [0, 0, 1])
        if np.linalg.norm(perp) < 1e-3:
            perp = np.cross(back, [0, 1, 0])
        perp /= np.linalg.norm(perp)
        # two H at +/- half the HOH angle around the O-ion axis
        h1 = o + dOH * (np.cos(ang / 2) * back + np.sin(ang / 2) * perp)
        h2 = o + dOH * (np.cos(ang / 2) * back - np.sin(ang / 2) * perp)
        atoms += Atoms("OHH", positions=[o, h1, h2])
    return atoms, 0  # anchor = ion index 0


def scan_permeant(membrane, permeant, anchor, z_rel, fmax=0.15, max_steps=60):
    """Drive `permeant` (prebuilt Atoms, anchor atom index) through the pore."""
    calc = rs.get_calc()
    cx, cy, z0 = rs.pore_centre(membrane)
    n_mem = len(membrane)
    energies = []
    for zr in z_rel:
        system = membrane.copy()
        mol = permeant.copy()
        mol.translate([cx, cy, z0 + zr] - mol.positions[anchor])
        system += mol
        system.calc = calc
        cons = [FixAtoms(indices=list(range(n_mem))),
                FixedPlane(n_mem + anchor, [0, 0, 1])]
        system.set_constraint(cons)
        try:
            BFGS(system, logfile=None).run(fmax=fmax, steps=max_steps)
            energies.append(system.get_potential_energy())
        except Exception as exc:  # pragma: no cover
            print(f"    [warn] z={zr:+.1f}: {exc}")
            energies.append(np.nan)
    e = np.array(energies)
    return e - e[0]


def run(material_name, membrane, z_rel, outdir, fmax, max_steps, do_ions):
    print(f"\n=== {material_name}: {len(membrane)} atoms ===")
    profiles, stats = {}, {}

    # lone water reference
    water, wa = rs.make_permeant("H2O")
    print("  H2O ...", flush=True)
    profiles["H2O"] = scan_permeant(membrane, water, wa, z_rel, fmax, max_steps)
    stats["H2O"] = rs.barrier_stats(z_rel, profiles["H2O"])
    print(f"    E_a(H2O) = {stats['H2O']['barrier_eV']} eV")

    for ion in do_ions:
        # bare ion
        bare, ba = rs.make_permeant(ion)
        print(f"  {ion} (bare) ...", flush=True)
        profiles[f"{ion}_bare"] = scan_permeant(membrane, bare, ba, z_rel,
                                                fmax, max_steps)
        stats[f"{ion}_bare"] = rs.barrier_stats(z_rel, profiles[f"{ion}_bare"])
        print(f"    E_a({ion} bare) = {stats[f'{ion}_bare']['barrier_eV']} eV")

        # hydrated ion
        hyd, ha = build_hydrated_ion(ion, ION_N[ion])
        print(f"  {ion}.(H2O){ION_N[ion]} (hydrated) ...", flush=True)
        profiles[f"{ion}_hyd"] = scan_permeant(membrane, hyd, ha, z_rel,
                                               fmax, max_steps)
        stats[f"{ion}_hyd"] = rs.barrier_stats(z_rel, profiles[f"{ion}_hyd"])
        print(f"    E_a({ion} hydrated) = {stats[f'{ion}_hyd']['barrier_eV']} eV")

    _plot(material_name, z_rel, profiles, outdir)
    return {"material": material_name, "natoms": len(membrane), "stats": stats,
            "profiles": {k: [round(float(x), 4) if not np.isnan(x) else None
                             for x in v] for k, v in profiles.items()}}


def _plot(name, z_rel, profiles, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    style = {"H2O": ("steelblue", "-"),
             "Na_bare": ("orange", "--"), "Na_hyd": ("darkorange", "-"),
             "Cl_bare": ("lightgreen", "--"), "Cl_hyd": ("seagreen", "-")}
    fig, ax = plt.subplots(figsize=(6.5, 4.4))
    for k, e in profiles.items():
        c, ls = style.get(k, ("gray", "-"))
        ax.plot(z_rel, e, "o-", color=c, ls=ls, label=k.replace("_", " "))
    ax.axvline(0, color="k", ls=":", lw=1, alpha=0.6)
    ax.set_xlabel("ion z relative to membrane plane (Å) — feed < 0 < permeate")
    ax.set_ylabel("energy relative to feed (eV)")
    ax.set_title(f"{name}: hydrated vs bare ion vs water (CHGNet surrogate)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    p = os.path.join(outdir, f"hydrated_{name}.png")
    fig.savefig(p, dpi=130)
    plt.close(fig)
    print(f"  wrote {p}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--pore-radius", type=float, default=2.5)
    ap.add_argument("--material", choices=["graphene", "MoS2", "both"],
                    default="graphene")
    ap.add_argument("--outdir", default=os.path.join(
        os.path.dirname(__file__), "results"))
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    if args.quick:
        z_rel = np.linspace(-5, 5, 7)
        do_ions = ["Na"]
        fmax, max_steps = 0.25, 25
    else:
        z_rel = np.linspace(-6, 6, 13)
        do_ions = ["Na", "Cl"]
        fmax, max_steps = 0.15, 60

    mats = []
    if args.material in ("graphene", "both"):
        g, _ = bgp.build_pore(6, 4, args.pore_radius, passivate="H")
        mats.append(("graphene", g))
    if args.material in ("MoS2", "both"):
        m, _ = bmp.build_pore(6, 6, args.pore_radius, rim="auto")
        mats.append(("MoS2", m))

    results = [run(name, mem, z_rel, args.outdir, fmax, max_steps, do_ions)
               for name, mem in mats]

    out = {"method": "CHGNet universal potential (DFT surrogate)",
           "pore_radius_A": args.pore_radius,
           "z_reaction_coord_A": [round(float(z), 2) for z in z_rel],
           "materials": results}
    jpath = os.path.join(args.outdir, "hydrated_pmf_results.json")
    with open(jpath, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\nWrote {jpath}")

    print("\n===== HYDRATED vs BARE SELECTIVITY (CHGNet surrogate) =====")
    print(f"{'material':9s} {'species':12s} {'E_a (eV)':>9s}")
    for m in results:
        for k, s in m["stats"].items():
            print(f"{m['material']:9s} {k:12s} {str(s['barrier_eV']):>9s}")


if __name__ == "__main__":
    main()
