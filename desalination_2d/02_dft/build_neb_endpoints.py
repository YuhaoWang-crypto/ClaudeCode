#!/usr/bin/env python3
"""
Milestone 2 -- build the two NEB endpoints for a Na+ (or Cl-/water) crossing a
graphene nanopore, and emit a ready-to-run neb.x input.

    python build_neb_endpoints.py --membrane ../figures/graphene_pore.xyz \
        --ion Na --out qe/neb_na_ready.in

Endpoints
---------
  initial : ion sitting in the feed reservoir, ~5 A above the pore centre
  final   : ion in the plane of the membrane, at the pore centre (the barrier
            top is usually near here -- CI-NEB refines the true saddle)

The barrier E_a from the resulting NEB is your FIRST physical selectivity
number and the ground-truth the MLP must later reproduce.
"""
import argparse
import numpy as np
from ase.io import read, write
from ase import Atom

CHARGE = {"Na": +1.0, "Cl": -1.0, "K": +1.0, "Li": +1.0}
PSEUDOS = {
    "C": "C.pbe-n-kjpaw_psl.1.0.0.UPF",
    "H": "H.pbe-kjpaw_psl.1.0.0.UPF",
    "Na": "Na.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Cl": "Cl.pbe-n-kjpaw_psl.1.0.0.UPF",
}


def pore_centre(membrane):
    """Return (cx, cy, z_plane) of the pore: in-plane cell centre, mean sheet z."""
    cx = membrane.cell[0, 0] / 2.0 + membrane.cell[1, 0] / 2.0
    cy = membrane.cell[1, 1] / 2.0
    z_plane = membrane.positions[membrane.numbers == 6][:, 2].mean()
    return cx, cy, z_plane


def make_endpoints(membrane_file, ion, height=5.0):
    membrane = read(membrane_file)
    cx, cy, z0 = pore_centre(membrane)

    initial = membrane.copy()
    initial.append(Atom(ion, position=(cx, cy, z0 + height)))

    final = membrane.copy()
    final.append(Atom(ion, position=(cx, cy, z0)))  # ion in the pore plane
    return initial, final


def write_neb_input(initial, final, out, ion, ecutwfc=60.0, num_images=7,
                    kpts=(3, 3, 1)):
    """Assemble a neb.x input by splicing ASE-formatted positions into the
    two-image PATH template."""
    def positions_block(atoms):
        # build ATOMIC_POSITIONS directly (robust, no dependence on the exact
        # formatting of ASE's espresso writer across versions)
        lines = ["ATOMIC_POSITIONS angstrom"]
        for sym, (x, y, z) in zip(atoms.get_chemical_symbols(),
                                  atoms.positions):
            lines.append(f"   {sym:<3s} {x:14.9f} {y:14.9f} {z:14.9f}")
        return "\n".join(lines)

    species = sorted(set(initial.get_chemical_symbols()))
    charge = CHARGE.get(ion, 0.0)
    cell = initial.cell
    species_block = "\n".join(
        f"   {s:<3s} {Atom(s).mass:8.3f}  {PSEUDOS[s]}" for s in species)

    inp = f"""BEGIN
BEGIN_PATH_INPUT
&PATH
    string_method   = 'neb'
    nstep_path      = 100
    opt_scheme      = 'broyden'
    num_of_images   = {num_images}
    k_max           = 0.3
    k_min           = 0.2
    CI_scheme       = 'auto'
    path_thr        = 0.05
/
END_PATH_INPUT
BEGIN_ENGINE_INPUT
&CONTROL
    prefix     = 'neb_{ion.lower()}'
    outdir     = './out'
    pseudo_dir = './pseudo'
    tprnfor    = .true.
/
&SYSTEM
    ibrav       = 0
    nat         = {len(initial)}
    ntyp        = {len(species)}
    ecutwfc     = {ecutwfc}
    ecutrho     = {ecutwfc * 8:.1f}
    occupations = 'smearing'
    smearing    = 'mv'
    degauss     = 0.01
    input_dft   = 'vdw-df2'
    tot_charge  = {charge}
/
&ELECTRONS
    conv_thr    = 1.0d-7
    mixing_beta = 0.3
/
CELL_PARAMETERS angstrom
   {cell[0,0]:14.9f} {cell[0,1]:14.9f} {cell[0,2]:14.9f}
   {cell[1,0]:14.9f} {cell[1,1]:14.9f} {cell[1,2]:14.9f}
   {cell[2,0]:14.9f} {cell[2,1]:14.9f} {cell[2,2]:14.9f}

ATOMIC_SPECIES
{species_block}

BEGIN_POSITIONS
FIRST_IMAGE
{positions_block(initial)}

LAST_IMAGE
{positions_block(final)}
END_POSITIONS

K_POINTS automatic
   {kpts[0]} {kpts[1]} {kpts[2]} 0 0 0
END_ENGINE_INPUT
END
"""
    with open(out, "w") as fh:
        fh.write(inp)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--membrane", required=True)
    ap.add_argument("--ion", default="Na", choices=list(PSEUDOS.keys() - {"C", "H"}) + ["Na", "Cl"])
    ap.add_argument("--height", type=float, default=5.0,
                    help="Initial ion height above the pore plane (A)")
    ap.add_argument("--num-images", type=int, default=7)
    ap.add_argument("--ecutwfc", type=float, default=60.0)
    ap.add_argument("--out", default="qe/neb_ready.in")
    args = ap.parse_args()

    initial, final = make_endpoints(args.membrane, args.ion, args.height)
    write("neb_initial.xyz", initial)
    write("neb_final.xyz", final)
    write_neb_input(initial, final, args.out, args.ion,
                    ecutwfc=args.ecutwfc, num_images=args.num_images)
    print(f"Wrote {args.out}")
    print(f"  ion={args.ion} (charge {CHARGE.get(args.ion, 0.0):+.0f}), "
          f"nat={len(initial)}, images={args.num_images}")
    print("  endpoints: neb_initial.xyz (feed) -> neb_final.xyz (pore plane)")
    print("  run: neb.x -inp", args.out, "> neb.out ; "
          "then parse with parse_neb_barrier.py")


if __name__ == "__main__":
    main()
