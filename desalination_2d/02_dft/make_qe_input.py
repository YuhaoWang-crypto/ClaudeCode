#!/usr/bin/env python3
"""
Milestone 1 -- turn a builder geometry into a READY-TO-RUN Quantum ESPRESSO
SCF input (real nat / cell / positions, no __PLACEHOLDER__ left).

    python make_qe_input.py --xyz ../figures/graphene_pore.xyz \
        --out qe/scf_graphene_ready.in --ecutwfc 60 --kpts 3 3 1

The output is exactly what `modal_run_dft.py` / `slurm/qe_scf.slurm` consume.
Uses ASE's espresso writer so the formatting is guaranteed valid; we only
supply the &control/&system/&electrons namelists and pseudopotential names.

Convergence note: `--ecutwfc` and `--kpts` are the two knobs to converge first.
A porous supercell is large in-plane, so 2x2x1 or 3x3x1 k-points usually
suffice; always converge the total energy w.r.t. both before trusting numbers.
"""
import argparse
from ase.io import read, write

# SSSP-efficiency PBE PAW pseudos (download to ./pseudo). Names match the
# scf_graphene.in template so the two stay consistent.
PSEUDOS = {
    "C": "C.pbe-n-kjpaw_psl.1.0.0.UPF",
    "H": "H.pbe-kjpaw_psl.1.0.0.UPF",
    "O": "O.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Na": "Na.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Cl": "Cl.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Mo": "Mo.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "S": "S.pbe-n-kjpaw_psl.1.0.0.UPF",
}


def build_input_data(ecutwfc, ecutrho_ratio, degauss, functional, tot_charge):
    return {
        "control": {
            "calculation": "scf",
            "prefix": "membrane",
            "outdir": "./out",
            "pseudo_dir": "./pseudo",
            "tprnfor": True,
            "tstress": True,
            "verbosity": "high",
        },
        "system": {
            "ecutwfc": ecutwfc,
            "ecutrho": ecutwfc * ecutrho_ratio,
            "occupations": "smearing",
            "smearing": "mv",
            "degauss": degauss,
            "input_dft": functional,
            "tot_charge": tot_charge,
        },
        "electrons": {
            "conv_thr": 1.0e-7,
            "mixing_beta": 0.3,
            "electron_maxstep": 200,
        },
    }


def make(xyz, out, ecutwfc=60.0, ecutrho_ratio=8.0, kpts=(3, 3, 1),
         degauss=0.01, functional="vdw-df2", tot_charge=0.0):
    atoms = read(xyz)
    species = sorted(set(atoms.get_chemical_symbols()))
    missing = [s for s in species if s not in PSEUDOS]
    if missing:
        raise SystemExit(f"No pseudopotential mapped for: {missing}. "
                         f"Add it to PSEUDOS in make_qe_input.py.")
    pseudopotentials = {s: PSEUDOS[s] for s in species}

    input_data = build_input_data(ecutwfc, ecutrho_ratio, degauss,
                                  functional, tot_charge)
    write(out, atoms, format="espresso-in",
          input_data=input_data,
          pseudopotentials=pseudopotentials,
          kpts=tuple(kpts))
    print(f"Wrote {out}")
    print(f"  species : {species}")
    print(f"  natoms  : {len(atoms)}")
    print(f"  cell(A) : {[round(x,3) for x in atoms.cell.lengths()]}")
    print(f"  ecutwfc : {ecutwfc} Ry   kpts: {tuple(kpts)}   dft: {functional}")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xyz", required=True)
    ap.add_argument("--out", default="qe/scf_ready.in")
    ap.add_argument("--ecutwfc", type=float, default=60.0)
    ap.add_argument("--ecutrho-ratio", type=float, default=8.0)
    ap.add_argument("--kpts", type=int, nargs=3, default=[3, 3, 1])
    ap.add_argument("--degauss", type=float, default=0.01)
    ap.add_argument("--functional", default="vdw-df2")
    ap.add_argument("--tot-charge", type=float, default=0.0)
    args = ap.parse_args()
    make(args.xyz, args.out, args.ecutwfc, args.ecutrho_ratio, tuple(args.kpts),
         args.degauss, args.functional, args.tot_charge)


if __name__ == "__main__":
    main()
