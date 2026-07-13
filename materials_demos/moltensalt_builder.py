"""
Client-ready initial-structure builder for real molten-salt systems.

This bridges demo #4 (LJ proxy) to a production run: it builds a physically
reasonable initial configuration for a REAL molten salt at a target density and
writes the initial structure for CP2K / VASP AIMD. Swap composition/density and
you have the starting point for the DeePMD data-generation campaign.

Supported presets (extend freely):
    LiCl-KCl  eutectic (58.5 : 41.5 mol%, ~1.6 g/cm^3 at 723 K)
    FLiNaK    (LiF-NaF-KF 46.5:11.5:42 mol%, ~2.02 g/cm^3 at 773 K)
    NaCl      (~1.54 g/cm^3 near melting)

Usage:
    python3 moltensalt_builder.py FLiNaK --natoms 100 --out flinak_init
"""

import argparse
import os
import numpy as np
from ase import Atoms
from ase.data import atomic_masses, atomic_numbers

# preset: ion -> (count ratio), density g/cm^3, reference temperature K
PRESETS = {
    "LiCl-KCl": dict(ions={"Li": 0.585, "K": 0.415, "Cl": 1.0}, rho=1.60, T=723),
    "FLiNaK":   dict(ions={"Li": 0.465, "Na": 0.115, "K": 0.42, "F": 1.0}, rho=2.02, T=773),
    "NaCl":     dict(ions={"Na": 1.0, "Cl": 1.0}, rho=1.54, T=1100),
}

AMU_TO_G = 1.66053906660e-24
ANG3_TO_CM3 = 1e-24


def build(preset, n_target, seed=0):
    spec = PRESETS[preset]
    ratios = spec["ions"]
    # cations sum must balance anions for charge neutrality (mono-valent ions here)
    cations = {k: v for k, v in ratios.items() if k in ("Li", "Na", "K")}
    anions = {k: v for k, v in ratios.items() if k in ("F", "Cl")}
    csum = sum(cations.values())
    # normalise so #cations == #anions
    unit = {**{k: v / csum for k, v in cations.items()},
            **{k: v / sum(anions.values()) for k, v in anions.items()}}

    # scale to ~n_target atoms (each "formula unit" here = 1 cation + 1 anion)
    n_pairs = max(1, round(n_target / 2))
    counts = {}
    for el, frac in unit.items():
        counts[el] = max(1, round(frac * n_pairs))
    symbols = []
    for el, c in counts.items():
        symbols += [el] * c
    n = len(symbols)

    # cell edge from target density
    mass_g = sum(atomic_masses[atomic_numbers[s]] for s in symbols) * AMU_TO_G
    vol_cm3 = mass_g / spec["rho"]
    vol_ang3 = vol_cm3 / ANG3_TO_CM3
    a = vol_ang3 ** (1.0 / 3.0)

    # random placement on a slightly jittered grid (avoids atomic overlap)
    rng = np.random.default_rng(seed)
    m = int(np.ceil(n ** (1.0 / 3.0)))
    grid = [(i, j, k) for i in range(m) for j in range(m) for k in range(m)]
    grid = np.array(grid[:n], dtype=float)
    pos = (grid + 0.5) / m * a + rng.normal(0, a / (m * 12), size=(n, 3))
    pos = np.mod(pos, a)
    rng.shuffle(symbols)

    atoms = Atoms(symbols=symbols, positions=pos, cell=[a, a, a], pbc=True)
    return atoms, spec, counts, a


def write_inputs(atoms, preset, spec, counts, a, outdir):
    os.makedirs(outdir, exist_ok=True)
    from ase.io import write
    # VASP POSCAR (sort=True groups atoms by element as VASP requires)
    write(os.path.join(outdir, "POSCAR"), atoms, format="vasp", direct=True, sort=True)
    # xyz for CP2K
    write(os.path.join(outdir, "init.xyz"), atoms, format="xyz")

    comp = "  ".join(f"{el}:{c}" for el, c in counts.items())
    with open(os.path.join(outdir, "SYSTEM.md"), "w") as fh:
        fh.write(f"# {preset} molten-salt initial structure\n\n")
        fh.write(f"- composition (atom counts): {comp}\n")
        fh.write(f"- total atoms: {len(atoms)}\n")
        fh.write(f"- target density: {spec['rho']} g/cm^3 (ref T = {spec['T']} K)\n")
        fh.write(f"- cubic cell edge: {a:.3f} Angstrom\n")
        fh.write(f"- files: POSCAR (VASP), init.xyz (CP2K), cell edge above\n\n")
        fh.write("Next: run NVT AIMD at the reference T (see hpc_templates/), "
                 "sample uncorrelated frames, label E/forces, package with dpdata "
                 "into deepmd/npy (see moltensalt_mlp_pipeline.py pack_deepmd()).\n")
    return outdir


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("preset", choices=list(PRESETS), nargs="?", default="FLiNaK")
    ap.add_argument("--natoms", type=int, default=100)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    atoms, spec, counts, a = build(args.preset, args.natoms)
    outdir = args.out or os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                      f"moltensalt_{args.preset}")
    write_inputs(atoms, args.preset, spec, counts, a, outdir)
    print(f"[{args.preset}] {len(atoms)} atoms  cell={a:.2f} A  "
          f"rho={spec['rho']} g/cm^3 @ {spec['T']} K")
    print(f"  composition: {counts}")
    print(f"  wrote POSCAR + init.xyz + SYSTEM.md -> {outdir}")


if __name__ == "__main__":
    main()
