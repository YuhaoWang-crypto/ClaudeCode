#!/usr/bin/env python3
"""
Milestone 3 (data side) -- turn DFT labels into a DeePMD training set, and
validate the data plumbing locally.

Two entry points:

  1. from_dft(...)      : the REAL path. Uses `dpdata` to read QE / VASP / CP2K
                          outputs and write DeePMD 'npy' systems. Run this on
                          your DFT/AIMD results (Modal or HPC).

  2. plumbing_demo(...) : a LOCAL, dependency-light check. It builds a handful
                          of perturbed graphene(+water) frames, labels them with
                          ASE's EMT calculator (a cheap stand-in for QE -- NOT
                          physically meaningful, purely to exercise the export),
                          writes them in the exact DeePMD raw/npy layout by hand,
                          and reloads them to assert the shapes are correct.
                          Swapping EMT -> QE is a one-line change; the plumbing
                          this validates is identical.

    python prepare_deepmd_data.py --demo
    python prepare_deepmd_data.py --from-dft ./qe_runs --out ../data/train/pore_crossing
"""
import argparse
import os
import pathlib
import numpy as np


# ---------------------------------------------------------------- real path --
def from_dft(input_dir, out_dir, fmt="qe/pw"):
    """Convert DFT outputs in `input_dir` to a DeePMD npy system in `out_dir`.

    fmt examples: 'qe/pw' (pw.x scf/relax), 'vasp/outcar', 'cp2k/output'.
    Requires `pip install dpdata`.
    """
    try:
        import dpdata
    except ImportError:
        raise SystemExit("pip install dpdata to use --from-dft")

    ms = dpdata.MultiSystems()
    found = 0
    for root, _, files in os.walk(input_dir):
        for f in files:
            path = os.path.join(root, f)
            try:
                sys = dpdata.LabeledSystem(path, fmt=fmt)
                if len(sys) > 0:
                    ms.append(sys)
                    found += 1
            except Exception:
                continue
    if found == 0:
        raise SystemExit(f"No parseable {fmt} outputs under {input_dir}")
    ms.to_deepmd_npy(out_dir)
    print(f"Wrote DeePMD dataset to {out_dir} from {found} DFT runs.")
    return out_dir


# ------------------------------------------------------- local plumbing demo --
def _write_deepmd_npy(out_dir, symbols, type_map, coords, energies, forces,
                      boxes):
    """Write a DeePMD 'npy' system by hand (no dpdata needed).

    Layout:
        out_dir/type.raw          (natoms ints, one per atom)
        out_dir/type_map.raw      (element names, one per line)
        out_dir/set.000/coord.npy (nframes, natoms*3)
        out_dir/set.000/box.npy   (nframes, 9)
        out_dir/set.000/energy.npy(nframes,)
        out_dir/set.000/force.npy (nframes, natoms*3)
    """
    os.makedirs(os.path.join(out_dir, "set.000"), exist_ok=True)
    type_idx = [type_map.index(s) for s in symbols]
    with open(os.path.join(out_dir, "type.raw"), "w") as fh:
        fh.write(" ".join(map(str, type_idx)) + "\n")
    with open(os.path.join(out_dir, "type_map.raw"), "w") as fh:
        fh.write("\n".join(type_map) + "\n")
    nf = len(energies)
    np.save(os.path.join(out_dir, "set.000", "coord.npy"),
            np.asarray(coords).reshape(nf, -1))
    np.save(os.path.join(out_dir, "set.000", "box.npy"),
            np.asarray(boxes).reshape(nf, 9))
    np.save(os.path.join(out_dir, "set.000", "energy.npy"),
            np.asarray(energies).reshape(nf))
    np.save(os.path.join(out_dir, "set.000", "force.npy"),
            np.asarray(forces).reshape(nf, -1))


def to_extxyz(system_dir, out_xyz, energy_key="REF_energy",
              forces_key="REF_forces"):
    """Convert a DeePMD 'npy' system to extended-XYZ for MACE / NequIP.

    DeePMD trains from the npy layout; MACE and NequIP read extended-XYZ with a
    per-frame energy in info[energy_key] and per-atom forces in
    arrays[forces_key]. This lets the SAME DFT data feed all three trainers for
    an apples-to-apples bake-off.
    """
    from ase import Atoms
    from ase.io import write

    type_map = pathlib.Path(
        os.path.join(system_dir, "type_map.raw")).read_text().split()
    type_idx = [int(x) for x in pathlib.Path(
        os.path.join(system_dir, "type.raw")).read_text().split()]
    symbols = [type_map[i] for i in type_idx]
    natoms = len(symbols)

    set_dirs = sorted(d for d in os.listdir(system_dir)
                      if d.startswith("set."))
    frames = []
    for sd in set_dirs:
        p = os.path.join(system_dir, sd)
        coord = np.load(os.path.join(p, "coord.npy")).reshape(-1, natoms, 3)
        box = np.load(os.path.join(p, "box.npy")).reshape(-1, 3, 3)
        energy = np.load(os.path.join(p, "energy.npy")).reshape(-1)
        force = np.load(os.path.join(p, "force.npy")).reshape(-1, natoms, 3)
        for k in range(coord.shape[0]):
            atoms = Atoms(symbols=symbols, positions=coord[k],
                          cell=box[k], pbc=True)
            atoms.info[energy_key] = float(energy[k])
            atoms.arrays[forces_key] = force[k]
            frames.append(atoms)

    os.makedirs(os.path.dirname(os.path.abspath(out_xyz)), exist_ok=True)
    write(out_xyz, frames, format="extxyz")
    print(f"Wrote {len(frames)} frames -> {out_xyz}")
    print(f"  keys: info['{energy_key}']  arrays['{forces_key}']  "
          f"({natoms} atoms, symbols {sorted(set(symbols))})")
    return out_xyz


def plumbing_demo(out_dir, n_frames=8, seed=0):
    """Generate EMT-labelled graphene frames and export as a DeePMD system."""
    from ase.build import graphene
    from ase.calculators.emt import EMT

    rng = np.random.default_rng(seed)
    base = graphene(size=(4, 4, 1), vacuum=8.0)
    symbols = base.get_chemical_symbols()
    type_map = sorted(set(symbols))

    coords, energies, forces, boxes = [], [], [], []
    for _ in range(n_frames):
        atoms = base.copy()
        atoms.positions += rng.normal(0.0, 0.05, size=atoms.positions.shape)
        atoms.calc = EMT()
        energies.append(atoms.get_potential_energy())
        forces.append(atoms.get_forces())
        coords.append(atoms.get_positions())
        boxes.append(np.array(atoms.cell))

    _write_deepmd_npy(out_dir, symbols, type_map, coords, energies, forces,
                      boxes)

    # ---- reload and assert shapes (the actual validation) ------------------
    natoms = len(symbols)
    c = np.load(os.path.join(out_dir, "set.000", "coord.npy"))
    e = np.load(os.path.join(out_dir, "set.000", "energy.npy"))
    f = np.load(os.path.join(out_dir, "set.000", "force.npy"))
    b = np.load(os.path.join(out_dir, "set.000", "box.npy"))
    assert c.shape == (n_frames, natoms * 3), c.shape
    assert e.shape == (n_frames,), e.shape
    assert f.shape == (n_frames, natoms * 3), f.shape
    assert b.shape == (n_frames, 9), b.shape
    print("Plumbing demo PASS:")
    print(f"  wrote {n_frames} frames of {natoms}-atom graphene to {out_dir}")
    print(f"  type_map={type_map}  coord{c.shape}  energy{e.shape}  "
          f"force{f.shape}  box{b.shape}")
    print("  -> swap EMT for QE (from_dft) and the same layout carries real "
          "DFT labels into DeePMD.")
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--demo", action="store_true",
                    help="Run the local EMT plumbing validation")
    ap.add_argument("--from-dft", metavar="DIR",
                    help="Directory of DFT outputs to convert via dpdata")
    ap.add_argument("--to-extxyz", metavar="SYSTEM_DIR",
                    help="Convert a DeePMD npy system to extended-XYZ "
                         "(for MACE/NequIP)")
    ap.add_argument("--fmt", default="qe/pw")
    ap.add_argument("--out", default="../data/train/demo_system")
    args = ap.parse_args()

    if args.from_dft:
        from_dft(args.from_dft, args.out, fmt=args.fmt)
    elif args.to_extxyz:
        to_extxyz(args.to_extxyz, args.out)
    else:
        plumbing_demo(args.out if args.demo else "../data/train/demo_system")


if __name__ == "__main__":
    main()
