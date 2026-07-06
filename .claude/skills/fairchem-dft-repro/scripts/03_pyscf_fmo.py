"""
03_pyscf_fmo.py
Reproduce the Frontier Molecular Orbital (FMO) deliverable for C13H10N4S.

fairchem MLIPs output energy/forces only and CANNOT produce molecular orbitals,
so this deliverable is reproduced with PySCF (open-source quantum chemistry),
matching the report's Gaussian level of theory as closely as possible:
    B3LYP / 6-31G(d)   (report used B3LYP/6-31G(d) for geometry; PBE1PBE/def2-TZVP
                        for single points -- we report B3LYP/6-31G(d) HOMO/LUMO).

Note on scope: the report's tabulated HOMO/LUMO (gap 0.065 eV) were computed on
the molecule *adsorbed on the Fe(110) metal slab*, so that tiny gap reflects
metallic states, not the molecule itself. A Gaussian-basis QM code cannot treat
that periodic metal. Here we compute the chemically meaningful FMO of the FREE
inhibitor molecule -- the standard corrosion-inhibitor descriptor (HOMO, LUMO,
gap, hardness, electronegativity, electrophilicity).

Geometry: UMA(omol)-optimized gas-phase C13H10N4S.
"""
import os, json, time
import numpy as np
from ase.io import read, write

HERE = os.path.dirname(os.path.abspath(__file__))
SDIR = os.path.join(HERE, "..", "structures")
RDIR = os.path.join(HERE, "..", "results")
HARTREE_EV = 27.211386245988


def uma_optimize_gasphase():
    from ase.optimize import LBFGS
    from fairchem.core import pretrained_mlip, FAIRChemCalculator
    atoms = read(os.path.join(SDIR, "inhibitor_C13H10N4S.xyz"))
    atoms.info.update({"charge": 0, "spin": 1})
    p = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cpu")
    atoms.calc = FAIRChemCalculator(p, task_name="omol")
    LBFGS(atoms, logfile=os.path.join(RDIR, "opt_inhibitor_gas.log")).run(
        fmax=0.03, steps=300)
    write(os.path.join(SDIR, "relaxed_inhibitor_gas.xyz"), atoms)
    return atoms


def pyscf_fmo(atoms):
    from pyscf import gto, dft
    mol = gto.Mole()
    mol.atom = [(s, tuple(p)) for s, p in
                zip(atoms.get_chemical_symbols(), atoms.get_positions())]
    mol.basis = "6-31g*"
    mol.charge = 0
    mol.spin = 0            # closed shell
    mol.build()
    mf = dft.RKS(mol)
    mf.xc = "b3lyp"
    mf.conv_tol = 1e-8
    t0 = time.time()
    e_tot = mf.kernel()
    dt = time.time() - t0

    mo_e = mf.mo_energy                    # Hartree
    occ = mf.mo_occ
    homo_idx = int(np.max(np.where(occ > 0)))
    lumo_idx = homo_idx + 1
    E_homo = float(mo_e[homo_idx] * HARTREE_EV)
    E_lumo = float(mo_e[lumo_idx] * HARTREE_EV)
    gap = E_lumo - E_homo

    # global reactivity descriptors (corrosion-inhibitor standard)
    I = -E_homo                # ionization potential (Koopmans)
    A = -E_lumo                # electron affinity
    chi = (I + A) / 2.0        # electronegativity
    eta = (I - A) / 2.0        # chemical hardness
    softness = 1.0 / eta if eta != 0 else float("inf")
    omega = chi ** 2 / (2.0 * eta) if eta != 0 else float("inf")
    mu = float(np.linalg.norm(mf.dip_moment(unit="Debye", verbose=0)))

    res = {
        "level_of_theory": "B3LYP/6-31G(d) (PySCF), free molecule",
        "note": ("report's tabulated gap 0.065 eV was for the molecule adsorbed "
                 "on Fe(110) metal; here is the free-molecule FMO"),
        "E_total_Ha": float(e_tot),
        "HOMO_index": homo_idx, "LUMO_index": lumo_idx,
        "E_HOMO_eV": E_homo, "E_LUMO_eV": E_lumo, "gap_eV": gap,
        "ionization_potential_eV": I, "electron_affinity_eV": A,
        "electronegativity_chi_eV": chi, "hardness_eta_eV": eta,
        "softness_eV^-1": softness, "electrophilicity_omega_eV": omega,
        "dipole_Debye": mu, "scf_seconds": dt, "n_basis": mol.nao,
    }
    print(f"  SCF done: E={e_tot:.6f} Ha in {dt:.0f}s, {mol.nao} basis fns")
    print(f"  HOMO={E_homo:.3f} eV  LUMO={E_lumo:.3f} eV  gap={gap:.3f} eV")
    print(f"  chi={chi:.2f}  eta={eta:.2f}  omega={omega:.2f} eV  mu={mu:.2f} D")
    return res, mf, mol, homo_idx, lumo_idx


def dump_cube(mf, mol, homo_idx, lumo_idx):
    """Write HOMO/LUMO cube files for visualization (optional)."""
    try:
        from pyscf.tools import cubegen
        cubegen.orbital(mol, os.path.join(RDIR, "HOMO.cube"),
                        mf.mo_coeff[:, homo_idx], nx=60, ny=60, nz=60)
        cubegen.orbital(mol, os.path.join(RDIR, "LUMO.cube"),
                        mf.mo_coeff[:, lumo_idx], nx=60, ny=60, nz=60)
        print("  HOMO.cube / LUMO.cube written")
    except Exception as e:
        print(f"  cube generation skipped: {e}")


def main():
    print("[FMO of C13H10N4S] optimizing gas-phase geometry with UMA(omol)...")
    atoms = uma_optimize_gasphase()
    print("[FMO] running B3LYP/6-31G(d) in PySCF...")
    res, mf, mol, hi, li = pyscf_fmo(atoms)
    dump_cube(mf, mol, hi, li)
    with open(os.path.join(RDIR, "inhibitor_fmo.json"), "w") as fh:
        json.dump(res, fh, indent=2)
    print("-> results/inhibitor_fmo.json")


if __name__ == "__main__":
    main()
