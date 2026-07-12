"""
Demo #6 — Why L- and D-amino acids give different SERS spectra on gold:
a molecular-level DFT demonstration.

Motivation (from the chat context)
-----------------------------------
Client (旺仔秋秋唐): wants to compare L- vs D-arginine adsorbed on gold
nanoparticles by SERS, via DFT: stable adsorption geometry, adsorption energy,
and vibrational (Raman) frequencies on the Au surface.

THE KEY PHYSICS THIS DEMO MAKES CONCRETE
----------------------------------------
Enantiomers are mirror images. In an ACHIRAL environment (isolated molecule, or
an achiral probe) their vibrational spectra are RIGOROUSLY identical - mirror
symmetry maps one Hessian onto the other. Therefore a measured L/D SERS
DIFFERENCE cannot come from the free molecule; it must originate in the
DIFFERENT ADSORPTION GEOMETRY the two enantiomers adopt on the (locally chiral)
gold surface. This demo proves the first half directly with DFT and thereby
shows exactly why the production job needs the periodic Au slab.

WHAT THIS IS (runnable here)
----------------------------
Using L-alanine as the smallest chiral amino-acid proxy for arginine:
    1. build L-alanine (RDKit) and optimise it (PySCF, PBE/6-31G),
    2. construct D-alanine as the exact mirror image,
    3. compute harmonic vibrational frequencies of both,
    4. show the two spectra are identical to numerical precision,
    5. compute an Au-adsorption energy (single Au-atom probe) to exercise the
       adsorption-energy workflow.

WHAT THIS IS *NOT*
------------------
Not the real arginine-on-Au(111) periodic calculation. Production needs: a
periodic Au slab (VASP/QE), the true adsorption footprint, dispersion (D3),
implicit/explicit solvent, and surface-enhanced Raman tensors. Alanine is a
stand-in and a single Au atom is only a binding probe. The point proven here -
free-molecule L=D, difference lives on the surface - is exactly the rationale
for doing the slab calculation.
"""

import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BASIS = "6-31g"
XC = "pbe"


def alanine_geometry():
    """Return (symbols, coords Angstrom) for L-alanine via RDKit embedding."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.AddHs(Chem.MolFromSmiles("C[C@@H](N)C(=O)O"))  # L-alanine
    AllChem.EmbedMolecule(mol, randomSeed=1)
    AllChem.MMFFOptimizeMolecule(mol)
    conf = mol.GetConformer()
    syms = [a.GetSymbol() for a in mol.GetAtoms()]
    xyz = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
    return syms, xyz


def build_mol(syms, xyz, charge=0, spin=0):
    from pyscf import gto
    atom = [[s, tuple(c)] for s, c in zip(syms, xyz)]
    return gto.M(atom=atom, basis=BASIS, charge=charge, spin=spin, verbose=0)


def optimize(syms, xyz):
    from pyscf import dft
    from pyscf.geomopt.berny_solver import optimize as berny_opt
    mol = build_mol(syms, xyz)
    mf = dft.RKS(mol); mf.xc = XC
    mol_eq = berny_opt(mf, maxsteps=80)
    return mol_eq


def frequencies(mol):
    from pyscf import dft
    from pyscf.hessian import thermo
    mf = dft.RKS(mol); mf.xc = XC
    e = mf.kernel()
    hess = mf.Hessian().kernel()
    res = thermo.harmonic_analysis(mol, hess)
    freqs = np.real(res["freq_wavenumber"])
    freqs = np.sort(freqs[freqs > 50.0])  # drop translations/rotations/imaginary
    return e, freqs


def mirror(mol):
    """Exact mirror image (x -> -x): builds the opposite enantiomer."""
    syms = [mol.atom_symbol(i) for i in range(mol.natm)]
    xyz = mol.atom_coords(unit="Angstrom").copy()
    xyz[:, 0] *= -1.0
    return build_mol(syms, xyz)


def au_adsorption_energy(mol_ala):
    """Single Au-atom binding probe: E_ads = E(ala-Au) - E(ala) - E(Au)."""
    from pyscf import dft
    syms = [mol_ala.atom_symbol(i) for i in range(mol_ala.natm)]
    xyz = mol_ala.atom_coords(unit="Angstrom")
    # place Au ~2.4 A from the carboxylate/amino region (near the N atom)
    n_idx = [i for i, s in enumerate(syms) if s == "N"][0]
    au_pos = xyz[n_idx] + np.array([0.0, 0.0, 2.4])

    def energy(atom, spin):
        m = mol_ala.copy()
        from pyscf import gto
        m = gto.M(atom=atom, basis={"default": BASIS, "Au": "lanl2dz"},
                  ecp={"Au": "lanl2dz"}, spin=spin, verbose=0)
        mf = dft.UKS(m); mf.xc = XC; mf.level_shift = 0.2; mf.max_cycle = 150
        return mf.kernel()

    ala_atoms = [[s, tuple(c)] for s, c in zip(syms, xyz)]
    complex_atoms = ala_atoms + [["Au", tuple(au_pos)]]
    au_atoms = [["Au", (0.0, 0.0, 0.0)]]

    e_complex = energy(complex_atoms, spin=1)   # ala(singlet)+Au(doublet) -> doublet
    e_ala = energy(ala_atoms, spin=0)
    e_au = energy(au_atoms, spin=1)
    HART2KCAL = 627.509
    return (e_complex - e_ala - e_au) * HART2KCAL


def main():
    print("Building and optimising L-alanine (PBE/6-31G)...")
    syms, xyz0 = alanine_geometry()
    mol_L = optimize(syms, xyz0)
    print("Computing L-alanine vibrational frequencies...")
    eL, freq_L = frequencies(mol_L)

    print("Constructing D-alanine (mirror image) and its frequencies...")
    mol_D = mirror(mol_L)
    eD, freq_D = frequencies(mol_D)

    n = min(len(freq_L), len(freq_D))
    max_diff = np.max(np.abs(freq_L[:n] - freq_D[:n]))
    print(f"\n  E(L) = {eL:.6f} Ha   E(D) = {eD:.6f} Ha   |dE| = {abs(eL-eD):.2e} Ha")
    print(f"  #modes: L={len(freq_L)} D={len(freq_D)}")
    print(f"  max |freq_L - freq_D| = {max_diff:.4f} cm^-1  "
          f"(-> enantiomers vibrate identically in vacuum)")

    print("\nComputing single-atom Au adsorption energy (workflow probe)...")
    try:
        e_ads = au_adsorption_energy(mol_L)
        print(f"  E_ads(L-alanine + Au probe) = {e_ads:+.1f} kcal/mol")
        print("  (by mirror symmetry the D value on this achiral probe is identical;")
        print("   the real L/D SERS split needs the extended chiral Au slab)")
    except Exception as exc:
        e_ads = None
        print(f"  (Au probe skipped: {type(exc).__name__}: {exc})")

    # ---- figure: overlaid broadened IR/Raman-style spectra ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def broaden(freqs, grid, width=15.0):
        s = np.zeros_like(grid)
        for f in freqs:
            s += np.exp(-0.5 * ((grid - f) / width) ** 2)
        return s

    grid = np.linspace(0, 3800, 2000)
    sL = broaden(freq_L, grid)
    sD = broaden(freq_D, grid)

    fig, ax = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    ax[0].plot(grid, sL, color="#2c6fbb", lw=2, label="L-alanine")
    ax[0].plot(grid, sD, color="#d1495b", lw=1.5, ls="--", label="D-alanine")
    ax[0].set_ylabel("intensity (a.u.)")
    ax[0].set_title(f"Free-molecule vibrational spectra overlie exactly  "
                    f"(max Δν = {max_diff:.3f} cm⁻¹)")
    ax[0].legend(); ax[0].grid(alpha=0.3)

    ax[1].plot(grid, sL - sD, color="#555", lw=1.5)
    ax[1].axhline(0, color="0.7", lw=0.8)
    ax[1].set_ylabel("L − D difference")
    ax[1].set_xlabel("wavenumber (cm⁻¹)")
    ax[1].set_title("Difference spectrum is ~zero: chirality is invisible in vacuum "
                    "→ the SERS difference must come from adsorption geometry on Au")
    ax[1].grid(alpha=0.3)
    ax[1].set_ylim(-1, 1)

    fig.suptitle("Demo #6 · Origin of L/D SERS difference on gold\n"
                 "DFT (PBE/6-31G) · alanine proxy  (free-molecule L≡D; difference lives on the surface)",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = os.path.join(HERE, "sers_chirality_demo.png")
    fig.savefig(out, dpi=150)
    print(f"\nsaved figure -> {out}")


if __name__ == "__main__":
    main()
