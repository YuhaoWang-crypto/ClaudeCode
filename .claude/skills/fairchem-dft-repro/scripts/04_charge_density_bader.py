"""
04_charge_density_bader.py
Open-source reproduction of the DIFFERENTIAL CHARGE DENSITY (Δρ) + BADER CHARGE
deliverable -- the part fairchem/UMA fundamentally cannot do (an MLIP has no
electron density).

Engine: PySCF (open-source DFT, PBE -- same functional family as the report's
VASP calculation) -> real-space electron density on a grid (Gaussian cube) ->
  * Δρ = ρ(complex) - ρ(cluster) - ρ(adsorbate)   [same grid]   (charge density
    redistribution on binding, i.e. the report's "differential charge density")
  * Bader partitioning of ρ(complex) -> per-atom charges -> net charge transfer
    between the molecule and the oxide (the analog of the report's "0.10 e").

NOTE ON SIZE. Full ab-initio DFT on the real 86-92 atom RSV-TiOx complexes (or
the 128-atom magnetic Fe(110) slab) does not finish on this CPU-only box. This
script therefore demonstrates the *complete, working pipeline* on a small,
closed-shell, chemically-faithful analog: formic acid (HCOOH) -- which mimics
the carboxyl/hydroxyl group through which rosuvastatin binds Ti -- adsorbed on a
Ti2O4 cluster (Ti(IV), d0, closed shell => robust SCF). The identical pipeline
(scripts 05_*) generates ready-to-run inputs for the full systems on GPU/HPC.
"""
import os, json, time
import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.optimize import LBFGS

HERE = os.path.dirname(os.path.abspath(__file__))
SDIR = os.path.join(HERE, "..", "structures")
RDIR = os.path.join(HERE, "..", "results")
BOHR = 0.52917721067

# ---------------------------------------------------------------------------
# 1. Build + UMA-relax the small model (HCOOH on Ti2O4)
# ---------------------------------------------------------------------------
def build_and_relax():
    from fairchem.core import pretrained_mlip, FAIRChemCalculator
    # Ti2O4 (anatase-like) rough geometry
    tio2 = Atoms("Ti2O4", positions=[
        (0.00, 0.00, 0.00), (2.80, 0.00, 0.00),
        (1.40, 1.20, 0.00), (1.40, -1.20, 0.00),
        (-1.10, 0.00, 1.00), (3.90, 0.00, -1.00)])
    # formic acid HCOOH, carboxyl O pointing toward a Ti
    hcooh = Atoms("CHOOH", positions=[
        (1.40, 0.00, 2.60),    # C
        (2.20, 0.00, 3.20),    # H (formyl)
        (0.30, 0.00, 3.10),    # O (=O, will point to Ti)
        (1.55, 0.00, 1.35),    # O (-OH, near Ti)
        (0.70, 0.00, 0.95)])   # H (hydroxyl)
    p = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cpu")
    def relax(a, spin=1):
        a = a.copy(); a.info.update({"charge": 0, "spin": spin})
        a.calc = FAIRChemCalculator(p, task_name="omol")
        LBFGS(a, logfile=None).run(fmax=0.05, steps=200)
        return a
    clu0 = relax(tio2)
    ad0 = relax(hcooh)
    n_clu = len(clu0)
    comp = relax(clu0 + ad0)
    # IMPORTANT: for a valid Δρ the fragments must be at their FROZEN in-complex
    # geometry (same grid positions), otherwise ∫Δρ ≠ 0 and Δρ is meaningless.
    clu = comp[:n_clu]           # cluster atoms, frozen at complex coords
    ad = comp[n_clu:]            # adsorbate atoms, frozen at complex coords
    write(os.path.join(SDIR, "demo_cluster_Ti2O4.xyz"), clu)
    write(os.path.join(SDIR, "demo_adsorbate_HCOOH.xyz"), ad)
    write(os.path.join(SDIR, "demo_complex_HCOOH_Ti2O4.xyz"), comp)
    return clu, ad, comp

# ---------------------------------------------------------------------------
# 2. PySCF PBE densities on a common grid  -> Δρ cube
# ---------------------------------------------------------------------------
def scf_density(atoms, label, basis="def2-svp"):
    from pyscf import gto, dft, scf
    mol = gto.M(atom=[(s, tuple(p)) for s, p in
                zip(atoms.get_chemical_symbols(), atoms.get_positions())],
                basis=basis, charge=0, spin=0, verbose=1)
    mf = dft.RKS(mol); mf.xc = "pbe"; mf = mf.density_fit()
    mf = scf.addons.smearing_(mf, sigma=0.01, method="fermi")
    mf.conv_tol = 1e-6; mf.max_cycle = 100
    t = time.time(); e = mf.kernel(); dt = time.time() - t
    print(f"  [{label}] E={e:.4f} Ha converged={mf.converged} "
          f"{dt:.0f}s ({mol.nao} bfn)", flush=True)
    return mol, mf

def density_on_coords(mol, mf, coords_bohr):
    from pyscf import dft as pdft
    ao = pdft.numint.eval_ao(mol, coords_bohr)     # coords already in Bohr
    return pdft.numint.eval_rho(mol, ao, mf.make_rdm1())

# ---------------------------------------------------------------------------
# 3. Bader charges on the complex density
# ---------------------------------------------------------------------------
def bader_charges(cube_path):
    """Run the Henkelman-group `bader` binary (the standard grid Bader tool used
    with VASP) on a Gaussian cube; parse ACF.dat -> electrons per atom."""
    import subprocess
    bexe = os.path.join(HERE, "..", "tools", "bader")
    workdir = os.path.dirname(cube_path)
    subprocess.run([bexe, os.path.basename(cube_path)], cwd=workdir,
                   check=True, capture_output=True)
    elec = []
    with open(os.path.join(workdir, "ACF.dat")) as f:
        for line in f:
            parts = line.split()
            # data rows: index X Y Z CHARGE MINDIST VOL  (>=7 cols, int index)
            if len(parts) >= 7 and parts[0].isdigit():
                elec.append(float(parts[4]))       # CHARGE column = electrons
    return np.asarray(elec, dtype=float)

def main():
    print("[demo] building + UMA-relaxing HCOOH/Ti2O4 model...", flush=True)
    clu, ad, comp = build_and_relax()
    print(f"  cluster {clu.get_chemical_formula()}, adsorbate "
          f"{ad.get_chemical_formula()}, complex {comp.get_chemical_formula()}",
          flush=True)

    print("[demo] PySCF PBE densities...", flush=True)
    from pyscf.tools import cubegen
    mc, mfc = scf_density(comp, "complex")
    ms, mfs = scf_density(clu, "cluster")
    ma, mfa = scf_density(ad, "adsorbate")

    # ONE common grid, defined by PySCF's Cube from the complex (correct header)
    cube = cubegen.Cube(mc, nx=96, ny=96, nz=96, margin=4.0)
    coords = cube.get_coords()            # Bohr
    dv = abs(np.linalg.det(cube.box)) / (cube.nx * cube.ny * cube.nz)  # Bohr^3/voxel
    rho_c = density_on_coords(mc, mfc, coords)
    rho_s = density_on_coords(ms, mfs, coords)
    rho_a = density_on_coords(ma, mfa, coords)
    drho = rho_c - rho_s - rho_a          # differential charge density
    print(f"  ∫ρ_complex dV = {rho_c.sum()*dv:.2f} e (N_elec={sum(comp.get_atomic_numbers())})",
          flush=True)
    print(f"  ∫Δρ dV = {drho.sum()*dv:+.3f} e (should be ~0)", flush=True)

    cube.write(drho.reshape(cube.nx, cube.ny, cube.nz),
               os.path.join(RDIR, "demo_delta_rho.cube"), comment="differential charge density")
    cube.write(rho_c.reshape(cube.nx, cube.ny, cube.nz),
               os.path.join(RDIR, "demo_rho_complex.cube"), comment="total density")
    print("  wrote demo_delta_rho.cube, demo_rho_complex.cube", flush=True)

    # Bader on complex total density
    res = {"model_system": "HCOOH on Ti2O4 (closed-shell demo)",
           "engine": "PySCF PBE/def2-SVP (RI-J, Fermi smearing)",
           "note": "tractable analog of RSV-TiOx; full systems via scripts/05_*"}
    try:
        elec = bader_charges(os.path.join(RDIR, "demo_rho_complex.cube"))
        Zval = comp.get_atomic_numbers().astype(float)
        q = Zval - elec                          # net charge = Z - Bader electrons
        n_clu = len(clu)
        q_mol = float(q[n_clu:].sum()); q_clu = float(q[:n_clu].sum())
        res["bader_total_electrons"] = float(elec.sum())
        res["bader_net_charge_adsorbate_e"] = q_mol
        res["bader_net_charge_cluster_e"] = q_clu
        res["charge_transfer_e"] = abs(q_mol)
        res["direction"] = ("cluster->adsorbate" if q_mol < 0
                            else "adsorbate->cluster")
        print(f"  Bader: adsorbate net charge = {q_mol:+.3f} e, "
              f"cluster = {q_clu:+.3f} e  ({res['direction']})", flush=True)
    except Exception as e:
        import traceback; traceback.print_exc()
        res["bader_error"] = str(e)
        print(f"  Bader step failed ({e}); Δρ cube still produced.", flush=True)

    json.dump(res, open(os.path.join(RDIR, "charge_density_bader_demo.json"), "w"),
              indent=2)
    print("-> results/charge_density_bader_demo.json", flush=True)

if __name__ == "__main__":
    main()
