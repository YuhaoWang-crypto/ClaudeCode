"""
05_full_systems_inputs.py
Production-ready drivers to compute the differential charge density (Δρ) and
Bader charges for the FULL project systems -- the same working pipeline validated
in 04_charge_density_bader.py, pointed at the real geometries. These do NOT finish
on a 4-core CPU; run them on a GPU node / HPC / workstation.

  (A) RSV-TiOx complexes  -> PySCF (Gaussian-basis DFT, like the report's Gaussian
                             cluster work). Call run_bader_pipeline() on the
                             UMA-relaxed complexes in ../structures/.
  (B) C13H10N4S / Fe(110) -> Quantum ESPRESSO (plane-wave PBE/PAW, the open
                             equivalent of the report's VASP). This script writes
                             ready-to-run pw.x inputs; post-process with pp.x +
                             the Henkelman `bader` code (see printed instructions).
"""
import os, subprocess, json
import numpy as np
from ase.io import read, write

HERE = os.path.dirname(os.path.abspath(__file__))
SDIR = os.path.join(HERE, "..", "structures")
RDIR = os.path.join(HERE, "..", "results")
BADER = os.path.join(HERE, "..", "tools", "bader")


# ============================================================ (A) RSV-TiOx =====
def run_bader_pipeline(complex_xyz, n_cluster_atoms, tag, basis="def2-svp",
                       spin=0, nx=120, ny=120, nz=120):
    """Δρ + Bader for a molecule/cluster complex with PySCF PBE. Fragments are
    frozen at their in-complex geometry (required for a valid Δρ)."""
    from pyscf import gto, dft, scf
    from pyscf.tools import cubegen

    comp = read(complex_xyz)
    clu = comp[:n_cluster_atoms]
    ad = comp[n_cluster_atoms:]

    def density(atoms, label):
        mol = gto.M(atom=[(s, tuple(p)) for s, p in
                    zip(atoms.get_chemical_symbols(), atoms.get_positions())],
                    basis=basis, charge=0, spin=spin, verbose=2)
        mf = dft.RKS(mol) if spin == 0 else dft.UKS(mol)
        mf.xc = "pbe"; mf = mf.density_fit()
        mf = scf.addons.smearing_(mf, sigma=0.01, method="fermi")
        mf.conv_tol = 1e-6; mf.max_cycle = 150
        mf.kernel()
        return mol, mf

    mc, mfc = density(comp, "complex")
    ms, mfs = density(clu, "cluster")
    ma, mfa = density(ad, "adsorbate")

    cube = cubegen.Cube(mc, nx=nx, ny=ny, nz=nz, margin=4.0)
    coords = cube.get_coords()
    from pyscf import dft as pdft
    def rho(mol, mf):
        ao = pdft.numint.eval_ao(mol, coords)
        return pdft.numint.eval_rho(mol, ao, mf.make_rdm1())
    rc, rs, ra = rho(mc, mfc), rho(ms, mfs), rho(ma, mfa)
    drho = rc - rs - ra
    cube.write(drho.reshape(nx, ny, nz), os.path.join(RDIR, f"{tag}_delta_rho.cube"),
               comment="differential charge density")
    cube.write(rc.reshape(nx, ny, nz), os.path.join(RDIR, f"{tag}_rho.cube"),
               comment="total density")

    subprocess.run([BADER, f"{tag}_rho.cube"], cwd=RDIR, check=True,
                   capture_output=True)
    elec = []
    for line in open(os.path.join(RDIR, "ACF.dat")):
        p = line.split()
        if len(p) >= 7 and p[0].isdigit():
            elec.append(float(p[4]))
    elec = np.array(elec)
    q = comp.get_atomic_numbers() - elec
    out = {"tag": tag, "charge_transfer_e": float(abs(q[n_cluster_atoms:].sum())),
           "adsorbate_net_e": float(q[n_cluster_atoms:].sum())}
    json.dump(out, open(os.path.join(RDIR, f"{tag}_bader.json"), "w"), indent=2)
    print(f"[{tag}] charge transfer = {out['charge_transfer_e']:.3f} e")
    return out


# ==================================================== (B) Fe(110) via QE =======
def write_qe_inputs():
    """Write pw.x SCF inputs for slab / molecule / complex (plane-wave PBE/PAW).
    You must supply PAW/USPP pseudopotentials (e.g. SSSP efficiency library from
    https://www.materialscloud.org/discover/sssp ) and set the pseudo_dir."""
    from ase.io.espresso import write_espresso_in
    qedir = os.path.join(HERE, "..", "qe_inputs"); os.makedirs(qedir, exist_ok=True)
    pseudos = {"Fe": "Fe.pbe-spn-kjpaw_psl.1.0.0.UPF",
               "C": "C.pbe-n-kjpaw_psl.1.0.0.UPF", "H": "H.pbe-kjpaw_psl.1.0.0.UPF",
               "N": "N.pbe-n-kjpaw_psl.1.0.0.UPF", "S": "S.pbe-n-kjpaw_psl.1.0.0.UPF"}
    inp = dict(ecutwfc=50, ecutrho=400, occupations="smearing", smearing="mv",
               degauss=0.02, nspin=2, starting_magnetization=0.3,
               conv_thr=1e-6, mixing_beta=0.3)
    systems = {"slab": "fe110_slab_relaxed_deep.xyz",
               "complex": "relaxed_fe110_best.xyz",
               "mol": "relaxed_inhibitor_in_slabcell.xyz"}
    for tag, fn in systems.items():
        atoms = read(os.path.join(SDIR, fn))
        os.makedirs(os.path.join(qedir, tag), exist_ok=True)
        pwi = os.path.join(qedir, tag, "espresso.pwi")
        idata = {"control": {"calculation": "scf", "prefix": tag,
                             "pseudo_dir": "../pseudos", "tprnfor": True},
                 "system": inp, "electrons": {"conv_thr": 1e-6, "mixing_beta": 0.3}}
        with open(pwi, "w") as fd:
            write_espresso_in(fd, atoms, input_data=idata,
                              pseudopotentials=pseudos, kpts=(3, 3, 1))
        print(f"  wrote qe_inputs/{tag}/espresso.pwi")
    # instructions
    readme = os.path.join(qedir, "README.md")
    open(readme, "w").write(
        "# Fe(110) Δρ + Bader with Quantum ESPRESSO\n\n"
        "1. Get PAW pseudos (SSSP efficiency) into qe_inputs/pseudos/.\n"
        "2. Run SCF for slab, mol, complex:\n"
        "   `mpirun pw.x -in espresso.pwi > scf.out` (each folder)\n"
        "3. Export each charge density to cube with pp.x (plot_num=0).\n"
        "4. Δρ = ρ(complex) − ρ(slab) − ρ(mol): subtract the three cubes on the\n"
        "   same FFT grid (use `pp.x` `plot_num=0` with identical `nr1,nr2,nr3`,\n"
        "   or ASE/pymatgen to combine).\n"
        "5. Bader: `bader complex_CHGCAR` (Henkelman code) → ACF.dat → sum the\n"
        "   molecule atoms vs slab atoms for the net charge transfer.\n")
    print(f"  wrote {readme}")


if __name__ == "__main__":
    print("(A) To compute the full RSV-TiOx Δρ+Bader (needs HPC/GPU), run e.g.:")
    print("    run_bader_pipeline('../structures/relaxed_complex_n-type_TiO2-RSV.xyz',"
          " n_cluster_atoms=31, tag='rsv_tio2')")
    print("    run_bader_pipeline('../structures/relaxed_complex_p-type_TiO-RSV.xyz',"
          " n_cluster_atoms=25, tag='rsv_tio')")
    print("(B) Writing Quantum ESPRESSO inputs for Fe(110)...")
    write_qe_inputs()
