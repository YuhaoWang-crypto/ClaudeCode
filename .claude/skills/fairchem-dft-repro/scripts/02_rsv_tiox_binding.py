"""
02_rsv_tiox_binding.py
Reproduce the RSV (rosuvastatin) - TiOx binding deliverable with fairchem UMA
(OMol task, closed-shell singlets).

  DeltaE_bind = E(complex) - E(cluster) - E(RSV)          [electronic]
  DeltaG_bind ~ DeltaE_bind + DeltaG_thermal(harmonic)    [optional, DO_VIB=1]
  Ka          = exp(-DeltaG_bind / RT)     (report's relationship)

Original report (Gaussian, cluster model):
  n-type TiO2-RSV : dGbind = -1.808 eV,  Ka = 3.8e30
  p-type TiO -RSV : dGbind = -3.372 eV,  Ka = 1.0e57
The report notes Ka magnitudes are only meaningful in RELATIVE terms.

UMA OMol is trained on wB97M-V DFT (not B3LYP), so absolute numbers differ from
the report; the physically meaningful output is the RELATIVE binding strength
(p-type >> n-type) and the sign/scale of the binding energy.
"""
import os, json, time
import numpy as np
from ase.io import read, write
from ase.optimize import LBFGS
from fairchem.core import pretrained_mlip, FAIRChemCalculator

HERE = os.path.dirname(os.path.abspath(__file__))
SDIR = os.path.join(HERE, "..", "structures")
RDIR = os.path.join(HERE, "..", "results")
os.makedirs(RDIR, exist_ok=True)

FMAX = 0.05
MAXSTEPS = 300
KB = 8.617333262e-5   # eV/K
T = 298.15
RT = KB * T
DO_VIB = os.environ.get("DO_VIB", "0") == "1"

print("Loading UMA (uma-s-1p1, omol task)...")
predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cpu")

def omol_calc():
    return FAIRChemCalculator(predictor, task_name="omol")

def set_cs(atoms, charge=0, spin=1):
    atoms.info.update({"charge": charge, "spin": spin})  # spin = multiplicity
    return atoms

def relax(atoms, label, charge=0, spin=1):
    atoms = atoms.copy(); set_cs(atoms, charge, spin)
    atoms.calc = omol_calc()
    t0 = time.time()
    opt = LBFGS(atoms, logfile=os.path.join(RDIR, f"opt_{label}.log"))
    opt.run(fmax=FMAX, steps=MAXSTEPS)
    e = float(atoms.get_potential_energy())
    fmax = float(np.sqrt((atoms.get_forces() ** 2).sum(axis=1)).max())
    write(os.path.join(SDIR, f"relaxed_{label}.xyz"), atoms)
    print(f"  [{label}] E={e:.3f} eV steps={opt.get_number_of_steps()} "
          f"fmax={fmax:.3f} {time.time()-t0:.0f}s ({len(atoms)} atoms)")
    return atoms, e

def make_complex(rsv, cluster, ti_o_dist=2.1):
    """Dock RSV onto the cluster: put an RSV carboxyl/hydroxyl O near an
    exposed surface Ti atom (Ti is the Lewis-acidic binding site)."""
    rsv = rsv.copy(); cluster = cluster.copy()
    # candidate binding O atoms on RSV (carboxylic acid + hydroxyls)
    o_idx = [i for i, s in enumerate(rsv.symbols) if s == "O"]
    # pick the most exposed Ti of the cluster (max distance from cluster COM)
    ti_idx = [i for i, s in enumerate(cluster.symbols) if s == "Ti"]
    ccom = cluster.get_center_of_mass()
    ti = ti_idx[int(np.argmax([np.linalg.norm(cluster.positions[i] - ccom)
                               for i in ti_idx]))]
    ti_pos = cluster.positions[ti]
    outward = (ti_pos - ccom); outward /= (np.linalg.norm(outward) + 1e-9)
    # translate RSV so its first carboxyl O sits ti_o_dist beyond Ti, outward
    target = ti_pos + outward * ti_o_dist
    o_anchor = o_idx[0]
    rsv.translate(target - rsv.positions[o_anchor])
    # nudge whole RSV further out along outward to avoid clashes
    rsv.translate(outward * 1.5)
    comp = cluster + rsv
    return comp

def harmonic_G(atoms, label, charge=0, spin=1, linear=False):
    """Harmonic Gibbs free energy (eV) at T via ASE Vibrations + IdealGasThermo.
    Expensive: 6N force evaluations. Used only when DO_VIB=1."""
    from ase.vibrations import Vibrations
    from ase.thermochemistry import IdealGasThermo
    atoms = atoms.copy(); set_cs(atoms, charge, spin)
    atoms.calc = omol_calc()
    e_pot = float(atoms.get_potential_energy())
    vibdir = os.path.join(RDIR, f"vib_{label}")
    vib = Vibrations(atoms, name=vibdir)
    t0 = time.time()
    vib.run()
    energies = vib.get_energies()  # complex, includes imaginary as 0j
    vib.clean()
    real = np.array([e.real for e in energies if abs(e.imag) < 1e-6 and e.real > 1e-3])
    geometry = "linear" if linear else "nonlinear"
    thermo = IdealGasThermo(vib_energies=real, geometry=geometry,
                            potentialenergy=e_pot, atoms=atoms,
                            symmetrynumber=1, spin=(spin - 1) / 2.0)
    G = thermo.get_gibbs_energy(temperature=T, pressure=101325.0, verbose=False)
    print(f"  [vib {label}] G={G:.3f} eV, {len(real)} real modes, {time.time()-t0:.0f}s")
    return float(G)

def ka_from_dG(dG):
    # cap the exponent for reporting (values are astronomically large, as the
    # report itself notes)
    x = -dG / RT
    return {"exponent_base_e": x, "log10_Ka": x / np.log(10)}

def run_system(name, cluster_file, report_dG, report_Ka):
    print(f"\n=== {name} ===")
    res = {"report_dGbind_eV": report_dG, "report_Ka": report_Ka}
    rsv = read(os.path.join(SDIR, "rosuvastatin_RSV.xyz"))
    clu = read(os.path.join(SDIR, cluster_file))
    rsv_r, E_rsv = relax(rsv, "rsv")
    clu_lbl = os.path.splitext(cluster_file)[0]
    clu_r, E_clu = relax(clu, clu_lbl)
    comp = make_complex(rsv_r, clu_r)
    comp_r, E_comp = relax(comp, f"complex_{name}")
    dE = E_comp - E_clu - E_rsv
    res.update({"E_rsv": E_rsv, "E_cluster": E_clu, "E_complex": E_comp,
                "dE_bind_eV": float(dE)})
    res["Ka_from_dE"] = ka_from_dG(dE)
    print(f"  dE_bind(UMA) = {dE:.3f} eV | report dGbind = {report_dG} eV")

    if DO_VIB:
        G_rsv = harmonic_G(rsv_r, "rsv")
        G_clu = harmonic_G(clu_r, clu_lbl)
        G_comp = harmonic_G(comp_r, f"complex_{name}")
        dG = G_comp - G_clu - G_rsv
        res.update({"G_rsv": G_rsv, "G_cluster": G_clu, "G_complex": G_comp,
                    "dG_bind_eV": float(dG), "Ka_from_dG": ka_from_dG(dG)})
        print(f"  dG_bind(UMA, harmonic) = {dG:.3f} eV")
    return res

def main():
    out = {"model": "uma-s-1p1", "task": "omol", "T_K": T, "RT_eV": RT,
           "fmax": FMAX, "do_vib": DO_VIB, "systems": {}}
    out["systems"]["n-type_TiO2-RSV"] = run_system(
        "n-type_TiO2-RSV", "tio2_cluster_ntype.xyz", -1.808251890, "3.8e30")
    out["systems"]["p-type_TiO-RSV"] = run_system(
        "p-type_TiO-RSV", "tio_cluster_ptype.xyz", -3.372142689, "1.0e57")

    # relative comparison (the physically meaningful part)
    dEn = out["systems"]["n-type_TiO2-RSV"]["dE_bind_eV"]
    dEp = out["systems"]["p-type_TiO-RSV"]["dE_bind_eV"]
    out["relative"] = {"dE_p_minus_n_eV": dEp - dEn,
                       "stronger_binder": "p-type_TiO" if dEp < dEn else "n-type_TiO2",
                       "report_says_stronger": "p-type_TiO"}
    with open(os.path.join(RDIR, "rsv_tiox_binding.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print("\n-> results/rsv_tiox_binding.json")
    print(f"Relative: p-type minus n-type dE = {dEp-dEn:.3f} eV "
          f"(report: p-type binds stronger)")

if __name__ == "__main__":
    main()
