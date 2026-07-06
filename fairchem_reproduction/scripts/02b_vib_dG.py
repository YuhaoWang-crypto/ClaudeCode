"""
02b_vib_dG.py
Approximate ΔG_bind for the RSV-TiOx complexes via harmonic vibrational free
energies from UMA (omol) forces, on the ALREADY-relaxed geometries.

  G(species) = E_elec + ZPE + G_vib(T) + G_trans+rot (ideal gas)
  ΔG_bind    = G(complex) - G(cluster) - G(RSV)
  Ka         = exp(-ΔG_bind / RT)

CAVEAT: MLIP-derived Hessians are noisier than energies; low-frequency and
spurious-imaginary modes make the absolute entropy approximate. Treat ΔG_bind
here as an *estimate*; ΔE_bind and the relative trend (script 02) are the robust
quantities, as the original report itself emphasises.
"""
import os, json, time
import numpy as np
from ase.io import read
from ase.vibrations import Vibrations
from ase.thermochemistry import IdealGasThermo
from fairchem.core import pretrained_mlip, FAIRChemCalculator

HERE = os.path.dirname(os.path.abspath(__file__))
SDIR = os.path.join(HERE, "..", "structures")
RDIR = os.path.join(HERE, "..", "results")
T, KB = 298.15, 8.617333262e-5
RT = KB * T

predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cpu")

def G_harmonic(fname, label):
    atoms = read(os.path.join(SDIR, fname))
    atoms.set_pbc(False)
    atoms.info.update({"charge": 0, "spin": 1})
    atoms.calc = FAIRChemCalculator(predictor, task_name="omol")
    e_pot = float(atoms.get_potential_energy())
    vibdir = os.path.join(RDIR, f"vib_{label}")
    for f in (list(__import__('glob').glob(vibdir + "*"))):
        try: os.remove(f)
        except OSError: pass
    vib = Vibrations(atoms, name=vibdir, delta=0.015)
    t0 = time.time()
    vib.run()
    energies = vib.get_energies()
    n_imag = int(sum(1 for e in energies if abs(e.imag) > 1e-6))
    real = np.array([e.real for e in energies if abs(e.imag) < 1e-6 and e.real > 6e-3])
    vib.clean()
    thermo = IdealGasThermo(vib_energies=real, geometry="nonlinear",
                            potentialenergy=e_pot, atoms=atoms,
                            symmetrynumber=1, spin=0)
    G = float(thermo.get_gibbs_energy(temperature=T, pressure=101325.0, verbose=False))
    print(f"  [{label}] E={e_pot:.3f}  G={G:.3f} eV  ({len(real)} real modes, "
          f"{n_imag} imag) {time.time()-t0:.0f}s  ({len(atoms)} atoms)")
    return {"E_elec": e_pot, "G": G, "n_real": len(real), "n_imag": n_imag}

def ka(dG):
    x = -dG / RT
    return {"dG_bind_eV": float(dG), "ln_Ka": float(x), "log10_Ka": float(x/np.log(10))}

def main():
    print("Computing harmonic Gibbs energies (this is the slow, ~1 hr step)...")
    G = {}
    G["rsv"] = G_harmonic("relaxed_rsv.xyz", "rsv")
    G["tio2"] = G_harmonic("relaxed_tio2_cluster_ntype.xyz", "tio2")
    G["tio"] = G_harmonic("relaxed_tio_cluster_ptype.xyz", "tio")
    G["comp_n"] = G_harmonic("relaxed_complex_n-type_TiO2-RSV.xyz", "comp_n")
    G["comp_p"] = G_harmonic("relaxed_complex_p-type_TiO-RSV.xyz", "comp_p")

    dG_n = G["comp_n"]["G"] - G["tio2"]["G"] - G["rsv"]["G"]
    dG_p = G["comp_p"]["G"] - G["tio"]["G"] - G["rsv"]["G"]
    out = {"model": "uma-s-1p1", "task": "omol", "method": "harmonic ideal-gas ΔG",
           "T_K": T, "RT_eV": RT, "G_components": G,
           "n-type_TiO2-RSV": {**ka(dG_n), "report_dG": -1.808251890,
                               "report_Ka": "3.8e30"},
           "p-type_TiO-RSV": {**ka(dG_p), "report_dG": -3.372142689,
                              "report_Ka": "1.0e57"},
           "relative": {"dG_p_minus_n_eV": float(dG_p - dG_n),
                        "stronger": "p-type" if dG_p < dG_n else "n-type"}}
    with open(os.path.join(RDIR, "rsv_tiox_dG.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\nΔG_bind: n-type={dG_n:.3f} eV, p-type={dG_p:.3f} eV")
    print(f"  -> results/rsv_tiox_dG.json")

if __name__ == "__main__":
    main()
