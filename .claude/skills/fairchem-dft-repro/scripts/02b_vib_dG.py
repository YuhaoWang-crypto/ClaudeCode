"""
02b_vib_dG.py
Approximate ΔG_bind for the RSV-TiOx complexes via harmonic vibrational free
energies from UMA (omol) forces, on the ALREADY-relaxed geometries.

  F(species) = E_elec + ZPE + F_vib(T)         (harmonic oscillator, ASE)
  ΔG_bind ≈ ΔF_bind = F(complex) - F(cluster) - F(RSV)
  Ka      = exp(-ΔG_bind / RT)

We use the HARMONIC-OSCILLATOR free energy (ASE HarmonicThermo) for every
species, i.e. we neglect the gas-phase translational/rotational partition
functions and take G ≈ F (the PV term is negligible for these bound/condensed
systems). This is the standard treatment for adsorption/binding free energies.
It slightly UNDERESTIMATES the entropic penalty of association (the free ligand
really does have trans/rot entropy), so |ΔG_bind| here is an upper bound.

CAVEAT: MLIP-derived Hessians are noisier than the energies; floppy 90-atom
complexes produce several small imaginary modes (reported below), which are
discarded. Treat ΔG_bind as an *estimate*; ΔE_bind and the relative trend
(script 02) are the robust quantities, as the report itself emphasises.
"""
import os, json, time, glob
import numpy as np
from ase.io import read
from ase.vibrations import Vibrations
from ase.thermochemistry import HarmonicThermo
from fairchem.core import pretrained_mlip, FAIRChemCalculator

HERE = os.path.dirname(os.path.abspath(__file__))
SDIR = os.path.join(HERE, "..", "structures")
RDIR = os.path.join(HERE, "..", "results")
T, KB = 298.15, 8.617333262e-5
RT = KB * T
IMAG_CUT = 1.5e-3   # eV (~12 cm^-1); modes below this magnitude are dropped

predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cpu")

def F_harmonic(fname, label):
    atoms = read(os.path.join(SDIR, fname))
    atoms.set_pbc(False)
    atoms.info.update({"charge": 0, "spin": 1})
    atoms.calc = FAIRChemCalculator(predictor, task_name="omol")
    e_pot = float(atoms.get_potential_energy())
    vibdir = os.path.join(RDIR, f"vib_{label}")
    # ASE Vibrations resumes from any cached displacements in vibdir
    vib = Vibrations(atoms, name=vibdir, delta=0.015)
    t0 = time.time()
    vib.run()                       # skips already-cached displacements
    energies = vib.get_energies()
    n_imag = int(sum(1 for e in energies if abs(e.imag) > 1e-9))
    real = np.array([e.real for e in energies
                     if abs(e.imag) < 1e-9 and e.real > IMAG_CUT])
    thermo = HarmonicThermo(vib_energies=real, potentialenergy=e_pot)
    F = float(thermo.get_helmholtz_energy(temperature=T, verbose=False))
    dt = time.time() - t0
    print(f"  [{label}] E={e_pot:.3f}  F={F:.3f} eV  "
          f"({len(real)} real modes kept, {n_imag} imag dropped) "
          f"{dt:.0f}s ({len(atoms)} atoms)", flush=True)
    return {"E_elec": e_pot, "F": F, "n_real": len(real), "n_imag": n_imag,
            "seconds": round(dt)}

def ka(dG):
    x = -dG / RT
    return {"dG_bind_eV": float(dG), "ln_Ka": float(x),
            "log10_Ka": float(x / np.log(10))}

def main():
    print("Computing harmonic free energies (slow: full Hessians on CPU)...",
          flush=True)
    species = [("relaxed_rsv.xyz", "rsv"),
               ("relaxed_tio2_cluster_ntype.xyz", "tio2"),
               ("relaxed_tio_cluster_ptype.xyz", "tio"),
               ("relaxed_complex_n-type_TiO2-RSV.xyz", "comp_n"),
               ("relaxed_complex_p-type_TiO-RSV.xyz", "comp_p")]
    F = {}
    partial = os.path.join(RDIR, "rsv_tiox_dG_partial.json")
    for fname, label in species:
        F[label] = F_harmonic(fname, label)
        json.dump(F, open(partial, "w"), indent=2)   # checkpoint after each

    dG_n = F["comp_n"]["F"] - F["tio2"]["F"] - F["rsv"]["F"]
    dG_p = F["comp_p"]["F"] - F["tio"]["F"] - F["rsv"]["F"]
    out = {"model": "uma-s-1p1", "task": "omol",
           "method": "harmonic-oscillator ΔF≈ΔG (trans/rot neglected)",
           "T_K": T, "RT_eV": RT, "F_components": F,
           "n-type_TiO2-RSV": {**ka(dG_n), "report_dG": -1.808251890,
                               "report_Ka": "3.8e30"},
           "p-type_TiO-RSV": {**ka(dG_p), "report_dG": -3.372142689,
                              "report_Ka": "1.0e57"},
           "relative": {"dG_p_minus_n_eV": float(dG_p - dG_n),
                        "stronger": "p-type" if dG_p < dG_n else "n-type",
                        "report_stronger": "p-type"}}
    json.dump(out, open(os.path.join(RDIR, "rsv_tiox_dG.json"), "w"), indent=2)
    print(f"\nΔG_bind(approx): n-type={dG_n:.3f} eV, p-type={dG_p:.3f} eV",
          flush=True)
    print("  -> results/rsv_tiox_dG.json", flush=True)

if __name__ == "__main__":
    main()
