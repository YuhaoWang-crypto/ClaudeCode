"""
06b_pcet_oxidation.py
Proton-coupled (PCET) oxidation potentials for DA / AA / UA.

The 1-electron radical-cation route (06_*) mis-orders these molecules because
their real oxidation is proton-coupled. Here we compute the physically correct
first oxidation as a 1H+/1e- PCET forming the NEUTRAL radical (antioxidant
mechanism):

    RedH(aq)  ->  Red•(aq) + H+(aq) + e-(vac)

    E°(vs SHE, pH 0) = [E(Red•,aq) - E(RedH,aq) + G(H+,aq)] - E_abs(SHE)
    E°(pH 7)         = E°(pH 0) - 0.05916 * 7          (Nernst, 1H+/1e-)

The molecule that most easily loses H (strongest antioxidant) oxidises at the
lowest potential -> ascorbic acid should now come out easiest, matching experiment.

For each analyte we remove every O-H / N-H hydrogen, relax the neutral radical
with UMA, evaluate a solvated (ddCOSMO water) B3LYP/6-31G* energy, and take the
MOST STABLE radical as the oxidation product (no hand-picking of the site).

Constants (Camaioni-Schwerdtfeger / Tissandier proton, IUPAC SHE):
    G(H+,aq)  = -11.72 eV   (= -270.3 kcal/mol: gas -6.28 + 1atm->1M 1.89 - 265.9 solv)
    E_abs(SHE) = 4.44 V     (shifting to 4.28 V moves ALL values by +0.16 V uniformly)
"""
import os, json, time
import numpy as np
from ase.io import read, write
from ase.optimize import LBFGS

HERE = os.path.dirname(os.path.abspath(__file__))
SDIR = os.path.join(HERE, "..", "structures")
RDIR = os.path.join(HERE, "..", "results")
HARTREE_EV = 27.211386245988
G_HPLUS_AQ = -11.72     # eV, absolute aqueous proton free energy
E_ABS_SHE = 4.44        # V
PH = 7.0

ANALYTES = ["dopamine", "ascorbic_acid", "uric_acid"]

def solvated_energy(atoms, charge, spin2s):
    from pyscf import gto, dft
    from pyscf.solvent import ddcosmo
    mol = gto.M(atom=[(s, tuple(p)) for s, p in
                zip(atoms.get_chemical_symbols(), atoms.get_positions())],
                basis="6-31g*", charge=charge, spin=spin2s, verbose=0)
    mf = dft.RKS(mol) if spin2s == 0 else dft.UKS(mol)
    mf.xc = "b3lyp"
    mf = ddcosmo.ddcosmo_for_scf(mf); mf.with_solvent.eps = 78.3553
    mf.conv_tol = 1e-7
    return float(mf.kernel()) * HARTREE_EV      # eV

def on_hydrogens(atoms):
    """indices of H atoms bonded to O or N (the redox-active, ionisable protons)."""
    pos = atoms.get_positions(); sym = atoms.get_chemical_symbols()
    hs = []
    for i, s in enumerate(sym):
        if s != "H":
            continue
        d = np.linalg.norm(pos - pos[i], axis=1)
        for j, sj in enumerate(sym):
            if j != i and sj in ("O", "N") and d[j] < 1.30:
                hs.append(i); break
    return hs

def best_radical(neutral, predictor, tag):
    """Remove each O-H/N-H, UMA-relax the neutral radical, return the lowest
    solvated-energy radical (E in eV) and its structure."""
    from fairchem.core import FAIRChemCalculator
    best = None
    for h in on_hydrogens(neutral):
        rad = neutral.copy(); del rad[h]
        rad.info.update({"charge": 0, "spin": 2})    # neutral doublet radical
        rad.calc = FAIRChemCalculator(predictor, task_name="omol")
        LBFGS(rad, logfile=None).run(fmax=0.05, steps=300)
        E = solvated_energy(rad, 0, 1)
        if best is None or E < best[0]:
            best = (E, rad, h)
    write(os.path.join(SDIR, f"radical_{tag}.xyz"), best[1])
    return best

def main():
    from fairchem.core import pretrained_mlip
    predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cpu")
    t0 = time.time()
    out = {"method": "1H+/1e- PCET (neutral radical), UMA geom + PySCF "
                     "B3LYP/6-31G* ddCOSMO(water)",
           "G_Hplus_aq_eV": G_HPLUS_AQ, "E_abs_SHE_V": E_ABS_SHE, "pH": PH,
           "results": {}}
    for name in ANALYTES:
        neu = read(os.path.join(SDIR, f"analyte_{name}.xyz"))
        E_neu = solvated_energy(neu, 0, 0)
        E_rad, rad, hsite = best_radical(neu, predictor, name)
        E0 = (E_rad - E_neu) + G_HPLUS_AQ - E_ABS_SHE     # V vs SHE at pH 0
        E7 = E0 - 0.05916 * PH
        out["results"][name] = {"E_ox_pH0_vs_SHE_V": E0, "E_ox_pH7_vs_SHE_V": E7,
                                "n_radical_sites_tried": len(on_hydrogens(neu))}
        print(f"[{name}] E_ox(pH0)={E0:+.2f} V  E_ox(pH7)={E7:+.2f} V", flush=True)
    order = sorted(out["results"], key=lambda k: out["results"][k]["E_ox_pH7_vs_SHE_V"])
    out["oxidation_order_easiest_first"] = order
    out["experimental_order_pH7"] = ["ascorbic_acid", "dopamine", "uric_acid"]
    json.dump(out, open(os.path.join(RDIR, "pcet_oxidation.json"), "w"), indent=2)
    print(f"\nPredicted easiest->hardest (pH7): {order}")
    print(f"Experimental (pH7):               {out['experimental_order_pH7']}")
    print(f"[{time.time()-t0:.0f}s] -> results/pcet_oxidation.json")

if __name__ == "__main__":
    main()
