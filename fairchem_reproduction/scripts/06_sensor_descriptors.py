"""
06_sensor_descriptors.py
Sensor-material screening demo: predict a MOLECULAR DESCRIPTOR FINGERPRINT for a
panel of analytes, the quantities that actually drive an electrochemical sensor's
signal (not just Ka).

Panel: dopamine (DA), ascorbic acid (AA), uric acid (UA) -- the classic
electrochemical-sensing trio whose oxidation peaks overlap at a bare electrode;
a good sensor material must SEPARATE them.

Per analyte we compute (UMA for geometries, PySCF for electronics):
  * HOMO, LUMO, gap                         (B3LYP/6-31G*)
  * global reactivity: chi, eta, omega       (electronegativity/hardness/electrophilicity)
  * dipole moment
  * adiabatic ionization potential (gas)     IP = E(cation) - E(neutral)
  * 1-electron oxidation potential vs SHE    E_ox = [E(cat,solv)-E(neu,solv)] - 4.44 V
                                             (ddCOSMO water)

CAVEAT (honest): DA/AA/UA oxidation in water is proton-coupled (PCET), so the
absolute experimental peak potentials also depend on pH; here we compute the pure
1-electron oxidation as a transferable REACTIVITY descriptor. The RELATIVE
ordering and separation -- what matters for sensor selectivity -- is the point,
not the absolute volt.
"""
import os, json, time
import numpy as np
from ase import Atoms
from ase.io import write
from ase.optimize import LBFGS

HERE = os.path.dirname(os.path.abspath(__file__))
SDIR = os.path.join(HERE, "..", "structures")
RDIR = os.path.join(HERE, "..", "results")
HARTREE_EV = 27.211386245988
SHE_ABS = 4.44        # V, absolute potential of the standard hydrogen electrode

ANALYTES = {
    "dopamine":       "NCCc1ccc(O)c(O)c1",
    "ascorbic_acid":  "OC[C@@H](O)[C@@H]1OC(=O)C(O)=C1O",
    "uric_acid":      "O=C1NC(=O)C2=C(N1)NC(=O)N2",
}

def mol_from_smiles(smiles, seed=1):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    p = AllChem.ETKDGv3(); p.randomSeed = seed
    AllChem.EmbedMolecule(m, p); AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
    conf = m.GetConformer()
    syms = [a.GetSymbol() for a in m.GetAtoms()]
    pos = np.array([list(conf.GetAtomPosition(i)) for i in range(m.GetNumAtoms())])
    return Atoms(symbols=syms, positions=pos)

def uma_relax(atoms, charge, spin, predictor):
    from fairchem.core import FAIRChemCalculator
    a = atoms.copy(); a.info.update({"charge": charge, "spin": spin})
    a.calc = FAIRChemCalculator(predictor, task_name="omol")
    LBFGS(a, logfile=None).run(fmax=0.05, steps=300)
    return a

def pyscf_energy(atoms, charge, spin2s, solvent=False, want_orbitals=False):
    """spin2s = number of unpaired electrons (0 singlet, 1 doublet)."""
    from pyscf import gto, dft
    mol = gto.M(atom=[(s, tuple(p)) for s, p in
                zip(atoms.get_chemical_symbols(), atoms.get_positions())],
                basis="6-31g*", charge=charge, spin=spin2s, verbose=0)
    mf = dft.RKS(mol) if spin2s == 0 else dft.UKS(mol)
    mf.xc = "b3lyp"
    if solvent:
        from pyscf.solvent import ddcosmo
        mf = ddcosmo.ddcosmo_for_scf(mf)
        mf.with_solvent.eps = 78.3553          # water
    mf.conv_tol = 1e-7
    e = mf.kernel()
    out = {"E": float(e)}
    if want_orbitals:
        mo_e, occ = mf.mo_energy, mf.mo_occ
        homo = int(np.max(np.where(occ > 0))); lumo = homo + 1
        out["HOMO"] = float(mo_e[homo] * HARTREE_EV)
        out["LUMO"] = float(mo_e[lumo] * HARTREE_EV)
        out["dipole"] = float(np.linalg.norm(mf.dip_moment(unit="Debye", verbose=0)))
    return out

def descriptors(name, smiles, predictor):
    print(f"[{name}]", flush=True)
    neu = mol_from_smiles(smiles)
    neu = uma_relax(neu, 0, 1, predictor)          # neutral singlet
    cat = uma_relax(neu, 1, 2, predictor)          # radical cation doublet
    write(os.path.join(SDIR, f"analyte_{name}.xyz"), neu)

    g_neu = pyscf_energy(neu, 0, 0, solvent=False, want_orbitals=True)
    g_cat = pyscf_energy(cat, 1, 1, solvent=False)
    s_neu = pyscf_energy(neu, 0, 0, solvent=True)
    s_cat = pyscf_energy(cat, 1, 1, solvent=True)

    HOMO, LUMO = g_neu["HOMO"], g_neu["LUMO"]
    gap = LUMO - HOMO
    I = -HOMO; A = -LUMO
    chi = (I + A) / 2; eta = (I - A) / 2
    omega = chi**2 / (2*eta) if eta else float("nan")
    IP_gas = (g_cat["E"] - g_neu["E"]) * HARTREE_EV
    dGox = (s_cat["E"] - s_neu["E"]) * HARTREE_EV          # eV (1 e-)
    Eox_SHE = dGox - SHE_ABS                                # V vs SHE

    d = {"HOMO_eV": HOMO, "LUMO_eV": LUMO, "gap_eV": gap,
         "chi_eV": chi, "eta_eV": eta, "omega_eV": omega,
         "dipole_D": g_neu["dipole"], "IP_gas_eV": IP_gas,
         "E_ox_vs_SHE_V": Eox_SHE}
    print(f"  HOMO={HOMO:.2f} gap={gap:.2f} eta={eta:.2f} IP={IP_gas:.2f} "
          f"E_ox(SHE)={Eox_SHE:+.2f} V", flush=True)
    return d

def main():
    from fairchem.core import pretrained_mlip
    predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cpu")
    t0 = time.time()
    table = {name: descriptors(name, smi, predictor)
             for name, smi in ANALYTES.items()}
    # rank by ease of oxidation (higher HOMO / lower E_ox = easier to oxidize)
    order = sorted(table, key=lambda k: table[k]["E_ox_vs_SHE_V"])
    out = {"panel": "DA / AA / UA electrochemical sensing benchmark",
           "level": "UMA(omol) geometries + PySCF B3LYP/6-31G*, ddCOSMO(water)",
           "SHE_abs_V": SHE_ABS, "descriptors": table,
           "oxidation_order_easiest_first": order,
           "note": ("real DA/AA/UA oxidation is proton-coupled (PCET); absolute "
                    "peak potentials are pH-dependent. Relative ordering / peak "
                    "separation is the sensor-selectivity signal.")}
    json.dump(out, open(os.path.join(RDIR, "sensor_descriptors.json"), "w"), indent=2)
    print(f"\nEasiest->hardest to oxidize: {order}  ({time.time()-t0:.0f}s)")
    print("-> results/sensor_descriptors.json")

if __name__ == "__main__":
    main()
