"""
sensor_screen.py -- parameterised electrochemical-sensor descriptor screening.

Given a panel of analytes (name -> SMILES), compute for each the molecular
descriptor fingerprint and the proton-coupled (PCET) oxidation potential, then
write a results table and figures. Swap the analytes to screen any drug /
metabolite / interferent set for sensor-material selection.

Usage:
    python sensor_screen.py --analytes analytes.json [--ph 7] [--outdir out]

analytes.json:  {"dopamine": "NCCc1ccc(O)c(O)c1", "uric_acid": "O=C1NC(=O)..."}

Descriptors (UMA geometries + PySCF B3LYP/6-31G*, ddCOSMO water):
  HOMO, LUMO, gap, chi, eta, omega, dipole, adiabatic IP,
  E_ox (1e- radical cation)  and  E_ox_PCET (1H+/1e- neutral radical, vs SHE at pH).

Requires: torch, fairchem-core, ase, pyscf, rdkit, matplotlib  (see setup.sh).
"""
import os, sys, json, argparse, time
import numpy as np

HARTREE_EV = 27.211386245988
G_HPLUS_AQ = -11.72     # eV, absolute aqueous proton free energy
E_ABS_SHE = 4.44        # V, absolute SHE potential

def mol_from_smiles(smiles, seed=1):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from ase import Atoms
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    p = AllChem.ETKDGv3(); p.randomSeed = seed
    if AllChem.EmbedMolecule(m, p) != 0:
        AllChem.EmbedMolecule(m, AllChem.ETKDGv2())
    AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
    conf = m.GetConformer()
    syms = [a.GetSymbol() for a in m.GetAtoms()]
    pos = np.array([list(conf.GetAtomPosition(i)) for i in range(m.GetNumAtoms())])
    return Atoms(symbols=syms, positions=pos)

def uma_relax(atoms, charge, spin, predictor):
    from fairchem.core import FAIRChemCalculator
    from ase.optimize import LBFGS
    a = atoms.copy(); a.info.update({"charge": charge, "spin": spin})
    a.calc = FAIRChemCalculator(predictor, task_name="omol")
    LBFGS(a, logfile=None).run(fmax=0.05, steps=300)
    return a

def pyscf_scf(atoms, charge, spin2s, solvent=False, orbitals=False):
    from pyscf import gto, dft
    mol = gto.M(atom=[(s, tuple(p)) for s, p in
                zip(atoms.get_chemical_symbols(), atoms.get_positions())],
                basis="6-31g*", charge=charge, spin=spin2s, verbose=0)
    mf = dft.RKS(mol) if spin2s == 0 else dft.UKS(mol)
    mf.xc = "b3lyp"
    if solvent:
        from pyscf.solvent import ddcosmo
        mf = ddcosmo.ddcosmo_for_scf(mf); mf.with_solvent.eps = 78.3553
    mf.conv_tol = 1e-7
    e = mf.kernel()
    out = {"E_eV": float(e) * HARTREE_EV}
    if orbitals:
        occ = mf.mo_occ; homo = int(np.max(np.where(occ > 0)))
        out["HOMO"] = float(mf.mo_energy[homo] * HARTREE_EV)
        out["LUMO"] = float(mf.mo_energy[homo + 1] * HARTREE_EV)
        out["dipole"] = float(np.linalg.norm(mf.dip_moment(unit="Debye", verbose=0)))
    return out

def on_hydrogens(atoms):
    pos = atoms.get_positions(); sym = atoms.get_chemical_symbols(); hs = []
    for i, s in enumerate(sym):
        if s != "H":
            continue
        d = np.linalg.norm(pos - pos[i], axis=1)
        if any(sym[j] in ("O", "N") and j != i and d[j] < 1.30 for j in range(len(sym))):
            hs.append(i)
    return hs

def analyte_descriptors(name, smiles, predictor, ph):
    neu0 = mol_from_smiles(smiles)
    neu = uma_relax(neu0, 0, 1, predictor)
    cat = uma_relax(neu, 1, 2, predictor)
    g_neu = pyscf_scf(neu, 0, 0, orbitals=True)
    g_cat = pyscf_scf(cat, 1, 1)
    s_neu = pyscf_scf(neu, 0, 0, solvent=True)
    s_cat = pyscf_scf(cat, 1, 1, solvent=True)

    HOMO, LUMO = g_neu["HOMO"], g_neu["LUMO"]
    I, A = -HOMO, -LUMO
    chi, eta = (I + A) / 2, (I - A) / 2
    omega = chi**2 / (2*eta) if eta else float("nan")
    IP = g_cat["E_eV"] - g_neu["E_eV"]
    Eox_cation = (s_cat["E_eV"] - s_neu["E_eV"]) - E_ABS_SHE

    # PCET: most stable neutral radical from removing an O-H/N-H
    best = None
    for h in on_hydrogens(neu):
        rad = neu.copy(); del rad[h]
        rad = uma_relax(rad, 0, 2, predictor)
        E = pyscf_scf(rad, 0, 1, solvent=True)["E_eV"]
        if best is None or E < best:
            best = E
    Eox_pcet_ph0 = (best - s_neu["E_eV"]) + G_HPLUS_AQ - E_ABS_SHE if best else float("nan")
    Eox_pcet = Eox_pcet_ph0 - 0.05916 * ph if best else float("nan")

    return {"HOMO_eV": HOMO, "LUMO_eV": LUMO, "gap_eV": LUMO - HOMO,
            "chi_eV": chi, "eta_eV": eta, "omega_eV": omega,
            "dipole_D": g_neu["dipole"], "IP_gas_eV": IP,
            "E_ox_cation_vs_SHE_V": Eox_cation,
            "E_ox_PCET_pH%.0f_vs_SHE_V" % ph: Eox_pcet}

def plot(table, ph, outdir):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    names = list(table); ekey = "E_ox_PCET_pH%.0f_vs_SHE_V" % ph
    fig, (a, b) = plt.subplots(1, 2, figsize=(12, 5))
    V = np.linspace(min(table[n][ekey] for n in names) - 0.2,
                    max(table[n][ekey] for n in names) + 0.2, 600)
    for n in names:
        e = table[n][ekey]
        a.plot(V, np.exp(-0.5*((V-e)/0.03)**2), lw=2, label=f"{n}: {e:+.2f} V")
        a.fill_between(V, np.exp(-0.5*((V-e)/0.03)**2), alpha=0.15)
    a.set_xlabel(f"PCET oxidation potential vs SHE (V), pH {ph:.0f}")
    a.set_ylabel("Response (a.u.)"); a.set_title("Predicted voltammogram"); a.legend(fontsize=8, frameon=False)
    keys = ["HOMO_eV","LUMO_eV","gap_eV","eta_eV","omega_eV","dipole_D","IP_gas_eV",ekey]
    M = np.array([[table[n][k] for k in keys] for n in names], float)
    Z = (M - M.mean(0)) / (M.std(0) + 1e-9)
    im = b.imshow(Z, aspect="auto", cmap="coolwarm", vmin=-1.5, vmax=1.5)
    b.set_xticks(range(len(keys))); b.set_xticklabels([k.split("_")[0] for k in keys], rotation=30, fontsize=8)
    b.set_yticks(range(len(names))); b.set_yticklabels(names)
    for i in range(len(names)):
        for j in range(len(keys)):
            b.text(j, i, f"{M[i,j]:.2f}", ha="center", va="center", fontsize=7)
    b.set_title("Descriptor fingerprint (z-score)"); fig.colorbar(im, ax=b, shrink=0.7)
    fig.tight_layout(); path = os.path.join(outdir, "sensor_screen.png")
    fig.savefig(path, dpi=150); print("wrote", path)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--analytes", required=True, help="JSON {name: SMILES}")
    ap.add_argument("--ph", type=float, default=7.0)
    ap.add_argument("--outdir", default="sensor_out")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    analytes = json.load(open(args.analytes))

    from fairchem.core import pretrained_mlip
    predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cpu")
    t0 = time.time(); table = {}
    for name, smi in analytes.items():
        print(f"[{name}] ...", flush=True)
        table[name] = analyte_descriptors(name, smi, predictor, args.ph)
        ek = "E_ox_PCET_pH%.0f_vs_SHE_V" % args.ph
        print(f"  HOMO={table[name]['HOMO_eV']:.2f}  E_ox(PCET,pH{args.ph:.0f})="
              f"{table[name][ek]:+.2f} V", flush=True)
    ekey = "E_ox_PCET_pH%.0f_vs_SHE_V" % args.ph
    order = sorted(table, key=lambda k: table[k][ekey])
    json.dump({"pH": args.ph, "descriptors": table,
               "oxidation_order_easiest_first": order},
              open(os.path.join(args.outdir, "sensor_screen.json"), "w"), indent=2)
    plot(table, args.ph, args.outdir)
    print(f"\nEasiest->hardest to oxidise: {order}  [{time.time()-t0:.0f}s]")
    print(f"-> {args.outdir}/sensor_screen.json")

if __name__ == "__main__":
    main()
