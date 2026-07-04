import modal, json

app = modal.App("is621-qm")
img = (modal.Image.micromamba(python_version="3.11")
       .micromamba_install("psi4", "xtb", "openbabel", "numpy",
                           channels=["conda-forge"]))


@app.function(image=img, cpu=8.0, memory=16384, timeout=3600)
def run_qm(cluster_pdb: str, charge: int = -1):
    import subprocess, os, re, numpy as np
    os.makedirs("/work", exist_ok=True); os.chdir("/work")
    open("qm_raw.pdb", "w").write(cluster_pdb)
    log = {}

    # 1) add hydrogens at pH 7 with openbabel, export xyz
    subprocess.run(["obabel", "qm_raw.pdb", "-O", "qm_h.xyz", "-p", "7"], check=True,
                   capture_output=True)
    xyz = open("qm_h.xyz").read()
    natom = int(xyz.splitlines()[0])
    log["n_atoms_after_H"] = natom

    # 2) GFN2-xTB geometry optimization (semi-empirical, fast)
    r = subprocess.run(["xtb", "qm_h.xyz", "--chrg", str(charge), "--gfn", "2", "--opt", "tight"],
                       capture_output=True, text=True)
    xtb_out = r.stdout
    e_xtb = None
    for line in xtb_out.splitlines():
        if "TOTAL ENERGY" in line:
            m = re.search(r"(-?\d+\.\d+)", line); e_xtb = float(m.group(1)) if m else None
    log["xtb_total_energy_Eh"] = e_xtb
    log["xtb_converged"] = os.path.exists("xtbopt.xyz")
    opt_xyz = open("xtbopt.xyz").read() if os.path.exists("xtbopt.xyz") else xyz

    # 3) Psi4 DFT single point (B3LYP-D3BJ/def2-SVP) on the xTB-optimized geometry
    import psi4
    psi4.set_memory("12 GB"); psi4.set_num_threads(8)
    lines = opt_xyz.splitlines()
    atom_lines = lines[2:2 + int(lines[0])]
    geom = "\n".join(atom_lines)
    Z = {"H":1,"C":6,"N":7,"O":8,"P":15,"MG":12,"Mg":12,"NA":11,"CL":17}
    tot_e = sum(Z.get(a.split()[0], 0) for a in atom_lines) - charge
    mult = 1 if tot_e % 2 == 0 else 2
    log["n_electrons"] = tot_e; log["multiplicity"] = mult
    mol = psi4.geometry(f"{geom}\nunits angstrom\nsymmetry c1\nno_reorient\nno_com\n")
    mol.set_molecular_charge(charge); mol.set_multiplicity(mult); mol.update_geometry()
    psi4.set_options({"basis": "def2-SVP", "scf_type": "df",
                      "reference": "rks" if mult == 1 else "uks"})
    e_dft = psi4.energy("b3lyp", molecule=mol)
    log["psi4_b3lyp_def2svp_Eh"] = float(e_dft)
    log["psi4_method"] = "B3LYP/def2-SVP single point on GFN2-xTB geometry (POC; D3 omitted)"
    log["opt_xyz"] = opt_xyz
    return log


@app.local_entrypoint()
def main():
    pdb = open("mdqm/qm_cluster_raw.pdb").read()
    out = run_qm.remote(pdb, charge=-2)
    out_small = {k: v for k, v in out.items() if k != "opt_xyz"}
    print(json.dumps(out_small, indent=2))
    json.dump(out, open("mdqm/qm_result.json", "w"))
    open("mdqm/qm_optimized.xyz", "w").write(out.get("opt_xyz", ""))
    print("WROTE mdqm/qm_result.json")
