import modal, json

app = modal.App("is621-qmmm-barrier")
img = (modal.Image.micromamba(python_version="3.11")
       .micromamba_install("psi4", "xtb", "openbabel", "numpy", channels=["conda-forge"]))


@app.function(image=img, cpu=8.0, memory=16384, timeout=5400)
def barrier(cluster_pdb: str, label: str, charge: int = -2,
            rc_start: float = 3.39, rc_end: float = 1.65, nsteps: int = 14):
    import subprocess, os, re, numpy as np
    wd = f"/work/{label}"; os.makedirs(wd, exist_ok=True); os.chdir(wd)
    open("raw.pdb", "w").write(cluster_pdb)
    log = {"label": label}

    # 1) add H at pH 7, get xyz
    subprocess.run(["obabel", "raw.pdb", "-O", "h.xyz", "-p", "7"], check=True, capture_output=True)
    lines = open("h.xyz").read().splitlines()
    n = int(lines[0]); els = []; xyz = []
    for l in lines[2:2 + n]:
        p = l.split(); els.append(p[0]); xyz.append([float(x) for x in p[1:4]])
    xyz = np.array(xyz); els = np.array(els)
    Z = {"H": 1, "C": 6, "N": 7, "O": 8, "P": 15, "Mg": 12, "MG": 12, "Na": 11, "Cl": 17}
    tot_e = sum(Z[e] for e in els) - charge
    if tot_e % 2 == 1:      # keep closed shell; nudge charge by 1 (report it)
        charge += 1; tot_e -= 1
    log["charge"] = charge; log["n_electrons"] = tot_e; log["n_atoms"] = n

    # 2) identify reaction-coordinate atoms and anchors (1-based for xtb)
    iP = int(np.where(els == "P")[0][0])
    P = xyz[iP]
    dP = np.linalg.norm(xyz - P, axis=1)
    # O3' = oxygen 2.5-3.6 A from P that is bonded to a carbon (leaving group), pick closest such
    o_cands = [i for i in range(n) if els[i] == "O" and 2.5 < dP[i] < 3.8
               and any(els[j] == "C" and np.linalg.norm(xyz[i]-xyz[j]) < 1.7 for j in range(n))]
    iO3 = min(o_cands, key=lambda i: dP[i])
    anchors = [i for i in range(n) if els[i] != "H" and dP[i] > 5.0]   # rigid scaffold caps
    log["iP"] = iP+1; log["iO3"] = iO3+1; log["d_P_O3_start"] = float(dP[iO3])
    log["n_anchors"] = len(anchors)

    fix_str = ",".join(str(a+1) for a in anchors)
    inp = f"""$fix
  atoms: {fix_str}
$end
$constrain
  force constant=1.5
  distance: {iP+1},{iO3+1},{dP[iO3]:.3f}
$end
$scan
  1: {rc_start:.3f},{rc_end:.3f},{nsteps}
$end
"""
    open("scan.inp", "w").write(inp)

    # 3) run xtb relaxed scan (GFN2)
    r = subprocess.run(["xtb", "h.xyz", "--input", "scan.inp", "--gfn", "2",
                        "--chrg", str(charge), "--opt", "tight"],
                       capture_output=True, text=True)
    log["xtb_ok"] = os.path.exists("xtbscan.log")
    if not log["xtb_ok"]:
        log["xtb_tail"] = r.stdout[-1500:] + r.stderr[-1500:]; return log

    # 4) parse scan profile (each frame comment holds an energy in Eh)
    prof = []
    sc = open("xtbscan.log").read().splitlines()
    i = 0
    while i < len(sc):
        na = int(sc[i]); comment = sc[i+1]
        m = re.search(r"(-?\d+\.\d+)", comment)
        e = float(m.group(1)) if m else None
        prof.append(e); i += na + 2
    prof = np.array(prof, float)
    coords = np.linspace(rc_start, rc_end, len(prof))
    kcal = 627.509
    e_rel = (prof - prof.min()) * kcal
    # reactant = bonded end (rc ~1.65, last frames), product = cleaved (rc ~3.39, first frames)
    react_i = int(np.argmin(np.abs(coords - rc_end)))
    ts_i = int(np.argmax(prof))
    barrier_fwd = float((prof[ts_i] - prof[react_i]) * kcal)   # cleavage from re-ligated reactant
    log["scan_coords"] = coords.tolist()
    log["scan_e_rel_kcal"] = e_rel.tolist()
    log["ts_coord"] = float(coords[ts_i]); log["barrier_kcal"] = barrier_fwd

    # 5) Psi4 DFT single points at reactant and TS geometries
    frames = []
    i = 0
    while i < len(sc):
        na = int(sc[i]); block = sc[i+2:i+2+na]; frames.append(block); i += na+2
    def dft(block):
        import psi4
        psi4.core.clean(); psi4.set_memory("12 GB"); psi4.set_num_threads(8)
        geom = "\n".join(block)
        mol = psi4.geometry(f"{geom}\nunits angstrom\nsymmetry c1\nno_reorient\nno_com\n")
        mol.set_molecular_charge(charge); mol.set_multiplicity(1); mol.update_geometry()
        psi4.set_options({"basis": "def2-SVP", "scf_type": "df", "reference": "rks"})
        return float(psi4.energy("b3lyp"))
    try:
        e_r = dft(frames[react_i]); e_ts = dft(frames[ts_i])
        log["dft_barrier_kcal"] = (e_ts - e_r) * 627.509
        log["dft_note"] = "B3LYP/def2-SVP single points at xTB reactant & TS"
    except Exception as ex:
        log["dft_error"] = str(ex)[:300]
    return log


@app.local_entrypoint()
def main():
    wt = open("mdqm/qm_a13_wt.pdb").read()
    ser = open("mdqm/qm_a13_ser.pdb").read()
    out = {}
    out["WT"] = barrier.remote(wt, "WT")
    out["A13S"] = barrier.remote(ser, "A13S")
    for k in out:
        print(k, {kk: out[k].get(kk) for kk in ("charge", "barrier_kcal", "dft_barrier_kcal", "ts_coord", "xtb_ok")})
    if out["WT"].get("barrier_kcal") and out["A13S"].get("barrier_kcal"):
        print("xTB  ddG-dagger (A13S - WT) =", round(out["A13S"]["barrier_kcal"]-out["WT"]["barrier_kcal"], 2), "kcal/mol")
        if out["WT"].get("dft_barrier_kcal") and out["A13S"].get("dft_barrier_kcal"):
            print("DFT  ddG-dagger (A13S - WT) =", round(out["A13S"]["dft_barrier_kcal"]-out["WT"]["dft_barrier_kcal"], 2), "kcal/mol")
    json.dump(out, open("mdqm/qmmm_barrier.json", "w"))
    print("WROTE mdqm/qmmm_barrier.json")
