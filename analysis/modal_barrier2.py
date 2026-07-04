import modal, json

app = modal.App("is621-barrier2")
img = (modal.Image.micromamba(python_version="3.11")
       .micromamba_install("psi4", "xtb", "openbabel", "numpy", channels=["conda-forge"]))


@app.function(image=img, cpu=8.0, memory=16384, timeout=7200)
def run(pdb: str, label: str, charge: int = 0):
    import subprocess, os, re, numpy as np
    wd = f"/work/{label}"; os.makedirs(wd, exist_ok=True); os.chdir(wd)
    open("raw.pdb", "w").write(pdb)
    log = {"label": label}
    subprocess.run(["obabel", "raw.pdb", "-O", "h.xyz", "-p", "7"], check=True, capture_output=True)

    def read_xyz(fn):
        L = open(fn).read().splitlines(); n = int(L[0])
        els = [l.split()[0] for l in L[2:2+n]]
        xyz = np.array([[float(x) for x in l.split()[1:4]] for l in L[2:2+n]])
        return els, xyz
    def write_xyz(fn, els, xyz, comment=""):
        with open(fn, "w") as f:
            f.write(f"{len(els)}\n{comment}\n")
            for e, c in zip(els, xyz): f.write(f"{e} {c[0]:.6f} {c[1]:.6f} {c[2]:.6f}\n")
    els, xyz = read_xyz("h.xyz"); n = len(els)
    Z = {"H":1,"C":6,"N":7,"O":8,"P":15,"Mg":12,"MG":12}
    Zsum = sum(Z[e] for e in els)
    # formal charge of this anionic two-metal cluster is ~0 to -2; pick even-electron candidates
    cand_charges = [c for c in (0, -1, -2) if (Zsum - c) % 2 == 0]
    log["n_atoms"] = n; log["Zsum"] = Zsum; log["charge_candidates"] = cand_charges

    iP = int([i for i in range(n) if els[i]=="P"][0]); P = xyz[iP]
    dP = np.linalg.norm(xyz-P, axis=1)
    # O3' leaving group = O bonded to a C, ~1.6 A from P (re-ligated diester)
    o_es = [i for i in range(n) if els[i]=="O" and 1.3<dP[i]<1.9
            and any(els[j]=="C" and np.linalg.norm(xyz[i]-xyz[j])<1.7 for j in range(n))]
    iO3 = max(o_es, key=lambda i: min(np.linalg.norm(xyz[i]-xyz[j]) for j in range(n) if els[j]=="C"))
    # anchors: heavy cap atoms far from P (scaffold), gentle
    anchors = [i for i in range(n) if els[i]!="H" and dP[i]>5.5]
    log["iP"]=iP+1; log["iO3"]=iO3+1; log["n_anchors"]=len(anchors)
    fixs = ",".join(str(a+1) for a in anchors)

    # ---- 1. relax reactant (cap anchors only); try candidate charges w/ SCF aids ----
    open("fix.inp","w").write(f"$fix\n  atoms: {fixs}\n$end\n")
    # GFN-FF pre-relaxation clears steric clashes in the hand-built cluster (robust, no SCF)
    subprocess.run(["xtb","h.xyz","--gfnff","--input","fix.inp","--opt","normal"],
                   capture_output=True, text=True)
    start = "xtbopt.xyz" if os.path.exists("xtbopt.xyz") else "h.xyz"
    if start == "xtbopt.xyz":
        os.rename("xtbopt.xyz", "ff.xyz"); start = "ff.xyz"
    log["gfnff_preopt"] = (start == "ff.xyz")
    charge = None; relax_tail = ""
    for cc in cand_charges:
        for f in ("xtbopt.xyz","xtblast.xyz"):
            if os.path.exists(f): os.remove(f)
        r = subprocess.run(["xtb",start,"--input","fix.inp","--gfn","2","--chrg",str(cc),
                            "--opt","normal","--etemp","1000","--iterations","500"],
                           capture_output=True, text=True)
        if os.path.exists("xtbopt.xyz"):
            charge = cc; break
        relax_tail = r.stdout[-800:] + r.stderr[-400:]
    if charge is None:
        log["error"]="reactant relax failed for all charge candidates"; log["tail"]=relax_tail; return log
    log["charge"] = charge; log["n_electrons"] = Zsum - charge
    els_r, xyz_r = read_xyz("xtbopt.xyz")
    d = lambda a,b: float(np.linalg.norm(xyz_r[a]-xyz_r[b]))
    imgs = [i for i in range(n) if els[i] in ("Mg","MG")]
    react_geom = {"P_O3": d(iP,iO3), "MgMg": d(imgs[0],imgs[1]) if len(imgs)>1 else None,
                  "Mg1_P": d(imgs[0],iP)}
    log["reactant_geom"]=react_geom
    write_xyz("react.xyz", els_r, xyz_r)

    # ---- 2. constrained relaxed scan: break P-O3' (cleavage) 1.62 -> 3.4 ----
    d0 = d(iP,iO3)
    scan = (f"$fix\n  atoms: {fixs}\n$end\n$constrain\n  force constant=1.0\n"
            f"  distance: {iP+1},{iO3+1},{d0:.3f}\n$end\n$scan\n  1: {d0:.3f},3.40,18\n$end\n")
    open("scan.inp","w").write(scan)
    r2 = subprocess.run(["xtb","react.xyz","--input","scan.inp","--gfn","2","--chrg",str(charge),
                        "--opt","normal","--etemp","400","--iterations","500"],
                        capture_output=True, text=True)
    if not os.path.exists("xtbscan.log"):
        log["error"]="scan failed"; log["tail"]=r2.stdout[-1200:]; return log
    sc = open("xtbscan.log").read().splitlines()
    E, frames = [], []; i=0
    while i < len(sc):
        na=int(sc[i]); m=re.search(r"(-?\d+\.\d+)", sc[i+1]); E.append(float(m.group(1)))
        frames.append(sc[i+2:i+2+na]); i+=na+2
    E=np.array(E); kcal=627.509
    coords=np.linspace(d0,3.40,len(E))
    e_rel=(E-E.min())*kcal
    # quality gate
    max_step=float(np.max(np.abs(np.diff(e_rel))))
    ts_i=int(np.argmax(E)); react_i=0
    barrier=float((E[ts_i]-E[react_i])*kcal)
    smooth = max_step < 60 and 3 < barrier < 80 and ts_i not in (0,len(E)-1)
    log.update(scan_coords=coords.tolist(), scan_e_rel_kcal=e_rel.tolist(),
               max_step_kcal=max_step, ts_coord=float(coords[ts_i]),
               barrier_kcal=barrier, smooth_gate_pass=bool(smooth))

    # ---- 3. DFT single points at reactant & TS (only if smooth) ----
    if smooth:
        import psi4
        def dft(block):
            psi4.core.clean(); psi4.set_memory("12 GB"); psi4.set_num_threads(8)
            mol=psi4.geometry("\n".join(block)+"\nunits angstrom\nsymmetry c1\nno_reorient\nno_com\n")
            mol.set_molecular_charge(charge); mol.set_multiplicity(1); mol.update_geometry()
            psi4.set_options({"basis":"def2-SVP","scf_type":"df","reference":"rks"})
            return float(psi4.energy("b3lyp"))
        try:
            er=dft(frames[react_i]); ets=dft(frames[ts_i])
            log["dft_barrier_kcal"]=(ets-er)*627.509
        except Exception as ex:
            log["dft_error"]=str(ex)[:200]
    return log


@app.local_entrypoint()
def main():
    wt=open("mdqm/reactant_wt.pdb").read(); ser=open("mdqm/reactant_ser.pdb").read()
    out={"WT":run.remote(wt,"WT"), "A13S":run.remote(ser,"A13S")}
    for k,v in out.items():
        print(k,{x:v.get(x) for x in ("charge","reactant_geom","barrier_kcal","dft_barrier_kcal","max_step_kcal","smooth_gate_pass","ts_coord","error")})
    if all(out[k].get("smooth_gate_pass") for k in out):
        dd=out["A13S"]["barrier_kcal"]-out["WT"]["barrier_kcal"]
        print(f"\nxTB ddG-dagger (A13S-WT) = {dd:+.2f} kcal/mol")
        if out["WT"].get("dft_barrier_kcal") and out["A13S"].get("dft_barrier_kcal"):
            print(f"DFT ddG-dagger (A13S-WT) = {out['A13S']['dft_barrier_kcal']-out['WT']['dft_barrier_kcal']:+.2f} kcal/mol")
    else:
        print("\nGATE NOT PASSED — barrier not trustworthy; not reporting ddG.")
    json.dump(out, open("mdqm/barrier2.json","w"))
    print("WROTE mdqm/barrier2.json")
