"""
硼化物 MB2 电子结构描述符 (GPAW/PBE): N(E_F) 与 B-p 在 E_F 的占比。
MgB2 的超导来自 B σ 带(px,py)被掏空穴 -> E_F 处 B-p 态占主导;
过渡金属硼化物 E_F 处是 d 态主导, B-p 占比低。此描述符应能挑出 MgB2。
"""
import modal

image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("gcc", "g++", "libxc-dev", "libblas-dev", "liblapack-dev",
                 "libopenmpi-dev", "openmpi-bin", "libfftw3-dev")
    .pip_install("gpaw", "ase", "numpy", "scipy")
    .run_commands("gpaw install-data --register /opt/gpaw-data")
    .env({"OMP_NUM_THREADS": "1"})
)
app = modal.App("gpaw-diboride", image=image)


@app.function(cpu=8.0, memory=16384, timeout=2400)
def descriptor(item: tuple, kpts: int = 10, ecut: float = 420.0) -> dict:
    cif, name = item
    import io, glob, os
    import numpy as np
    from ase.io import read as ase_read
    from gpaw import GPAW, PW, FermiDirac

    cand = glob.glob("/opt/gpaw-data/gpaw-setups-*")
    if cand:
        os.environ["GPAW_SETUP_PATH"] = cand[0]

    atoms = ase_read(io.StringIO(cif), format="cif")
    calc = GPAW(mode=PW(ecut), xc="PBE", kpts=(kpts, kpts, kpts),
                nbands="150%", occupations=FermiDirac(0.05),
                convergence={"density": 1e-6}, txt="/tmp/g.txt")
    atoms.calc = calc
    atoms.get_potential_energy()
    ef = float(calc.get_fermi_level())

    # 稳健的 N(E_F): 本征值高斯展开 (对金属可靠)
    sigma = 0.20
    wk = calc.get_k_point_weights()
    nsp = calc.get_number_of_spins()
    spin_deg = 2.0 if nsp == 1 else 1.0
    ntot = 0.0
    for s in range(nsp):
        for k in range(len(wk)):
            eps = np.asarray(calc.get_eigenvalues(kpt=k, spin=s))
            ntot += wk[k] * spin_deg * np.sum(
                np.exp(-0.5 * ((eps - ef) / sigma) ** 2)) / (sigma * np.sqrt(2 * np.pi))
    ntot = float(ntot)

    # B-p 占比: 用网格 raw_pdos / raw_dos 的比值 (比值对单点bug不敏感)
    frac_Bp = 0.0; n_Bp = 0.0
    try:
        dosc = calc.dos()
        eg = np.linspace(ef - 0.4, ef + 0.4, 9)
        tot_g = np.asarray(dosc.raw_dos(eg, spin=0, width=0.2)).mean()
        syms = atoms.get_chemical_symbols()
        bp = np.zeros_like(eg)
        for a in [i for i, s in enumerate(syms) if s == "B"]:
            bp += np.asarray(dosc.raw_pdos(eg, a, 1, None, spin=0, width=0.2))
        n_Bp = float(bp.mean())
        frac_Bp = n_Bp / tot_g if tot_g > 1e-6 else 0.0
    except Exception as e:
        frac_Bp = -1.0
    res = {"name": str(name), "formula": str(atoms.get_chemical_formula()),
           "N_Ef_total": float(round(ntot, 3)), "Bp_fraction": float(round(frac_Bp, 3))}
    print(f"[{name}] N(Ef)={ntot:.3f} states/eV  B-p占比={frac_Bp:.2f}")
    return res


@app.local_entrypoint()
def main(candidates: str, out: str):
    import json
    cifs = json.load(open(candidates))
    items = [(cif, name) for name, cif in cifs.items()]
    print(f"并行提交 {len(items)} 个硼化物 DFT...")
    results = list(descriptor.map(items))   # 并行执行
    json.dump(results, open(out, "w"))
    for r in sorted(results, key=lambda z: -z["Bp_fraction"]):
        print(f"{r['name'].split('|')[1]:18s} N(Ef)={r['N_Ef_total']:.2f}  B-p占比={r['Bp_fraction']:.2f}")
    print("已保存:", out)
