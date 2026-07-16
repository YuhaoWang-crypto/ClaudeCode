"""
真实 DFT (GPAW/PBE) 态密度 —— 纯 SrTiO3 vs Nb 掺杂, 看 n 型能级。
返回 DOS 能量网格、总 DOS、Nb-4d/Ti-3d/O-2p 投影 DOS、费米能、HOMO-LUMO 间隙。
"""
import modal

image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("gcc", "g++", "libxc-dev", "libblas-dev", "liblapack-dev",
                 "libopenmpi-dev", "openmpi-bin", "libfftw3-dev")
    .pip_install("gpaw", "ase", "numpy", "scipy")
    .run_commands("gpaw install-data --register /opt/gpaw-data")
    .env({"GPAW_SETUP_PATH": "/opt/gpaw-data/gpaw-setups-24.11.0",
          "OMP_NUM_THREADS": "1"})
)
app = modal.App("gpaw-sto-dos", image=image)


@app.function(cpu=16.0, memory=32768, timeout=5400)
def dos(cif: str, name: str, kpts: int = 4, ecut: float = 400.0) -> dict:
    import io, glob, os
    import numpy as np
    from ase.io import read as ase_read
    from gpaw import GPAW, PW, FermiDirac

    # 定位 PAW 数据集目录
    cand = glob.glob("/opt/gpaw-data/gpaw-setups-*")
    if cand:
        os.environ["GPAW_SETUP_PATH"] = cand[0]

    atoms = ase_read(io.StringIO(cif), format="cif")
    calc = GPAW(mode=PW(ecut), xc="PBE", kpts=(kpts, kpts, kpts),
                nbands="200%",  # 足够空带以覆盖导带(否则金属体系 E_F 浮在已算带之上)
                convergence={"bands": -10},
                occupations=FermiDirac(0.05), txt="/tmp/gpaw.txt")
    atoms.calc = calc
    atoms.get_potential_energy()
    ef = float(calc.get_fermi_level())

    # 新版 GPAW DOS 接口
    egrid = np.linspace(ef - 9.0, ef + 7.0, 1200)
    dosc = calc.dos()
    total = np.asarray(dosc.raw_dos(egrid, spin=0, width=0.12))
    syms = atoms.get_chemical_symbols()
    L = {"s": 0, "p": 1, "d": 2}
    def pdos(elem, ang):
        idxs = [i for i, s in enumerate(syms) if s == elem]
        if not idxs:
            return None
        acc = np.zeros_like(egrid)
        for a in idxs:
            acc += np.asarray(dosc.raw_pdos(egrid, a, L[ang], None, spin=0, width=0.12))
        return [round(float(x), 4) for x in acc]

    # 间隙估计: 找 Ef 附近 DOS≈0 的区间宽度
    near = (egrid > ef - 4) & (egrid < ef + 4)
    tol = 0.02 * total.max()
    gap = float(np.sum((total < tol) & near) * (egrid[1] - egrid[0]))

    res = {"name": name, "formula": atoms.get_chemical_formula(),
           "E_fermi": ef, "dos_gap_est": round(gap, 2),
           "energy_rel_Ef": [round(float(x - ef), 4) for x in egrid],
           "dos_total": [round(float(x), 4) for x in total],
           "pdos_Ti_d": pdos("Ti", "d"), "pdos_Nb_d": pdos("Nb", "d"),
           "pdos_O_p": pdos("O", "p")}
    print(f"[{name}] Ef={ef:.3f} eV  gap_est~{gap:.2f} eV")
    return res


@app.local_entrypoint()
def main(candidates: str, out: str):
    import json
    cifs = json.load(open(candidates))
    results = []
    for name, cif in cifs.items():
        kp = 6 if "pure" in name else 4
        print(f"跑 DFT: {name} (kpts={kp})...")
        results.append(dos.remote(cif, name, kpts=kp))
    json.dump(results, open(out, "w"))
    for r in results:
        print(f"{r['name']}: Ef={r['E_fermi']:.3f}  gap_est={r.get('dos_gap_est')} eV")
    print("已保存:", out)
