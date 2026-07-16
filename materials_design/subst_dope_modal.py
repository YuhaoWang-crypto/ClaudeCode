"""
替换/掺杂 -> 性质变化, FairChem UMA 模拟 (Modal GPU)。
对每个结构: 变胞弛豫 -> 能量/晶格; 可选 Γ 最大声子; 可选 Birch-Murnaghan 体模量(EOS)。
"""
import os
import modal

image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install("fairchem-core>=2.0", "ase", "pymatgen", "huggingface-hub")
)
app = modal.App("uma-subst-dope", image=image)


@app.function(gpu="A10G", timeout=3600,
              secrets=[modal.Secret.from_dict({"HF_TOKEN": os.environ.get("HF_TOKEN", "")})])
def analyze(cifs: dict, model: str = "uma-s-1p1", task: str = "omat",
            fmax: float = 0.03, do_phonon: bool = True, do_eos: bool = True) -> list:
    import io, os, tempfile
    import numpy as np
    from ase.io import read as ase_read
    from ase.optimize import FIRE
    from ase.filters import FrechetCellFilter
    from ase.vibrations import Vibrations
    from ase.units import invcm, GPa
    from ase.eos import EquationOfState
    from fairchem.core import pretrained_mlip, FAIRChemCalculator

    predictor = pretrained_mlip.get_predict_unit(model, device="cuda")
    calc = FAIRChemCalculator(predictor, task_name=task)

    out = []
    for name, cif in cifs.items():
        try:
            atoms = ase_read(io.StringIO(cif), format="cif")
            atoms.calc = calc
            FIRE(FrechetCellFilter(atoms), logfile=None).run(fmax=fmax, steps=400)
            n = len(atoms)
            rec = {"name": name, "formula": atoms.get_chemical_formula(), "n_atoms": n,
                   "E_per_atom": round(float(atoms.get_potential_energy()) / n, 4),
                   "lattice_abc": [round(float(x), 4) for x in atoms.cell.cellpar()[:3]],
                   "volume_per_atom": round(float(atoms.get_volume()) / n, 4)}
            # 最近邻金属-配体键长 (取第一个替换/掺杂原子周围)
            try:
                from ase.geometry import get_distances
                D = atoms.get_all_distances(mic=True)
                np.fill_diagonal(D, 1e9)
                rec["min_bond"] = round(float(D.min()), 4)
            except Exception:
                pass
            # Γ 最大声子
            if do_phonon and n <= 14:
                try:
                    with tempfile.TemporaryDirectory() as td:
                        vib = Vibrations(atoms, name=os.path.join(td, "v"))
                        vib.run()
                        en = vib.get_energies()
                        real = np.array([e.real for e in en if abs(e.imag) < 1e-6])
                        rec["phonon_max_cm-1"] = round(float(real.max() / invcm), 1) if real.size else None
                except Exception as pe:
                    rec["phonon_error"] = repr(pe)[:80]
            # 体模量 (EOS): 等比例缩放体积, 内部坐标不动, 拟合 Birch-Murnaghan
            if do_eos and n <= 14:
                try:
                    cell0 = atoms.cell.copy(); vols = []; ens = []
                    for s in np.linspace(0.94, 1.06, 7):
                        a2 = atoms.copy(); a2.calc = calc
                        a2.set_cell(cell0 * s ** (1/3), scale_atoms=True)
                        vols.append(a2.get_volume()); ens.append(a2.get_potential_energy())
                    eos = EquationOfState(vols, ens, eos="birchmurnaghan")
                    v0, e0, B = eos.fit()
                    rec["bulk_modulus_GPa"] = round(float(B / GPa), 1)
                except Exception as be:
                    rec["eos_error"] = repr(be)[:80]
            out.append(rec)
            print(f"[OK] {name}: {rec['formula']} a={rec['lattice_abc']} "
                  f"phonon={rec.get('phonon_max_cm-1')} B={rec.get('bulk_modulus_GPa')}")
        except Exception as e:
            out.append({"name": name, "error": repr(e)[:120]})
            print(f"[ERR] {name}: {e!r}")
    return out


@app.local_entrypoint()
def main(candidates: str, out: str):
    import json
    cifs = json.load(open(candidates))
    print(f"提交 {len(cifs)} 个结构...")
    res = analyze.remote(cifs)
    for r in res:
        print(json.dumps(r, ensure_ascii=False))
    json.dump(res, open(out, "w"), ensure_ascii=False, indent=2)
    print("已保存:", out)
