"""
FairChem UMA structure relaxer on Modal (GPU).

给定一组结构(CIF 字符串),用 Meta FairChem 的 UMA 通用机器学习势(MLIP)做几何弛豫,
返回弛豫后能量/原子、最大受力、体积变化等,用于对手工设计的候选材料做稳定性 sanity check。

需要环境变量 HF_TOKEN (UMA 权重在 HuggingFace 上是 gated 的)。

用法:
    modal run fairchem_relax_modal.py            # 跑内置示例
"""
import os

import modal

image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install(
        "fairchem-core>=2.0",
        "ase",
        "pymatgen",
        "huggingface-hub",
    )
)

app = modal.App("fairchem-uma-relax", image=image)


@app.function(gpu="A10G", timeout=3600, secrets=[modal.Secret.from_dict({"HF_TOKEN": os.environ.get("HF_TOKEN", "")})])
def relax(cifs: dict, model: str = "uma-s-1p1", task: str = "omat", fmax: float = 0.05,
          steps: int = 300, phonon_max_atoms: int = 12) -> list:
    """cifs: {name: cif_string}. 弛豫每个结构; 对小胞额外算 Γ 点最大声子频率。"""
    import io
    import os
    import tempfile

    import numpy as np
    from ase.io import read as ase_read
    from ase.optimize import FIRE
    from ase.filters import FrechetCellFilter
    from ase.vibrations import Vibrations
    from ase.units import invcm

    from fairchem.core import pretrained_mlip, FAIRChemCalculator

    predictor = pretrained_mlip.get_predict_unit(model, device="cuda")
    calc = FAIRChemCalculator(predictor, task_name=task)

    out = []
    for name, cif in cifs.items():
        try:
            atoms = ase_read(io.StringIO(cif), format="cif")
            atoms.calc = calc
            e0 = float(atoms.get_potential_energy())
            v0 = float(atoms.get_volume())
            opt = FIRE(FrechetCellFilter(atoms), logfile=None)
            opt.run(fmax=fmax, steps=steps)
            e1 = float(atoms.get_potential_energy())
            v1 = float(atoms.get_volume())
            fmax_final = float(np.abs(atoms.get_forces()).max())
            cell = atoms.cell.cellpar()
            rec = {
                "name": name,
                "formula": atoms.get_chemical_formula(),
                "n_atoms": len(atoms),
                "E_relaxed_per_atom": round(e1 / len(atoms), 4),
                "volume_change_pct": round(100.0 * (v1 - v0) / v0, 2),
                "lattice_abc": [round(float(x), 4) for x in cell[:3]],
                "fmax_final": round(fmax_final, 4),
                "converged": bool(fmax_final <= fmax),
            }
            # Γ 点声子: 最大声子频率 (多声子弛豫/电声耦合的关键量)
            if len(atoms) <= phonon_max_atoms:
                try:
                    with tempfile.TemporaryDirectory() as td:
                        vib = Vibrations(atoms, name=os.path.join(td, "vib"))
                        vib.run()
                        energies = vib.get_energies()  # eV, 可能含虚频(复数)
                        real = np.array([e.real for e in energies if abs(e.imag) < 1e-6])
                        max_cm = float(real.max() / invcm) if real.size else None
                        rec["phonon_max_cm-1"] = round(max_cm, 1) if max_cm else None
                except Exception as pe:
                    rec["phonon_error"] = repr(pe)
            out.append(rec)
            print(f"[OK] {name}: {rec['formula']} E={rec['E_relaxed_per_atom']} eV/at "
                  f"abc={rec['lattice_abc']} phonon_max={rec.get('phonon_max_cm-1')}")
        except Exception as e:
            out.append({"name": name, "error": repr(e)})
            print(f"[ERR] {name}: {e!r}")
    return out


@app.local_entrypoint()
def main(candidates: str = "/home/user/candidates3.json"):
    import json
    cifs = json.load(open(candidates))
    print(f"提交 {len(cifs)} 个候选结构到 FairChem UMA 弛豫...")
    res = relax.remote(cifs)
    print("\n==== FairChem UMA 弛豫结果 ====")
    for r in res:
        print(json.dumps(r, ensure_ascii=False))
    json.dump(res, open("/home/user/fairchem_results.json", "w"), ensure_ascii=False, indent=2)
