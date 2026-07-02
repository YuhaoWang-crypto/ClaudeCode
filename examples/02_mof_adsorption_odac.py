"""
示例 2：MOF 吸水 / 吸附 —— 用 UMA 的 `odac` 任务（Open Direct Air Capture，
专为 MOF + CO2/H2O 吸附训练）算吸附能。

吸附能 E_ads = E(MOF+吸附质) - E(MOF) - E(吸附质)

重要限制：
  - 单点吸附能 ≠ 吸水量/等温线。要得到 water uptake、工作容量，需要
    巨正则蒙特卡洛（GCMC，如 RASPA）或热力学积分，MLIP 只是提供能量的一环。
  - 需要一个真实 MOF 结构（CIF）。下面用占位说明流程。

模板，API 以安装版本为准。
"""
from ase.io import read
from ase.build import molecule
from fairchem.core import pretrained_mlip, FAIRChemCalculator


def energy(atoms, calc):
    atoms = atoms.copy()
    atoms.calc = calc
    return atoms.get_potential_energy()


def main():
    predictor = pretrained_mlip.get_predict_unit("uma-s-1", device="cpu")
    calc = FAIRChemCalculator(predictor, task_name="odac")  # MOF 吸附任务

    # 1) 读入 MOF 空胞（换成你的 CIF）
    mof = read("your_mof.cif")

    # 2) 吸附质：水分子
    h2o = molecule("H2O")

    # 3) 把水放进孔道的某个位置（真实工作流应枚举/采样多个吸附位点，取最低能）
    guest = mof.copy()
    h2o.translate(mof.get_center_of_mass() - h2o.get_center_of_mass())
    guest += h2o

    e_mof = energy(mof, calc)
    e_h2o = energy(h2o, calc)
    e_complex = energy(guest, calc)
    print("吸附能 E_ads (eV):", e_complex - e_mof - e_h2o)


if __name__ == "__main__":
    main()
