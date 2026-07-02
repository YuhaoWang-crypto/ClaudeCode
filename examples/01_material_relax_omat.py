"""
示例 1：无机材料 —— 用 UMA 的 `omat` 任务做能量/结构弛豫。
适用：晶体稳定性、形成能、弹性/力学、相对能量排序。

注意：这是模板，API 以你安装的 fairchem-core 版本为准（本文件基于 v2.x UMA 接口）。
运行前需能访问 HF 且已 login。
"""
from ase.build import bulk
from ase.optimize import BFGS
from fairchem.core import pretrained_mlip, FAIRChemCalculator


def main():
    # uma-s-1 = small 通用模型；也可用 uma-m-1（更准更慢）
    predictor = pretrained_mlip.get_predict_unit("uma-s-1", device="cpu")
    calc = FAIRChemCalculator(predictor, task_name="omat")  # 无机材料任务

    atoms = bulk("Cu", "fcc", a=3.6)
    atoms.calc = calc
    print("初始能量 (eV):", atoms.get_potential_energy())

    # 弛豫晶胞内原子（如需变胞用 ase.filters.FrechetCellFilter）
    BFGS(atoms).run(fmax=0.05, steps=100)
    print("弛豫后能量 (eV):", atoms.get_potential_energy())


if __name__ == "__main__":
    main()
