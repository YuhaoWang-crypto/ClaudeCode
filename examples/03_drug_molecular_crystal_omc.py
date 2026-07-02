"""
示例 3：药物结晶 / 有机分子晶体 —— 用 UMA 的 `omc` 任务
（Open Molecular Crystals，专为有机分子晶体训练）算晶格能、比较多晶型。

为什么不用 OMat24：OMat24 训练在无机晶体上，不覆盖靠范德华力+氢键
维系的有机分子晶体。药物晶体应走 `omc`（晶体）/`omol`（孤立分子）任务。

典型用途：
  - 给定候选晶体结构，算相对晶格能，为多晶型（polymorph）排序
  - 结合 ASE 变胞弛豫

重要限制：
  - 完整的晶体结构预测（CSP）是"生成候选结构 + 能量排序"两步；这里只做
    第二步的能量排序。生成候选需要 CSP 工具（如 Genarris/CSP 流程）。
  - 弱相互作用对色散很敏感，注意模型是否已含色散校正，必要时叠加 D3。

模板，API 以安装版本为准。
"""
from ase.io import read
from ase.filters import FrechetCellFilter
from ase.optimize import BFGS
from fairchem.core import pretrained_mlip, FAIRChemCalculator


def lattice_energy_per_molecule(crystal, n_molecules, e_gas_per_molecule, calc):
    c = crystal.copy()
    c.calc = calc
    BFGS(FrechetCellFilter(c)).run(fmax=0.03, steps=200)  # 变胞弛豫
    e_crystal = c.get_potential_energy()
    return e_crystal / n_molecules - e_gas_per_molecule


def main():
    predictor = pretrained_mlip.get_predict_unit("uma-s-1", device="cpu")
    crystal_calc = FAIRChemCalculator(predictor, task_name="omc")   # 分子晶体
    mol_calc = FAIRChemCalculator(predictor, task_name="omol")      # 孤立分子

    # 孤立分子参考能（同一分子在真空）
    mol = read("drug_molecule.xyz")
    mol.calc = mol_calc
    e_gas = mol.get_potential_energy()

    # 遍历候选多晶型 CIF，按晶格能排序
    for cif, nmol in [("polymorph_I.cif", 4), ("polymorph_II.cif", 4)]:
        crystal = read(cif)
        e_latt = lattice_energy_per_molecule(crystal, nmol, e_gas, crystal_calc)
        print(f"{cif}: 每分子晶格能 = {e_latt:.4f} eV")


if __name__ == "__main__":
    main()
