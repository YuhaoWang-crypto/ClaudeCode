"""Write analysis_report.md from the computed result tables."""

from __future__ import annotations
import pandas as pd
import report_data as R

OUT = R.OUT


def md_table(hdr, rows):
    out = ["| " + " | ".join(hdr) + " |",
           "|" + "|".join(["---"] * len(hdr)) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(x) for x in r) + " |")
    return "\n".join(out)


def main():
    d = R.load()
    v = R.validation_stats()
    n = len(d)
    parts = []
    A = parts.append

    A("# ADMET 与类药性分析报告 — 30 个 FDA 已批准药物面板\n")
    A(f"生成日期: 2026-09-16  ·  化合物数: {n}  ·  预测终点: 41 (ADMET-AI) "
      f"+ 17 理化/规则描述符 + 5 个结构警示目录\n")

    A("## 1. 方法概要\n")
    A("- **标准化**: RDKit `MolStandardize` (Cleanup → 最大片段 → 去电荷) → 规范 SMILES\n"
      "- **理化描述符**: MW, cLogP (Crippen), TPSA, HBD/HBA, 可旋转键, 芳环数, "
      "sp³ 碳比例, 摩尔折射率, QED\n"
      "- **类药性规则**: Lipinski 五规则、Veber、Ghose、Egan\n"
      "- **结构警示**: RDKit `FilterCatalog` — PAINS A/B/C、Brenk、NIH\n"
      "- **ADMET 预测**: ADMET-AI (Swanson et al. 2024), Chemprop-RDKit 集成模型, "
      "在 TDC ADMET Benchmark Group 41 个数据集上训练\n"
      "- **参考背景**: ADMET-AI 内置 DrugBank 已批准药物参考集 (2,579 个药物) 的百分位数\n")

    A("## 2. 理化性质与类药性\n")
    A(md_table(*R.table_physchem(d)) + "\n")
    A(md_table(*R.table_rules(d)) + "\n")

    A("## 3. 关键 ADMET 终点\n")
    A(md_table(*R.table_endpoints(d)) + "\n")

    A("## 4. 连续型 (回归) 终点\n")
    A(md_table(*R.table_regression(d)) + "\n")

    A("## 5. 高风险化合物\n")
    A(md_table(*R.table_flagged(d)) + "\n")

    A("## 6. DrugBank 已批准药物参考百分位数\n")
    A(md_table(*R.table_percentiles(d)) + "\n")

    A("## 7. Demo 验证 — 已知阳性/阴性对照药物\n")
    A(md_table(*R.table_validation()) + "\n")
    A(f"DMPK / hERG / BBB 类终点共 {v['dmpk_tot']} 个阳性对照, "
      f"命中 {v['dmpk_rec']} 个;\n"
      f"Tox21 核受体类终点共 {v['nr_tot']} 个阳性对照, 仅命中 {v['nr_rec']} 个。\n")

    A("## 8. 主要局限\n")
    A("- ADMET-AI 的核受体 (NR-*) 与应激反应 (SR-*) 终点来自 Tox21 qHTS 激动模式"
      "定量高通量筛选, 数据高度不平衡; 本面板的验证显示它们**无法**重现已批准药物"
      "的既定药理活性 (他莫昔芬/ER、比卡鲁胺/AR、罗格列酮/PPARγ、阿那曲唑/芳香化酶"
      "均低于 0.5)。\n"
      "- hERG 模型在本面板上敏感度高但特异性低 (30 个药物中 21 个 > 0.5), "
      "在全为上市药物的集合上会高估风险。\n"
      "- 个别回归终点超出物理范围 (他莫昔芬血浆蛋白结合率预测 105.5%), "
      "提示适用域边界。\n"
      "- 所有预测为计算筛选辅助, 不能替代 GLP 体外/体内安全药理学试验。\n")

    txt = "\n".join(parts)
    with open(f"{OUT}/analysis_report.md", "w") as fh:
        fh.write(txt)
    print(f"✓ wrote {OUT}/analysis_report.md ({len(txt)} chars)")


if __name__ == "__main__":
    main()
