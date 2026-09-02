"""六个客户场景 —— 同一批 kernel, 不同 spec, 不同输出。

这是「模型产品 vs 数据库」的实测: 数据库对同一个基因永远返回同一行;
力场对同一个化合物在不同浓度/组合/上下文下返回不同的结论、不同的装配层数、
甚至不同的**可回答性** (有时它必须拒绝出数)。

A vs B  : 同一个药, 只改浓度和时长 -> 输出不同
A vs A2 : 同一个浓度, 换一个药 -> 一个有选择性窗口, 一个是广谱毒性
C vs C2 : 同一个组合结构, 一个命中实测上位性, 一个没有 -> 证据等级不同
D       : 覆盖不全的真实客户输入 -> 看引擎在哪里降级、在哪里拒答
E       : 换细胞上下文 -> 引擎必须说清主轴不覆盖该系
"""

from __future__ import annotations

from .spec import AssaySpec, Compound, ReadoutRequest

SCENARIOS: dict[str, AssaySpec] = {}


# --- A · 肝毒性早筛, 有选择性窗口 ----------------------------------------
# Pemetrexed → TYMS。TYMS 在 HepG2 (−0.818) 比 RPE1 (−0.171) 必需得多,
# 这是一个真正的肝癌选择性脆弱点, 引擎应该给出一个 >1 的窗口。
SCENARIOS["A_selective_window"] = AssaySpec(
    customer="CRO-A · 肝毒性早筛",
    note="单药, 靶点在肝癌系选择性必需, 问还有多少剂量余量",
    context="hepg2", normal_reference="rpe1", exposure_h=72.0,
    compounds=[Compound(name="Pemetrexed", conc_uM=0.5, ic50_uM=0.1, hill=1.0)],
)

# --- A2 · 同样浓度, 换成广谱必需靶点 -------------------------------------
# Bortezomib → PSMB5。HepG2 −1.847 / RPE1 −1.873, 两系几乎一样必需。
# 同一个引擎、同一套参数, 结论应该完全相反: 没有选择性区间。
SCENARIOS["A2_broad_toxicity"] = AssaySpec(
    customer="CRO-A · 同条件对照药",
    note="与 A 同浓度同时长, 换成两系共同必需的靶点",
    context="hepg2", normal_reference="rpe1", exposure_h=72.0,
    compounds=[Compound(name="Bortezomib", conc_uM=0.5, ic50_uM=0.1, hill=1.0)],
)

# --- B · 同一个药, 客户换了浓度和时长 ------------------------------------
# 与 A 完全同药同靶点。数据库会返回同一行; 力场返回不同的占据率、
# 不同的窗口位置、不同的通路幅度。
SCENARIOS["B_same_drug_higher_dose"] = AssaySpec(
    customer="CRO-A · 同药加压条件",
    note="与 A 同药, 浓度 ×20、暴露 ×2",
    context="hepg2", normal_reference="rpe1", exposure_h=144.0,
    compounds=[Compound(name="Pemetrexed", conc_uM=10.0, ic50_uM=0.1, hill=1.0)],
)

# --- C · 联用方案, 命中实测上位性 ----------------------------------------
# PLK4+STIL 在 Norman 2019 的 126 对里有实测, 且两个基因在 Chronos 上
# 都真的必需 (PLK4 −1.18/−1.15, STIL −0.70/−0.54)。交叉项会被启用。
SCENARIOS["C_combo_measured"] = AssaySpec(
    customer="Biotech-B · 联用方案排序",
    note="双靶组合, 靶点对在 ComboMap 中有实测上位性 → 启用交叉项",
    context="hepg2", normal_reference="rpe1", exposure_h=72.0,
    compounds=[
        Compound(name="cmpd-PLK4i", target="PLK4", conc_uM=1.0, ic50_uM=0.2),
        Compound(name="cmpd-STILi", target="STIL", conc_uM=1.0, ic50_uM=0.3),
    ],
)

# --- C2 · 同样结构的组合, 但没测过 ---------------------------------------
# TOP2A + TYMS 都必需, 但这一对不在 ComboMap 的 126 对里。
# 引擎应该明确回落到 Bliss 加和, 并用经验分布给出区间, 而不是硬编一个协同数。
SCENARIOS["C2_combo_unmeasured"] = AssaySpec(
    customer="Biotech-B · 无实测证据的组合",
    note="双靶组合, 靶点对不在 ComboMap 里 → 回落加和 + 经验区间",
    context="hepg2", normal_reference="rpe1", exposure_h=72.0,
    compounds=[
        Compound(name="Doxorubicin", conc_uM=0.5, ic50_uM=0.05),
        Compound(name="Pemetrexed", conc_uM=0.5, ic50_uM=0.1),
    ],
)

# --- D · 三药, 覆盖不全 ---------------------------------------------------
# 真实客户输入更常见的样子: 一个有全层覆盖, 一个缺 IC50 且是多靶点药,
# 一个完全查不到。看引擎逐条降级并在该拒答的地方拒答。
SCENARIOS["D_partial_coverage"] = AssaySpec(
    customer="Pharma-C · 三药组合可行性",
    note="覆盖不全: 缺 IC50、多靶点、完全未命中",
    context="hepg2", normal_reference="rpe1", exposure_h=72.0,
    compounds=[
        Compound(name="Paclitaxel", conc_uM=0.1, ic50_uM=0.005),
        Compound(name="Sorafenib", conc_uM=5.0),
        Compound(name="内部代号-X271", conc_uM=1.0, ic50_uM=1.0),
    ],
)

# --- E · 换细胞上下文 ----------------------------------------------------
# 同一个组合换到 T-ALL。引擎必须说清: Chronos 只覆盖 HepG2/RPE1,
# 换上下文只影响通路层与置信度加权, 不改变毒性主轴的数值。
SCENARIOS["E_context_switch"] = AssaySpec(
    customer="Biotech-B · 换适应症上下文",
    note="与 C 同组合, 上下文换成 Jurkat (T-ALL)",
    context="jurkat", normal_reference="rpe1", exposure_h=72.0,
    compounds=[
        Compound(name="cmpd-PLK4i", target="PLK4", conc_uM=1.0, ic50_uM=0.2),
        Compound(name="cmpd-STILi", target="STIL", conc_uM=1.0, ic50_uM=0.3),
    ],
    readout=ReadoutRequest(viability_curve=True, pathway_profile=True,
                           combination=True, top_genes=False),
)
