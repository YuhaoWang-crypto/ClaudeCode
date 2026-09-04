"""
案例 — ASFV p72 (B646L) 走完整条链路

把 assaysim 的每个模型逐段套到一个真实抗原上，并**先判定适用性再算**。
结论提前说：链路只有一部分适用，不适用的那几段本案例明确不出数。

真实输入（全部有出处，无一是编的）
----------------------------------
UniProt P22776 (CAPSH_ASFB7)   646 aa, 73,179 Da        主衣壳蛋白 p72
Wang N et al. 2019 Science     衣壳共 17,280 个蛋白      PMID 31624094
Liu S et al. 2019 Cell Res     T=277；p72 8,280 拷贝     PMID 31649031
                               以三聚体呈假六聚体排列
                               顶点 60 个五邻体蛋白
ChEMBL_37                      ASFV 靶点 pchembl 活性 0 条

运行:  python3 -m assaysim.case_asfv_p72
"""

from __future__ import annotations

import dataclasses
import math

import numpy as np

from . import agglutination as ag
from . import agid
from . import neutralization as neu
from . import viral_dynamics as vd

# --- 真实常量 ------------------------------------------------------------
P72_MW = 73179.0  # Da, UniProt P22776
P72_LEN = 646  # aa
P72_COPIES = 8280  # 每个病毒粒子，Liu 2019 Cell Res
P72_TRIMERS = P72_COPIES // 3  # 2760 个假六聚体壳粒
CAPSID_T = 277
CAPSID_TOTAL_PROTEINS = 17280  # Wang 2019 Science
IGG_MW = 150000.0

RULE = "=" * 78


def h(title: str) -> None:
    print(f"\n{RULE}\n{title}\n{RULE}")


def applicability() -> None:
    h("第 0 步 — 适用性判定（先做这个，再谈算什么）")
    rows = [
        ("M1 近邻热力学", "✅ 适用",
         "B646L 是 WOAH 指定的 ASFV PCR 靶基因，引物设计与包容性是纯序列问题"),
        ("M5 AGID 反应-扩散", "✅ 直接适用",
         "p72 是 ASFV 血清学的主抗原；这正是 AGID/免疫扩散试剂盒卖的东西"),
        ("M7 凝集聚集", "✅ 直接适用",
         "p72 包被乳胶的平板凝集；桥联占据模型不依赖病毒生物学"),
        ("M3 占据模型 Kd→NT50", "⚠️ 数学能跑，生物学不成立",
         "抗 p72 抗体不是中和抗体 —— 见下方结构性理由"),
        ("M4 单孔 ODE", "⚠️ 结构适用，参数缺",
         "ASFV 长在原代猪肺泡巨噬细胞，速率参数必须自行标定"),
        ("M6 非细胞→细胞桥接", "❌ 不适用",
         "ChEMBL_37 里 ASFV 靶点的 pchembl 活性为 0 条，无锚点可标定也无法验证"),
    ]
    for name, verdict, why in rows:
        print(f"  {name:22s} {verdict:16s} {why}")

    print(f"""
  为什么抗 p72 抗体不中和（这是个结构性理由，不是经验之谈）:
    ASFV 的成熟胞外病毒粒子在衣壳之外还有一层来自宿主细胞膜的囊膜。
    p72 位于囊膜**之内**的衣壳上，抗体在胞外根本够不到它。
    p72 因此是极好的**诊断**抗原（感染后抗体滴度高且持久），
    却是无效的**中和**靶点 —— 这也是 ASFV 至今没有有效疫苗的原因之一。

    模型不会替你拒绝这个输入。把 Kd 喂进占据模型，它照样吐一个 NT50。
    判断适用性是使用者的责任，本步骤就是干这个的。""")


def step_diffusion() -> dict:
    h("第 1 步 — 从真实分子量推扩散系数（M5 的输入）")
    forms = {
        "p72 单体 (E. coli His-tag)": P72_MW,
        "p72 三聚体 (天然壳粒)": P72_MW * 3,
        "猪 IgG (抗体侧)": IGG_MW,
    }
    out = {}
    for name, mw in forms.items():
        sf = 1.50 if "IgG" in name else 1.25
        rh = agid.hydrodynamic_radius_from_MW(mw, sf)
        D_free = agid.D_from_MW(mw, sf)
        D_gel = agid.D_from_MW(mw, sf, gel_factor=0.60)
        out[name] = {"mw": mw, "rh": rh, "D_free": D_free, "D_gel": D_gel}
        print(f"  {name:28s} MW {mw:9,.0f} Da  R_h {rh:5.2f} nm  "
              f"D(水) {D_free:.2e}  D(1%琼脂糖,f=0.6) {D_gel:.2e} cm²/s")
    print("""
  ⚠️ 凝胶阻滞因子 0.60 是**假定值**，不是本案例算出来的。琼脂糖对 70–220 kDa
     蛋白的阻滞随凝胶浓度变化很大，必须用自家凝胶实测一次（一次标定即可复用）。
     下面沉淀线位置的绝对值随这个因子移动，但浓度窗口的**形状**不随它变。""")
    return out


def step_agid(D: dict) -> None:
    h("第 2 步 — AGID：p72 抗原的沉淀线与浓度窗口")
    D_ag = D["p72 三聚体 (天然壳粒)"]["D_gel"]
    D_ab = D["猪 IgG (抗体侧)"]["D_gel"]
    print(f"  抗原 D = {D_ag:.2e} cm²/s（p72 三聚体）   抗体 D = {D_ab:.2e} cm²/s（猪 IgG）")

    base = agid.AgidSetup(L_cm=0.6, n_grid=161, D_ag=D_ag, D_ab=D_ab,
                          A_well=1.0, B_well=1.0, nu=2.0, k_p=5.0, log_width=0.9)
    r = agid.simulate_agid(base, t_end_s=24 * 3600.0, n_out=25)
    print(f"\n  孔间距 {base.L_cm} cm，24 h 后：")
    for i in (6, 12, 18, 24):
        print(f"    {r.t[i] / 3600:5.1f} h   带心 {r.band_position(i):.4f} cm   "
              f"带宽 SD {r.band_sharpness(i):.4f} cm")

    print("\n  沉淀线位置随抗原孔浓度的移动（诊断上的意义：线偏了说明比例不对）:")
    for a in (0.1, 0.3, 1.0, 3.0, 10.0):
        rr = agid.simulate_agid(dataclasses.replace(base, A_well=a),
                                t_end_s=24 * 3600.0, n_out=4)
        print(f"    A/B = {a:5.1f}   带心 {rr.band_position():.4f} cm   "
              f"带宽 SD {rr.band_sharpness():.4f} cm")

    h("第 3 步 — 试管沉淀：p72 抗原的可用浓度窗口（替代棋盘滴定）")
    print("  固定抗体量 = 1，扫抗原量。ν = 每个 p72 三聚体上被同时结合的 IgG 数。")
    for nu in (2.0, 4.0):
        # 网格必须细到能分辨窗口本身，否则报出来的宽度是分辨率假象而非物理
        hk = agid.heidelberger_kendall_curve(
            [10 ** x for x in np.linspace(-2.5, 2.0, 181)], B_total=1.0, nu=nu)
        grid_log = 4.5 / 180
        peak = max(hk, key=lambda d: d["precipitate"])
        ok = [d for d in hk if d["precipitate"] > 0.5 * peak["precipitate"]]
        lo, hi = ok[0]["ratio"], ok[-1]["ratio"]
        width = math.log10(hi / lo)
        print(f"\n    ν = {nu:.0f}:  峰在 A/B = {peak['ratio']:.3g}（理论 1/ν = {1 / nu:.3g}）")
        print(f"            半高浓度窗口 A/B ∈ [{lo:.3g}, {hi:.3g}]，"
              f"跨 {width:.2f} 个 log（网格分辨率 {grid_log:.3f} log，"
              f"窗口是它的 {width / grid_log:.0f} 倍，非分辨率假象）")
        print(f"            超出上界后抗体进入可溶复合物：A/B={hk[-1]['ratio']:.3g} 时 "
              f"可溶 {hk[-1]['soluble']:.3g} vs 沉淀 {hk[-1]['precipitate']:.2e}")
    print("""
  这是本案例里**唯一能直接省钱**的输出：不做棋盘滴定，先算出该配多少抗原。
  ⚠️ 但窗口的绝对位置取决于 ν 与等价带宽度，二者都需要一次实测标定；
     模型给的是"窗口有多宽、峰在哪个比例"，不是"配 3.7 µg/mL"。""")


def step_agglutination() -> None:
    h("第 4 步 — p72 包被乳胶的平板凝集")
    kd = 5.0  # nM，占位：抗 p72 单抗的典型量级，必须实测
    print(f"  假定抗 p72 单抗 Kd = {kd} nM（⚠️ 占位值，必须 SPR/BLI 实测）")
    print(f"  最优抗体浓度 = Kd = {ag.optimal_ab_concentration(kd)} nM（θ=0.5 时桥联效率取 1）\n")
    concs = [kd * 10 ** x for x in (-2, -1, -0.5, 0, 0.5, 1, 2)]
    pc = ag.prozone_curve(concs, kd=kd, K0=2e-3, n0=1e3, t_read_s=6.0, k_max=40)
    print(f"  {'[Ab]/Kd':>9} {'θ':>7} {'桥联效率':>9} {'可见凝集':>9}")
    for d in pc:
        print(f"  {d['ab'] / kd:9.3g} {d['theta']:7.3f} "
              f"{ag.bridging_efficiency(d['theta']):9.3f} {d['visible']:9.3f}")
    lo = min(d["ab"] / kd for d in pc if d["visible"] > 0.5 * max(x["visible"] for x in pc))
    hi = max(d["ab"] / kd for d in pc if d["visible"] > 0.5 * max(x["visible"] for x in pc))
    print(f"\n  半高可见窗口: [Ab]/Kd ∈ [{lo:.3g}, {hi:.3g}]，跨 {math.log10(hi / lo):.2f} 个 log")
    print("""  诊断上的含义：单一抗体浓度出阴性无法区分"抗体太少"和"抗体太多"，
  所以 p72 乳胶凝集必须做系列稀释 —— 这是模型的推论，不是操作惯例。""")


def step_occupancy() -> None:
    h("第 5 步 — 占据模型：数学跑得通，生物学不成立")
    print(f"  真实化学计量（Liu 2019 Cell Res, PMID 31649031）:")
    print(f"    T = {CAPSID_T}，衣壳共 {CAPSID_TOTAL_PROTEINS:,} 个蛋白")
    print(f"    p72 拷贝数 n = {P72_COPIES:,}，以 {P72_TRIMERS:,} 个三聚体壳粒排列\n")
    for label, n in (("按 p72 单体计", P72_COPIES), ("按三聚体壳粒计", P72_TRIMERS)):
        v = neu.Virion("ASFV p72", n_spikes=n, k_hits=1)
        amp = neu.amplification_factor(v)
        print(f"    {label:16s} n={n:6,d}   NT50/Kd = {amp:.3e}   即 Kd/{1 / amp:,.0f}")
    print(f"""
  这个数**不要用**。理由在第 0 步：p72 在囊膜之内，胞外抗体够不到。
  占据模型假设"抗体能结合到位点上"，这个前提对 p72 在胞外病毒粒子上不成立。
  模型算出 Kd/{1 / neu.amplification_factor(neu.Virion('x', P72_COPIES, 1)):,.0f} 这样一个漂亮的放大倍数，
  但它描述的是一个不存在的中和过程。

  这一步的价值恰恰在于：它演示了**模型不会拒绝无效输入**。""")


def step_ode() -> None:
    h("第 6 步 — 单孔 ODE：结构适用，参数必须标定")
    print("""  ASFV 长在原代猪肺泡巨噬细胞 (PAM)，不是传代细胞系。文献里的一步生长曲线
  潜隐期约 8–12 h、24–48 h 达峰，但 β/δ/p/c 四个速率常数没有可直接引用的公认值，
  且不可从别的病毒×细胞系迁移过来。

  下面用一组**量级示意**参数说明标定后能拿到什么，绝不是 ASFV 的真实参数：""")
    par = vd.CellVirusParams(
        name="ASFV / PAM（量级示意，未标定）",
        beta=2e-7, k_eclipse=1 / 10.0, delta=1 / 20.0, p=3.0, c=1 / 10.0, T0=1e5,
        source="⚠️ 占位参数，仅用于演示标定后的产出形态")
    print(f"\n    R₀ = {par.R0:.2f}    渐近增长率 λ = {par.growth_rate():.4f} /h "
          f"（倍增时间 {math.log(2) / par.growth_rate():.1f} h）")
    tr = vd.simulate(par, moi=0.01, t_end_h=120)
    for t in (24, 48, 72, 96):
        print(f"    {t:3d} h   CPE {tr.cpe_at(t):.3f}")

    drug = vd.Drug("假想复制抑制剂", moa="replication", ec50_mol=100.0, hill=1.0)
    print("\n  同一分子层效力 (EC50 = 100 nM) 在不同实验设计下测到的孔水平 EC50:")
    for moi in (0.001, 0.01, 0.1):
        for rt in (48.0, 72.0):
            res = vd.apparent_ec50(par, drug, moi=moi, readout_h=rt)
            print(f"    MOI {moi:<6g} 读板 {rt:>4.0f} h  ->  表观 EC50 {res['ec50_apparent']:8.2f} nM "
                  f"（位移 {res['shift_log10']:+.2f} log，表观 Hill {res['hill_apparent']:.2f}）")
    print("""
  标定之后这张表就是**实验设计工具**：先在模型里选 MOI 和读板点，再上手做。
  未标定时它只说明一件事 —— 孔水平 EC50 有很大一部分是设计的产物。""")


def step_bridge() -> None:
    h("第 7 步 — 非细胞→细胞桥接：明确不出预测")
    print("""  ChEMBL_37 检索结果：ASFV 靶点 (CHEMBL613714, ORGANISM) 带 pchembl_value
  的活性记录 **0 条**。

  这意味着三件事同时成立：
    1. 无法为 ASFV 标定桥接的斜率与截距
    2. 无法用留出集验证任何 ASFV 预测
    3. 拿别的病毒（HIV / 流感）标定好的桥接迁移过来 —— 我们已经测过代价：
       RMSE 从 1.01 涨到 1.49，外加 +0.84 log（约 7 倍）的系统性偏差

  所以本案例对 ASFV 的小分子细胞效力**不出任何预测**。
  这不是模型跑不动，是没有可验证的依据。""")


def main() -> None:
    print(RULE)
    print("assaysim 案例：ASFV p72 (B646L) 全链路")
    print(f"UniProt P22776 · {P72_LEN} aa · {P72_MW:,.0f} Da · "
          f"衣壳 {P72_COPIES:,} 拷贝 (T={CAPSID_T})")
    print(RULE)
    applicability()
    D = step_diffusion()
    step_agid(D)
    step_agglutination()
    step_occupancy()
    step_ode()
    step_bridge()
    h("小结")
    print("""  6 段链路里：2 段直接适用并给出可用输出（AGID 浓度窗口、凝集稀释窗口），
  2 段结构适用但参数缺（占据模型、单孔 ODE），1 段明确不出预测（桥接），
  1 段属诊断侧未在本案例展开（引物设计）。

  对 p72 这样一个**诊断抗原**，链路真正的落点是免疫沉淀/凝集那一侧，
  不是抗病毒那一侧 —— 这与它在 ASFV 生物学里的角色一致。""")


if __name__ == "__main__":
    main()
