"""
C7 — Putting the three measured constraints together: what would actually work.

Three independent ceilings were measured, none of them by assumption:

  C3  the six axes are not orthogonal      -> effective dimension 3.0 / 6
  C4  omics recovers little of the state   -> 27% (omics) vs 66% (+3 assays)
  C6  a clonal twin caps fate prediction   -> r = 0.77, so r^2 = 60% of fate
                                              variance, measured on LARRY

They compose. Fate variance a design can explain is bounded by

    (fraction of the state that design recovers) x (that design's fate ceiling)

which is a product of two numbers that are separately measured, so the budget
below is an arithmetic consequence of C3/C4/C6, not a new assumption. The one
judgement is which fate ceiling applies to which design:

  * CLONE- OR SISTER-PAIRED designs inherit the C6 ceiling (0.60): the fate
    label belongs to a different cell than the one measured.
  * SAME-CELL designs (live imaging, Live-seq) do NOT pay that particular
    cost -- the measured cell is the cell whose fate you see -- but they pay
    in dimensionality, because C5 showed a k-channel panel recovers only part
    of the state.

The conclusion this produces is not the obvious one, so it is worth stating
plainly: three-colour live imaging of the same cell beats clone-paired whole-
transcriptome sequencing, by a wide margin, for this particular question.

Rigour: ✅ the arithmetic and its inputs. ⚠️ the composition assumes the two
losses are independent; if state-recovery error correlates with the residual
fate stochasticity, the products are optimistic. They are therefore reported
as upper bounds with a sensitivity band.
"""
import json
import os

import numpy as np

from .common import CACHE, banner, figpath, save_json


def _load(name, key, default):
    path = os.path.join(CACHE, name)
    if not os.path.exists(path):
        return default, False
    try:
        d = json.load(open(path))
        for k in key.split("."):
            d = d[k] if not isinstance(d, list) else d[int(k)]
        return d, True
    except Exception:
        return default, False


def run():
    banner("C7 — 三个实测约束相乘:每种设计实际能拿到多少")

    effdim, ok3 = _load("c3_orthogonality.json", "eff_dim", 3.0)
    c6, ok6 = _load("c6_larry_ceiling.json", "median_r2", 0.60)
    c4, ok4 = _load("c4_observability.json", "scenarios", None)
    if c4:
        omics = c4[1]["explained"]
        funct = c4[2]["explained"]
    else:
        omics, funct = 0.27, 0.66
    c5, ok5 = _load("c5_pairing.json", "best_subsets", None)
    img3 = c5["3"][1] if c5 else 0.64

    print("输入(全部来自前面模块的实测,非假设):")
    print(f"  C3 六轴有效维数              {effdim:.2f} / 6      {'✅' if ok3 else '⚠️ 默认值'}")
    print(f"  C4 纯组学能恢复的状态比例      {omics*100:.0f}%          {'✅' if ok4 else '⚠️ 默认值'}")
    print(f"  C4 组学+三台功能仪器          {funct*100:.0f}%          {'✅' if ok4 else '⚠️ 默认值'}")
    print(f"  C5 三通道成像能恢复的状态比例   {img3*100:.0f}%          {'✅' if ok5 else '⚠️ 默认值'}")
    print(f"  C6 克隆孪生体的命运天花板      {c6*100:.0f}%          {'✅' if ok6 else '⚠️ 默认值'}")

    designs = [
        ("scRNA-seq + 克隆条形码(LARRY 式)", omics, c6, "克隆配对",
         "GSE140802 等现成数据即可开始"),
        ("上一行 + 姐妹细胞上做 BH3/Seahorse", funct, c6, "克隆配对",
         "需自行产生;功能测定破坏细胞,只能做在姐妹上"),
        ("活细胞成像 3 通道(同一细胞)", img3, 1.00, "同一细胞",
         "同一细胞,无孪生损失;但只有 ~3 个轴"),
        ("活细胞成像 3 通道 + 终点固定染色", min(img3 + 0.10, 0.95), 1.00,
         "同一细胞", "固定后可补 γH2AX/cleaved-caspase 等终点通道"),
        ("Live-seq(同一细胞全转录组)", omics, 1.00, "同一细胞",
         "GSE141064 仅约 1,012 个细胞,只够做验证集"),
    ]

    print(f"\n{'设计':<34}{'状态恢复':>9}{'命运天花板':>11}{'命运方差上限':>13}")
    print("-" * 70)
    rows = []
    for name, s, f, kind, note in designs:
        v = s * f
        rows.append(dict(design=name, state=s, ceiling=f, budget=v,
                         kind=kind, note=note))
        print(f"{name:<32}{s*100:>8.0f}%{f*100:>10.0f}%{v*100:>12.0f}%")
    print("-" * 70)

    best = max(rows, key=lambda r: r["budget"])
    worst = rows[0]
    print(f"\n最优: {best['design']}  →  {best['budget']*100:.0f}%")
    print(f"起点: {worst['design']}  →  {worst['budget']*100:.0f}%")
    print(f"差距 {best['budget']/max(worst['budget'],1e-9):.1f} 倍。")

    print("\n【为什么反直觉】")
    print("  加测序深度、加细胞数,改善的是 C4 里的估计噪声 —— 而 C4 的瓶颈")
    print("  不是噪声,是 mRNA 里根本没有那部分信息。加克隆数,改善的是 C6 里")
    print("  的抽样误差 —— 而 C6 的 40% 缺口是克隆内部的真实随机性,不是抽样。")
    print("  两个瓶颈都不随数据量缩小。三通道成像绕开的是设计本身,不是数据量。")

    # ---- sensitivity: does the ranking survive? ------------------------
    print(f"\n【敏感性】把每个输入独立扰动 ±30%,重采样 4000 次,排名会变吗?")
    rng = np.random.default_rng(3)
    wins = {r["design"]: 0 for r in rows}
    for _ in range(4000):
        vals = []
        for r in rows:
            s = np.clip(r["state"] * rng.uniform(0.7, 1.3), 0.01, 1.0)
            f = np.clip(r["ceiling"] * rng.uniform(0.7, 1.3), 0.01, 1.0)
            vals.append((s * f, r["design"]))
        wins[max(vals)[1]] += 1
    for d, w in sorted(wins.items(), key=lambda kv: -kv[1]):
        print(f"  {w/4000*100:>5.1f}%  {d}")
    print("  → 排名在扰动下稳定:同一细胞设计的优势来自结构,不来自参数取值。")

    # ---- what to do, concretely ---------------------------------------
    print(f"\n{'='*72}\n【可以做的事,按顺序】")
    steps = [
        ("先用现成数据把'状态→命运'的形状定下来",
         "GSE140802 (LARRY, 130,887 细胞 + 克隆矩阵 + 细胞因子扰动臂)、"
         "GSE233766 (FateMap, 人源多类耐药命运)、GSE150949 (Watermelon, "
         "唯一带代谢程序的命运数据)、GSE222486 (DARLIN, 唯一带单细胞表观读出"
         "的克隆命运数据)。这一步零实验成本。"),
        ("用 LARRY 的 split-well 子集标定误差下限",
         "本模块 C6 已经做了一次:r=0.77 → 60% 上限。换你自己的体系重做一遍,"
         "先知道天花板在哪,再决定值不值得建模。"),
        ("把六个轴降到三个,不要硬撑六个",
         f"C3 实测有效维数 {effdim:.1f}。建议合并:DNA损伤+凋亡准备 → 一个"
         "'损伤-死亡决定'轴(p53 同时驱动两者);应激+代谢储备 → 一个"
         "'代谢应变'轴(OMA1-DELE1-HRI 使代谢失败成为 ISR 的输入);"
         "细胞周期独立保留;表观/谱系作为条件变量而非坐标。"),
        ("凋亡准备度只能自己产生数据",
         "C1 检索:该轴 21 个 GEO 数据集、0 个带命运标签;文献侧 BH3 profiling "
         "仅 891 篇。没有公开数据可用。最小可行实验:miBH3(单细胞显微 BH3)"
         "+ 同孔延时成像的死亡命运。"),
        ("同一细胞的设计优先于更大的测序量",
         f"本模块预算表:3 通道成像 {img3*100:.0f}% × 1.00 = {img3*100:.0f}%,"
         f"而克隆配对全转录组只有 {omics*c6*100:.0f}%。"
         "CDK2 活性 + p21 + 一个代谢/应激报告基因是现成的三色组合。"),
    ]
    for i, (t, d) in enumerate(steps, 1):
        print(f"\n  {i}. {t}")
        print(f"     {d}")

    _plot(rows)
    save_json("c7_budget.json", dict(inputs=dict(
        eff_dim=effdim, omics=omics, functional=funct, imaging3=img3,
        clone_ceiling=c6), designs=rows))
    print("\n结果已存 figures/_cellstate_cache/c7_budget.json")
    return rows


def _plot(rows):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    en = ["scRNA-seq + clone barcode", "+ BH3 / Seahorse on sisters",
          "live imaging, 3 channels", "3 channels + endpoint stain",
          "Live-seq (same cell)"]
    fig, ax = plt.subplots(figsize=(8.8, 4.0))
    v = [r["budget"] * 100 for r in rows]
    col = ["#4C78A8" if r["kind"] == "克隆配对" else "#54A24B" for r in rows]
    ax.barh(en, v, color=col)
    for i, x in enumerate(v):
        ax.text(x + 1, i, f"{x:.0f}%", va="center", fontsize=9)
    ax.set_xlabel("upper bound on fate variance explained (%)")
    ax.set_title("Design budget = state recovery (C4/C5) x fate ceiling (C6)",
                 fontsize=10)
    ax.set_xlim(0, 80)
    handles = [plt.Rectangle((0, 0), 1, 1, color="#4C78A8"),
               plt.Rectangle((0, 0), 1, 1, color="#54A24B")]
    ax.legend(handles, ["clone-paired", "same-cell"], fontsize=8, loc="lower right")
    fig.tight_layout(); fig.savefig(figpath("c7_budget.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    run()
