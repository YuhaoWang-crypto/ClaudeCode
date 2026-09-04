"""
M8 — 选择性指数 SI = CC50/EC50 为什么预测不了：噪声传播的定量说明

背景
----
两份独立工作对同一个端点给出了相反的结论：
  · 一份把 SI 直接建模，三个模型 R² 全部 ≤ 0，判定"明确不出预测"
  · 另一份把 CC50 单独建模，骨架 CV Spearman 达到 0.535

作者把差异归因于"数据拼装方式不同"。本模块检验另一个更简单的解释：

    log SI = log CC50 − log EC50

    · **噪声叠加**：Var(noise_SI) = Var(noise_CC50) + Var(noise_EC50)
    · **信号可能反而收缩**：若 CC50 与 EC50 跨化合物正相关（强效化合物
      往往也更毒），相减会抵消掉共同变异，Var(signal_SI) 变小

两头一起挤，SI 的信噪比可以远低于 CC50，即使二者用完全相同的建模流程。

可达 R² 上限
------------
任何预测器都不可能优于目标量自身的重现性：

    R²_max = 1 − Var(noise) / Var(total)

Var(noise) 由**同一化合物在不同论文间的重复测量**估计。这是跨论文重现性，
包含实验室间差异，正是任何用文献数据训练的模型所面对的噪声。

⚠️ 本模块给的是**上限**，不是某个具体模型的性能。R²_max ≤ 0 意味着
   "该端点在这批数据上不可预测"；R²_max 高不保证某个模型一定能达到。

实测结论（ChEMBL_37，HIV 体系 6,983 个配对化合物）
--------------------------------------------------
    端点     总 SD    噪声 SD   R²上限
    CC50     0.620    0.185    0.911
    EC50     1.273    0.616    0.766
    SI       1.207    0.643    0.716

假说**方向成立但量级不足**：噪声传播确实把上限从 0.91 压到 0.72
（CC50 与 EC50 跨化合物 r=+0.347，把 SI 的总 SD 从独立情形下的 1.416
压到 1.207，收缩 15%），但 0.72 远不足以解释一个实测 R² ≤ 0 的模型。

=> 那个负结果的主因**不在**噪声传播，得回到数据拼装口径、留出划分或模型本身。
   本模块因此是一次**否定自己假说**的记录，不是佐证。

运行:  python3 -m assaysim.si_noise
"""

from __future__ import annotations

import json
import math
import re
import statistics as st
from pathlib import Path

import numpy as np

from .chembl_data import CELLULAR_FORMATS, _usable, collapse_per_molecule, fetch_target_activities

CC50_CACHE = Path("/tmp/cc50_rows.json")

_CELL_PATTERNS = {
    "MT-4": r"\bMT-?4\b", "CEM": r"\bCEM\b", "Vero E6": r"\bVero[- ]?E6\b",
    "Huh-7": r"\bHuh-?7\b", "HepG2": r"\bHep\s?G-?2\b", "MDCK": r"\bMDCK\b",
    "HeLa": r"\bHeLa\b", "A549": r"\bA-?549\b",
}
_CELL_RE = {k: re.compile(v, re.I) for k, v in _CELL_PATTERNS.items()}
_no_match = re.compile(r"(?!x)x")

# 极差 -> SD 的换算因子（n≈5 时极差的期望约为 2.5σ）
RANGE_TO_SD = 2.5


def load_cc50_by_molecule(min_n: int = 1, cell: str | None = None) -> dict[str, dict]:
    """从 CC50 原始记录聚合到每个分子的 log10(CC50 nM) 中位数与离散度。

    cell=None 时把**所有细胞系混在一起**（对应"泛细胞毒性"模型的设定）。
    ⚠️ 混合时，细胞系之间的真实差异会被算进"噪声"，因为对一个不区分细胞系的
       模型而言它确实不可约。要得到干净的**测量**噪声，需指定单一细胞系。
    """
    if not CC50_CACHE.exists():
        raise FileNotFoundError(
            f"缺少 {CC50_CACHE}；先运行抓取脚本把 ChEMBL CC50 记录缓存下来")
    rows = json.loads(CC50_CACHE.read_text())

    by_mol: dict[str, list[float]] = {}
    for a in rows:
        if a.get("standard_units") != "nM":
            continue
        if cell is not None:
            # ⚠️ activity 端点**没有** assay_cell_type 字段（细胞系在 assay 端点的
            #    cell_chembl_id）。早先版本按该字段过滤，结果 0/107782 全被滤掉。
            #    这里按 assay_description 文本匹配，与被审计方的策展口径一致。
            if not _CELL_RE.get(cell, _no_match).search(a.get("assay_description") or ""):
                continue
        if a.get("standard_relation") not in ("=", None):
            continue
        v = a.get("standard_value")
        mol = a.get("molecule_chembl_id")
        if not mol or v in (None, ""):
            continue
        try:
            f = float(v)
        except (TypeError, ValueError):
            continue
        if f <= 0:
            continue
        by_mol.setdefault(mol, []).append(math.log10(f))

    out = {}
    for mol, vals in by_mol.items():
        if len(vals) < min_n:
            continue
        out[mol] = {"log": st.median(vals), "n": len(vals),
                    "spread": (max(vals) - min(vals)) if len(vals) > 1 else 0.0}
    return out


def noise_sd(entries: dict, min_n: int = 3) -> float:
    """由重复测量的极差中位数估计噪声 SD（log10 单位）。"""
    sp = [e["spread"] for e in entries.values() if e["n"] >= min_n]
    if not sp:
        return float("nan")
    return st.median(sp) / RANGE_TO_SD


def ceiling_r2(total_sd: float, noise_sd_: float) -> float:
    """R²_max = 1 − Var(noise)/Var(total)。"""
    if total_sd <= 0:
        return float("nan")
    return 1.0 - (noise_sd_ ** 2) / (total_sd ** 2)


def analyse(verbose: bool = True, cell: str | None = None) -> dict:
    """对 HIV-1 这套体系做 CC50 / EC50 / SI 三个端点的信噪比对比。

    cell 指定后只用该细胞系的 CC50，去掉细胞系间差异这一项。
    """
    cc = load_cc50_by_molecule(cell=cell)
    cell_acts = fetch_target_activities("CHEMBL378", verbose=False)
    ec = collapse_per_molecule(cell_acts, {"EC50", "IC50"}, CELLULAR_FORMATS)

    shared = sorted(set(cc) & set(ec))
    if verbose:
        print(f"CC50 有记录的分子      : {len(cc):,}")
        print(f"HIV 细胞法 EC50 的分子 : {len(ec):,}")
        print(f"两者都有（可算 SI）    : {len(shared):,}\n")

    # 端点数值（log10 nM 尺度；EC50 侧的 pchembl = 9 − log10(nM)，先换回来）
    log_cc = np.array([cc[m]["log"] for m in shared])
    log_ec = np.array([9.0 - ec[m]["pchembl"] for m in shared])
    log_si = log_cc - log_ec  # SI = CC50/EC50

    cc_sub = {m: cc[m] for m in shared}
    ec_sub = {m: {"spread": ec[m]["spread"], "n": ec[m]["n"]} for m in shared}

    n_cc = noise_sd(cc_sub)
    n_ec = noise_sd(ec_sub)
    n_si = math.sqrt(n_cc ** 2 + n_ec ** 2)  # 独立误差叠加

    r = float(np.corrcoef(log_cc, log_ec)[0, 1])
    res = {
        "n_pairs": len(shared),
        "corr_cc_ec": r,
        "endpoints": {
            "CC50": {"total_sd": float(np.std(log_cc)), "noise_sd": n_cc},
            "EC50": {"total_sd": float(np.std(log_ec)), "noise_sd": n_ec},
            "SI":   {"total_sd": float(np.std(log_si)), "noise_sd": n_si},
        },
    }
    for k, v in res["endpoints"].items():
        v["r2_max"] = ceiling_r2(v["total_sd"], v["noise_sd"])
        v["snr"] = (v["total_sd"] ** 2 - v["noise_sd"] ** 2) / v["noise_sd"] ** 2

    if verbose:
        print(f"log10 CC50 与 log10 EC50 跨化合物相关系数 r = {r:+.3f}")
        print("  （正相关 => 相减会抵消共同变异，SI 的信号被压缩）\n")
        print(f"  {'端点':<6} {'总 SD':>8} {'噪声 SD':>9} {'信噪比':>8} {'R²上限':>9}")
        for k, v in res["endpoints"].items():
            print(f"  {k:<6} {v['total_sd']:8.3f} {v['noise_sd']:9.3f} "
                  f"{v['snr']:8.2f} {v['r2_max']:9.3f}")
        # 若 CC50 与 EC50 完全独立，SI 的总方差会是多少
        indep = math.sqrt(np.var(log_cc) + np.var(log_ec))
        print(f"\n  若二者独立，SI 的总 SD 应为 {indep:.3f}；"
              f"实测 {res['endpoints']['SI']['total_sd']:.3f}"
              f"（被相关性压缩了 {100 * (1 - res['endpoints']['SI']['total_sd'] / indep):.0f}%）")
    return res


def main() -> None:
    print("=" * 74)
    print("SI = CC50/EC50 的可预测性上限 — 噪声传播分析（ChEMBL_37 真实数据）")
    print("=" * 74 + "\n")
    print("--- 混合全部细胞系（对应『泛细胞毒性』模型的设定） ---")
    res = analyse()
    e = res["endpoints"]
    for c in ("MT-4", "CEM"):
        try:
            sub = analyse(verbose=False, cell=c)
        except Exception:
            continue
        if sub["n_pairs"] < 50:
            continue
        se = sub["endpoints"]
        print(f"\n--- 仅 {c}（去掉细胞系间差异） N={sub['n_pairs']:,} ---")
        print(f"  {'端点':<6} {'总 SD':>8} {'噪声 SD':>9} {'R²上限':>9}")
        for k, v in se.items():
            print(f"  {k:<6} {v['total_sd']:8.3f} {v['noise_sd']:9.3f} {v['r2_max']:9.3f}")
    print("\n" + "=" * 74)
    print("结论")
    print("=" * 74)
    si_ceil = e["SI"]["r2_max"]
    cc_ceil = e["CC50"]["r2_max"]
    lower = si_ceil < cc_ceil
    # 被审计方报告的 SI 实测 R² <= 0。噪声传播要能解释它，上限本身必须逼近 0。
    sufficient = si_ceil < 0.1
    print(f"  CC50 R²上限 {e['CC50']['r2_max']:.3f}   SI R²上限 {e['SI']['r2_max']:.3f}")
    print(f"  方向：SI 上限{'低于' if lower else '不低于'} CC50 —— 假说方向"
          f"{'成立' if lower else '不成立'}")
    print(f"  充分性：SI 上限 {si_ceil:.3f}，而被审计方报告的实测 R² ≤ 0。")
    if not sufficient:
        print(f"""
  => **噪声传播不足以解释那个负结果。**
     它把上限从 {cc_ceil:.2f} 压到 {si_ceil:.2f}，但一个报 R² ≤ 0 的模型
     离这个上限还差得远。所以差异的主因在别处 —— 数据拼装口径、
     留出划分方式，或该模型本身。我原先"两个噪声量相除所以测不准"
     的推测**方向对、量级不够**，不能作为解释。""")
    else:
        print("\n  => 噪声传播足以解释：SI 端点本身就不可预测。")
    print("""
  ⚠️ 这是**上限**，不是任何具体模型的实测性能。上限低说明端点本身不可预测；
     上限高不保证某个模型能达到。噪声由跨论文重复测量估计，含实验室间差异。""")


if __name__ == "__main__":
    main()
