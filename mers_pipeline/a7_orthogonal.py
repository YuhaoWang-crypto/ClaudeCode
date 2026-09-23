"""A7 —— 正交细胞层检查 与 药物重定位模块的可行性评估。

方案原设想
----------
A5「正交细胞水平层」：对同时具有蛋白酶 IC50 和细胞抗病毒 EC50 的化合物，
检验预测蛋白酶效力与细胞抗病毒活性的相关性。
A4「药物重定位筛选」：用 A3 选出的模型 + 对接一致性评分扫描已批准/临床药物库，
输出三层 applicability-domain 分级。

数据现实
--------
1) 正交层：MERS-CoV 的 132 个细胞水平活性化合物与 19 个 Mpro 活性化合物
   **只有 1 个交集**（PLpro 为 0）。n=1 无法计算相关性 —— 该分析不可执行。
2) 重定位：其前提是"已选出一个经验证的模型"。但 A6 明确拒绝了模型选择
   （两种策略的 Spearman 95% CI 都跨 0 且互相重叠），A3 的 redock 也未通过门槛。
   在没有经验证打分函数的情况下输出排序清单，等于给出无依据的预测。

本模块因此改为做三件**可执行且有信息量**的事：
  (1) 如实量化正交层的重叠情况，说明为什么相关性算不出来；
  (2) 用 applicability domain 分析**定量说明**重定位预测会有多不可靠 ——
      把候选药物库与 MERS 训练数据的化学空间距离算出来；
  (3) 给出一个机制性问题的可检验形式：细胞水平抗 MERS 活性能否由蛋白酶抑制解释。

产出：results/trackA/a7_orthogonal.json, figures/mers/a7_applicability.png
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator

from .config import FIGS, RES_A, log, setup_cjk_fonts

RDLogger.DisableLog("rdApp.*")
_MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def fps(smiles: list[str]):
    out = []
    for s in smiles:
        m = Chem.MolFromSmiles(s)
        out.append(_MORGAN.GetFingerprint(m) if m is not None else None)
    return [f for f in out if f is not None]


def nearest_neighbour_similarity(query_fps, ref_fps) -> np.ndarray:
    """每个查询分子到参考集的最大 Tanimoto 相似度。"""
    return np.array([max(DataStructs.BulkTanimotoSimilarity(q, ref_fps)) for q in query_fps])


def ad_tier(sim: float) -> str:
    """三层 applicability domain 分级（方案要求的分层口径）。"""
    if sim >= 0.50:
        return "high"
    if sim >= 0.30:
        return "borderline"
    return "out-of-domain"


def main() -> dict:
    mpro = pd.read_csv(RES_A / "mers_mpro_ic50.csv")
    plpro = pd.read_csv(RES_A / "mers_plpro_ic50.csv")
    cell = pd.read_csv(RES_A / "mers_cell_antiviral.csv")
    bench = pd.read_csv(RES_A / "benchmark_pooled.csv")

    # ---------- (1) 正交层重叠 ----------
    ov_mpro = set(mpro["smiles_std"]) & set(cell["smiles_std"])
    ov_plpro = set(plpro["smiles_std"]) & set(cell["smiles_std"])
    ov_bench = set(bench["smiles_std"]) & set(cell["smiles_std"])
    log(f"正交层重叠: Mpro活性∩细胞 {len(ov_mpro)} | PLpro∩细胞 {len(ov_plpro)} | "
        f"基准集(含非活性)∩细胞 {len(ov_bench)}")

    pairs = []
    mm, cc = mpro.set_index("smiles_std"), cell.set_index("smiles_std")
    for s in ov_mpro:
        pairs.append({
            "molecule_chembl_id": mm.loc[s, "molecule_chembl_id"],
            "p_protease": round(float(mm.loc[s, "pActivity"]), 2),
            "p_cell": round(float(cc.loc[s, "pActivity"]), 2),
            "cell_endpoint": cc.loc[s, "standard_type"],
        })
    orthogonal = {
        "n_protease_compounds": int(len(mpro)),
        "n_cell_compounds": int(len(cell)),
        "n_overlap_mpro": len(ov_mpro),
        "n_overlap_plpro": len(ov_plpro),
        "paired_compounds": pairs,
        "correlation_computable": len(ov_mpro) >= 5,
        "verdict": (
            f"只有 {len(ov_mpro)} 个化合物同时具备 MERS-CoV Mpro 生化 IC50 与细胞抗病毒效力"
            f"（PLpro 为 {len(ov_plpro)} 个）。方案设想的'蛋白酶效力 vs 细胞活性相关性'"
            f"在公开数据上**无法计算**。这本身是有信息量的结果：MERS-CoV 的生化数据与"
            f"细胞数据几乎来自完全不同的化合物集合，机制一致性在现有公开数据下无法闭环。"
        ),
    }

    # ---------- (2) 重定位的 applicability domain ----------
    # 候选库：MERS 细胞水平活性化合物（有实测抗病毒活性的真实药物样分子）
    # 参考集：MERS Mpro 基准集（模型若存在，只能在这个化学空间里可信）
    ref = fps(bench["smiles_std"].tolist())
    lib_smiles = cell["smiles_std"].tolist()
    lib = fps(lib_smiles)
    sims = nearest_neighbour_similarity(lib, ref)
    tiers = [ad_tier(s) for s in sims]
    tier_counts = pd.Series(tiers).value_counts().to_dict()
    log(f"AD 分层（{len(sims)} 个候选 vs MERS Mpro 基准集 {len(ref)} 个）: {tier_counts}")
    log(f"  到训练集最近邻 Tanimoto: 中位 {np.median(sims):.3f}  "
        f"90分位 {np.percentile(sims,90):.3f}  最大 {sims.max():.3f}")

    # 训练集自身的内部相似度，作为"什么算近"的参照尺度
    self_sims = []
    for i, f in enumerate(ref):
        others = ref[:i] + ref[i + 1:]
        if others:
            self_sims.append(max(DataStructs.BulkTanimotoSimilarity(f, others)))
    self_sims = np.array(self_sims)
    log(f"  基准集内部最近邻 Tanimoto: 中位 {np.median(self_sims):.3f}")

    repurposing = {
        "prerequisite_met": False,
        "why": (
            "重定位筛选依赖'A3 选出的模型'，但 A6 已按预先声明的口径拒绝模型选择"
            "（两种数据策略的 scaffold-split Spearman 95% CI 都跨过 0 且互相重叠），"
            "且 A5 的 redock 验证未通过 2 Å 门槛。没有经验证的打分函数，"
            "输出排序候选清单就是给出无依据的预测，因此本模块**不产出候选药物清单**。"
        ),
        "applicability_domain_analysis": {
            "library": "MERS-CoV 细胞水平活性化合物（ChEMBL CHEMBL4296578）",
            "n_library": int(len(lib)),
            "reference_set": "MERS-CoV Mpro 基准集",
            "n_reference": int(len(ref)),
            "tier_counts": tier_counts,
            "tier_thresholds": {"high": ">=0.50", "borderline": "0.30-0.50",
                                "out_of_domain": "<0.30"},
            "nn_tanimoto": {
                "median": round(float(np.median(sims)), 3),
                "p90": round(float(np.percentile(sims, 90)), 3),
                "max": round(float(sims.max()), 3),
            },
            "reference_internal_nn_tanimoto_median": round(float(np.median(self_sims)), 3),
            "interpretation": (
                f"候选库中 {tier_counts.get('out-of-domain', 0)}/{len(sims)} 个分子落在"
                f"域外（到训练集最近邻 Tanimoto < 0.30）。即便存在一个有效模型，"
                f"对这部分分子的预测也只是外推。这为'为什么不输出候选清单'提供了定量依据。"
            ),
        },
    }

    # ---------- (3) 可检验的机制性问题 ----------
    mechanism = {
        "question": "MERS-CoV 细胞水平抗病毒活性能否由 Mpro/PLpro 抑制解释？",
        "testable_form": (
            "把 132 个细胞水平活性化合物与 60 个**实测蛋白酶非活性**化合物一起对接到 Mpro，"
            "检验前者的对接分数是否系统性优于后者。若不优于，说明这些化合物的细胞活性"
            "主要不是通过 Mpro 抑制实现的。"),
        "status": "已给出可执行形式；其结论受限于 A5 中 redock 未通过这一前提，"
                  "即便出现差异也只能作探索性解读。",
        "n_cell_active": int(len(cell)),
        "n_measured_protease_inactive": int((bench["active"] == 0).sum()),
    }

    out = {"orthogonal_cell_layer": orthogonal,
           "drug_repurposing": repurposing,
           "mechanism_question": mechanism}
    (RES_A / "a7_orthogonal.json").write_text(json.dumps(out, indent=2, ensure_ascii=False))
    _plot(sims, self_sims, tier_counts)
    log(f"A7 完成 -> {RES_A/'a7_orthogonal.json'}")
    return out


def _plot(sims: np.ndarray, self_sims: np.ndarray, tier_counts: dict) -> None:
    setup_cjk_fonts()
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
    ax = axes[0]
    ax.hist(sims, bins=24, range=(0, 1), color="#2b6cb0", alpha=0.85,
            label=f"重定位候选库 (n={len(sims)})")
    ax.hist(self_sims, bins=24, range=(0, 1), color="#dd6b20", alpha=0.6,
            label=f"基准集内部 (n={len(self_sims)})")
    for x, c, lab in [(0.30, "#e53e3e", "域外门槛 0.30"), (0.50, "#38a169", "high 门槛 0.50")]:
        ax.axvline(x, color=c, ls="--", lw=1.5, label=lab)
    ax.set_xlabel("到 MERS Mpro 基准集的最近邻 Tanimoto")
    ax.set_ylabel("化合物数")
    ax.set_title("化学空间距离")
    ax.legend(frameon=False, fontsize=8.5)
    ax.grid(alpha=0.25, axis="y")

    ax = axes[1]
    order = ["high", "borderline", "out-of-domain"]
    vals = [tier_counts.get(o, 0) for o in order]
    cols = ["#38a169", "#dd6b20", "#e53e3e"]
    bars = ax.bar(["high", "borderline", "域外"], vals, color=cols, width=0.6)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 1, str(v), ha="center", fontsize=11)
    ax.set_ylabel("化合物数")
    ax.set_title("applicability domain 三层分级")
    ax.grid(alpha=0.25, axis="y")

    fig.suptitle("A7：重定位预测的适用域分析（说明为何不输出候选清单）", fontsize=12.5)
    fig.tight_layout()
    out = FIGS / "a7_applicability.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log(f"图已保存 -> {out}")


if __name__ == "__main__":
    main()
