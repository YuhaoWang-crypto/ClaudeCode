"""A8 —— PLpro 对接与文献 IC50 的定性对照。

方案对 PLpro 的定位：「数据量过小（~20 IC50），不训练 ML 模型，仅做对接 +
与文献 IC50 的定性对照」。本模块执行这一项。

必须同时声明的两个前提
----------------------
1. **PLpro 没有 redock 验证**（见 A3）：4RF1/4RF0 是 PLpro–泛素复合物，
   唯一非溶剂小分子 3CN 只有 4 个重原子；4RNA 是 apo。因此位姿与打分的
   可信度显著低于 Mpro（而 Mpro 本身的 redock 也未通过 2 Å 门槛）。
2. **n = 9**：任何相关系数的置信区间都会宽到几乎不含信息。本模块照样把
   CI 算出来，正是为了让这一点显式可见，而不是只报一个点估计。

产出：results/trackA/a8_plpro.json, figures/mers/a8_plpro.png
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy import stats

from .a5_docking import screen
from .config import FIGS, RANDOM_SEED, RES_A, log, setup_cjk_fonts


def spearman_ci(x, y, n_boot=10000, seed=RANDOM_SEED):
    rng = np.random.default_rng(seed)
    x, y = np.asarray(x, float), np.asarray(y, float)
    pt = stats.spearmanr(x, y).statistic
    vals = []
    for _ in range(n_boot):
        i = rng.integers(0, len(x), len(x))
        if len(np.unique(x[i])) < 3 or len(np.unique(y[i])) < 3:
            continue
        v = stats.spearmanr(x[i], y[i]).statistic
        if np.isfinite(v):
            vals.append(v)
    return float(pt), float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def main(exhaustiveness: int = 16) -> dict:
    d = pd.read_csv(RES_A / "mers_plpro_ic50.csv")
    log(f"PLpro 化合物 {len(d)} 个（IC50/Kd），对接到 4RNA (apo)")

    dk = screen(d, "PLPRO_4RNA", exhaustiveness)
    dk.to_csv(RES_A / "docking_scores_plpro_4rna.csv", index=False)

    ok = dk[dk["best_affinity"].notna()].copy()
    log(f"成功对接 {len(ok)}/{len(dk)}")
    if len(ok) < 4:
        out = {"error": f"成功对接的化合物只有 {len(ok)} 个，无法做对照"}
        (RES_A / "a8_plpro.json").write_text(json.dumps(out, indent=2, ensure_ascii=False))
        return out

    # Vina 亲和力越负越好；pActivity 越大越好 -> 预期负相关
    sp, lo, hi = spearman_ci(ok["best_affinity"], ok["pActivity"])
    pe = stats.pearsonr(ok["best_affinity"], ok["pActivity"])
    log(f"对接亲和力 vs 实测 pActivity: Spearman={sp:+.3f} 95%CI [{lo:+.3f}, {hi:+.3f}]")
    log(f"  Pearson={pe.statistic:+.3f} p={pe.pvalue:.3f}")
    log(f"  亲和力范围 {ok['best_affinity'].min():.2f} 到 {ok['best_affinity'].max():.2f} kcal/mol"
        f"（跨度仅 {ok['best_affinity'].max()-ok['best_affinity'].min():.2f}）")
    log(f"  实测 pActivity 范围 {ok['pActivity'].min():.2f} 到 {ok['pActivity'].max():.2f}")

    out = {
        "receptor": "PLPRO_4RNA (4RNA, 1.79 Å, apo)",
        "n_compounds": int(len(ok)),
        "docking_affinity_kcal_mol": {
            "min": round(float(ok["best_affinity"].min()), 2),
            "max": round(float(ok["best_affinity"].max()), 2),
            "median": round(float(ok["best_affinity"].median()), 2),
            "span": round(float(ok["best_affinity"].max() - ok["best_affinity"].min()), 2),
        },
        "measured_pActivity": {
            "min": round(float(ok["pActivity"].min()), 2),
            "max": round(float(ok["pActivity"].max()), 2),
            "span": round(float(ok["pActivity"].max() - ok["pActivity"].min()), 2),
        },
        "docking_vs_measured": {
            "spearman": round(sp, 3), "spearman_ci95": [round(lo, 3), round(hi, 3)],
            "pearson": round(float(pe.statistic), 3), "pearson_p": round(float(pe.pvalue), 4),
            "expected_sign": "负（亲和力越负 = 结合越强 = pActivity 越大）",
        },
        "per_compound": ok[["molecule_chembl_id", "pActivity", "value_nM",
                            "best_affinity", "standard_type"]]
                        .round(3).to_dict("records"),
        "caveats": [
            "PLpro 没有 redock 验证（4RF1/4RF0 为 PLpro–泛素复合物，3CN 仅 4 个重原子；"
            "4RNA 为 apo），位姿与打分可信度显著低于 Mpro。",
            f"n={len(ok)}，Spearman 的 95% CI 宽度为 {hi-lo:.2f} —— 这个区间几乎覆盖了"
            f"相关系数的全部可能取值，因此该点估计不构成任何证据。",
            "对接亲和力的动态范围本身很窄，而 Vina 打分误差约 2–3 kcal/mol，"
            "二者同量级时打分无法区分化合物。",
            "受体为 apo 结构，未建模配体诱导的 BL2 环构象变化。",
        ],
        "verdict": (
            f"对接分数与实测效力的 Spearman = {sp:+.3f}，95% CI [{lo:+.3f}, {hi:+.3f}]。"
            f"区间跨过 0 且几乎覆盖全域，**无法据此判断对接是否捕捉到了 PLpro 的构效关系**。"
            f"这与方案对 PLpro 的预期一致（数据量过小、富集统计功效弱），"
            f"结果按定性/探索性呈现。"
        ),
    }
    (RES_A / "a8_plpro.json").write_text(json.dumps(out, indent=2, ensure_ascii=False))
    _plot(ok, sp, lo, hi)
    log(f"A8 完成 -> {RES_A/'a8_plpro.json'}")
    log(out["verdict"])
    return out


def _plot(ok: pd.DataFrame, sp: float, lo: float, hi: float) -> None:
    setup_cjk_fonts()
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    ax.scatter(ok["best_affinity"], ok["pActivity"], s=95, c="#805ad5",
               edgecolor="white", lw=1.3, zorder=3)
    for _, r in ok.iterrows():
        ax.annotate(str(r["molecule_chembl_id"]).replace("CHEMBL", ""),
                    (r["best_affinity"], r["pActivity"]),
                    fontsize=7, xytext=(4, 4), textcoords="offset points", color="#4a5568")
    ax.set_xlabel("Vina 对接亲和力 (kcal/mol，越负越强)")
    ax.set_ylabel("实测 pIC50 / pKd (MERS-CoV PLpro)")
    ax.set_title(f"A8：PLpro 对接 vs 文献实测效力（n={len(ok)}）\n"
                 f"Spearman = {sp:+.2f}  95%CI [{lo:+.2f}, {hi:+.2f}] —— 区间几乎覆盖全域",
                 fontsize=11)
    ax.grid(alpha=0.25)
    fig.text(0.5, -0.02,
             "受体 4RNA 为 apo 结构，无 redock 验证 —— 结果仅作定性/探索性呈现",
             ha="center", fontsize=9, color="#c53030")
    fig.tight_layout()
    o = FIGS / "a8_plpro.png"
    fig.savefig(o, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log(f"图已保存 -> {o}")


if __name__ == "__main__":
    import sys
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 16)
