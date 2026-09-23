"""A2 —— 跨冠状病毒迁移前提的直接检验（方案变更的依据）。

原方案设想："SARS-CoV-2 + SARS-CoV Mpro 合并训练，MERS-CoV Mpro 作为外部独立测试集"。
A1 的数据把这个设想证伪了两次：

  1) 去重后 MERS-CoV Mpro 只有 19 个带精确 IC50 的化合物，其中 17 个的
     标准化结构**完全相同**地出现在 SARS 训练集里。它们不是同系物，是同一分子。
     用它做"外部测试集"等同于在训练集上评估。
  2) 正因为是同一分子在两种酶上的配对测量，我们反而得到了一个更有价值的
     实验设计：直接测量"SARS Mpro 效力"能否预测"MERS Mpro 效力"——
     这正是任何跨冠状病毒迁移模型所依赖的前提。

本模块不训练模型，只做这个配对分析，并给出 bootstrap 置信区间与效力等价性检验。
产出：results/trackA/a2_transferability.json / a2_paired_potency.csv
      figures/mers/a2_transferability.png
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy import stats

from .config import FIGS, RANDOM_SEED, RES_A, log, setup_cjk_fonts


def bootstrap_ci(x, y, fn, n_boot: int = 10000, seed: int = RANDOM_SEED):
    """对配对样本做 bootstrap，返回 (点估计, 2.5%, 97.5%)。"""
    rng = np.random.default_rng(seed)
    x, y = np.asarray(x, float), np.asarray(y, float)
    point = fn(x, y)
    n = len(x)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        if len(np.unique(x[idx])) < 3 or len(np.unique(y[idx])) < 3:
            continue
        try:
            vals.append(fn(x[idx], y[idx]))
        except Exception:  # noqa: BLE001
            continue
    vals = np.asarray([v for v in vals if np.isfinite(v)])
    return float(point), float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def main() -> dict:
    mers = pd.read_csv(RES_A / "mers_mpro_ic50.csv")
    sars = pd.read_csv(RES_A / "sars_mpro_train.csv")

    paired = mers.merge(
        sars, on="smiles_std", suffixes=("_mers", "_sars"), how="inner"
    )[
        ["smiles_std", "molecule_chembl_id_mers", "pActivity_mers", "pActivity_sars",
         "value_nM_mers", "value_nM_sars", "n_measurements_sars", "scaffold_mers"]
    ].rename(columns={"molecule_chembl_id_mers": "molecule_chembl_id",
                      "scaffold_mers": "scaffold"})
    paired["delta_p"] = paired["pActivity_mers"] - paired["pActivity_sars"]
    paired = paired.sort_values("pActivity_mers", ascending=False).reset_index(drop=True)
    paired.to_csv(RES_A / "a2_paired_potency.csv", index=False)

    n = len(paired)
    x = paired["pActivity_sars"].to_numpy()
    y = paired["pActivity_mers"].to_numpy()
    log(f"配对化合物 n={n}（MERS Mpro 全部 {len(mers)} 个中的 {n} 个）")

    sp, sp_lo, sp_hi = bootstrap_ci(x, y, lambda a, b: stats.spearmanr(a, b).statistic)
    pe, pe_lo, pe_hi = bootstrap_ci(x, y, lambda a, b: stats.pearsonr(a, b).statistic)
    sp_p = stats.spearmanr(x, y).pvalue
    pe_p = stats.pearsonr(x, y).pvalue

    # 单调性是迁移模型真正需要的：排序能否迁移
    # Kendall tau 对小样本更稳健
    kt = stats.kendalltau(x, y)

    # 效力平移：MERS 上普遍更弱/更强？
    dp = paired["delta_p"].to_numpy()
    wilcox = stats.wilcoxon(dp)
    dp_mean, dp_lo, dp_hi = bootstrap_ci(dp, dp, lambda a, _b: float(np.mean(a)))

    # 用"SARS 效力直接当作 MERS 预测"能得到多少 R²（迁移模型的上限基准）
    ss_res = float(np.sum((y - x) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2_identity = 1 - ss_res / ss_tot
    mae_identity = float(np.mean(np.abs(y - x)))
    # 对照：只预测 MERS 均值（no-skill 基线）
    mae_mean = float(np.mean(np.abs(y - np.mean(y))))

    log(f"Spearman = {sp:+.3f}  95%CI [{sp_lo:+.3f}, {sp_hi:+.3f}]  p={sp_p:.3f}")
    log(f"Pearson  = {pe:+.3f}  95%CI [{pe_lo:+.3f}, {pe_hi:+.3f}]  p={pe_p:.3f}")
    log(f"Kendall  = {kt.statistic:+.3f}  p={kt.pvalue:.3f}")
    log(f"Δp(MERS-SARS) 均值 = {dp_mean:+.3f} 95%CI [{dp_lo:+.3f}, {dp_hi:+.3f}]  "
        f"Wilcoxon p={wilcox.pvalue:.3f}")
    log(f"直接迁移 R² = {r2_identity:+.3f} | MAE = {mae_identity:.3f} log "
        f"(vs 预测均值基线 MAE = {mae_mean:.3f})")

    # MERS 数据本身的动态范围 —— 决定任何回归任务是否可行
    span = float(mers["pActivity"].max() - mers["pActivity"].min())
    log(f"MERS Mpro pActivity 动态范围 = {span:.2f} log；标准差 = {mers['pActivity'].std():.2f}")

    verdict_transfer = (
        "不成立" if (sp_lo < 0 < sp_hi or sp <= 0.3) else "部分成立"
    )
    log(f"迁移前提判定：{verdict_transfer}")

    out = {
        "n_paired": int(n),
        "n_mers_total": int(len(mers)),
        "leakage_fraction": round(n / len(mers), 3),
        "spearman": {"point": sp, "ci95": [sp_lo, sp_hi], "p": float(sp_p)},
        "pearson": {"point": pe, "ci95": [pe_lo, pe_hi], "p": float(pe_p)},
        "kendall": {"point": float(kt.statistic), "p": float(kt.pvalue)},
        "delta_pActivity_mers_minus_sars": {
            "mean": dp_mean, "ci95": [dp_lo, dp_hi],
            "wilcoxon_p": float(wilcox.pvalue),
            "median": float(np.median(dp)),
        },
        "direct_transfer_baseline": {
            "r2": r2_identity, "mae_log": mae_identity,
            "mae_predict_mean_baseline": mae_mean,
        },
        "mers_dynamic_range_log": span,
        "mers_pActivity_sd": float(mers["pActivity"].std()),
        "verdict": {
            "external_test_set_valid": False,
            "reason": f"{n}/{len(mers)} 的 MERS Mpro 化合物结构完全出现在 SARS 训练集中",
            "potency_transfer_supported": verdict_transfer,
        },
    }
    (RES_A / "a2_transferability.json").write_text(json.dumps(out, indent=2, ensure_ascii=False))
    _plot(paired, sp, sp_lo, sp_hi, n)
    return out


def _plot(paired: pd.DataFrame, sp: float, lo: float, hi: float, n: int) -> None:
    setup_cjk_fonts()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    x = paired["pActivity_sars"]
    y = paired["pActivity_mers"]

    ax = axes[0]
    ax.scatter(x, y, s=70, c="#2b6cb0", edgecolor="white", linewidth=1.2, zorder=3)
    lim = [min(x.min(), y.min()) - 0.3, max(x.max(), y.max()) + 0.3]
    ax.plot(lim, lim, "--", c="#a0aec0", lw=1.4, label="完全迁移 (y = x)", zorder=2)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("SARS-CoV / SARS-CoV-2 Mpro  pIC50")
    ax.set_ylabel("MERS-CoV Mpro  pIC50")
    ax.set_title(f"同一分子的配对效力 (n={n})\nSpearman = {sp:+.2f}  95%CI [{lo:+.2f}, {hi:+.2f}]")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=0.25, zorder=1)

    ax = axes[1]
    dp = paired["delta_p"]
    ax.hist(dp, bins=9, color="#dd6b20", edgecolor="white")
    ax.axvline(0, color="#2d3748", ls="--", lw=1.4)
    ax.axvline(dp.mean(), color="#c53030", lw=2, label=f"均值 {dp.mean():+.2f}")
    ax.set_xlabel("Δ pIC50  (MERS - SARS)")
    ax.set_ylabel("化合物数")
    ax.set_title("效力在两种酶之间的偏移")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=0.25, axis="y")

    fig.suptitle("A2：跨冠状病毒迁移前提的直接检验", fontsize=13, y=1.00)
    fig.tight_layout()
    out = FIGS / "a2_transferability.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log(f"图已保存 -> {out}")


if __name__ == "__main__":
    main()
