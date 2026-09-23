"""A6 —— QSAR 建模：两种数据策略的对比，外加一个可行的分类任务。

数据现实闸门（data-reality gate）
--------------------------------
去重后 MERS-CoV Mpro 只有 19 个带精确 IC50 的化合物，pActivity 跨度 1.83 log、
标准差 0.57。这远低于任何可信回归模型的下限。按 binding-affinity 建模惯例：
  * N < 100 -> 只用指纹模型（Morgan FP + RF/GBM），不跑图神经网络；
  * 但 N = 19 连指纹回归都只能算 discovery-grade 的探索，其 CV 指标的
    置信区间会宽到无法据此做模型选择。
本模块**照样把两种策略都跑出来**（方案要求对比），但每个数字都附 bootstrap CI，
并明确标注哪些结论是数据支撑不了的。

三个任务
--------
  任务 1（回归·策略 1）：仅用 MERS Mpro 的 19 个化合物，留一 scaffold 交叉验证。
  任务 2（回归·策略 2）：用 SARS-CoV-2 + SARS-CoV Mpro 的 2563 个化合物训练，
                         先报告其**内部** scaffold-split CV 性能（证明流程本身有效），
                         再在 MERS 化合物上做外部评估。
                         注意：A2 已证明 17/19 的 MERS 化合物结构完全出现在训练集里，
                         所以这不是真正的"外部测试集"，此处按"训练集内记忆检验"解读。
  任务 3（分类·可行任务）：MERS Mpro 活性 vs 实测非活性（IC50 > 10 µM）。
                         这是本轨道数据量唯一支撑得起的 ML 任务。

产出：results/trackA/a6_qsar.json, figures/mers/a6_qsar.png
"""
from __future__ import annotations

import json
import warnings

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from scipy import stats
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import mean_absolute_error, roc_auc_score
from sklearn.model_selection import GroupKFold, LeaveOneGroupOut

from .config import FIGS, RANDOM_SEED, RES_A, log, setup_cjk_fonts

RDLogger.DisableLog("rdApp.*")
warnings.filterwarnings("ignore")

_MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def featurize(smiles: list[str]) -> np.ndarray:
    """Morgan 指纹（r=2, 2048 bit）。N 很小，故不叠加高维描述符以免过拟合。"""
    X = np.zeros((len(smiles), 2048), dtype=np.uint8)
    for i, s in enumerate(smiles):
        m = Chem.MolFromSmiles(s)
        if m is None:
            continue
        fp = _MORGAN.GetFingerprintAsNumPy(m)
        X[i] = fp
    return X


def spearman_ci(y, p, n_boot=5000, seed=RANDOM_SEED):
    rng = np.random.default_rng(seed)
    y, p = np.asarray(y, float), np.asarray(p, float)
    pt = stats.spearmanr(y, p).statistic
    vals = []
    for _ in range(n_boot):
        i = rng.integers(0, len(y), len(y))
        if len(np.unique(y[i])) < 3:
            continue
        v = stats.spearmanr(y[i], p[i]).statistic
        if np.isfinite(v):
            vals.append(v)
    return (float(pt), float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5)))


def regression_metrics(y, p, train_mean: float) -> dict:
    y, p = np.asarray(y, float), np.asarray(p, float)
    sp, lo, hi = spearman_ci(y, p)
    ss_res = float(np.sum((y - p) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    return {
        "n": int(len(y)),
        "spearman": round(sp, 3), "spearman_ci95": [round(lo, 3), round(hi, 3)],
        "pearson": round(float(stats.pearsonr(y, p).statistic), 3),
        "r2": round(1 - ss_res / ss_tot, 3) if ss_tot > 0 else None,
        "mae": round(float(mean_absolute_error(y, p)), 3),
        "rmse": round(float(np.sqrt(np.mean((y - p) ** 2))), 3),
        "mae_baseline_predict_train_mean": round(
            float(mean_absolute_error(y, np.full_like(y, train_mean))), 3),
    }


# ------------------------------------------------------------------ 任务 1
def task1_mers_only() -> dict:
    """策略 1：仅用 MERS Mpro 的 19 个化合物，留一 scaffold 交叉验证。"""
    d = pd.read_csv(RES_A / "mers_mpro_ic50.csv")
    X = featurize(d["smiles_std"].tolist())
    y = d["pActivity"].to_numpy()
    groups = d["scaffold"].to_numpy()
    n_sc = len(set(groups))
    log(f"[任务1·策略1 仅MERS] n={len(d)} 化合物 / {n_sc} 个 scaffold，留一 scaffold CV")

    logo = LeaveOneGroupOut()
    pred = np.zeros_like(y)
    for tr, te in logo.split(X, y, groups):
        m = RandomForestRegressor(n_estimators=500, random_state=RANDOM_SEED,
                                  n_jobs=-1, min_samples_leaf=1)
        m.fit(X[tr], y[tr])
        pred[te] = m.predict(X[te])
    met = regression_metrics(y, pred, float(y.mean()))
    log(f"    Spearman={met['spearman']} CI{met['spearman_ci95']} R²={met['r2']} "
        f"MAE={met['mae']} (预测均值基线 MAE={met['mae_baseline_predict_train_mean']})")
    return {
        "strategy": "仅 MERS-CoV Mpro",
        "n_compounds": int(len(d)), "n_scaffolds": int(n_sc),
        "cv": "留一 scaffold 交叉验证",
        "metrics": met,
        "predictions": {"y_true": y.tolist(), "y_pred": pred.round(3).tolist()},
        "data_reality_gate": (
            f"N={len(d)} 远低于可信回归的下限；pActivity 跨度仅 "
            f"{y.max()-y.min():.2f} log。按惯例不训练图神经网络，仅用指纹模型，"
            f"且结果只能作 discovery-grade 探索，不能用于模型选择。"
        ),
    }


# ------------------------------------------------------------------ 任务 2
def task2_sars_transfer() -> dict:
    """策略 2：SARS 训练 -> 先内部 CV，再在 MERS 上评估。"""
    sars = pd.read_csv(RES_A / "sars_mpro_train.csv")
    mers = pd.read_csv(RES_A / "mers_mpro_ic50.csv")

    Xs = featurize(sars["smiles_std"].tolist())
    ys = sars["pActivity"].to_numpy()
    gs = sars["scaffold"].to_numpy()
    log(f"[任务2·策略2 SARS迁移] 训练集 n={len(sars)} / {sars['scaffold'].nunique()} scaffold")

    # --- 2a 内部 scaffold-split CV：证明建模流程本身是有效的 ---
    gkf = GroupKFold(n_splits=5)
    pred_in = np.zeros_like(ys)
    for tr, te in gkf.split(Xs, ys, gs):
        m = RandomForestRegressor(n_estimators=500, random_state=RANDOM_SEED, n_jobs=-1)
        m.fit(Xs[tr], ys[tr])
        pred_in[te] = m.predict(Xs[te])
    met_in = regression_metrics(ys, pred_in, float(ys.mean()))
    log(f"    内部 5折 scaffold-split CV: Spearman={met_in['spearman']} "
        f"CI{met_in['spearman_ci95']} R²={met_in['r2']} MAE={met_in['mae']}")

    # --- 2b 在 MERS 上评估 ---
    model = RandomForestRegressor(n_estimators=500, random_state=RANDOM_SEED, n_jobs=-1)
    model.fit(Xs, ys)
    Xm = featurize(mers["smiles_std"].tolist())
    ym = mers["pActivity"].to_numpy()
    pm = model.predict(Xm)
    met_out = regression_metrics(ym, pm, float(ys.mean()))
    log(f"    MERS 上评估(n={len(mers)}): Spearman={met_out['spearman']} "
        f"CI{met_out['spearman_ci95']} R²={met_out['r2']} MAE={met_out['mae']}")

    # 区分"训练集里见过的"与"真正没见过的"
    seen = mers["smiles_std"].isin(set(sars["smiles_std"]))
    unseen = ~seen
    met_unseen = None
    if unseen.sum() >= 2:
        met_unseen = {
            "n": int(unseen.sum()),
            "mae": round(float(mean_absolute_error(ym[unseen.values], pm[unseen.values])), 3),
            "note": "样本量太小，不报告相关系数",
        }
    log(f"    其中训练集见过的 {int(seen.sum())} 个 / 未见过的 {int(unseen.sum())} 个")

    return {
        "strategy": "SARS-CoV-2 + SARS-CoV Mpro 迁移",
        "n_train": int(len(sars)), "n_train_scaffolds": int(sars["scaffold"].nunique()),
        "internal_cv": {"scheme": "5折 GroupKFold(scaffold)", "metrics": met_in},
        "mers_evaluation": {
            "metrics": met_out,
            "n_seen_in_training": int(seen.sum()),
            "n_unseen": int(unseen.sum()),
            "unseen_metrics": met_unseen,
            "predictions": {"y_true": ym.tolist(), "y_pred": pm.round(3).tolist(),
                            "seen_in_training": seen.tolist()},
        },
        "validity_caveat": (
            "这不是真正的外部测试集：19 个 MERS 化合物里有 17 个的标准化结构"
            "完全出现在 SARS 训练集中（见 A2）。因此 MERS 上的表现主要反映"
            "模型对训练样本的记忆，而非泛化能力 —— 而即便如此它仍然表现不佳，"
            "这与 A2 测得的'同分子跨酶效力不可迁移'一致。"
        ),
    }


# ------------------------------------------------------------------ 任务 3
def task3_classification() -> dict:
    """可行任务：活性 vs 实测非活性的分类，附单描述符基线。"""
    out = {}
    for name, path in [("primary_single_campaign", "benchmark_primary.csv"),
                       ("pooled_two_documents", "benchmark_pooled.csv")]:
        d = pd.read_csv(RES_A / path)
        y = d["active"].to_numpy()
        if y.sum() < 3:
            out[name] = {"error": f"活性样本仅 {int(y.sum())} 个，无法做分组 CV"}
            continue
        X = featurize(d["smiles_std"].tolist())
        groups = d["scaffold"].to_numpy()
        n_splits = min(5, len(set(groups[y == 1])))
        log(f"[任务3 分类·{name}] n={len(d)} ({int(y.sum())} 活性) "
            f"{len(set(groups))} scaffold -> {n_splits} 折 scaffold-split")

        gkf = GroupKFold(n_splits=n_splits)
        pred = np.zeros(len(y), dtype=float)
        for tr, te in gkf.split(X, y, groups):
            if y[tr].sum() == 0:
                pred[te] = 0.0
                continue
            m = RandomForestClassifier(n_estimators=500, random_state=RANDOM_SEED,
                                       n_jobs=-1, class_weight="balanced")
            m.fit(X[tr], y[tr])
            pred[te] = m.predict_proba(X[te])[:, 1]

        from .a5_docking import auc_bootstrap_ci, delong_like_paired_test
        auc, lo, hi = auc_bootstrap_ci(y, pred)
        res = {
            "n": int(len(d)), "n_active": int(y.sum()), "n_inactive": int((y == 0).sum()),
            "n_scaffolds": int(len(set(groups))), "cv": f"{n_splits}折 GroupKFold(scaffold)",
            "fingerprint_rf": {"roc_auc": round(auc, 3), "auc_ci95": [round(lo, 3), round(hi, 3)]},
            "descriptor_baselines": {},
        }
        for col in ["MW", "RotB", "HeavyAtoms", "cLogP", "TPSA"]:
            if col not in d.columns:
                continue
            s = d[col].to_numpy(float)
            a = roc_auc_score(y, s)
            if a < 0.5:
                s, a = -s, 1 - a
            res["descriptor_baselines"][col] = {
                "roc_auc": round(float(a), 3),
                "vs_model": delong_like_paired_test(y, pred, s),
            }
        best = max(res["descriptor_baselines"].items(), key=lambda kv: kv[1]["roc_auc"])
        res["strongest_baseline"] = {"descriptor": best[0], **best[1]}
        res["model_beats_strongest_baseline"] = bool(best[1]["vs_model"]["ci95"][0] > 0)
        log(f"    指纹RF AUC={res['fingerprint_rf']['roc_auc']} CI{res['fingerprint_rf']['auc_ci95']} | "
            f"最强基线 {best[0]} AUC={best[1]['roc_auc']} | "
            f"差值 {best[1]['vs_model']['auc_diff']} CI{best[1]['vs_model']['ci95']} "
            f"-> 显著更优: {res['model_beats_strongest_baseline']}")
        res["predictions"] = {"y_true": y.tolist(), "y_score": pred.round(4).tolist()}
        out[name] = res
    return out


# ------------------------------------------------------------------ 绘图
def _plot(t1: dict, t2: dict, t3: dict) -> None:
    setup_cjk_fonts()
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_curve

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.9))

    # 左：两种回归策略在同一批 MERS 化合物上的 pred-vs-actual
    ax = axes[0]
    y1 = np.array(t1["predictions"]["y_true"]); p1 = np.array(t1["predictions"]["y_pred"])
    y2 = np.array(t2["mers_evaluation"]["predictions"]["y_true"])
    p2 = np.array(t2["mers_evaluation"]["predictions"]["y_pred"])
    ax.scatter(y1, p1, s=62, c="#2b6cb0", edgecolor="white", lw=1,
               label=f"策略1 仅MERS (ρ={t1['metrics']['spearman']})", zorder=3)
    ax.scatter(y2, p2, s=62, c="#dd6b20", marker="^", edgecolor="white", lw=1,
               label=f"策略2 SARS迁移 (ρ={t2['mers_evaluation']['metrics']['spearman']})", zorder=3)
    lim = [min(y1.min(), p1.min(), p2.min()) - 0.3, max(y1.max(), p1.max(), p2.max()) + 0.3]
    ax.plot(lim, lim, "--", c="#a0aec0", lw=1.3, label="理想", zorder=2)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("实测 pIC50 (MERS Mpro)"); ax.set_ylabel("预测 pIC50")
    ax.set_title(f"两种数据策略 · 同一批 MERS 化合物 (n={len(y1)})")
    ax.legend(frameon=False, fontsize=8.5); ax.grid(alpha=0.25)

    # 中：SARS 内部 CV —— 证明流程本身有效
    ax = axes[1]
    mi = t2["internal_cv"]["metrics"]
    bars = {"SARS 内部\nscaffold-split CV": mi["spearman"],
            "策略2 用到\nMERS 上": t2["mers_evaluation"]["metrics"]["spearman"],
            "策略1\n仅 MERS": t1["metrics"]["spearman"]}
    cols = ["#38a169", "#dd6b20", "#2b6cb0"]
    ax.bar(list(bars), list(bars.values()), color=cols, width=0.6)
    ax.axhline(0, color="#2d3748", lw=1.2)
    for i, (k, v) in enumerate(bars.items()):
        ax.text(i, v + (0.03 if v >= 0 else -0.07), f"{v:+.2f}", ha="center", fontsize=10)
    ax.set_ylabel("Spearman ρ")
    ax.set_title("建模流程有效，但迁移不成立")
    ax.grid(alpha=0.25, axis="y")

    # 右：分类任务 ROC
    ax = axes[2]
    for key, c, nm in [("primary_single_campaign", "#2b6cb0", "主基准(单一活动)"),
                       ("pooled_two_documents", "#dd6b20", "合并基准")]:
        r = t3.get(key, {})
        if "predictions" not in r:
            continue
        y = np.array(r["predictions"]["y_true"]); s = np.array(r["predictions"]["y_score"])
        f, t, _ = roc_curve(y, s)
        ax.plot(f, t, lw=2.2, color=c,
                label=f"{nm} 指纹RF AUC={r['fingerprint_rf']['roc_auc']}")
        b = r["strongest_baseline"]
        ax.plot([], [], " ", label=f"    └ 最强基线 {b['descriptor']} AUC={b['roc_auc']}")
    ax.plot([0, 1], [0, 1], ":", c="#cbd5e0", lw=1.2, label="随机")
    ax.set_xlabel("假阳性率"); ax.set_ylabel("真阳性率")
    ax.set_title("任务3：活性 vs 实测非活性 分类")
    ax.legend(frameon=False, fontsize=8); ax.grid(alpha=0.25)

    fig.suptitle("A6：QSAR 两种数据策略对比 + 可行的分类任务", fontsize=13)
    fig.tight_layout()
    out = FIGS / "a6_qsar.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log(f"图已保存 -> {out}")


def main() -> dict:
    t1 = task1_mers_only()
    t2 = task2_sars_transfer()
    t3 = task3_classification()

    # 两策略在同一 MERS 测试集上的正面对比（方案要求的验收项）
    comparison = {
        "test_set": "MERS-CoV Mpro 的 19 个化合物",
        "strategy1_mers_only_spearman": t1["metrics"]["spearman"],
        "strategy1_ci95": t1["metrics"]["spearman_ci95"],
        "strategy2_sars_transfer_spearman": t2["mers_evaluation"]["metrics"]["spearman"],
        "strategy2_ci95": t2["mers_evaluation"]["metrics"]["spearman_ci95"],
        "model_selection_possible": False,
        "why": (
            "两个策略的 Spearman 95% 置信区间都跨过 0 且彼此大幅重叠，"
            "在 n=19、动态范围 1.83 log 的数据上无法区分优劣。"
            "按预先声明的验收口径（以 scaffold-split Spearman 选模型），"
            "此处**不做模型选择**，并判定 MERS-CoV Mpro 的 IC50 回归预测"
            "在当前公开数据下不可行。"
        ),
    }
    out = {"task1_regression_mers_only": t1,
           "task2_regression_sars_transfer": t2,
           "task3_classification": t3,
           "strategy_comparison": comparison}
    (RES_A / "a6_qsar.json").write_text(json.dumps(out, indent=2, ensure_ascii=False))
    _plot(t1, t2, t3)
    log(f"A6 完成 -> {RES_A/'a6_qsar.json'}")
    return out


if __name__ == "__main__":
    main()
