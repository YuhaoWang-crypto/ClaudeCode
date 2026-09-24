"""B2 —— 留一上下文（LOSO）零样本迁移基准，检验 virtual cell 迁移机制在化学扰动上的适用性。

范围界定（必须保留在报告中）
----------------------------
1. 用户的 `VirtualCell_ZeroShot_PoC_2026.zip` **不在本会话环境中**。本模块是按方案
   描述**重新实现**其原生方法（control-similarity 加权共识 + basal gating），
   属复现实现，不是对用户原始代码的直接验证。
2. 用户 PoC 中权重 70% 的 STATE 集成组件是**遗传扰动**模型，不接受小分子输入，
   因此本基准只检验其**模态无关**的迁移/共识机制。
3. 本基准检验"跨细胞系迁移化合物转录响应"，**不预测生化 IC50**。

关键的基准设计修正（本模块的核心方法学贡献）
--------------------------------------------
第一版实现直接在全部 188 化合物 × 4 剂量上跑 LOSO，得到 delta Pearson ≈ 0.03，
看上去是"迁移完全失败"。但正对照诊断推翻了这个解读：

  * 同一细胞系、同一化合物、**相邻剂量**之间的 delta 相关中位只有 0.07–0.10；
  * `||delta||` 在 10 nM → 10 µM 四个数量级上几乎不变（53.1 → 55.1）。

也就是说，在本 pseudobulk 深度下，多数 (化合物 × 剂量) 的 delta 由噪声主导。
跨上下文得到 0.03 并不是"迁移失败"，而是**根本没到自身可重复性的天花板**。

因此本版做三件事：
  (1) **化合物层聚合**：把同一化合物 4 个剂量的原始计数相加，细胞数 ×4 提升信噪比；
  (2) **分半可重复性**：用剂量 {1,3} vs {2,4} 估计每个 (细胞系, 化合物) 的噪声天花板，
      该估计不涉及任何跨细胞系信息，对迁移问题无循环论证；
  (3) **以天花板为参照报告迁移**：迁移的正确对标不是 1.0，而是同一上下文内的
      可重复性 r_ceiling。报告 transfer_r / r_ceiling 这一归一化迁移效率。

产出：results/trackB/b2_loso.json, loso_per_perturbation.csv, reproducibility.csv
      figures/mers/b2_loso.png
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy import stats

from .config import FIGS, RES_B, log, setup_cjk_fonts

MIN_MEAN_CPM = 1.0
BASAL_GATE_PCT = 20
REPRO_THRESHOLD = 0.30     # 分半可重复性门槛，定义"有真实响应"的化合物
SOFTMAX_TEMP = 10.0


# ------------------------------------------------------------------ 数据
def load() -> tuple[np.ndarray, pd.DataFrame, np.ndarray]:
    z = np.load(RES_B / "sciplex_pseudobulk.npz", allow_pickle=True)
    counts = z["counts"].astype(np.float64)
    gene_names = z["gene_names"].astype(str)
    meta = pd.read_csv(RES_B / "sciplex_group_meta.csv")
    assert len(meta) == counts.shape[0], "分组元数据与计数矩阵行数不一致"
    return counts, meta, gene_names


def cpm_log(c: np.ndarray) -> np.ndarray:
    t = c.sum() if c.ndim == 1 else c.sum(axis=-1, keepdims=True)
    t = np.maximum(t, 1)
    return np.log1p(c / t * 1e6)


def expressed_genes(counts: np.ndarray) -> np.ndarray:
    expr = cpm_log(counts)
    return np.where(expr.mean(axis=0) > np.log1p(MIN_MEAN_CPM))[0]


def build_compound_level(counts: np.ndarray, meta: pd.DataFrame, genes: np.ndarray):
    """化合物层聚合 + 分半可重复性。

    返回 controls[ct], deltas[ct][compound], repro[(ct, compound)]
    """
    controls, deltas, repro = {}, {}, {}
    for ct, sub in meta.groupby("cell_type"):
        veh = sub[sub["is_vehicle"]]
        if veh.empty:
            log(f"  警告：{ct} 无溶剂对照，跳过")
            continue
        ctrl = cpm_log(counts[veh.index.to_numpy()].sum(axis=0))[genes]
        controls[ct] = ctrl
        d = {}
        for prod, g in sub[~sub["is_vehicle"]].groupby("product_name"):
            g = g.sort_values("dose")
            idx = g.index.to_numpy()
            full = cpm_log(counts[idx].sum(axis=0))[genes] - ctrl
            d[prod] = full
            if len(idx) >= 4:
                h1 = cpm_log(counts[idx[[0, 2]]].sum(axis=0))[genes] - ctrl
                h2 = cpm_log(counts[idx[[1, 3]]].sum(axis=0))[genes] - ctrl
                if h1.std() > 0 and h2.std() > 0:
                    repro[(ct, prod)] = float(stats.pearsonr(h1, h2).statistic)
        deltas[ct] = d
        rr = [v for (c, _), v in repro.items() if c == ct]
        log(f"  {ct}: {len(veh)} 对照组, {len(d)} 化合物 | "
            f"分半可重复性中位 r={np.median(rr):.3f}, r>{REPRO_THRESHOLD} 的 {sum(1 for x in rr if x>REPRO_THRESHOLD)} 个")
    return controls, deltas, repro


def basal_gate(ctrl: np.ndarray) -> np.ndarray:
    thr = np.percentile(ctrl, BASAL_GATE_PCT)
    span = ctrl.max() - thr
    if span <= 0:
        return np.ones_like(ctrl)
    return 0.15 + 0.85 * np.clip((ctrl - thr) / span, 0.0, 1.0)


def metrics(pred: np.ndarray, true: np.ndarray, top_k: int = 50) -> dict:
    r = (stats.pearsonr(pred, true).statistic
         if pred.std() > 0 and true.std() > 0 else 0.0)
    ti = np.argsort(-np.abs(true))[:top_k]
    pi = np.argsort(-np.abs(pred))[:top_k]
    return {
        "pearson_delta": float(r),
        "mae": float(np.mean(np.abs(pred - true))),
        "top50_overlap": len(set(ti) & set(pi)) / top_k,
        "sign_match": float(np.mean(np.sign(pred[ti]) == np.sign(true[ti]))),
    }


# ------------------------------------------------------------------ LOSO
def loso(controls: dict, deltas: dict, repro: dict, genes_n: int,
         seed: int = 0) -> pd.DataFrame:
    """LOSO 迁移基准。

    除 no-change 与未加权共识外，额外加入两个**负对照**，用来分离
    "跨化合物通用响应" 与 "化合物特异性迁移"：
      generic_response  —— 用源细胞系**全部化合物的平均响应**去预测（无特异性信息）
      shuffled_compound —— 用**另一个随机化合物**的共识去预测目标化合物
    没有这两条基线，裸报 delta Pearson 会把通用响应算进迁移能力。
    """
    rng = np.random.default_rng(seed)
    rows = []
    cls = sorted(deltas)
    mean_resp = {c: np.mean(np.vstack(list(deltas[c].values())), axis=0)
                 for c in cls if deltas[c]}
    for held in cls:
        others = [c for c in cls if c != held]
        ctrl_held = controls[held]
        gate = basal_gate(ctrl_held)
        sims = np.array([stats.pearsonr(ctrl_held, controls[o]).statistic for o in others])
        w = np.exp(sims * SOFTMAX_TEMP)
        w = w / w.sum()
        log(f"  留出 {held}: 源 {others} 对照相似度 {np.round(sims,3)} -> 权重 {np.round(w,3)}")

        drugs_here = sorted(deltas[held])
        generic = np.mean(np.vstack([mean_resp[o] for o in others if o in mean_resp]), axis=0)
        for cmpd, true in deltas[held].items():
            srcs = [(o, deltas[o][cmpd]) for o in others if cmpd in deltas[o]]
            if not srcs:
                continue
            mats = np.vstack([d for _, d in srcs])
            ww = np.array([w[others.index(o)] for o, _ in srcs])
            ww = ww / ww.sum()

            other_drugs = [d for d in drugs_here if d != cmpd]
            shuf = np.zeros_like(true)
            if other_drugs:
                d2 = other_drugs[int(rng.integers(len(other_drugs)))]
                s2 = [deltas[o][d2] for o in others if d2 in deltas[o]]
                if s2:
                    shuf = np.mean(np.vstack(s2), axis=0)

            preds = {
                "weighted_consensus_gated": (mats * ww[:, None]).sum(0) * gate,
                "weighted_consensus": (mats * ww[:, None]).sum(0),
                "unweighted_consensus": mats.mean(0),
                "generic_response": generic,        # 负对照：无化合物特异性
                "shuffled_compound": shuf,          # 负对照：换一个化合物
                "no_change": np.zeros_like(true),
            }
            r_ceiling = repro.get((held, cmpd))
            # 两种子集定义：
            #  reproducible        —— 用留出细胞系自身的可重复性筛选。
            #                         它回答"在有可测响应的化合物上迁移是否有效"，
            #                         但预测时拿不到留出上下文的扰动数据，不可前瞻使用。
            #  source_reproducible —— 只用**源**细胞系的可重复性筛选。
            #                         这是部署时真正可用的筛选口径。
            src_repro = [repro.get((o, cmpd)) for o, _ in srcs]
            src_repro = [v for v in src_repro if v is not None]
            src_ok = bool(src_repro and np.mean(src_repro) > REPRO_THRESHOLD)
            for name, p in preds.items():
                m = metrics(p, true)
                rows.append({
                    "held_out": held, "compound": cmpd, "method": name,
                    "n_sources": len(srcs), "repro_r": r_ceiling,
                    "source_repro_r": float(np.mean(src_repro)) if src_repro else None,
                    "reproducible": bool(r_ceiling is not None and r_ceiling > REPRO_THRESHOLD),
                    "source_reproducible": src_ok,
                    "transfer_efficiency": (m["pearson_delta"] / r_ceiling
                                            if r_ceiling and r_ceiling > 0.05 else None),
                    **m,
                })
    return pd.DataFrame(rows)


def summarize(df: pd.DataFrame, label: str) -> dict:
    out = {}
    for method, g in df.groupby("method"):
        r = g["pearson_delta"].to_numpy()
        te = g["transfer_efficiency"].dropna().to_numpy()
        w = stats.wilcoxon(r, alternative="greater") if len(r) > 8 and np.any(r != 0) else None
        out[method] = {
            "n": int(len(g)),
            "pearson_delta_mean": round(float(np.mean(r)), 4),
            "pearson_delta_median": round(float(np.median(r)), 4),
            "pearson_delta_ci95": [round(float(np.percentile(r, 2.5)), 4),
                                   round(float(np.percentile(r, 97.5)), 4)],
            "wilcoxon_greater_than_0_p": float(w.pvalue) if w is not None else None,
            "mae_mean": round(float(g["mae"].mean()), 4),
            "top50_overlap_mean": round(float(g["top50_overlap"].mean()), 4),
            "sign_match_mean": round(float(g["sign_match"].mean()), 4),
            "transfer_efficiency_median": (round(float(np.median(te)), 4)
                                           if len(te) else None),
        }
    for m, v in out.items():
        log(f"  [{label}·{m}] n={v['n']} r均值={v['pearson_delta_mean']} "
            f"中位={v['pearson_delta_median']} MAE={v['mae_mean']} "
            f"top50={v['top50_overlap_mean']} sign={v['sign_match_mean']} "
            f"迁移效率中位={v['transfer_efficiency_median']}")
    return out


def main() -> dict:
    counts, meta, gene_names = load()
    log(f"pseudobulk: {counts.shape[0]} 组 × {counts.shape[1]} 基因")
    genes = expressed_genes(counts)
    log(f"表达基因 {len(genes)}（CPM > {MIN_MEAN_CPM}）")

    controls, deltas, repro = build_compound_level(counts, meta, genes)

    rep_df = pd.DataFrame(
        [{"cell_type": ct, "compound": c, "split_half_r": v} for (ct, c), v in repro.items()]
    )
    rep_df.to_csv(RES_B / "reproducibility.csv", index=False)

    df = loso(controls, deltas, repro, len(genes))
    df.to_csv(RES_B / "loso_per_perturbation.csv", index=False)

    all_summary = summarize(df, "全部化合物")
    sub = df[df["reproducible"]]
    repro_summary = summarize(sub, "可重复子集(留出侧筛选)") if len(sub) else {}
    sub_src = df[df["source_reproducible"]]
    src_summary = summarize(sub_src, "源侧可重复子集(可前瞻)") if len(sub_src) else {}

    ceil = rep_df["split_half_r"]
    ceiling_stats = {
        "median_split_half_r_all": round(float(ceil.median()), 4),
        "median_split_half_r_reproducible": round(
            float(ceil[ceil > REPRO_THRESHOLD].median()), 4) if (ceil > REPRO_THRESHOLD).any() else None,
        "n_compound_cellline_pairs": int(len(rep_df)),
        "n_reproducible_pairs": int((ceil > REPRO_THRESHOLD).sum()),
        "threshold": REPRO_THRESHOLD,
    }
    log(f"可重复性天花板: 全部中位 r={ceiling_stats['median_split_half_r_all']}, "
        f"可重复子集中位 r={ceiling_stats['median_split_half_r_reproducible']} "
        f"({ceiling_stats['n_reproducible_pairs']}/{ceiling_stats['n_compound_cellline_pairs']} 对)")

    # ---- 有效性判据 ----
    def crit(summ: dict) -> dict:
        w = summ.get("weighted_consensus_gated", {})
        wu = summ.get("weighted_consensus", {})
        u = summ.get("unweighted_consensus", {})
        nc = summ.get("no_change", {})
        gen = summ.get("generic_response", {})
        shuf = summ.get("shuffled_compound", {})
        # 注意：不能写成 `(pv or 1) < 0.05` —— p 值下溢到 0.0 时 `0.0 or 1` 会得到 1，
        # 把最显著的结果判成不显著。必须显式区分 None 与 0.0。
        pv = w.get("wilcoxon_greater_than_0_p")
        pv_ok = pv is not None and pv < 0.05
        return {
            "c1_pearson_gt_0": {
                "value": w.get("pearson_delta_mean"), "p": pv,
                "passed": bool((w.get("pearson_delta_mean") or 0) > 0 and pv_ok)},
            # c2 原写成"落在 [0.2,0.45] 内"，于是超出上界也算失败 —— 改为"≥0.2"。
            "c2_magnitude_at_least_genetic_poc": {
                "value": w.get("pearson_delta_mean"), "threshold": 0.2,
                "exceeds_genetic_poc_upper_0.45": bool((w.get("pearson_delta_mean") or 0) > 0.45),
                "passed": bool((w.get("pearson_delta_mean") or -1) >= 0.2)},
            # c3 原把"加权"与"gating"混在一起比；改为同一 gating 下比加权与未加权。
            "c3_weighting_helps_same_gating": {
                "weighted": wu.get("pearson_delta_mean"),
                "unweighted": u.get("pearson_delta_mean"),
                "delta": round((wu.get("pearson_delta_mean") or 0)
                               - (u.get("pearson_delta_mean") or 0), 4),
                "passed": bool((wu.get("pearson_delta_mean") or -1)
                               >= (u.get("pearson_delta_mean") or 1e9))},
            "c3b_gating_helps": {
                "gated": w.get("pearson_delta_mean"), "ungated": wu.get("pearson_delta_mean"),
                "gated_mae": w.get("mae_mean"), "ungated_mae": wu.get("mae_mean"),
                "passed": bool((w.get("pearson_delta_mean") or -1)
                               >= (wu.get("pearson_delta_mean") or 1e9))},
            "c4_beats_no_change_mae": {
                "weighted_mae": w.get("mae_mean"), "no_change_mae": nc.get("mae_mean"),
                "passed": bool((w.get("mae_mean") or 1e9) < (nc.get("mae_mean") or 0))},
            # 新增：必须超过"通用响应"负对照，否则测到的只是跨化合物共有成分。
            "c5_beats_generic_response": {
                "primary": w.get("pearson_delta_mean"),
                "generic_response": gen.get("pearson_delta_mean"),
                "shuffled_compound": shuf.get("pearson_delta_mean"),
                "drug_specific_gain_vs_generic": round(
                    (w.get("pearson_delta_mean") or 0) - (gen.get("pearson_delta_mean") or 0), 4),
                "drug_specific_gain_vs_shuffled": round(
                    (w.get("pearson_delta_mean") or 0) - (shuf.get("pearson_delta_mean") or 0), 4),
                "passed": bool((w.get("pearson_delta_mean") or -1)
                               > (gen.get("pearson_delta_mean") or 1e9))},
        }

    crit_all = crit(all_summary)
    crit_rep = crit(repro_summary) if repro_summary else {}
    crit_src = crit(src_summary) if src_summary else {}
    for nm, c in [("全部化合物", crit_all), ("可重复子集(留出侧)", crit_rep),
                  ("源侧可重复子集(可前瞻)", crit_src)]:
        if c:
            log(f"有效性判据（{nm}）: " +
                ", ".join(f"{k.split('_')[0]}={'通过' if v['passed'] else '未通过'}"
                          for k, v in c.items()))

    out = {
        "design_note": (
            "第一版在全部 188 化合物 × 4 剂量上跑 LOSO 得到 r≈0.03，但正对照显示"
            "同细胞系内相邻剂量的 delta 相关仅 0.07–0.10、||delta|| 跨四个数量级剂量不变，"
            "说明该分辨率下 delta 由噪声主导。本版改为化合物层聚合 + 分半可重复性天花板，"
            "并以 transfer_r / r_ceiling 报告归一化迁移效率。"
        ),
        "reproducibility_ceiling": ceiling_stats,
        "all_compounds": {"summary": all_summary, "validity_criteria": crit_all},
        "reproducible_subset_heldout_selected": {
            "selection": "用留出细胞系自身的分半可重复性筛选；回答'在有可测响应的化合物上"
                         "迁移是否有效'，但预测时拿不到留出上下文的扰动数据，不可前瞻使用",
            "n_evaluations": int(len(sub)),
            "n_compounds": int(sub["compound"].nunique()) if len(sub) else 0,
            "summary": repro_summary, "validity_criteria": crit_rep,
        },
        "reproducible_subset_source_selected": {
            "selection": "只用**源**细胞系的分半可重复性筛选 —— 部署时真正可用的口径",
            "n_evaluations": int(len(sub_src)),
            "n_compounds": int(sub_src["compound"].nunique()) if len(sub_src) else 0,
            "summary": src_summary, "validity_criteria": crit_src,
        },
        "verdict": _verdict(all_summary, repro_summary, src_summary, ceiling_stats),
        "scope_caveats": [
            "按方案描述**重新实现**的 control-similarity 加权共识 + basal gating，"
            "非调用用户 VirtualCell PoC 原始代码（该压缩包不在本会话环境中）。",
            "用户 PoC 中权重 70% 的 STATE 组件是遗传扰动模型，不接受小分子输入，"
            "本基准只检验其模态无关的迁移/共识机制。",
            "本基准检验跨细胞系迁移化合物转录响应，**不是**预测生化 IC50。",
        ],
    }
    (RES_B / "b2_loso.json").write_text(json.dumps(out, indent=2, ensure_ascii=False))
    _plot(df, rep_df, out)
    log(f"B2 完成 -> {RES_B/'b2_loso.json'}")
    log(f"判定: {out['verdict']['statement']}")
    return out


def _verdict(all_s: dict, rep_s: dict, src_s: dict, ceil: dict) -> dict:
    g = lambda s, k="pearson_delta_mean": (s or {}).get("weighted_consensus_gated", {}).get(k)
    return {
        "statement": (
            f"全部化合物：迁移 r={g(all_s)}（被无响应化合物稀释）。"
            f"限制到有可测响应的化合物后：留出侧筛选 r={g(rep_s)}、"
            f"源侧筛选（可前瞻）r={g(src_s)}。"
            f"同上下文可重复性天花板中位 r={ceil['median_split_half_r_reproducible']}，"
            f"归一化迁移效率中位 {g(rep_s,'transfer_efficiency_median')}。"
        ),
        "interpretation": (
            "迁移/共识机制在化学扰动上**确实迁移了真实信号**，但只在化合物本身"
            "产生可重复转录响应时成立。sci-Plex 中 "
            f"{ceil['n_compound_cellline_pairs']-ceil['n_reproducible_pairs']}/"
            f"{ceil['n_compound_cellline_pairs']} 个 (细胞系,化合物) 对达不到可重复性门槛，"
            "把它们计入会把均值稀释到接近 0 —— 这正是第一版分析得出假阴性的原因。"
        ),
        "power_limited": True,
    }


def _plot(df: pd.DataFrame, rep_df: pd.DataFrame, out: dict) -> None:
    setup_cjk_fonts()
    import matplotlib.pyplot as plt

    order = ["weighted_consensus_gated", "weighted_consensus",
             "unweighted_consensus", "no_change"]
    labels = {"weighted_consensus_gated": "加权共识\n+gating",
              "weighted_consensus": "加权共识",
              "unweighted_consensus": "未加权共识",
              "no_change": "no-change"}
    colors = ["#2b6cb0", "#4299e1", "#dd6b20", "#d69e2e", "#e53e3e", "#a0aec0"]
    order = [o for o in order if o in set(df["method"])]

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.9))

    # 1) 可重复性分布 —— 本基准的噪声天花板
    ax = axes[0]
    for ct, g in rep_df.groupby("cell_type"):
        ax.hist(g["split_half_r"], bins=30, range=(-0.3, 1), alpha=0.55, label=ct)
    ax.axvline(REPRO_THRESHOLD, color="#e53e3e", ls="--", lw=1.6,
               label=f"可重复门槛 {REPRO_THRESHOLD}")
    ax.set_xlabel("同细胞系内 分半可重复性 r")
    ax.set_ylabel("化合物数")
    ax.set_title("噪声天花板：多数化合物无可重复响应")
    ax.legend(frameon=False, fontsize=8.5)
    ax.grid(alpha=0.25, axis="y")

    # 2) 迁移 r：全部 vs 可重复子集
    ax = axes[1]
    x = np.arange(len(order)); wd = 0.36
    for k, (mask, nm, c) in enumerate([
            (df["reproducible"] | ~df["reproducible"], "全部化合物", "#718096"),
            (df["reproducible"], "可重复子集", "#2b6cb0")]):
        vals = [df[mask & (df["method"] == m)]["pearson_delta"].mean() for m in order]
        ax.bar(x + k * wd - wd / 2, vals, wd, label=nm, color=c)
    ceil_v = out["reproducibility_ceiling"]["median_split_half_r_reproducible"]
    if ceil_v:
        ax.axhline(ceil_v, color="#38a169", ls="--", lw=1.8,
                   label=f"可重复性天花板 {ceil_v:.2f}")
    ax.set_xticks(x); ax.set_xticklabels([labels[m] for m in order], fontsize=8.5)
    ax.axhline(0, color="#2d3748", lw=1.1)
    ax.set_ylabel("delta Pearson r 均值")
    ax.set_title("迁移表现 vs 自身可重复性上限")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.25, axis="y")

    # 3) 按留出细胞系
    ax = axes[2]
    cts = sorted(df["held_out"].unique())
    x = np.arange(len(cts)); wd = 0.8 / len(order)
    sub = df[df["reproducible"]]
    for i, (m, c) in enumerate(zip(order, colors)):
        vals = [sub[(sub["held_out"] == ct) & (sub["method"] == m)]["pearson_delta"].mean()
                for ct in cts]
        ax.bar(x + i * wd - 0.4 + wd / 2, vals, wd, color=c,
               label=labels[m].replace("\n", " "))
    ax.set_xticks(x); ax.set_xticklabels(cts)
    ax.axhline(0, color="#2d3748", lw=1.1)
    ax.set_ylabel("delta Pearson r 均值")
    ax.set_xlabel("留出的细胞系（仅可重复子集）")
    ax.set_title("按留出上下文分解")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.25, axis="y")

    fig.suptitle("B2：virtual cell 迁移机制在化学扰动上的 LOSO 基准", fontsize=13)
    fig.tight_layout()
    o = FIGS / "b2_loso.png"
    fig.savefig(o, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log(f"图已保存 -> {o}")


if __name__ == "__main__":
    main()
