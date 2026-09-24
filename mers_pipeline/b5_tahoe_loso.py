"""B5 —— Tahoe-100M 上的 LOSO 跨上下文迁移基准。

与 B2（sci-Plex）完全同一口径，便于直接对比：
  * 化合物层聚合（Tahoe：对 3 个浓度的 log2FC 取均值）
  * 分半可重复性作为噪声天花板（Tahoe：0.5 µM 与 5 µM 两个浓度互为半份）
  * control-similarity 加权共识 + basal gating
  * 同一套指标与四条有效性判据
  * 两种子集口径：留出侧筛选（回答"有可测响应时是否迁移"）
                  与源侧筛选（部署时真正可用）

Tahoe 相对 sci-Plex 的关键差异
------------------------------
  * **50 个细胞系**（sci-Plex 只有 3 个）→ 每折有 49 个源上下文，
    共识平均的收益不再被"只有 2 个源"限制。这是本轮扩展的核心价值。
  * delta 由 DESeq2 在完整 pseudobulk 上给出（相对同板 DMSO 对照），
    而非本流程自己从计数算出。
  * 基础表达用 baseMean 代理（见 b4 的说明）。

产出：results/trackB/b5_tahoe_loso.json, tahoe_loso_per_compound.csv
      figures/mers/b5_tahoe_loso.png
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy import stats

from .config import FIGS, RES_B, log, setup_cjk_fonts

REPRO_THRESHOLD = 0.30
SOFTMAX_TEMP = 10.0
MIN_NONNAN_FRAC = 0.98      # 基因需在这么高比例的条件里非 NaN 才纳入
MIN_GENES_PER_COND = 5000   # 条件级 QC：DESeq2 可检验基因过少的条件视为失败
BASAL_GATE_PCT = 20


def load(name: str = "tahoe_delta"):
    z = np.load(RES_B / f"{name}.npz", allow_pickle=True)
    # 50 药物版写的是 tahoe_meta.csv，全量版写的是 <name>_meta.csv
    for cand in (RES_B / f"{name}_meta.csv", RES_B / "tahoe_meta.csv"):
        if cand.exists():
            meta = pd.read_csv(cand)
            break
    else:
        raise FileNotFoundError(f"找不到 {name} 对应的条件元数据")
    return (z["delta"], z["genes"].astype(str), z["basal"],
            z["basal_cells"].astype(str), meta)


def qc_conditions(delta: np.ndarray, meta: pd.DataFrame) -> np.ndarray:
    """条件级 QC：DESeq2 能检验的基因数过少的条件是失败的比较，不是"无响应"。

    这些条件填 0 后会伪装成"没有变化"，从而稀释迁移指标 —— 必须先剔除。
    """
    n_ok = (~np.isnan(delta)).sum(axis=1)
    keep = n_ok >= MIN_GENES_PER_COND
    dropped = meta.loc[~keep, ["cell_line", "drug", "conc"]].assign(n_genes=n_ok[~keep])
    log(f"条件级 QC：{len(meta)} -> {int(keep.sum())}"
        f"（剔除可检验基因 < {MIN_GENES_PER_COND} 的 {int((~keep).sum())} 个条件）")
    for _, r in dropped.iterrows():
        log(f"    剔除 {r['cell_line']} / {r['drug']} @ {r['conc']} µM（{r['n_genes']} 个基因）")
    return keep


def select_genes(delta: np.ndarray, genes: np.ndarray) -> np.ndarray:
    ok = (~np.isnan(delta)).mean(axis=0) >= MIN_NONNAN_FRAC
    idx = np.where(ok)[0]
    log(f"基因筛选：{len(genes)} -> {len(idx)}（要求在 ≥{MIN_NONNAN_FRAC:.0%} 的条件里非 NaN）")
    return idx


def build(delta: np.ndarray, meta: pd.DataFrame, gidx: np.ndarray):
    """化合物层聚合 + 分半可重复性（两个高浓度互为半份）。"""
    d = delta[:, gidx]
    d = np.nan_to_num(d, nan=0.0)   # 余下的少量 NaN 视作"无变化"
    meta = meta.reset_index(drop=True)

    agg: dict[str, dict[str, np.ndarray]] = {}
    repro: dict[tuple[str, str], float] = {}
    for (cl, drug), g in meta.groupby(["cell_line", "drug"]):
        rows = g.index.to_numpy()
        agg.setdefault(cl, {})[drug] = d[rows].mean(axis=0)
        # 分半：0.5 µM vs 5 µM
        hi = g.index[np.isclose(g["conc"], 5.0)].to_numpy()
        mid = g.index[np.isclose(g["conc"], 0.5)].to_numpy()
        if len(hi) and len(mid):
            a, b = d[hi].mean(axis=0), d[mid].mean(axis=0)
            if a.std() > 0 and b.std() > 0:
                repro[(cl, drug)] = float(stats.pearsonr(a, b).statistic)
    for cl in sorted(agg):
        rr = [v for (c, _), v in repro.items() if c == cl]
        if rr:
            log(f"  {cl}: {len(agg[cl])} 化合物 | 分半可重复性中位 r={np.median(rr):.3f} "
                f"| r>{REPRO_THRESHOLD} 的 {sum(1 for x in rr if x > REPRO_THRESHOLD)} 个")
    return agg, repro


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
    return {"pearson_delta": float(r),
            "mae": float(np.mean(np.abs(pred - true))),
            "top50_overlap": len(set(ti) & set(pi)) / top_k,
            "sign_match": float(np.mean(np.sign(pred[ti]) == np.sign(true[ti])))}


def loso(agg: dict, repro: dict, controls: dict, seed: int = 0) -> pd.DataFrame:
    """LOSO 迁移基准。

    除 no-change 与未加权共识外，额外加入两个**负对照**，用来分离
    "跨药物通用响应" 与 "药物特异性迁移"：

      generic_response —— 用源细胞系**全部药物的平均响应**去预测。
                          这个预测里没有任何药物特异性信息。
      shuffled_compound —— 用**另一个随机药物**的共识去预测目标药物。

    没有这两条基线，裸报 delta Pearson 会把通用响应算进迁移能力。
    """
    rng = np.random.default_rng(seed)
    rows = []
    cls = sorted(agg)
    log(f"LOSO：{len(cls)} 个细胞系，每折 {len(cls)-1} 个源上下文")
    # 每个细胞系在全部药物上的平均响应 —— 通用响应成分
    mean_resp = {c: np.mean(np.vstack(list(agg[c].values())), axis=0) for c in cls}
    for held in cls:
        others = [c for c in cls if c != held and c in controls]
        if not others or held not in controls:
            continue
        ctrl_held = controls[held]
        gate = basal_gate(ctrl_held)
        sims = np.array([stats.pearsonr(ctrl_held, controls[o]).statistic for o in others])
        w = np.exp(sims * SOFTMAX_TEMP)
        w = w / w.sum()

        drugs_here = sorted(agg[held])
        generic = np.mean(np.vstack([mean_resp[o] for o in others]), axis=0)
        for cmpd, true in agg[held].items():
            srcs = [(o, agg[o][cmpd]) for o in others if cmpd in agg[o]]
            if not srcs:
                continue
            mats = np.vstack([v for _, v in srcs])
            ww = np.array([w[others.index(o)] for o, _ in srcs])
            ww = ww / ww.sum()

            # 负对照：随机换一个别的药物的共识
            other_drugs = [d for d in drugs_here if d != cmpd]
            shuf = np.zeros_like(true)
            if other_drugs:
                d2 = other_drugs[int(rng.integers(len(other_drugs)))]
                srcs2 = [agg[o][d2] for o in others if d2 in agg[o]]
                if srcs2:
                    shuf = np.mean(np.vstack(srcs2), axis=0)

            preds = {
                "weighted_consensus_gated": (mats * ww[:, None]).sum(0) * gate,
                "weighted_consensus": (mats * ww[:, None]).sum(0),
                "unweighted_consensus": mats.mean(0),
                "generic_response": generic,        # 负对照：无药物特异性
                "shuffled_compound": shuf,          # 负对照：换一个药物
                "no_change": np.zeros_like(true),
            }
            rc = repro.get((held, cmpd))
            src_rep = [repro.get((o, cmpd)) for o, _ in srcs]
            src_rep = [v for v in src_rep if v is not None]
            src_ok = bool(src_rep and np.mean(src_rep) > REPRO_THRESHOLD)
            for name, p in preds.items():
                m = metrics(p, true)
                rows.append({
                    "held_out": held, "compound": cmpd, "method": name,
                    "n_sources": len(srcs), "repro_r": rc,
                    "source_repro_r": float(np.mean(src_rep)) if src_rep else None,
                    "reproducible": bool(rc is not None and rc > REPRO_THRESHOLD),
                    "source_reproducible": src_ok,
                    "transfer_efficiency": (m["pearson_delta"] / rc
                                            if rc and rc > 0.05 else None),
                    **m})
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
            "transfer_efficiency_median": (round(float(np.median(te)), 4) if len(te) else None),
        }
    for m, v in out.items():
        log(f"  [{label}·{m}] n={v['n']} r均值={v['pearson_delta_mean']} "
            f"中位={v['pearson_delta_median']} MAE={v['mae_mean']} "
            f"top50={v['top50_overlap_mean']} sign={v['sign_match_mean']} "
            f"迁移效率中位={v['transfer_efficiency_median']}")
    return out


def criteria(summ: dict, primary: str = "weighted_consensus_gated") -> dict:
    """有效性判据。

    相对第一版的三处修正：
      * c2 原写成"落在 [0.2, 0.45] 区间内"，于是**超出上界也算失败** —— 这是错的，
        超过遗传扰动 PoC 的水平是更好而不是更差。改为"≥ 0.2"，并单独标注是否超出上界。
      * c3 原拿 weighted_gated 与 unweighted 比，把"加权"和"gating"两件事混在一起。
        改为在**同一 gating 条件下**比较加权与未加权，另设 c3b 单独评估 gating。
      * 新增 c5：必须超过"通用响应"负对照，否则测到的只是跨药物共有成分，
        不是药物特异性迁移。
    """
    p = summ.get(primary, {})
    wg = summ.get("weighted_consensus_gated", {})
    w = summ.get("weighted_consensus", {})
    u = summ.get("unweighted_consensus", {})
    nc = summ.get("no_change", {})
    gen = summ.get("generic_response", {})
    shuf = summ.get("shuffled_compound", {})
    v = lambda d: d.get("pearson_delta_mean")

    # 注意：不能写成 `(pv or 1) < 0.05` —— p 值下溢到 0.0 时 `0.0 or 1` 会得到 1，
    # 把最显著的结果判成不显著。必须显式区分 None 与 0.0。
    pv = p.get("wilcoxon_greater_than_0_p")
    pv_ok = pv is not None and pv < 0.05

    return {
        "c1_pearson_gt_0": {
            "value": v(p), "p": pv,
            "passed": bool((v(p) or 0) > 0 and pv_ok)},
        "c2_magnitude_at_least_genetic_poc": {
            "value": v(p), "threshold": 0.2,
            "exceeds_genetic_poc_upper_0.45": bool((v(p) or 0) > 0.45),
            "passed": bool((v(p) or -1) >= 0.2)},
        "c3_weighting_helps_same_gating": {
            "weighted": v(w), "unweighted": v(u),
            "delta": round((v(w) or 0) - (v(u) or 0), 4),
            "passed": bool((v(w) or -1) >= (v(u) or 1e9))},
        "c3b_gating_helps": {
            "gated": v(wg), "ungated": v(w),
            "gated_mae": wg.get("mae_mean"), "ungated_mae": w.get("mae_mean"),
            "passed": bool((v(wg) or -1) >= (v(w) or 1e9))},
        "c4_beats_no_change_mae": {
            "primary_mae": p.get("mae_mean"), "no_change_mae": nc.get("mae_mean"),
            "passed": bool((p.get("mae_mean") or 1e9) < (nc.get("mae_mean") or 0))},
        "c5_beats_generic_response": {
            "primary": v(p), "generic_response": v(gen), "shuffled_compound": v(shuf),
            "drug_specific_gain_vs_generic": round((v(p) or 0) - (v(gen) or 0), 4),
            "drug_specific_gain_vs_shuffled": round((v(p) or 0) - (v(shuf) or 0), 4),
            "passed": bool((v(p) or -1) > (v(gen) or 1e9))},
    }


def main(name: str = "tahoe_delta", tag: str = "") -> dict:
    delta, genes, basal, basal_cells, meta = load(name)
    log(f"Tahoe delta: {delta.shape[0]} 条件 × {delta.shape[1]} 基因；"
        f"{meta['cell_line'].nunique()} 细胞系 × {meta['drug'].nunique()} 药物")
    keep = qc_conditions(delta, meta)
    n_dropped = int((~keep).sum())
    delta, meta = delta[keep], meta.loc[keep].reset_index(drop=True)
    gidx = select_genes(delta, genes)

    controls = {}
    for i, c in enumerate(basal_cells):
        v = np.nan_to_num(basal[i][gidx], nan=0.0)
        controls[c] = np.log1p(np.clip(v, 0, None))
    log(f"基础表达代理：{len(controls)} 个细胞系")

    agg, repro = build(delta, meta, gidx)
    df = loso(agg, repro, controls)
    df.to_csv(RES_B / f"tahoe_loso_per_compound{tag}.csv", index=False)

    rep_vals = pd.Series([v for v in repro.values()])
    ceiling = {
        "median_split_half_r_all": round(float(rep_vals.median()), 4),
        "median_split_half_r_reproducible": (
            round(float(rep_vals[rep_vals > REPRO_THRESHOLD].median()), 4)
            if (rep_vals > REPRO_THRESHOLD).any() else None),
        "n_compound_cellline_pairs": int(len(rep_vals)),
        "n_reproducible_pairs": int((rep_vals > REPRO_THRESHOLD).sum()),
        "threshold": REPRO_THRESHOLD,
        "split_half_definition": "同细胞系同化合物的 0.5 µM 与 5 µM 两个浓度互为半份",
        "caveat": ("这不是严格的重复测量分半：两个浓度的响应本身就可能不同，"
                   "因此它**低估**了真实可重复性，使 transfer_efficiency 可能 > 1。"
                   "sci-Plex 侧用的是交错剂量 {1,3} vs {2,4}，更接近重复分半。"
                   "两个数据集的 transfer_efficiency 不可直接互比。"),
    }
    log(f"可重复性天花板：全部中位 r={ceiling['median_split_half_r_all']}，"
        f"可重复子集中位 r={ceiling['median_split_half_r_reproducible']}"
        f"（{ceiling['n_reproducible_pairs']}/{ceiling['n_compound_cellline_pairs']} 对）")

    s_all = summarize(df, "全部化合物")
    s_held = summarize(df[df["reproducible"]], "可重复子集(留出侧)") if df["reproducible"].any() else {}
    s_src = summarize(df[df["source_reproducible"]], "源侧可重复子集(可前瞻)") if df["source_reproducible"].any() else {}

    crits = {"all": criteria(s_all),
             "reproducible_heldout_selected": criteria(s_held) if s_held else {},
             "reproducible_source_selected": criteria(s_src) if s_src else {}}
    for k, c in crits.items():
        if c:
            log(f"有效性判据（{k}）: " + ", ".join(
                f"{n.split('_')[0]}={'通过' if v['passed'] else '未通过'}" for n, v in c.items()))

    # 与 sci-Plex 结果对比
    comp = None
    p = RES_B / "b2_loso.json"
    if p.exists():
        b2 = json.loads(p.read_text())
        g = lambda d, k="pearson_delta_mean": (d or {}).get("weighted_consensus_gated", {}).get(k)
        comp = {
            "sciplex": {
                "n_cell_lines": 3,
                "n_source_contexts_per_fold": 2,
                "all_compounds_r": g(b2["all_compounds"]["summary"]),
                "source_selected_r": g(b2["reproducible_subset_source_selected"]["summary"]),
                "ceiling_median_reproducible": b2["reproducibility_ceiling"]["median_split_half_r_reproducible"],
                "reproducible_fraction": round(
                    b2["reproducibility_ceiling"]["n_reproducible_pairs"]
                    / b2["reproducibility_ceiling"]["n_compound_cellline_pairs"], 3),
            },
            "tahoe": {
                "n_cell_lines": int(len(agg)),
                "n_source_contexts_per_fold": int(len(agg) - 1),
                "all_compounds_r": g(s_all),
                "source_selected_r": g(s_src),
                "ceiling_median_reproducible": ceiling["median_split_half_r_reproducible"],
                "reproducible_fraction": round(
                    ceiling["n_reproducible_pairs"] / ceiling["n_compound_cellline_pairs"], 3),
            },
        }

    out = {
        "dataset": "Tahoe-100M pseudobulk DE 子集（见 b4_tahoe.json）",
        "n_cell_lines": int(len(agg)),
        "n_compounds": int(meta["drug"].nunique()),
        "n_genes_evaluated": int(len(gidx)),
        "condition_qc": {"min_testable_genes": MIN_GENES_PER_COND,
                         "n_conditions_dropped": n_dropped,
                         "n_conditions_used": int(len(meta)),
                         "rationale": "DESeq2 可检验基因过少的条件是失败的比较，"
                                      "填 0 后会伪装成\"无响应\"并稀释迁移指标"},
        "reproducibility_ceiling": ceiling,
        "all_compounds": {"summary": s_all, "validity_criteria": crits["all"]},
        "reproducible_subset_heldout_selected": {
            "n_evaluations": int(df["reproducible"].sum()),
            "summary": s_held, "validity_criteria": crits["reproducible_heldout_selected"]},
        "reproducible_subset_source_selected": {
            "n_evaluations": int(df["source_reproducible"].sum()),
            "summary": s_src, "validity_criteria": crits["reproducible_source_selected"]},
        "comparison_with_sciplex": comp,
        "scope_caveats": [
            "按方案描述**重新实现**的 control-similarity 加权共识 + basal gating，"
            "非调用用户 VirtualCell PoC 原始代码（该压缩包不在本会话环境中）。",
            "用户 PoC 中权重 70% 的 STATE 组件是遗传扰动模型，不接受小分子输入，"
            "本基准只检验其模态无关的迁移/共识机制。",
            "基础表达用 baseMean 代理（含少量处理细胞，见 b4 说明）。",
            "本基准检验跨细胞系迁移化合物转录响应，**不是**预测生化 IC50。",
        ],
    }
    (RES_B / f"b5_tahoe_loso{tag}.json").write_text(json.dumps(out, indent=2, ensure_ascii=False))
    _plot(df, rep_vals, out, comp, tag)
    log(f"B5 完成 -> {RES_B/f'b5_tahoe_loso{tag}.json'}")
    return out


def _plot(df: pd.DataFrame, rep_vals: pd.Series, out: dict, comp: dict | None,
          tag: str = "") -> None:
    setup_cjk_fonts()
    import matplotlib.pyplot as plt

    order = ["weighted_consensus_gated", "weighted_consensus", "unweighted_consensus",
             "generic_response", "shuffled_compound", "no_change"]
    labels = {"weighted_consensus_gated": "加权共识\n+gating", "weighted_consensus": "加权共识",
              "unweighted_consensus": "未加权共识", "generic_response": "通用响应\n(负对照)",
              "shuffled_compound": "打乱药物\n(负对照)", "no_change": "no-change"}
    order = [o for o in order if o in set(df["method"])]
    colors = ["#2b6cb0", "#4299e1", "#dd6b20", "#d69e2e", "#e53e3e", "#a0aec0"]

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.9))

    ax = axes[0]
    ax.hist(rep_vals, bins=40, range=(-0.3, 1), color="#805ad5", alpha=0.8)
    ax.axvline(REPRO_THRESHOLD, color="#e53e3e", ls="--", lw=1.6,
               label=f"可重复门槛 {REPRO_THRESHOLD}")
    ax.set_xlabel("同细胞系内 分半可重复性 r（0.5 vs 5 µM）")
    ax.set_ylabel("(细胞系,化合物) 对数")
    ax.set_title(f"噪声天花板（{out['reproducibility_ceiling']['n_reproducible_pairs']}"
                 f"/{out['reproducibility_ceiling']['n_compound_cellline_pairs']} 对达标）")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=0.25, axis="y")

    ax = axes[1]
    x = np.arange(len(order)); wd = 0.27
    for k, (mask, nm, c) in enumerate([
            (pd.Series(True, index=df.index), "全部化合物", "#718096"),
            (df["reproducible"], "可重复(留出侧)", "#2b6cb0"),
            (df["source_reproducible"], "可重复(源侧)", "#38a169")]):
        vals = [df[mask & (df["method"] == m)]["pearson_delta"].mean() for m in order]
        ax.bar(x + (k - 1) * wd, vals, wd, label=nm, color=c)
    ceil_v = out["reproducibility_ceiling"]["median_split_half_r_reproducible"]
    if ceil_v:
        ax.axhline(ceil_v, color="#e53e3e", ls="--", lw=1.6, label=f"天花板 {ceil_v:.2f}")
    ax.set_xticks(x); ax.set_xticklabels([labels[m] for m in order], fontsize=8.5)
    ax.axhline(0, color="#2d3748", lw=1.1)
    ax.set_ylabel("delta Pearson r 均值")
    ax.set_title(f"Tahoe LOSO（{out['n_cell_lines']} 细胞系，每折 {out['n_cell_lines']-1} 个源）")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.25, axis="y")

    ax = axes[2]
    if comp:
        keys = ["all_compounds_r", "source_selected_r", "ceiling_median_reproducible"]
        names = ["全部化合物", "源侧可重复子集", "可重复性天花板"]
        x = np.arange(len(keys)); wd = 0.36
        for k, (ds, c) in enumerate([("sciplex", "#dd6b20"), ("tahoe", "#2b6cb0")]):
            vals = [comp[ds][kk] or 0 for kk in keys]
            lbl = (f"sci-Plex（3 细胞系）" if ds == "sciplex"
                   else f"Tahoe（{comp['tahoe']['n_cell_lines']} 细胞系）")
            b = ax.bar(x + (k - 0.5) * wd, vals, wd, label=lbl, color=c)
            for bb, v in zip(b, vals):
                ax.text(bb.get_x() + bb.get_width() / 2, v + 0.008, f"{v:.2f}",
                        ha="center", fontsize=8)
        ax.set_xticks(x); ax.set_xticklabels(names, fontsize=9)
        ax.set_ylabel("delta Pearson r")
        ax.set_title("两个数据集对比")
        ax.legend(frameon=False, fontsize=8.5)
        ax.grid(alpha=0.25, axis="y")

    fig.suptitle("B5：Tahoe-100M 上的跨上下文迁移基准（50 个细胞系）", fontsize=13)
    fig.tight_layout()
    o = FIGS / f"b5_tahoe_loso{tag}.png"
    fig.savefig(o, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log(f"图已保存 -> {o}")


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "full":
        main(name="tahoe_delta_full", tag="_full")
    else:
        main()
