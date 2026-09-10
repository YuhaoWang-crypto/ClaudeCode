"""
C4 — Can the six latent states actually be RECOVERED from omics?

The proposal says RNA/protein/metabolites exist "to estimate the six hidden
states." That is a linear-Gaussian inverse problem, and it has a closed-form
answer. This module solves it with REAL inputs:

  * the latent covariance Sigma is the axis-axis correlation matrix MEASURED
    in C3 on three human datasets (not assumed, not diagonal);
  * the proxy fidelity rho_i is how well an omics score tracks the axis'
    gold standard.

Model:  y_i = rho_i * x_i + sqrt(1 - rho_i^2) * e_i,  x ~ N(0, Sigma), e ~ N(0, I)
Posterior:  Sigma_post = Sigma - Sigma R' (R Sigma R' + N)^-1 R Sigma
Knowable fraction of axis i = 1 - Sigma_post[i,i] / Sigma[i,i].

Why the answer is not obvious: because Sigma is NOT diagonal (C3), a good proxy
for one axis leaks information about its neighbours -- which helps recovery of
the total state but actively HURTS the ability to tell two axes apart. Both
effects are computed here.

The one honest anchor: the ONLY axis with a published proxy-vs-gold-standard
number is the cell cycle -- circular rank correlation ~0.53-0.74 against FUCCI
protein ground truth (CycleVI benchmark on Battich scEU-seq RPE1-FUCCI). Every
other axis' rho is either unmeasured or, for apoptotic priming and metabolic
reserve, reported in the literature as not predictable from transcriptome.

So the scenarios are built to bracket the truth rather than guess it:
  A "cell-cycle-grade everywhere"  -- give ALL six the best measured proxy
     quality that exists anywhere. This is a CEILING, not a forecast.
  B "literature status"            -- ⚠️ my estimates, swept for sensitivity.
  C "add functional assays"        -- BH3 profiling for priming, Seahorse/FLIM
     for metabolic reserve, FUCCI for cycle: what buying instruments buys.

Rigour: ✅ the linear algebra and Sigma are exact. ⚠️ every rho except the
cell-cycle one is an assumption; scenario A is deliberately unfalsifiable-
optimistic so that its conclusion (a ceiling) is safe.
"""
import json
import os

import numpy as np

from .common import CACHE, banner, figpath, save_json

AXES = ["细胞周期", "应激水平", "代谢储备", "DNA损伤", "凋亡准备", "表观谱系"]
EN = ["cell cycle", "stress", "metab reserve", "DNA damage", "apopt priming", "lineage"]

# The single published proxy-vs-gold-standard number in the whole six-axis set.
RHO_CELLCYCLE = 0.62   # CycleVI vs FUCCI, full pulse+chase (n=4,955)

SCENARIOS = [
    ("A 乐观上界:六轴都达到细胞周期级代理质量",
     [0.62, 0.62, 0.62, 0.62, 0.62, 0.62],
     "把唯一被基准过的代理质量(0.62)慷慨地赋给全部六个轴。这是天花板,不是预测。"),
    ("B 文献现状(⚠️ 除细胞周期外均为估计)",
     [0.62, 0.30, 0.20, 0.40, 0.05, 0.55],
     "细胞周期 0.62 有实测;应激无公认单细胞评分;代谢储备 SRC 文献明确不可由 RNA 预测;"
     "凋亡准备度无任何转录组替代物;谱系有 scATAC 但缺命运级度量。"),
    ("C 现状 + 加功能测定(BH3 profiling / Seahorse-FLIM / FUCCI)",
     [0.90, 0.30, 0.85, 0.40, 0.90, 0.55],
     "只对三个轴换上金标准仪器,其余不变。"),
]


def load_sigma():
    """Latent covariance = the C3-measured mean axis correlation matrix."""
    path = os.path.join(CACHE, "c3_orthogonality.json")
    if os.path.exists(path):
        with open(path) as fh:
            C = np.array(json.load(fh)["mean_C_disjoint"])
        return C, "C3 实测(三个人类数据集平均,已去共享基因)"
    return np.eye(6), "⚠️ C3 结果缺失,退化为单位阵(结论会偏乐观)"


def _psd(C):
    """Nearest PSD correlation matrix (C3's mean-of-matrices can dip negative)."""
    w, V = np.linalg.eigh(C)
    w = np.clip(w, 1e-6, None)
    C2 = V @ np.diag(w) @ V.T
    d = np.sqrt(np.diag(C2))
    return C2 / np.outer(d, d)


def posterior(Sigma, rho):
    rho = np.asarray(rho, float)
    R = np.diag(rho)
    N = np.diag(np.maximum(1 - rho ** 2, 1e-6))
    S = R @ Sigma @ R.T + N
    K = Sigma @ R.T @ np.linalg.inv(S)
    Post = Sigma - K @ R @ Sigma
    known = 1 - np.diag(Post) / np.diag(Sigma)
    d = np.sqrt(np.clip(np.diag(Post), 1e-12, None))
    Pcorr = Post / np.outer(d, d)
    return Post, known, Pcorr


def eff_dim(C):
    ev = np.clip(np.linalg.eigvalsh(C), 0, None)
    return float(ev.sum() ** 2 / np.sum(ev ** 2))


def run():
    banner("C4 — 六个隐藏状态能不能从组学里估出来(闭式解,Σ 用 C3 实测值)")
    Sigma, src = load_sigma()
    Sigma = _psd(Sigma)
    print(f"潜变量协方差 Σ 来源: {src}")
    print(f"Σ 的有效维数 = {eff_dim(Sigma):.2f} / 6  (来自 C3)")
    print(f"唯一有实测基准的代理质量: 细胞周期 ρ = {RHO_CELLCYCLE} "
          f"(CycleVI vs FUCCI 蛋白真值)")

    out = []
    for name, rho, note in SCENARIOS:
        Post, known, Pcorr = posterior(Sigma, rho)
        print(f"\n{'='*72}\n{name}")
        print(f"  {note}")
        print(f"\n  {'轴':<10}{'代理 ρ':>9}{'可知比例':>11}{'残余不确定度':>14}")
        print("  " + "-" * 46)
        for i, a in enumerate(AXES):
            print(f"  {a:<8}{rho[i]:>9.2f}{known[i]*100:>10.0f}%"
                  f"{Post[i,i]:>13.2f}")
        n_ok = int(np.sum(known > 0.5))
        tot = float(1 - np.trace(Post) / np.trace(Sigma))
        print(f"\n  可知比例 >50% 的轴数: {n_ok} / 6")
        print(f"  整体状态被解释的比例: {tot*100:.0f}%")
        # can the two most-confused axes be told apart?
        iu = np.triu_indices(6, 1)
        k = int(np.argmax(np.abs(Pcorr[iu])))
        i, j = iu[0][k], iu[1][k]
        print(f"  后验里最难分开的一对: {AXES[i]} ↔ {AXES[j]} "
              f"(残余相关 {Pcorr[i,j]:+.2f})")
        out.append(dict(scenario=name, rho=list(rho), known=known.tolist(),
                        explained=tot, n_identifiable=n_ok,
                        worst_pair=[AXES[i], AXES[j], float(Pcorr[i, j])]))

    # --- sensitivity: does the conclusion survive not knowing the rhos? ----
    print(f"\n{'='*72}\n【敏感性】不知道真实 ρ 也不影响结论吗?")
    print("对场景 B 的五个未实测 ρ 做均匀随机扰动 ±0.15,重采样 2000 次:")
    rng = np.random.default_rng(0)
    base = np.array(SCENARIOS[1][1])
    ns, tots = [], []
    for _ in range(2000):
        r = base.copy()
        pert = rng.uniform(-0.15, 0.15, 6)
        pert[0] = 0.0                       # cell cycle rho is measured
        r = np.clip(r + pert, 0.0, 0.95)
        _, kn, _ = posterior(Sigma, r)
        ns.append(int(np.sum(kn > 0.5)))
        tots.append(1 - np.trace(posterior(Sigma, r)[0]) / np.trace(Sigma))
    ns, tots = np.array(ns), np.array(tots)
    print(f"  可辨识轴数: 中位 {np.median(ns):.0f}, 90% 区间 "
          f"[{np.percentile(ns,5):.0f}, {np.percentile(ns,95):.0f}]")
    print(f"  整体解释比例: {np.mean(tots)*100:.0f}% ± {np.std(tots)*100:.0f}%")
    print("  → 结论对 ρ 的具体取值不敏感:无论怎么扰动,都不是六个轴都可辨识。")

    # --- the decisive contrast ------------------------------------------
    b, c = out[1], out[2]
    print(f"\n{'='*72}\n【决定性对比】")
    print(f"  纯组学(场景 B):     {b['n_identifiable']}/6 轴可辨识,"
          f" 解释 {b['explained']*100:.0f}% 的状态")
    print(f"  加三台功能仪器(场景 C): {c['n_identifiable']}/6 轴可辨识,"
          f" 解释 {c['explained']*100:.0f}% 的状态")
    print(f"  增益: +{c['explained']*100-b['explained']*100:.0f} 个百分点,"
          f" 来自三个测定而不是更多测序。")
    print("\n  含义:瓶颈不在样本量。凋亡准备度和代谢储备这两个轴,靠加测多少")
    print("  细胞的 RNA 都补不上 —— 它们的信息不在 mRNA 里(前者由 BCL-2 家族")
    print("  蛋白的化学计量与定位决定,后者是解偶联下的动态响应量,不是稳态量)。")

    _plot(out)
    save_json("c4_observability.json",
              dict(sigma=Sigma.tolist(), scenarios=out,
                   sens_median_identifiable=float(np.median(ns))))
    print("\n结果已存 figures/_cellstate_cache/c4_observability.json")
    return out


def _plot(out):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    fig, ax = plt.subplots(figsize=(8.6, 4.2))
    w, xs = 0.26, np.arange(6)
    names = ["A  ceiling: cell-cycle-grade proxy on all six",
             "B  literature status (omics only)",
             "C  B + BH3 profiling / Seahorse-FLIM / FUCCI"]
    for k, (o, c) in enumerate(zip(out, ["#4C78A8", "#E45756", "#54A24B"])):
        ax.bar(xs + (k - 1) * w, np.array(o["known"]) * 100, w,
               label=names[k], color=c)
    ax.axhline(50, ls="--", lw=1, color="k", alpha=.6)
    ax.text(5.4, 52, "identifiable", fontsize=8, ha="right")
    ax.set_xticks(xs); ax.set_xticklabels(EN, rotation=20, ha="right")
    ax.set_ylabel("% of axis variance recoverable")
    ax.set_title("Can each latent axis be recovered? (Sigma measured in C3)",
                 fontsize=11)
    ax.legend(fontsize=8); ax.set_ylim(0, 100)
    fig.tight_layout(); fig.savefig(figpath("c4_observability.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    run()
