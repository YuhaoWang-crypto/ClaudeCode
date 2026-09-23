"""A5 —— 分子对接筛选与富集评估。

流程：RDKit ETKDG 生成 3D 构象 -> MMFF 优化 -> obabel 转 PDBQT -> AutoDock Vina 对接。
每个配体的结果单独缓存，支持中断续跑。

富集评估的关键对照：
  对接分数天然与分子大小相关。因此每个富集指标都与"仅用单个理化描述符"的
  同口径指标并列报告 —— 如果对接 AUC 不显著高于分子量 AUC，那么对接
  没有提供超出分子大小的信息。这是本模块的核心验收逻辑。

产出：results/trackA/docking_scores.csv
      results/trackA/a5_enrichment.json
      figures/mers/a5_enrichment.png
"""
from __future__ import annotations

import json
import subprocess
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from sklearn.metrics import roc_auc_score

from .config import DATA, FIGS, RANDOM_SEED, RES_A, log, setup_cjk_fonts

RDLogger.DisableLog("rdApp.*")
warnings.filterwarnings("ignore")

LIG_DIR = DATA / "ligands"
DOCK_DIR = DATA / "docking"
RECEPTOR_DIR = DATA / "receptors"
for _d in (LIG_DIR, DOCK_DIR):
    _d.mkdir(parents=True, exist_ok=True)


# ------------------------------------------------------------------ 配体准备
def prepare_ligand(smiles: str, lig_id: str) -> Path | None:
    """SMILES -> 3D -> MMFF 优化 -> PDBQT。已存在则跳过。"""
    pdbqt = LIG_DIR / f"{lig_id}.pdbqt"
    if pdbqt.exists() and pdbqt.stat().st_size > 0:
        return pdbqt
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = RANDOM_SEED
    if AllChem.EmbedMolecule(mol, params) != 0:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) != 0:
            return None
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    except Exception:  # noqa: BLE001
        pass
    sdf = LIG_DIR / f"{lig_id}.sdf"
    w = Chem.SDWriter(str(sdf))
    w.write(mol)
    w.close()
    r = subprocess.run(
        ["obabel", str(sdf), "-O", str(pdbqt), "-p", "7.4", "--partialcharge", "gasteiger"],
        capture_output=True, text=True,
    )
    if r.returncode != 0 or not pdbqt.exists() or pdbqt.stat().st_size == 0:
        return None
    return pdbqt


# ------------------------------------------------------------------ 对接
def dock_one(receptor_pdbqt: Path, lig_pdbqt: Path, box: dict, tag: str,
             exhaustiveness: int = 16, cpu: int = 0) -> dict | None:
    """单个配体对接，结果缓存为 JSON。"""
    cache = DOCK_DIR / f"{tag}.json"
    if cache.exists():
        try:
            return json.loads(cache.read_text())
        except Exception:  # noqa: BLE001
            pass
    from vina import Vina

    try:
        v = Vina(sf_name="vina", cpu=cpu, seed=RANDOM_SEED, verbosity=0)
        v.set_receptor(str(receptor_pdbqt))
        v.set_ligand_from_file(str(lig_pdbqt))
        v.compute_vina_maps(center=box["center"], box_size=box["size"])
        v.dock(exhaustiveness=exhaustiveness, n_poses=5)
        e = v.energies(n_poses=5)
        out = {
            "best_affinity": round(float(e[0][0]), 3),
            "mean_top3": round(float(np.mean([x[0] for x in e[:3]])), 3),
            "n_poses": len(e),
        }
        poses = DOCK_DIR / f"{tag}.poses.pdbqt"
        v.write_poses(str(poses), n_poses=3, overwrite=True)
    except Exception as exc:  # noqa: BLE001
        out = {"error": f"{exc.__class__.__name__}: {exc}"}
    cache.write_text(json.dumps(out))
    return out


def _dock_job(args: tuple) -> dict:
    """多进程 worker：每个进程用 1 个 CPU 跑一个配体。

    Vina 的内部多线程收益递减，把并行放在进程层面（N 个进程各占 1 核）
    比单进程开 N 线程快得多。
    """
    rec_q, smiles, lig_id, box, tag, exh = args
    lp = prepare_ligand(smiles, lig_id)
    if lp is None:
        return {"molecule_chembl_id": lig_id, "best_affinity": None,
                "error": "ligand_prep_failed"}
    res = dock_one(Path(rec_q), lp, box, tag, exh, cpu=1)
    return {"molecule_chembl_id": lig_id, **(res or {})}


def screen(df: pd.DataFrame, receptor_tag: str, exhaustiveness: int = 16,
           n_workers: int | None = None) -> pd.DataFrame:
    import multiprocessing as mp

    rec_q = RECEPTOR_DIR / f"{receptor_tag}.pdbqt"
    box = json.loads((RECEPTOR_DIR / f"{receptor_tag}.box.json").read_text())
    d = df.reset_index(drop=True)
    n_workers = n_workers or max(1, mp.cpu_count() - 1)

    jobs = []
    for i, r in d.iterrows():
        lig_id = r["molecule_chembl_id"] if isinstance(r["molecule_chembl_id"], str) else f"cpd{i}"
        jobs.append((str(rec_q), r["smiles_std"], lig_id, box,
                     f"{receptor_tag}__{lig_id}", exhaustiveness))

    log(f"  {receptor_tag}: {len(jobs)} 个配体，{n_workers} 个进程并行（每进程 1 核）")
    rows = []
    with mp.Pool(n_workers) as pool:
        for k, res in enumerate(pool.imap_unordered(_dock_job, jobs), 1):
            rows.append(res)
            if k % 10 == 0 or k == len(jobs):
                log(f"    {receptor_tag}: {k}/{len(jobs)}")
    out = pd.DataFrame(rows)
    return d.merge(out, on="molecule_chembl_id", how="left")


# ------------------------------------------------------------------ 富集指标
MIN_EF_K = 3   # 前 k 个里少于这么多化合物时，EF 退化到没有意义


def enrichment_factor(y_true: np.ndarray, score: np.ndarray, frac: float) -> float:
    """score 越大越好。EF = (前 frac 中的活性率) / (总体活性率)。

    注意：小数据集上高百分位的 EF 会退化 —— n=79 时 EF1% 只看 1 个化合物，
    取值只能是 0 或 1/Ra。此处对 k < MIN_EF_K 的情形返回 NaN 而不是给出
    一个看起来像指标的数字。
    """
    n = len(y_true)
    k = int(round(n * frac))
    if k < MIN_EF_K:
        return float("nan")
    order = np.argsort(-score)
    rate_top = y_true[order][:k].sum() / k
    rate_all = y_true.sum() / n
    return float(rate_top / rate_all) if rate_all > 0 else float("nan")


def bedroc(y_true: np.ndarray, score: np.ndarray, alpha: float = 20.0) -> float:
    """BEDROC（Truchon & Bayly 2007, eq. 36），强调排序靠前的富集，取值 [0, 1]。

        BEDROC = RIE · [Ra·sinh(α/2)] / [cosh(α/2) − cosh(α/2 − α·Ra)]
                 + 1 / [1 − e^(α(1−Ra))]
    """
    n = len(y_true)
    order = np.argsort(-score)
    y = y_true[order]
    ranks = np.where(y == 1)[0] + 1
    N_a = len(ranks)
    if N_a == 0 or N_a == n:
        return float("nan")
    ra = N_a / n
    rie = (np.sum(np.exp(-alpha * ranks / n))
           / (ra * (1 - np.exp(-alpha)) / (np.exp(alpha / n) - 1)))
    fac = ra * np.sinh(alpha / 2) / (np.cosh(alpha / 2) - np.cosh(alpha / 2 - alpha * ra))
    return float(rie * fac + 1 / (1 - np.exp(alpha * (1 - ra))))


def auc_bootstrap_ci(y: np.ndarray, s: np.ndarray, n_boot: int = 5000,
                     seed: int = RANDOM_SEED) -> tuple[float, float, float]:
    rng = np.random.default_rng(seed)
    point = roc_auc_score(y, s)
    vals = []
    idx_pos, idx_neg = np.where(y == 1)[0], np.where(y == 0)[0]
    for _ in range(n_boot):
        p = rng.choice(idx_pos, len(idx_pos), replace=True)
        q = rng.choice(idx_neg, len(idx_neg), replace=True)
        ii = np.concatenate([p, q])
        vals.append(roc_auc_score(y[ii], s[ii]))
    return float(point), float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def delong_like_paired_test(y: np.ndarray, s1: np.ndarray, s2: np.ndarray,
                            n_boot: int = 5000, seed: int = RANDOM_SEED) -> dict:
    """配对 bootstrap 检验两个打分函数的 AUC 差异（对接 vs 单描述符基线）。"""
    rng = np.random.default_rng(seed)
    idx_pos, idx_neg = np.where(y == 1)[0], np.where(y == 0)[0]
    diffs = []
    for _ in range(n_boot):
        p = rng.choice(idx_pos, len(idx_pos), replace=True)
        q = rng.choice(idx_neg, len(idx_neg), replace=True)
        ii = np.concatenate([p, q])
        diffs.append(roc_auc_score(y[ii], s1[ii]) - roc_auc_score(y[ii], s2[ii]))
    diffs = np.asarray(diffs)
    d = roc_auc_score(y, s1) - roc_auc_score(y, s2)
    # 双侧经验 p 值
    p = 2 * min((diffs <= 0).mean(), (diffs >= 0).mean())
    return {"auc_diff": round(float(d), 3),
            "ci95": [round(float(np.percentile(diffs, 2.5)), 3),
                     round(float(np.percentile(diffs, 97.5)), 3)],
            "p_empirical": round(float(min(1.0, p)), 4)}


def evaluate(df: pd.DataFrame, label: str) -> dict:
    d = df[df["best_affinity"].notna()].copy()
    y = d["active"].to_numpy()
    # Vina 亲和力越负越好 -> 取负号使"越大越好"
    score = -d["best_affinity"].to_numpy()

    if y.sum() == 0 or y.sum() == len(y):
        return {"label": label, "error": "缺少活性或非活性样本"}

    auc, lo, hi = auc_bootstrap_ci(y, score)
    res = {
        "label": label,
        "n": int(len(d)), "n_active": int(y.sum()), "n_inactive": int((y == 0).sum()),
        "docking": {
            "roc_auc": round(auc, 3), "auc_ci95": [round(lo, 3), round(hi, 3)],
            "EF1pct": round(enrichment_factor(y, score, 0.01), 2),
            "EF5pct": round(enrichment_factor(y, score, 0.05), 2),
            "EF10pct": round(enrichment_factor(y, score, 0.10), 2),
            "BEDROC_a20": round(bedroc(y, score, 20.0), 3),
            "mean_affinity_active": round(float(d.loc[y == 1, "best_affinity"].mean()), 2),
            "mean_affinity_inactive": round(float(d.loc[y == 0, "best_affinity"].mean()), 2),
        },
        "descriptor_baselines": {},
    }
    for col in ["MW", "HeavyAtoms", "RotB", "cLogP", "TPSA"]:
        if col not in d.columns:
            continue
        s = d[col].to_numpy(dtype=float)
        a = roc_auc_score(y, s)
        # 取方向更有利的一侧作为该描述符的最强基线（对基线宽容 = 对对接严格）
        if a < 0.5:
            s, a = -s, 1 - a
        res["descriptor_baselines"][col] = {
            "roc_auc": round(float(a), 3),
            "EF10pct": round(enrichment_factor(y, s, 0.10), 2),
            "BEDROC_a20": round(bedroc(y, s, 20.0), 3),
            "vs_docking": delong_like_paired_test(y, score, s),
        }
    best_base = max(res["descriptor_baselines"].items(),
                    key=lambda kv: kv[1]["roc_auc"], default=(None, None))
    if best_base[0]:
        res["strongest_baseline"] = best_base[0]
        res["docking_beats_strongest_baseline"] = bool(
            best_base[1]["vs_docking"]["ci95"][0] > 0
        )
    return res


def _plot(results: dict, scores: dict) -> None:
    setup_cjk_fonts()
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_curve

    labels = [k for k in scores if "active" in scores[k].columns]
    fig, axes = plt.subplots(1, len(labels), figsize=(5.4 * len(labels), 4.8), squeeze=False)
    for ax, key in zip(axes[0], labels):
        d = scores[key]
        d = d[d["best_affinity"].notna()]
        y = d["active"].to_numpy()
        if y.sum() == 0:
            continue
        s = -d["best_affinity"].to_numpy()
        fpr, tpr, _ = roc_curve(y, s)
        r = results[key]
        ax.plot(fpr, tpr, lw=2.4, color="#2b6cb0",
                label=f"Vina 对接  AUC={r['docking']['roc_auc']}")
        for col, c in [("MW", "#dd6b20"), ("RotB", "#718096")]:
            if col in r["descriptor_baselines"]:
                sv = d[col].to_numpy(float)
                if roc_auc_score(y, sv) < 0.5:
                    sv = -sv
                f2, t2, _ = roc_curve(y, sv)
                ax.plot(f2, t2, lw=1.6, ls="--", color=c,
                        label=f"仅 {col}  AUC={r['descriptor_baselines'][col]['roc_auc']}")
        ax.plot([0, 1], [0, 1], ":", color="#cbd5e0", lw=1.2, label="随机")
        ax.set_xlabel("假阳性率")
        ax.set_ylabel("真阳性率")
        ax.set_title(f"{key}\n{r['n_active']} 活性 / {r['n_inactive']} 非活性")
        ax.legend(frameon=False, fontsize=8.5, loc="lower right")
        ax.grid(alpha=0.25)
    fig.suptitle("A5：对接富集 vs 理化描述符基线", fontsize=13)
    fig.tight_layout()
    out = FIGS / "a5_enrichment.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log(f"图已保存 -> {out}")


def main(exhaustiveness: int = 16) -> dict:
    primary = pd.read_csv(RES_A / "benchmark_primary.csv")
    pooled = pd.read_csv(RES_A / "benchmark_pooled.csv")

    results, scores = {}, {}
    # Mpro：主受体 5WKK；用合并集跑（含主集全部化合物），再按文献切出主集
    log(f"Mpro 5WKK 对接筛选（{len(pooled)} 个化合物，exhaustiveness={exhaustiveness}）")
    dk = screen(pooled, "MPRO_5WKK", exhaustiveness)
    dk.to_csv(RES_A / "docking_scores_mpro_5wkk.csv", index=False)

    sub = dk[dk["document_chembl_id"] == "CHEMBL4495564"]
    scores["Mpro 5WKK · 主基准（单一筛选活动）"] = sub
    scores["Mpro 5WKK · 合并基准"] = dk
    results["Mpro 5WKK · 主基准（单一筛选活动）"] = evaluate(sub, "primary")
    results["Mpro 5WKK · 合并基准"] = evaluate(dk, "pooled")

    for k, r in results.items():
        if "error" in r:
            log(f"[{k}] {r['error']}")
            continue
        dd = r["docking"]
        log(f"[{k}] n={r['n']} ({r['n_active']}/{r['n_inactive']}) "
            f"对接 AUC={dd['roc_auc']} CI{dd['auc_ci95']} EF10%={dd['EF10pct']} "
            f"BEDROC={dd['BEDROC_a20']}")
        if "strongest_baseline" in r:
            b = r["descriptor_baselines"][r["strongest_baseline"]]
            log(f"    最强描述符基线 = {r['strongest_baseline']} AUC={b['roc_auc']} | "
                f"对接−基线 = {b['vs_docking']['auc_diff']} "
                f"CI{b['vs_docking']['ci95']} p={b['vs_docking']['p_empirical']} "
                f"-> 对接显著更优: {r['docking_beats_strongest_baseline']}")

    (RES_A / "a5_enrichment.json").write_text(json.dumps(results, indent=2, ensure_ascii=False))
    _plot(results, scores)
    log(f"A5 完成 -> {RES_A/'a5_enrichment.json'}")
    return results


if __name__ == "__main__":
    import sys
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 16)
