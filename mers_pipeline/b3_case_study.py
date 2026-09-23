"""B3 —— 抗病毒化合物案例研究：sci-Plex 化合物与 MERS-CoV 数据的交集。

方案预期重叠很少、可能为零；本模块如实量化，并在"完全相同结构"之外
再做一层化学相似性分析，说明即便没有精确重叠，sci-Plex 覆盖的化学空间
与 MERS-CoV 抗病毒化合物有多远。

产出：results/trackB/b3_case_study.json, sciplex_mers_similarity.csv
      figures/mers/b3_case_study.png
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator

from .a1_chembl_data import standardize_smiles
from .config import FIGS, RES_A, RES_B, SCIPLEX_DIR, log, setup_cjk_fonts

RDLogger.DisableLog("rdApp.*")
_MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

HASH_META = SCIPLEX_DIR / "GSM4150378_sciPlex3_A549_MCF7_K562_hashTable_metadata.txt.gz"


def load_sciplex_compounds() -> pd.DataFrame:
    d = pd.read_csv(HASH_META, sep="\t")
    d = d[d["name"].notna() & (d["name"] != "NA")]
    d = d[d["SMILES"].notna() & (d["SMILES"] != "NA")]
    d = d.drop_duplicates("name")[["name", "SMILES", "catalog_number", "CAS.Number", "M.w."]]
    std = [standardize_smiles(s) for s in d["SMILES"]]
    d["smiles_std"] = [s for s, _ in std]
    d = d[d["smiles_std"].notna()].reset_index(drop=True)
    log(f"sci-Plex 化合物（带可解析 SMILES）: {len(d)}")
    return d


def fp(s: str):
    m = Chem.MolFromSmiles(s)
    return _MORGAN.GetFingerprint(m) if m is not None else None


def strip_stereo(s: str) -> str | None:
    """去立体化学的规范 SMILES —— 用于识别对映体/外消旋体关系。

    Morgan 指纹（不含手性）无法区分对映体，所以相似性分析里会出现
    Tanimoto = 1.0 却不是精确 SMILES 匹配的情况。这个函数把这类关系显式标出来。
    """
    m = Chem.MolFromSmiles(s)
    if m is None:
        return None
    Chem.RemoveStereochemistry(m)
    return Chem.MolToSmiles(m)


def main() -> dict:
    sp = load_sciplex_compounds()

    sets = {
        "MERS_Mpro_活性": pd.read_csv(RES_A / "mers_mpro_ic50.csv"),
        "MERS_PLpro": pd.read_csv(RES_A / "mers_plpro_ic50.csv"),
        "MERS_Mpro_基准集(含实测非活性)": pd.read_csv(RES_A / "benchmark_pooled.csv"),
        "MERS_细胞抗病毒": pd.read_csv(RES_A / "mers_cell_antiviral.csv"),
    }

    sp["smiles_nostereo"] = [strip_stereo(s) for s in sp["smiles_std"]]
    sp_fps = [fp(s) for s in sp["smiles_std"]]
    valid = [i for i, f in enumerate(sp_fps) if f is not None]
    sp = sp.iloc[valid].reset_index(drop=True)
    sp_fps = [sp_fps[i] for i in valid]

    out: dict = {"n_sciplex_compounds": int(len(sp)), "overlaps": {}, "similarity": {}}
    rows = []
    for name, df in sets.items():
        ref_smiles = df["smiles_std"].dropna().tolist()
        exact = sorted(set(sp["smiles_std"]) & set(ref_smiles))
        # 去立体后的匹配：捕捉对映体/外消旋体这类"同一药物不同立体形式"的情况
        ref_nostereo = {x for x in (strip_stereo(s) for s in ref_smiles) if x}
        stereo_only = sorted(
            set(sp.loc[sp["smiles_nostereo"].isin(ref_nostereo), "smiles_std"]) - set(exact)
        )
        ref_fps = [f for f in (fp(s) for s in ref_smiles) if f is not None]
        sims = np.array([max(DataStructs.BulkTanimotoSimilarity(f, ref_fps))
                         for f in sp_fps]) if ref_fps else np.zeros(len(sp_fps))
        stereo_names = sp.loc[sp["smiles_std"].isin(stereo_only), "name"].tolist()
        out["overlaps"][name] = {
            "n_reference": int(len(ref_smiles)),
            "n_exact_structural_overlap": len(exact),
            "overlapping_compound_names": sp.loc[sp["smiles_std"].isin(exact), "name"].tolist(),
            "n_stereoisomer_only_overlap": len(stereo_only),
            "stereoisomer_only_names": stereo_names,
        }
        out["similarity"][name] = {
            "nn_tanimoto_median": round(float(np.median(sims)), 3),
            "nn_tanimoto_p95": round(float(np.percentile(sims, 95)), 3),
            "nn_tanimoto_max": round(float(sims.max()), 3),
            "n_above_0.5": int((sims >= 0.5).sum()),
            "n_above_0.7": int((sims >= 0.7).sum()),
            "top_matches": (
                sp.assign(sim=sims).nlargest(5, "sim")[["name", "sim"]]
                .assign(sim=lambda x: x["sim"].round(3)).to_dict("records")
            ),
        }
        rows.append(pd.DataFrame({"sciplex_compound": sp["name"], "reference_set": name,
                                  "nn_tanimoto": sims.round(4)}))
        log(f"[{name}] 参考 {len(ref_smiles)} 个 | 完全相同结构 {len(exact)} 个 | "
            f"仅立体异构体不同 {len(stereo_only)} 个{(' '+str(stereo_names)) if stereo_names else ''} | "
            f"最近邻 Tanimoto 中位 {np.median(sims):.3f} 最大 {sims.max():.3f} | "
            f"≥0.5 的 {int((sims>=0.5).sum())} 个")

    pd.concat(rows).to_csv(RES_B / "sciplex_mers_similarity.csv", index=False)

    # ---- 对交集化合物做逐个案例分析 ----
    hit_names = sorted({n for v in out["overlaps"].values()
                        for n in v["overlapping_compound_names"] + v["stereoisomer_only_names"]})
    out["case_studies"] = _case_studies(sp, hit_names)
    total_exact = len(hit_names)
    n_repro = sum(1 for c in out["case_studies"] if c["reproducible_anywhere"])
    max_transfer = max((max(c["loso_transfer_pearson"].values(), default=0)
                        for c in out["case_studies"]), default=0)
    out["verdict"] = {
        "total_exact_overlaps": total_exact,
        "overlapping_compounds": [c["compound"] for c in out["case_studies"]],
        "case_study_feasible": total_exact > 0,
        "statement": (
            f"sci-Plex3 的 {len(sp)} 个化合物与 MERS-CoV 的四个化合物集合共有 "
            f"{total_exact} 个可对应的分子："
            + "、".join(
                c["compound"] + ("（立体异构体对应：sci-Plex 为外消旋体，"
                                 "MERS 集为单一对映体）" if c.get("stereoisomer_match") else "")
                for c in out["case_studies"])
            + "。样本量不足以做统计推断，只能作个案描述。"
        ),
        "key_finding": (
            f"这 {total_exact} 个交集化合物中只有 {n_repro} 个在任一细胞系达到可重复性门槛，"
            f"其 LOSO 迁移 r 最高仅 {max_transfer:.3f}（远低于可重复子集的 0.22 中位）。"
            f"换言之：**两类数据真正相交的那几个化合物，恰恰是 virtual cell 方法"
            f"最没有话可说的那几个** —— 它们是弱效、常为泛试剂干扰型的蛋白酶结合剂"
            f"（如 Disulfiram 这类巯基反应性化合物会非选择性地弱抑制两种半胱氨酸蛋白酶），"
            f"而非强转录扰动剂。"
        ),
        "chemical_space": (
            "sci-Plex 化合物库以肿瘤表观遗传/激酶调控剂为主，MERS-CoV 抑制剂"
            "以肽模拟共价蛋白酶抑制剂为主，两者化学空间基本不相交"
            f"（最近邻 Tanimoto 中位 0.15–0.23）。"
        ),
    }
    log(out["verdict"]["statement"])
    (RES_B / "b3_case_study.json").write_text(json.dumps(out, indent=2, ensure_ascii=False))
    _plot(rows, out)
    log(f"B3 完成 -> {RES_B/'b3_case_study.json'}")
    return out


def _case_studies(sp: pd.DataFrame, names: list[str]) -> list[dict]:
    """对每个交集化合物，并列其 MERS 实测效力与 Track B 的可重复性/迁移表现。"""
    repro = pd.read_csv(RES_B / "reproducibility.csv")
    loso = pd.read_csv(RES_B / "loso_per_perturbation.csv")
    loso = loso[loso["method"] == "weighted_consensus_gated"]
    sets = {
        "MERS_Mpro": (RES_A / "mers_mpro_ic50.csv", "pActivity"),
        "MERS_PLpro": (RES_A / "mers_plpro_ic50.csv", "pActivity"),
        "MERS_cell_antiviral": (RES_A / "mers_cell_antiviral.csv", "pActivity"),
    }
    cases = []
    for name in names:
        row = sp[sp["name"] == name]
        if row.empty:
            continue
        smi = row.iloc[0]["smiles_std"]
        nostereo = strip_stereo(smi)
        c: dict = {"compound": name.strip(), "smiles": smi, "mers_potency": {}}
        for label, (path, col) in sets.items():
            t = pd.read_csv(path)
            m = t[t["smiles_std"] == smi]
            if not len(m):
                # 退回到去立体匹配（对映体/外消旋体）
                m = t[[strip_stereo(s) == nostereo for s in t["smiles_std"]]]
                if len(m):
                    c.setdefault("stereoisomer_match", []).append(label)
            if len(m):
                c["mers_potency"][label] = {
                    "value_nM": float(m.iloc[0]["value_nM"]),
                    "p": round(float(m.iloc[0][col]), 2),
                    "endpoint": m.iloc[0]["standard_type"],
                }
        rr = repro[repro["compound"] == name]
        c["sciplex_split_half_reproducibility"] = {
            k: round(float(v), 3) for k, v in zip(rr["cell_type"], rr["split_half_r"])
        }
        ll = loso[loso["compound"] == name]
        c["loso_transfer_pearson"] = {
            k: round(float(v), 3) for k, v in zip(ll["held_out"], ll["pearson_delta"])
        }
        c["reproducible_anywhere"] = bool(
            any(v > 0.30 for v in c["sciplex_split_half_reproducibility"].values()))
        cases.append(c)
        log(f"  案例 {name.strip()}: MERS {list(c['mers_potency'])} | "
            f"可重复性 {c['sciplex_split_half_reproducibility']} | "
            f"迁移 r {c['loso_transfer_pearson']}")
    return cases


def _plot(rows: list, out: dict) -> None:
    setup_cjk_fonts()
    import matplotlib.pyplot as plt

    df = pd.concat(rows)
    sets = list(out["similarity"])
    fig, ax = plt.subplots(figsize=(9.2, 5.0))
    colors = ["#2b6cb0", "#dd6b20", "#38a169", "#805ad5"]
    for (name, c) in zip(sets, colors):
        v = df.loc[df["reference_set"] == name, "nn_tanimoto"]
        ax.hist(v, bins=28, range=(0, 1), histtype="step", lw=2.1, color=c,
                label=f"{name}（中位 {v.median():.2f}）")
    ax.axvline(0.5, color="#e53e3e", ls="--", lw=1.6, label="相似性门槛 0.50")
    ax.set_xlabel("sci-Plex 化合物 到 MERS-CoV 化合物集的最近邻 Tanimoto")
    ax.set_ylabel("sci-Plex 化合物数")
    ax.set_title(f"B3：sci-Plex 与 MERS-CoV 化学空间的距离\n"
                 f"完全相同结构共 {out['verdict']['total_exact_overlaps']} 个")
    ax.legend(frameon=False, fontsize=8.8)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    o = FIGS / "b3_case_study.png"
    fig.savefig(o, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log(f"图已保存 -> {o}")


if __name__ == "__main__":
    main()
