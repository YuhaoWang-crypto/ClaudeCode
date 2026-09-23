"""A4 —— 构建 MERS-CoV Mpro 的活性/非活性基准集，并诊断混杂因素。

为什么不用合成诱饵（decoys）：
  ChEMBL 里 MERS Mpro 有 74 条 `IC50 > 10 µM` 的记录，来自同一个 FRET 筛选活动
  （文献 CHEMBL4495564）。这些是**实测的非活性化合物**，比性质匹配的合成诱饵
  严格得多 —— 合成诱饵只是"假定"不活性，而这些是测过的。

为什么要分两个基准集：
  活性物来自两篇文献。CHEMBL4196085 是一篇药化论文（12 个化合物，IC50 400-7500 nM，
  分子更大）；CHEMBL4495564 是筛选活动本身（10 个活性 + 74 个非活性）。
  把两者合并会引入**分子量混杂**：合并集上"仅用分子量"就能拿到 AUC 0.79，
  而单一活动集上仅用分子量只有 0.52（≈随机）。

  => 主基准 = 单一活动集（无混杂，但只有 8 个活性，统计功效低）
     次要基准 = 合并集（活性物更多，但必须与分子量基线并列报告）

产出：results/trackA/benchmark_primary.csv / benchmark_pooled.csv
      results/trackA/a4_benchmark.json
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score

from .a1_chembl_data import classify_enzyme, murcko_scaffold, standardize_smiles
from .config import CHEMBL_DIR, MAX_MW, MIN_MW, RES_A, log

RDLogger.DisableLog("rdApp.*")

PRIMARY_DOC = "CHEMBL4495564"   # 单一 FRET 筛选活动：活性与非活性同源


def descriptors(smi: str) -> dict:
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return {}
    return {
        "MW": Descriptors.MolWt(m),
        "cLogP": Crippen.MolLogP(m),
        "TPSA": rdMolDescriptors.CalcTPSA(m),
        "HBD": rdMolDescriptors.CalcNumHBD(m),
        "HBA": rdMolDescriptors.CalcNumHBA(m),
        "RotB": rdMolDescriptors.CalcNumRotatableBonds(m),
        "HeavyAtoms": m.GetNumHeavyAtoms(),
        "RingCount": rdMolDescriptors.CalcNumRings(m),
        "FracCSP3": rdMolDescriptors.CalcFractionCSP3(m),
    }


def build() -> pd.DataFrame:
    """把 MERS Mpro 的全部 IC50 记录整理成带 active 标签的化合物表。"""
    d = pd.read_csv(CHEMBL_DIR / "mers_pp1ab.raw.csv")
    d["enzyme"] = [classify_enzyme(x) for x in d["assay_description"]]
    g = d[(d["enzyme"] == "Mpro") & (d["standard_type"] == "IC50")].copy()
    g = g[g["standard_units"] == "nM"]
    g = g[g["data_validity_comment"].isna()]

    std = [standardize_smiles(s) for s in g["canonical_smiles"]]
    g["smiles_std"] = [s for s, _ in std]
    g = g[g["smiles_std"].notna()]
    g["value_nM"] = g["standard_value"].astype(float)
    g["active"] = (g["standard_relation"] == "=").astype(int)

    desc = pd.DataFrame([descriptors(s) for s in g["smiles_std"]], index=g.index)
    g = pd.concat([g, desc], axis=1)
    g = g[g["MW"].between(MIN_MW, MAX_MW)]

    # 每个结构一行；若既有活性又有非活性记录，以活性为准（更保守的标签）
    g = g.sort_values("active", ascending=False).drop_duplicates("smiles_std")
    g["scaffold"] = [murcko_scaffold(s) for s in g["smiles_std"]]
    g["pActivity"] = np.where(g["active"] == 1, -np.log10(g["value_nM"] * 1e-9), np.nan)
    cols = ["smiles_std", "molecule_chembl_id", "active", "value_nM", "pActivity",
            "standard_relation", "document_chembl_id", "assay_chembl_id", "scaffold"] + list(desc.columns)
    return g[cols].reset_index(drop=True)


def confound_report(df: pd.DataFrame, name: str) -> dict:
    """对每个描述符检查"单变量就能分类"的程度 —— 任何模型都必须显著超过它。"""
    out = {"name": name, "n_active": int(df["active"].sum()),
           "n_inactive": int((df["active"] == 0).sum()), "univariate_baselines": {}}
    for col in ["MW", "cLogP", "TPSA", "HBD", "HBA", "RotB", "HeavyAtoms",
                "RingCount", "FracCSP3"]:
        a, b = df.loc[df["active"] == 1, col], df.loc[df["active"] == 0, col]
        if a.std() == 0 and b.std() == 0:
            continue
        auc = roc_auc_score(df["active"], df[col])
        # AUC < 0.5 说明该描述符反向有判别力，取对称值看绝对判别强度
        auc_sym = max(auc, 1 - auc)
        p = mannwhitneyu(a, b).pvalue
        out["univariate_baselines"][col] = {
            "auc": round(float(auc), 3), "auc_symmetric": round(float(auc_sym), 3),
            "mannwhitney_p": float(p),
            "active_mean": round(float(a.mean()), 2), "inactive_mean": round(float(b.mean()), 2),
        }
    strongest = max(out["univariate_baselines"].items(), key=lambda kv: kv[1]["auc_symmetric"])
    out["strongest_single_descriptor"] = {"descriptor": strongest[0], **strongest[1]}
    return out


def main() -> dict:
    all_df = build()
    log(f"MERS Mpro IC50 化合物总表: {len(all_df)} 个（活性 {int(all_df['active'].sum())}，"
        f"非活性 {int((all_df['active']==0).sum())}）")

    primary = all_df[all_df["document_chembl_id"] == PRIMARY_DOC].reset_index(drop=True)
    pooled = all_df.reset_index(drop=True)

    primary.to_csv(RES_A / "benchmark_primary.csv", index=False)
    pooled.to_csv(RES_A / "benchmark_pooled.csv", index=False)

    rep_p = confound_report(primary, "primary_single_campaign")
    rep_o = confound_report(pooled, "pooled_two_documents")

    for rep in (rep_p, rep_o):
        s = rep["strongest_single_descriptor"]
        log(f"[{rep['name']}] {rep['n_active']} 活性 / {rep['n_inactive']} 非活性 | "
            f"最强单描述符 = {s['descriptor']} (AUC {s['auc']}, 对称 {s['auc_symmetric']}, "
            f"p={s['mannwhitney_p']:.3g})")

    out = {
        "primary_benchmark": {
            "document": PRIMARY_DOC,
            "rationale": "活性与非活性来自同一个 FRET 筛选活动，无跨文献批次混杂",
            **rep_p,
            "n_scaffolds": int(primary["scaffold"].nunique()),
            "caveat": "只有 8 个活性物，统计功效很低；AUC 的 95% CI 会很宽",
        },
        "pooled_benchmark": {
            "documents": sorted(pooled["document_chembl_id"].dropna().unique().tolist()),
            **rep_o,
            "n_scaffolds": int(pooled["scaffold"].nunique()),
            "caveat": (
                f"活性物混合了药化论文（分子更大）与筛选活动，引入分子量混杂："
                f"仅用分子量即可达到 AUC {rep_o['univariate_baselines']['MW']['auc']}。"
                f"任何模型在此集上的表现必须与该基线并列解读。"
            ),
        },
        "decoy_policy": (
            "使用 ChEMBL 中 IC50 > 10 µM 的**实测非活性**化合物，而非性质匹配的合成诱饵。"
            "合成诱饵只是假定不活性；实测非活性是同一实验测过的阴性结果。"
        ),
    }
    (RES_A / "a4_benchmark.json").write_text(json.dumps(out, indent=2, ensure_ascii=False))
    log(f"A4 完成 -> {RES_A/'a4_benchmark.json'}")
    return out


if __name__ == "__main__":
    main()
