"""A1 —— ChEMBL 数据收集与整理。

产出（results/trackA/）：
  mers_mpro_ic50.csv      MERS-CoV Mpro (3CLpro) 生化 IC50
  mers_plpro_ic50.csv     MERS-CoV PLpro 生化 IC50 / Kd
  mers_cell_antiviral.csv MERS-CoV 细胞水平 EC50/IC50
  sars_mpro_train.csv     SARS-CoV-2 + SARS-CoV Mpro 迁移训练集
  a1_summary.json         各数据集条数与整理口径

整理规则
  * 仅保留 standard_relation == '=' 且单位可换算到 nM 的记录
  * RDKit 标准化（去盐、中和、规范化 SMILES）
  * drug-like 过滤 150 <= MW <= 650
  * 同一化合物多次测量取中位数 pChEMBL
  * 标注 Bemis-Murcko scaffold（供 scaffold-split 使用）
  * 记录 ChEMBL 版本（数据为 CC BY-SA 3.0，报告中须署名并相同方式共享）
"""
from __future__ import annotations

import json
import math
import time
from typing import Iterable

import pandas as pd
import requests
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

from .config import CHEMBL_API, CHEMBL_DIR, MAX_MW, MIN_MW, MPRO_KEYS, PLPRO_KEYS, RES_A, TARGETS, log

RDLogger.DisableLog("rdApp.*")

SESSION = requests.Session()
SESSION.headers.update({"Accept": "application/json"})

# 单位换算到 nM
UNIT_TO_NM = {"nM": 1.0, "uM": 1e3, "µM": 1e3, "mM": 1e6, "M": 1e9, "pM": 1e-3}


def _get(url: str, params: dict, retries: int = 4) -> dict:
    for attempt in range(retries):
        try:
            r = SESSION.get(url, params=params, timeout=120)
            r.raise_for_status()
            return r.json()
        except Exception as exc:  # noqa: BLE001
            if attempt == retries - 1:
                raise
            wait = 2 ** (attempt + 1)
            log(f"  请求失败({exc.__class__.__name__})，{wait}s 后重试")
            time.sleep(wait)
    raise RuntimeError("unreachable")


def chembl_version() -> str:
    return _get(f"{CHEMBL_API}/status.json", {}).get("chembl_db_version", "unknown")


def fetch_activities(target_chembl_id: str, cache_name: str) -> pd.DataFrame:
    """分页抓取某靶标的全部活性记录，带本地缓存。"""
    cache = CHEMBL_DIR / f"{cache_name}.raw.csv"
    if cache.exists():
        log(f"  使用缓存 {cache.name}")
        return pd.read_csv(cache)

    rows: list[dict] = []
    offset, limit = 0, 1000
    while True:
        d = _get(
            f"{CHEMBL_API}/activity.json",
            {"target_chembl_id": target_chembl_id, "limit": limit, "offset": offset},
        )
        acts = d["activities"]
        rows.extend(acts)
        total = d["page_meta"]["total_count"]
        offset += limit
        log(f"  {target_chembl_id}: {min(offset, total)}/{total}")
        if offset >= total or not acts:
            break

    keep = [
        "molecule_chembl_id", "canonical_smiles", "standard_type", "standard_relation",
        "standard_value", "standard_units", "pchembl_value", "assay_chembl_id",
        "assay_description", "assay_type", "document_chembl_id", "target_chembl_id",
        "activity_comment", "data_validity_comment",
    ]
    df = pd.DataFrame(rows)
    for c in keep:
        if c not in df.columns:
            df[c] = None
    df = df[keep]
    df.to_csv(cache, index=False)
    return df


# ------------------------------------------------------------------ 标准化
def standardize_smiles(smi: str) -> tuple[str | None, float | None]:
    """去盐 -> 取最大片段 -> 规范化 SMILES；返回 (smiles, MW)。"""
    if not isinstance(smi, str) or not smi:
        return None, None
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None, None
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not frags:
        return None, None
    mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    try:
        Chem.SanitizeMol(mol)
    except Exception:  # noqa: BLE001
        return None, None
    return Chem.MolToSmiles(mol), Descriptors.MolWt(mol)


def murcko_scaffold(smi: str) -> str:
    try:
        sc = MurckoScaffold.MurckoScaffoldSmiles(smiles=smi, includeChirality=False)
        return sc or "ACYCLIC"
    except Exception:  # noqa: BLE001
        return "ACYCLIC"


def classify_enzyme(desc: str) -> str:
    d = (desc or "").lower()
    if any(k in d for k in PLPRO_KEYS):
        return "PLpro"
    if any(k in d for k in MPRO_KEYS):
        return "Mpro"
    return "other"


def to_nm(value, units) -> float | None:
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    f = UNIT_TO_NM.get(str(units))
    if f is None or v <= 0:
        return None
    return v * f


def clean(
    df: pd.DataFrame,
    types: Iterable[str],
    label: str,
    drug_like: bool = True,
) -> pd.DataFrame:
    """应用整理规则，返回每化合物一行的数据集。"""
    types = set(types)
    d = df[df["standard_type"].isin(types)].copy()
    n0 = len(d)
    # 只要精确关系
    d = d[d["standard_relation"] == "="]
    n1 = len(d)
    # 丢弃 ChEMBL 标记的可疑数据
    d = d[d["data_validity_comment"].isna()]
    # 单位换算
    d["value_nM"] = [to_nm(v, u) for v, u in zip(d["standard_value"], d["standard_units"])]
    d = d[d["value_nM"].notna() & (d["value_nM"] > 0)]
    n2 = len(d)

    # 结构标准化
    std = [standardize_smiles(s) for s in d["canonical_smiles"]]
    d["smiles_std"] = [s for s, _ in std]
    d["MW"] = [m for _, m in std]
    d = d[d["smiles_std"].notna()]
    n3 = len(d)

    if drug_like:
        d = d[(d["MW"] >= MIN_MW) & (d["MW"] <= MAX_MW)]
    n4 = len(d)

    d["pActivity"] = [-math.log10(v * 1e-9) for v in d["value_nM"]]

    # 同一化合物多次测量 -> 中位数
    agg = (
        d.groupby("smiles_std")
        .agg(
            molecule_chembl_id=("molecule_chembl_id", "first"),
            standard_type=("standard_type", lambda s: "/".join(sorted(set(s)))),
            value_nM=("value_nM", "median"),
            pActivity=("pActivity", "median"),
            n_measurements=("pActivity", "size"),
            MW=("MW", "first"),
            n_documents=("document_chembl_id", lambda s: s.nunique()),
        )
        .reset_index()
    )
    agg["scaffold"] = [murcko_scaffold(s) for s in agg["smiles_std"]]
    agg["dataset"] = label
    agg = agg.sort_values("pActivity", ascending=False).reset_index(drop=True)

    log(
        f"  [{label}] 类型过滤 {n0} -> '=' {n1} -> 单位 {n2} -> 结构 {n3} "
        f"-> drug-like {n4} -> 去重后化合物 {len(agg)}（scaffold {agg['scaffold'].nunique()}）"
    )
    return agg


def main() -> dict:
    ver = chembl_version()
    log(f"ChEMBL 版本: {ver}（数据许可 CC BY-SA 3.0）")
    summary: dict = {"chembl_version": ver, "datasets": {}}

    # ---- MERS 蛋白酶 ----
    log("抓取 MERS-CoV replicase pp1ab 活性…")
    mers = fetch_activities(TARGETS["MERS_PP1AB"], "mers_pp1ab")
    mers["enzyme"] = [classify_enzyme(x) for x in mers["assay_description"]]
    log(f"  酶归类: {mers['enzyme'].value_counts().to_dict()}")

    mpro = clean(mers[mers["enzyme"] == "Mpro"], {"IC50"}, "MERS_Mpro_IC50")
    plpro = clean(mers[mers["enzyme"] == "PLpro"], {"IC50", "Kd", "Ki"}, "MERS_PLpro_IC50_Kd")

    mpro.to_csv(RES_A / "mers_mpro_ic50.csv", index=False)
    plpro.to_csv(RES_A / "mers_plpro_ic50.csv", index=False)

    # ---- MERS 细胞水平 ----
    log("抓取 MERS-CoV 细胞水平抗病毒活性…")
    cell_raw = fetch_activities(TARGETS["MERS_CELL"], "mers_cell")
    cell = clean(cell_raw, {"EC50", "IC50", "EC90"}, "MERS_cell_antiviral")
    cell.to_csv(RES_A / "mers_cell_antiviral.csv", index=False)

    # ---- 跨冠状病毒迁移训练集 ----
    log("抓取 SARS-CoV-2 / SARS-CoV Mpro 迁移训练集…")
    frames = []
    for key, name in [("SARS2_MPRO", "sars2_mpro"), ("SARS_PP1AB_1", "sars_pp1ab_1"),
                      ("SARS_PP1AB_2", "sars_pp1ab_2")]:
        raw = fetch_activities(TARGETS[key], name)
        raw["enzyme"] = [classify_enzyme(x) for x in raw["assay_description"]]
        sub = raw[raw["enzyme"] == "Mpro"]
        log(f"  {TARGETS[key]}: 总 {len(raw)} 条，其中 Mpro 相关 {len(sub)} 条")
        if len(sub):
            sub = sub.copy()
            sub["source_target"] = TARGETS[key]
            frames.append(sub)
    sars_all = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    sars = clean(sars_all, {"IC50", "Ki"}, "SARS_Mpro_train") if len(sars_all) else pd.DataFrame()
    if len(sars):
        sars.to_csv(RES_A / "sars_mpro_train.csv", index=False)

    # ---- 训练/测试泄漏检查 ----
    overlap = set(sars["smiles_std"]) & set(mpro["smiles_std"]) if len(sars) else set()
    sc_overlap = set(sars["scaffold"]) & set(mpro["scaffold"]) if len(sars) else set()
    log(f"迁移训练集与 MERS Mpro 测试集: 完全相同结构 {len(overlap)} 个，共享 scaffold {len(sc_overlap)} 个")

    for name, d in [
        ("MERS_Mpro_IC50", mpro), ("MERS_PLpro", plpro),
        ("MERS_cell_antiviral", cell), ("SARS_Mpro_train", sars),
    ]:
        summary["datasets"][name] = {
            "n_compounds": int(len(d)),
            "n_scaffolds": int(d["scaffold"].nunique()) if len(d) else 0,
            "pActivity_min": float(d["pActivity"].min()) if len(d) else None,
            "pActivity_max": float(d["pActivity"].max()) if len(d) else None,
            "pActivity_median": float(d["pActivity"].median()) if len(d) else None,
        }
    summary["leakage_check"] = {
        "identical_structures_train_vs_MERS_test": len(overlap),
        "shared_scaffolds_train_vs_MERS_test": len(sc_overlap),
        "overlapping_smiles": sorted(overlap),
    }
    summary["rules"] = {
        "relation": "= only",
        "units": "converted to nM",
        "drug_like": f"{MIN_MW} <= MW <= {MAX_MW}",
        "duplicates": "median pActivity per standardized SMILES",
        "scaffold": "Bemis-Murcko (no chirality)",
    }
    (RES_A / "a1_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    log(f"A1 完成 -> {RES_A/'a1_summary.json'}")
    return summary


if __name__ == "__main__":
    main()
