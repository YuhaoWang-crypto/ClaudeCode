"""
真实数据层 — 从 ChEMBL 拉取配对的"非细胞 (生化)"与"细胞法"抗病毒活性。

核心思想
--------
ChEMBL 的 BAO (BioAssay Ontology) `bao_format` 字段天然区分测定形式:

    BAO_0000357  single protein format   -> 纯酶/纯蛋白，**非细胞**
    BAO_0000219  cell-based format       -> **细胞法**
    BAO_0000218  organism-based format   -> 整体动物
    BAO_0000221  tissue-based format

于是同一个分子若同时出现在两类测定里，就得到一对
(生化 pIC50, 细胞 pEC50) —— 这正是"能否用非细胞数据预测细胞数据"所需的真值。

这比用 assay_description 做文本匹配可靠得多。

数据来源
--------
ChEMBL_37 (2026-05-01), https://www.ebi.ac.uk/chembl/api/data
所有数据实时抓取并缓存到磁盘；本模块不含任何硬编码的活性数值。
"""

from __future__ import annotations

import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

BASE = "https://www.ebi.ac.uk/chembl/api/data"

CACHE_DIR = Path(__file__).resolve().parent.parent / "data" / "chembl"

# BAO 测定形式
BAO_PROTEIN = "BAO_0000357"  # 非细胞：纯蛋白
BAO_CELL = "BAO_0000219"  # 细胞法 (宿主细胞读数)
BAO_ORGANISM = "BAO_0000218"  # "生物体" —— 对病毒靶点即细胞培养内的抗病毒实验
BAO_TISSUE = "BAO_0000221"

# ⚠️ 关键领域细节：ChEMBL 把针对**病毒**靶点的细胞培养抗病毒实验标注为
# BAO_0000218 (organism-based)，因为这里的"生物体"指病毒本身。HIV / 流感
# 都无法脱离细胞培养，所以这些记录本质上就是细胞法。
# 若只认 BAO_0000219，会漏掉 HIV-1 organism 靶点 42470 条里的 40331 条。
CELLULAR_FORMATS = {BAO_CELL, BAO_ORGANISM}

# 我们关心的字段（用 only= 缩减 payload）
FIELDS = [
    "molecule_chembl_id",
    "canonical_smiles",
    "standard_type",
    "standard_value",
    "standard_units",
    "standard_relation",
    "pchembl_value",
    "bao_format",
    "bao_label",
    "assay_chembl_id",
    "assay_description",
    "assay_type",
    "target_chembl_id",
    "target_organism",
    "document_chembl_id",
    "document_year",
    "data_validity_comment",
    "potential_duplicate",
]

# 研究体系：每个 = (生化靶点 protein target, 细胞法靶点 organism target)
SYSTEMS = {
    "HIV1_RT": {
        "label": "HIV-1 逆转录酶",
        "biochem_target": "CHEMBL247",  # HIV-1 RT, SINGLE PROTEIN
        "cell_target": "CHEMBL378",  # HIV-1, ORGANISM
        "note": "NRTI/NNRTI 文献惯例同时报告酶抑制 IC50 与 MT-4/CEM 细胞 EC50",
    },
    "FLU_NA": {
        "label": "流感 A 神经氨酸酶",
        "biochem_target": "CHEMBL2051",  # Influenza A NA (A/PR/8/34), SINGLE PROTEIN
        "cell_target": "CHEMBL613740",  # Influenza A virus, ORGANISM
        "note": "NA 酶抑制 IC50 (MUNANA) 与 MDCK 细胞 CPE/空斑 EC50",
    },
}


def _get(path: str, retries: int = 4, **params) -> dict:
    url = f"{BASE}/{path}?{urllib.parse.urlencode(params)}"
    delay = 2.0
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=120) as fh:
                return json.load(fh)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            if attempt == retries - 1:
                raise RuntimeError(f"ChEMBL 请求失败: {url}") from exc
            time.sleep(delay)
            delay *= 2
    raise AssertionError("unreachable")


def fetch_target_activities(target_chembl_id: str, *, force: bool = False,
                            page_size: int = 1000, verbose: bool = True) -> list[dict]:
    """拉取某靶点全部带 pchembl_value 的活性记录，缓存到磁盘。"""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache = CACHE_DIR / f"{target_chembl_id}.json"
    if cache.exists() and not force:
        return json.loads(cache.read_text())

    out: list[dict] = []
    offset = 0
    total = None
    while True:
        r = _get(
            "activity.json",
            target_chembl_id=target_chembl_id,
            pchembl_value__isnull="false",
            only=",".join(FIELDS),
            limit=page_size,
            offset=offset,
        )
        if total is None:
            total = r["page_meta"]["total_count"]
            if verbose:
                print(f"  {target_chembl_id}: {total} 条活性记录，开始下载…")
        acts = r["activities"]
        out.extend(acts)
        offset += len(acts)
        if not acts or offset >= total:
            break
        if verbose and offset % 5000 == 0:
            print(f"    …{offset}/{total}")

    cache.write_text(json.dumps(out))
    if verbose:
        print(f"  {target_chembl_id}: 已缓存 {len(out)} 条 -> {cache.name}")
    return out


# --- 清洗 ----------------------------------------------------------------

# 只接受这些终点。IC50/EC50 是主力；Ki/Kd 也保留（生化侧）
BIOCHEM_TYPES = {"IC50", "Ki", "Kd"}
CELL_TYPES = {"EC50", "IC50"}


def _usable(a: dict) -> bool:
    """基本质量门：精确值、有 pchembl、无数据校验警告。"""
    if a.get("data_validity_comment"):
        return False
    if a.get("standard_relation") not in ("=", None):
        return False
    if a.get("pchembl_value") in (None, ""):
        return False
    try:
        float(a["pchembl_value"])
    except (TypeError, ValueError):
        return False
    return True


def split_by_format(activities: list[dict]) -> dict[str, list[dict]]:
    """按 BAO 测定形式分组。"""
    groups: dict[str, list[dict]] = {}
    for a in activities:
        groups.setdefault(a.get("bao_format") or "UNKNOWN", []).append(a)
    return groups


def median(xs: list[float]) -> float:
    s = sorted(xs)
    n = len(s)
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def collapse_per_molecule(activities: list[dict], allowed_types: set[str],
                          bao_format: str | set[str]) -> dict[str, dict]:
    """把同一分子的多条记录聚合成一个中位数 pchembl。

    同一化合物在不同论文里被反复测量，取中位数比取任意一条稳健，
    并且保留 n 与离散度，用于后面判断"真值本身有多不确定"。
    """
    accept = {bao_format} if isinstance(bao_format, str) else set(bao_format)
    by_mol: dict[str, list[dict]] = {}
    for a in activities:
        if a.get("bao_format") not in accept:
            continue
        if a.get("standard_type") not in allowed_types:
            continue
        if not _usable(a):
            continue
        by_mol.setdefault(a["molecule_chembl_id"], []).append(a)

    out: dict[str, dict] = {}
    for mol, recs in by_mol.items():
        vals = [float(r["pchembl_value"]) for r in recs]
        out[mol] = {
            "pchembl": median(vals),
            "n": len(vals),
            "spread": (max(vals) - min(vals)) if len(vals) > 1 else 0.0,
            "values": vals,
            "smiles": recs[0].get("canonical_smiles"),
            "types": sorted({r["standard_type"] for r in recs}),
        }
    return out


def build_pairs(system: str, *, force: bool = False, verbose: bool = True) -> list[dict]:
    """构建 (生化 pIC50, 细胞 pEC50) 配对表。

    返回每个元素:
        {molecule, p_biochem, p_cell, n_biochem, n_cell,
         spread_biochem, spread_cell, smiles, delta}
    delta = p_biochem - p_cell  (正值 = 细胞里效力更弱，即存在"细胞惩罚")
    """
    spec = SYSTEMS[system]
    if verbose:
        print(f"[{system}] {spec['label']}")

    bio_acts = fetch_target_activities(spec["biochem_target"], force=force, verbose=verbose)
    cell_acts = fetch_target_activities(spec["cell_target"], force=force, verbose=verbose)

    bio = collapse_per_molecule(bio_acts, BIOCHEM_TYPES, BAO_PROTEIN)
    cell = collapse_per_molecule(cell_acts, CELL_TYPES, CELLULAR_FORMATS)

    shared = sorted(set(bio) & set(cell))
    pairs = []
    for m in shared:
        b, c = bio[m], cell[m]
        pairs.append({
            "molecule": m,
            "p_biochem": b["pchembl"],
            "p_cell": c["pchembl"],
            "delta": b["pchembl"] - c["pchembl"],
            "n_biochem": b["n"],
            "n_cell": c["n"],
            "spread_biochem": b["spread"],
            "spread_cell": c["spread"],
            "smiles": b["smiles"] or c["smiles"],
        })

    if verbose:
        print(f"  非细胞 (纯蛋白) 有效分子: {len(bio)}")
        print(f"  细胞法      有效分子: {len(cell)}")
        print(f"  两者都有 -> 配对: {len(pairs)}")
    return pairs
