"""
M9 — 模式细胞的 CC50 数据覆盖：记录数 vs 去重化合物数

要检验的两件事
--------------
1. VC-CELL 护照声称各细胞系的 CC50 记录数是"页面构建时统计，非估计"。
   本模块从 ChEMBL_37 独立重数一遍，比对**排序**（绝对值会因策展口径而异，
   排序不会）。

2. 护照用**记录数**衡量覆盖度。原先怀疑这会高估窄化学空间体系（MT-4 的
   细胞毒性几乎只在 HIV 药物实验里测），应改用去重化合物数。

   **实测否定了这个怀疑**：全部 11 个细胞系的记录/化合物重复度都在
   1.05–1.38 之间，记录数与去重化合物数几乎等价（81.6% 的化合物只测过一次）。
   用记录数排序没有问题。

   但低重复度本身有更重要的后果：**几乎没有重复测量可用来估计噪声**，
   只有 6.8% 的化合物有 ≥3 次测量，噪声估计只能建立在这个被选择的子集上
   （被反复测的多半是参比化合物）。si_noise 模块受此限制。

运行:  python3 -m assaysim.cell_coverage
"""

from __future__ import annotations

import collections
import json
import math
from pathlib import Path

CC50_CACHE = Path("/tmp/cc50_rows.json")

# VC-CELL 护照上公布的记录数（vc_cell.html 内嵌的 CELLS 数组）
PASSPORT = {
    "MT-4": 6499, "Vero E6": 4549, "Huh-7": 3610, "HepG2": 2601,
    "CEM": 1939, "MDCK": 1476, "HeLa": 1415, "HEK293T": 1202,
    "A549": 752, "PBMC": 417, "MRC-5": 167,
}

# ChEMBL 的 assay_cell_type 写法与护照名称的对应（大小写/连字符不敏感匹配）
ALIASES = {
    "MT-4": ["MT-4", "MT4"],
    "Vero E6": ["Vero E6", "VeroE6", "Vero-E6"],
    "Huh-7": ["Huh-7", "Huh7", "HuH-7"],
    "HepG2": ["HepG2", "Hep G2", "HEPG2"],
    "CEM": ["CEM", "CCRF-CEM", "CEM-SS"],
    "MDCK": ["MDCK"],
    "HeLa": ["HeLa"],
    "HEK293T": ["HEK293T", "293T", "HEK293"],
    "A549": ["A549"],
    "PBMC": ["PBMC", "PBMCs"],
    "MRC-5": ["MRC-5", "MRC5"],
}


# --- 护照缺的那一列：动力学参数可得性 --------------------------------------
#
# Tier-C 数字孪生（VC-VK / assaysim M4）唯一的卡点不是 CC50，也不是组学，
# 而是 **每个 病毒×细胞系 组合的 β/k/δ/p/c 是否有人拟合并发表过**。
# 这些参数不可跨组合迁移（本仓库已量化：迁移代价 +0.84 log 系统性偏差）。
#
# 下表是 Europe PMC 摘要级检索的命中数（2026-09-04 执行，查询串附后）。
# ⚠️ 摘要级检索会漏掉正文里报了生长曲线但摘要没提建模的论文。可以断言的是
#    "没有可直接引用的已发表参数集"，不能断言"没有可用于拟合的数据"。
KINETIC_LIT = {
    "流感 (influenza)": {
        "hits": 14, "query": 'ABSTRACT:"target cell limited" AND ABSTRACT:"influenza"',
        "note": "Baccam 2006 等，参数可直接引用（注意其方程不含吸附项）"},
    "SARS-CoV-2": {
        "hits": 10, "query": 'ABSTRACT:"viral kinetic" AND ABSTRACT:"SARS-CoV-2"',
        "note": "多株系拟合已发表"},
    "HIV": {
        "hits": 9, "query": 'ABSTRACT:"target cell limited" AND ABSTRACT:"HIV"',
        "note": "该模型族的起源体系"},
    "RSV": {
        "hits": 3, "query": 'ABSTRACT:"viral kinetic" AND ABSTRACT:"respiratory syncytial"',
        "note": "偏体内"},
    "ASFV": {
        "hits": 0, "query": '"African swine fever" + eclipse phase / burst size / '
                            'mathematical model+replication（三项均 0）',
        "note": "有生长曲线的基因缺失株论文 7 篇，但无人发表拟合参数 —— 需自行标定"},
    "PRRSV": {
        "hits": 0, "query": 'ABSTRACT:"viral kinetic" AND ABSTRACT:"PRRSV"',
        "note": "同上"},
}


def norm(s: str) -> str:
    return s.strip().upper().replace("-", "").replace(" ", "")


def load_rows() -> list[dict]:
    if not CC50_CACHE.exists():
        raise FileNotFoundError(f"缺少 {CC50_CACHE}")
    return json.loads(CC50_CACHE.read_text())


def tally(rows: list[dict]) -> dict[str, dict]:
    """按护照里的细胞系统计记录数与去重化合物数。"""
    alias_map = {}
    for canon, alts in ALIASES.items():
        for a in alts:
            alias_map[norm(a)] = canon

    rec = collections.Counter()
    mols: dict[str, set] = collections.defaultdict(set)
    for a in rows:
        ct = a.get("assay_cell_type")
        if not ct:
            continue
        canon = alias_map.get(norm(ct))
        if canon is None:
            continue
        rec[canon] += 1
        m = a.get("molecule_chembl_id")
        if m:
            mols[canon].add(m)

    return {k: {"records": rec[k], "compounds": len(mols[k]),
                "reps": rec[k] / max(len(mols[k]), 1)}
            for k in PASSPORT}


def spearman(a: list[float], b: list[float]) -> float:
    def rank(x):
        order = sorted(range(len(x)), key=lambda i: x[i])
        r = [0.0] * len(x)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and x[order[j + 1]] == x[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r

    ra, rb = rank(a), rank(b)
    n = len(a)
    ma, mb = sum(ra) / n, sum(rb) / n
    num = sum((ra[i] - ma) * (rb[i] - mb) for i in range(n))
    da = math.sqrt(sum((r - ma) ** 2 for r in ra))
    db = math.sqrt(sum((r - mb) ** 2 for r in rb))
    return num / (da * db) if da and db else float("nan")


def print_kinetic_gap() -> None:
    print("=" * 84)
    print("护照缺的那一列：动力学参数可得性（Europe PMC 摘要级检索，2026-09-04）")
    print("=" * 84)
    print(f"\n{'病毒':<20} {'建模论文':>8}  说明")
    for k, v in KINETIC_LIT.items():
        print(f"{k:<20} {v['hits']:>8}  {v['note']}")
    print("""
CC50 覆盖和组学覆盖回答的是"这个细胞系能不能算毒性"；
动力学参数可得性回答的是"这个 病毒×细胞系 组合能不能立起一个数字孪生"。
后者才是 Tier-C 模型的真实卡点，而它对兽医病毒是**零**。

⚠️ 摘要级检索的边界：能断言"无可直接引用的已发表参数集"，
   不能断言"无可用于拟合的数据"——那些基因缺失株论文里的一步生长曲线就能拟合。""")


def main() -> None:
    try:
        rows = load_rows()
    except FileNotFoundError as e:
        print(f"⚠️ {e}\n   跳过 CC50 覆盖度重数，仅输出动力学参数可得性。\n")
        print_kinetic_gap()
        return
    t = tally(rows)
    print("=" * 84)
    print("模式细胞 CC50 覆盖度：独立重数 vs VC-CELL 护照")
    print(f"ChEMBL_37 全部 CC50 记录 {len(rows):,} 条")
    print("=" * 84)
    print(f"\n{'细胞系':<10} {'护照记录数':>10} {'本次记录数':>10} {'去重化合物':>10} "
          f"{'重复度':>7}  {'覆盖排名变化':>12}")

    by_rec = sorted(PASSPORT, key=lambda k: -t[k]["records"])
    by_cmp = sorted(PASSPORT, key=lambda k: -t[k]["compounds"])
    for k in sorted(PASSPORT, key=lambda x: -PASSPORT[x]):
        d = t[k]
        r1, r2 = by_rec.index(k) + 1, by_cmp.index(k) + 1
        shift = "—" if r1 == r2 else f"{r1} → {r2}"
        print(f"{k:<10} {PASSPORT[k]:10,d} {d['records']:10,d} {d['compounds']:10,d} "
              f"{d['reps']:7.2f}  {shift:>12}")

    ks = list(PASSPORT)
    sp_rec = spearman([PASSPORT[k] for k in ks], [t[k]["records"] for k in ks])
    sp_cmp = spearman([t[k]["records"] for k in ks], [t[k]["compounds"] for k in ks])
    print(f"\n护照记录数 vs 本次记录数，Spearman = {sp_rec:+.3f}")
    print(f"记录数 vs 去重化合物数，Spearman = {sp_cmp:+.3f}")
    print("""
⚠️ 绝对数不该完全一致：护照用的是本套件策展子集（按 assay 描述关键词筛过），
   本次是 ChEMBL 全部 standard_type=CC50 记录按 assay_cell_type 归类。
   可比的是**排序**。

覆盖度该用哪个数：记录/化合物的重复度差异越大，用记录数排序就越容易高估
那些"测得多但化学空间窄"的体系。""")

    print_kinetic_gap()


if __name__ == "__main__":
    main()
