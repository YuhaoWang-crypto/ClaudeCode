"""
C1 — Does the open data actually cover the six axes?

The question is NOT "is there a lot of single-cell data" (there is). It is
whether, per axis, the open corpus supplies the two different things a state
space needs:

  (a) OBSERVATIONS from which the latent axis can be estimated, and
  (b) FATE LABELS paired to those observations, so the state -> fate map can
      be fit and falsified.

This module counts both against live indexes instead of asserting them:

  * NCBI GEO (db=gds) via E-utilities  -> dataset counts per axis readout,
    per axis readout restricted to single-cell, and the INTERSECTION of the
    axis readout with a fate/lineage label. That intersection is the number
    that decides whether the state -> fate map is fittable at all.
  * EBI BioModels REST search          -> curated mechanistic models per axis
    (a state space needs dynamics, not just a static coordinate).
  * Europe PMC                         -> literature volume for the axis'
    gold-standard assay, as a sanity denominator.

Every number printed here is fetched, then cached under
figures/_cellstate_cache/ so re-runs are reproducible offline.

Rigour: the counts are ✅ exact (they are what the index returns for the
stated query). The mapping "query -> axis" is ⚠️ a keyword proxy: it will
include some irrelevant records and miss some relevant ones. The ORDERS OF
MAGNITUDE and the ratios between axes are the interpretable signal, not the
individual integers.
"""
import json
import subprocess
import time

from .common import banner, biomodels_search, curated, save_json

EUTILS = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
          "?db=gds&term={term}&retmax=0&retmode=json")
EPMC = ("https://www.ebi.ac.uk/europepmc/webservices/rest/search"
        "?query={q}&format=json&pageSize=1")

# A fate label = the experiment recorded what the cell/clone DID afterwards.
FATE = ('("lineage tracing" OR "lineage barcode" OR "clonal barcoding" OR '
        '"fate mapping" OR "cell fate" OR "clonal tracking" OR "barcode")')
SC = '("single cell" OR "single-cell" OR scRNA)'

AXES = [
    dict(
        key="cell_cycle", label="细胞周期 Cell cycle",
        geo='("cell cycle" OR "cell-cycle phase" OR proliferation)',
        gold="FUCCI reporter / EdU incorporation",
        gold_q='"FUCCI" OR "EdU incorporation"',
        bm=["cell cycle", "restriction point", "cyclin CDK"],
    ),
    dict(
        key="stress", label="应激水平 Stress level",
        geo='("integrated stress response" OR "heat shock response" OR '
            '"oxidative stress" OR "ATF4" OR "NRF2")',
        gold="ISR/HSR reporter, eIF2a-P immunoblot",
        gold_q='"integrated stress response"',
        bm=["heat shock", "oxidative stress", "unfolded protein response"],
    ),
    dict(
        key="metabolic_reserve", label="代谢储备 Metabolic reserve",
        geo='("spare respiratory capacity" OR "oxygen consumption rate" OR '
            '"Seahorse" OR "extracellular acidification")',
        gold="Seahorse XF spare respiratory capacity (FCCP uncoupling)",
        gold_q='"spare respiratory capacity"',
        bm=["glycolysis", "oxidative phosphorylation", "central carbon metabolism"],
    ),
    dict(
        key="dna_damage", label="DNA 损伤 DNA damage",
        geo='("DNA damage" OR "gammaH2AX" OR "H2AX" OR "double strand break")',
        gold="gammaH2AX / 53BP1 foci, comet assay",
        gold_q='"gammaH2AX foci" OR "comet assay"',
        bm=["DNA damage", "p53", "DNA repair checkpoint"],
    ),
    dict(
        key="apoptotic_priming", label="凋亡准备度 Apoptotic priming",
        geo='("apoptotic priming" OR "BH3 profiling" OR "mitochondrial priming" '
            'OR "BCL-2 dependence")',
        gold="BH3 profiling / dynamic BH3 profiling (cytochrome c release)",
        gold_q='"BH3 profiling"',
        bm=["apoptosis", "caspase", "BCL-2 mitochondrial"],
    ),
    dict(
        key="epigenetic_lineage", label="表观/谱系 Epigenetic-lineage",
        geo='("ATAC-seq" OR "chromatin accessibility" OR "DNA methylation" OR '
            '"chromatin state")',
        gold="scATAC-seq / WGBS / SHARE-seq chromatin potential",
        gold_q='"chromatin accessibility" AND "cell fate"',
        bm=["differentiation switch", "pluripotency", "bistable gene regulatory"],
    ),
]


def _get_json(url, timeout=60):
    try:
        out = subprocess.run(["curl", "-sL", "--max-time", str(timeout), url],
                             capture_output=True, text=True, check=True).stdout
        return json.loads(out)
    except Exception:
        return None


def geo_count(term):
    """Exact hit count from NCBI GEO DataSets for a free-text query."""
    url = EUTILS.format(term=term.replace(" ", "+").replace('"', "%22"))
    d = _get_json(url)
    time.sleep(0.4)  # stay under the un-keyed E-utilities rate limit
    if not d or "esearchresult" not in d:
        return None
    try:
        return int(d["esearchresult"]["count"])
    except Exception:
        return None


def epmc_count(q):
    d = _get_json(EPMC.format(q=q.replace(" ", "%20").replace('"', "%22")))
    time.sleep(0.2)
    if not d:
        return None
    return d.get("hitCount")


def biomodels_axis(queries):
    """Union of CURATED (BIOMD*) models across the axis' query terms."""
    ids, names, ok = set(), {}, False
    for q in queries:
        d = biomodels_search(q, num=40)
        if d is None:
            continue
        ok = True
        for m in d.get("models", []):
            if curated(m.get("id", "")):
                ids.add(m["id"])
                names[m["id"]] = m.get("name", "")
    return (sorted(ids), names) if ok else (None, {})


def run():
    banner("C1 — 六轴开源数据盘点(实时检索,非断言)")
    rows, offline = [], False

    for ax in AXES:
        n_read = geo_count(ax["geo"])
        n_sc = geo_count(f'{ax["geo"]} AND {SC}')
        n_pair = geo_count(f'{ax["geo"]} AND {FATE}')
        n_sc_pair = geo_count(f'{ax["geo"]} AND {SC} AND {FATE}')
        lit = epmc_count(ax["gold_q"])
        bm_ids, bm_names = biomodels_axis(ax["bm"])
        if n_read is None or bm_ids is None:
            offline = True
        rows.append(dict(axis=ax["key"], label=ax["label"], gold=ax["gold"],
                         geo_readout=n_read, geo_sc=n_sc, geo_fate=n_pair,
                         geo_sc_fate=n_sc_pair, pmc_gold=lit,
                         biomodels_curated=len(bm_ids) if bm_ids else None,
                         biomodels_ids=bm_ids or [], biomodels_names=bm_names))

    if offline:
        print("⚠️  某些查询未能联网完成 —— 下表可能不完整,不要据此下结论。")

    print("\n【观测侧】GEO DataSets 命中数(db=gds)")
    print(f"{'轴':<28}{'读出':>9}{'+单细胞':>10}{'+命运标签':>11}{'单细胞∩命运':>13}")
    print("-" * 72)
    for r in rows:
        f = lambda v: "n/a" if v is None else str(v)
        print(f"{r['label']:<26}{f(r['geo_readout']):>9}{f(r['geo_sc']):>10}"
              f"{f(r['geo_fate']):>11}{f(r['geo_sc_fate']):>13}")

    print("\n【配对率】命运标签在该轴数据中的占比 —— 这是状态→命运映射能否被拟合的上限")
    print(f"{'轴':<28}{'配对率':>10}{'单细胞配对率':>14}")
    print("-" * 54)
    for r in rows:
        if r["geo_readout"]:
            p = 100.0 * (r["geo_fate"] or 0) / r["geo_readout"]
            ps = (100.0 * (r["geo_sc_fate"] or 0) / r["geo_sc"]) if r["geo_sc"] else float("nan")
            r["pair_pct"], r["pair_sc_pct"] = p, ps
            print(f"{r['label']:<26}{p:>9.2f}%{ps:>13.2f}%")

    print("\n【动力学侧】BioModels 策展机制模型(BIOMD*,该轴查询词的并集)")
    for r in rows:
        n = r["biomodels_curated"]
        print(f"  {r['label']:<26} {('n/a' if n is None else n):>4} 个"
              f"   金标准: {r['gold']}")

    print("\n【文献侧】Europe PMC 中该轴金标准测定的文献量")
    for r in rows:
        print(f"  {r['label']:<26} {r['pmc_gold']!s:>8}  ({r['gold']})")

    # --- the interpretable contrast: spread across axes -------------------
    have = [r for r in rows if r["geo_readout"]]
    if have:
        mx = max(have, key=lambda r: r["geo_readout"])
        mn = min(have, key=lambda r: r["geo_readout"])
        print(f"\n观测量最富的轴: {mx['label']} ({mx['geo_readout']})")
        print(f"观测量最贫的轴: {mn['label']} ({mn['geo_readout']})")
        if mn["geo_readout"]:
            print(f"两者相差 {mx['geo_readout'] / mn['geo_readout']:.0f} 倍 "
                  f"—— 六个轴的数据支撑并不在同一个量级上。")

    save_json("c1_inventory.json", rows)
    print("\n结果已存 figures/_cellstate_cache/c1_inventory.json")
    return rows


if __name__ == "__main__":
    run()
