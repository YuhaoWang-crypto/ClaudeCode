"""Derive every number and table used in the report from the result CSVs.

Nothing here is hard-coded: each table is computed from
/home/user/results/admet_fda_panel/*.csv so the report and the CSV exports
cannot drift apart.
"""

from __future__ import annotations
import pandas as pd
import numpy as np

OUT = "/home/user/results/admet_fda_panel"


def load():
    return pd.read_csv(f"{OUT}/all_properties.csv")


# --------------------------------------------------------------------------- #
# Table 1 — physicochemical / drug-likeness summary
# --------------------------------------------------------------------------- #
PHYS_ROWS = [
    ("MW", "分子量 (Da)", "≤ 500", 1),
    ("cLogP", "cLogP", "≤ 5", 2),
    ("TPSA", "拓扑极性表面积 (Å²)", "≤ 140", 1),
    ("HBD", "氢键供体", "≤ 5", 1),
    ("HBA", "氢键受体", "≤ 10", 1),
    ("RotB", "可旋转键", "≤ 10", 1),
    ("FractionCSP3", "sp³ 碳比例", "—", 2),
    ("QED", "QED 类药性", "> 0.5 为佳", 2),
]


def table_physchem(d):
    rows = []
    for col, label, rule, nd in PHYS_ROWS:
        v = d[col]
        rows.append([label,
                     f"{v.median():.{nd}f}",
                     f"{v.mean():.{nd}f} ± {v.std():.{nd}f}",
                     f"{v.min():.{nd}f} – {v.max():.{nd}f}",
                     rule])
    return ["性质", "中位数", "均值 ± SD", "范围", "参考阈值"], rows


def table_rules(d):
    n = len(d)
    rows = []
    for col, label, note in [
        ("Lipinski_pass", "Lipinski 五规则 (≤1 项违反)", "口服吸收"),
        ("Veber_pass", "Veber 规则 (RotB≤10, TPSA≤140)", "口服生物利用度"),
        ("Ghose_pass", "Ghose 过滤", "类药化学空间"),
        ("Egan_pass", "Egan 卵 (TPSA/logP)", "被动跨膜吸收"),
    ]:
        k = int(d[col].sum())
        fails = ", ".join(d.loc[~d[col], "name"].tolist()[:6]) or "—"
        rows.append([label, f"{k}/{n}", f"{100*k/n:.0f}%", fails])
    return ["规则", "通过数", "通过率", "未通过化合物"], rows


# --------------------------------------------------------------------------- #
# Table 2 — key ADMET endpoint prevalence
# --------------------------------------------------------------------------- #
KEY_ENDPOINTS = [
    ("hERG", "hERG 钾通道阻断", "心脏 QT 延长风险", True),
    ("DILI", "药物性肝损伤 (DILI)", "肝毒性", True),
    ("AMES", "Ames 致突变", "遗传毒性", True),
    ("ClinTox", "临床试验毒性失败", "综合临床毒性", True),
    ("Carcinogens_Lagunin", "致癌性", "长期毒性", True),
    ("CYP1A2_Veith", "CYP1A2 抑制", "药物相互作用", True),
    ("CYP2C9_Veith", "CYP2C9 抑制", "药物相互作用", True),
    ("CYP2C19_Veith", "CYP2C19 抑制", "药物相互作用", True),
    ("CYP2D6_Veith", "CYP2D6 抑制", "药物相互作用", True),
    ("CYP3A4_Veith", "CYP3A4 抑制", "药物相互作用", True),
    ("Pgp_Broccatelli", "P-糖蛋白抑制", "外排转运/DDI", True),
    ("BBB_Martins", "血脑屏障穿透", "中枢暴露", True),
    ("HIA_Hou", "人肠道吸收", "口服吸收", True),
    ("Bioavailability_Ma", "口服生物利用度", "口服吸收", True),
]


def table_endpoints(d):
    n = len(d)
    rows = []
    for col, label, cat, _ in KEY_ENDPOINTS:
        v = d[col]
        k = int((v > 0.5).sum())
        top = d.loc[v.idxmax(), "name"]
        rows.append([label, cat, f"{v.median():.2f}",
                     f"{k}/{n} ({100*k/n:.0f}%)",
                     f"{top} ({v.max():.2f})"])
    return ["ADMET 终点", "关注类别", "中位预测值", "> 0.5 的化合物", "最高预测化合物"], rows


# --------------------------------------------------------------------------- #
# Table 3 — regression endpoints
# --------------------------------------------------------------------------- #
REG_ENDPOINTS = [
    ("Solubility_AqSolDB", "水溶解度 log S (mol/L)", 2),
    ("Lipophilicity_AstraZeneca", "亲脂性 logD7.4", 2),
    ("Caco2_Wang", "Caco-2 渗透性 log cm/s", 2),
    ("PPBR_AZ", "血浆蛋白结合率 (%)", 1),
    ("VDss_Lombardo", "稳态分布容积 log L/kg", 2),
    ("Half_Life_Obach", "半衰期 (h)", 1),
    ("Clearance_Hepatocyte_AZ", "肝细胞清除率 (µL/min/10⁶)", 1),
    ("Clearance_Microsome_AZ", "微粒体清除率 (mL/min/g)", 2),
    ("LD50_Zhu", "急性毒性 LD50 (-log mol/kg)", 2),
]


def table_regression(d):
    rows = []
    for col, label, nd in REG_ENDPOINTS:
        v = d[col]
        rows.append([label, f"{v.median():.{nd}f}",
                     f"{v.quantile(.25):.{nd}f} – {v.quantile(.75):.{nd}f}",
                     f"{d.loc[v.idxmin(),'name']} ({v.min():.{nd}f})",
                     f"{d.loc[v.idxmax(),'name']} ({v.max():.{nd}f})"])
    return ["回归终点", "中位数", "四分位距", "最低", "最高"], rows


# --------------------------------------------------------------------------- #
# flagged compounds
# --------------------------------------------------------------------------- #
def table_flagged(d, top=12):
    risk = ["hERG", "DILI", "AMES", "ClinTox", "Carcinogens_Lagunin"]
    sc = d[risk].gt(0.5).sum(axis=1)
    t = d.assign(n_risk=sc).sort_values(
        ["n_risk", "hERG"], ascending=False).head(top)
    rows = []
    for _, r in t.iterrows():
        flags = [lab for c, lab in
                 [("hERG", "hERG"), ("DILI", "DILI"), ("AMES", "Ames"),
                  ("ClinTox", "ClinTox"), ("Carcinogens_Lagunin", "致癌")]
                 if r[c] > 0.5]
        rows.append([r["name"], r["drug_class"],
                     f"{r['hERG']:.2f}", f"{r['DILI']:.2f}",
                     f"{r['AMES']:.2f}", int(r["alert_total"]),
                     ", ".join(flags) or "—"])
    return ["药物", "药理类别", "hERG", "DILI", "Ames", "结构警示", "触发标记"], rows


# --------------------------------------------------------------------------- #
# percentile context
# --------------------------------------------------------------------------- #
PCT_ENDPOINTS = ["hERG", "DILI", "AMES", "CYP3A4_Veith", "BBB_Martins",
                 "Solubility_AqSolDB"]
PCT_LABELS = {"hERG": "hERG 阻断", "DILI": "肝损伤", "AMES": "Ames 致突变",
              "CYP3A4_Veith": "CYP3A4 抑制", "BBB_Martins": "血脑屏障穿透",
              "Solubility_AqSolDB": "水溶解度"}


def table_percentiles(d, drugs=None):
    drugs = drugs or ["Ketoconazole", "Amiodarone", "Terfenadine", "Cisapride",
                      "Imatinib", "Atorvastatin", "Metformin", "Aspirin",
                      "Caffeine", "Ciprofloxacin"]
    hdr = ["药物"] + [PCT_LABELS[e] for e in PCT_ENDPOINTS]
    rows = []
    idx = d.set_index("name")
    for nm in drugs:
        r = idx.loc[nm]
        cells = []
        for e in PCT_ENDPOINTS:
            pc = r.get(f"{e}_drugbank_approved_percentile", np.nan)
            cells.append(f"{r[e]:.2f} ({pc:.0f}%)" if pd.notna(pc)
                         else f"{r[e]:.2f}")
        rows.append([nm] + cells)
    return hdr, rows


# --------------------------------------------------------------------------- #
# demo validation table
# --------------------------------------------------------------------------- #
VAL_LABELS = {
    "hERG": "hERG 阻断 / TdP 风险",
    "CYP3A4_Veith": "CYP3A4 抑制",
    "CYP2D6_Veith": "CYP2D6 抑制",
    "CYP1A2_Veith": "CYP1A2 抑制",
    "Pgp_Broccatelli": "P-糖蛋白抑制",
    "BBB_Martins": "血脑屏障穿透",
    "NR-Aromatase": "芳香化酶 (CYP19A1)",
    "NR-ER": "雌激素受体",
    "NR-PPAR-gamma": "PPAR-γ",
    "NR-AR": "雄激素受体",
}


def table_validation():
    v = pd.read_csv(f"{OUT}/demo_validation.csv")
    rows = []
    for ep, lab in VAL_LABELS.items():
        sub = v[v.endpoint == ep]
        pos = sub[sub.role == "positive control"]
        neg = sub[sub.role == "negative control"]
        rec = int(pos["above_0.5"].sum())
        verdict = ("全部命中" if rec == len(pos) else
                   f"部分命中" if rec else "未命中")
        pos_s = ", ".join(f"{r.drug} {r.score:.2f}" for r in pos.itertuples())
        rows.append([lab, pos_s, f"{neg.score.mean():.2f}",
                     f"{rec}/{len(pos)}", verdict])
    return ["终点 (阳性对照药物)", "阳性对照预测值", "阴性对照均值",
            "命中", "结论"], rows


def validation_stats():
    v = pd.read_csv(f"{OUT}/demo_validation.csv")
    dmpk = ["hERG", "CYP3A4_Veith", "CYP2D6_Veith", "Pgp_Broccatelli",
            "BBB_Martins", "CYP1A2_Veith"]
    nr = ["NR-ER", "NR-AR", "NR-Aromatase", "NR-PPAR-gamma"]
    p = v[v.role == "positive control"]
    return {
        "dmpk_rec": int(p[p.endpoint.isin(dmpk)]["above_0.5"].sum()),
        "dmpk_tot": int((p.endpoint.isin(dmpk)).sum()),
        "nr_rec": int(p[p.endpoint.isin(nr)]["above_0.5"].sum()),
        "nr_tot": int((p.endpoint.isin(nr)).sum()),
    }


if __name__ == "__main__":
    d = load()
    for fn in (table_physchem, table_rules, table_endpoints, table_regression,
               table_flagged, table_percentiles):
        hdr, rows = fn(d)
        print("\n" + " | ".join(hdr))
        for r in rows[:4]:
            print("  " + " | ".join(str(x) for x in r))
    print("\n", table_validation()[1])
    print(validation_stats())
