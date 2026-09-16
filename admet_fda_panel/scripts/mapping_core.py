"""Data-driven half of the safety-panel coverage mapping.

For every ADMET-AI endpoint that plausibly corresponds to a secondary-
pharmacology / safety-panel target, pull the positive drugs (p > 0.5) out of
the 30-drug FDA panel actually computed in this session. Nothing here is
hand-written: the example drugs in the mapping table come from the CSV.
"""

from __future__ import annotations
import pandas as pd

OUT = "/home/user/results/admet_fda_panel"

# ADMET-AI endpoint -> (Chinese target label, panel-relevance note)
ENDPOINT_TARGETS = {
    "hERG": "hERG (Kv11.1) 钾通道",
    "CYP1A2_Veith": "CYP1A2",
    "CYP2C9_Veith": "CYP2C9",
    "CYP2C19_Veith": "CYP2C19",
    "CYP2D6_Veith": "CYP2D6",
    "CYP3A4_Veith": "CYP3A4",
    "Pgp_Broccatelli": "P-糖蛋白 (ABCB1/MDR1)",
    "NR-AR": "雄激素受体 (AR)",
    "NR-AR-LBD": "雄激素受体 配体结合域",
    "NR-ER": "雌激素受体 (ER)",
    "NR-ER-LBD": "雌激素受体 配体结合域",
    "NR-Aromatase": "芳香化酶 (CYP19A1)",
    "NR-PPAR-gamma": "PPAR-γ",
    "NR-AhR": "芳香烃受体 (AhR)",
    "SR-MMP": "线粒体膜电位 (MMP)",
    "SR-ARE": "氧化应激 ARE/Nrf2 通路",
    "SR-p53": "p53 DNA 损伤通路",
    "SR-ATAD5": "ATAD5 基因毒性通路",
    "SR-HSE": "热休克反应 HSE",
}


def positives(d, endpoint, thr=0.5, top=3):
    """Return 'Drug 0.99, Drug 0.97' for the top predicted-positive drugs."""
    s = d.set_index("name")[endpoint].sort_values(ascending=False)
    hit = s[s > thr]
    if len(hit) == 0:
        best = s.head(1)
        return (f"无 (最高: {best.index[0]} {best.iloc[0]:.2f})", 0)
    txt = ", ".join(f"{k} {v:.2f}" for k, v in hit.head(top).items())
    if len(hit) > top:
        txt += f" 等 {len(hit)} 个"
    return txt, len(hit)


def endpoint_stats():
    d = pd.read_csv(f"{OUT}/all_properties.csv")
    rows = {}
    for ep, label in ENDPOINT_TARGETS.items():
        txt, n = positives(d, ep)
        rows[ep] = {"target": label, "positives": txt, "n_positive": n,
                    "median": float(d[ep].median()),
                    "max": float(d[ep].max())}
    return rows


def validation_lookup():
    """endpoint -> (recovered, total, verdict) from the demo validation."""
    v = pd.read_csv(f"{OUT}/demo_validation.csv")
    p = v[v.role == "positive control"]
    out = {}
    for ep, g in p.groupby("endpoint"):
        rec, tot = int(g["above_0.5"].sum()), len(g)
        out[ep] = (rec, tot,
                   "已验证" if rec == tot else "部分" if rec else "验证失败")
    return out


if __name__ == "__main__":
    st = endpoint_stats()
    vl = validation_lookup()
    for ep, r in st.items():
        v = vl.get(ep)
        vtxt = f"  [对照 {v[0]}/{v[1]} {v[2]}]" if v else ""
        print(f"{ep:22s} {r['target']:26s} n+={r['n_positive']:2d}  "
              f"{r['positives']}{vtxt}")
