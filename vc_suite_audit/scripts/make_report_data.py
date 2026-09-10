"""Consolidate every measured result into one JSON blob for the report page."""
import glob, json, os
import pandas as pd
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
R = HERE + "/results"

splits = pd.concat([pd.read_csv(f) for f in glob.glob(R + "/splits_*.csv")], ignore_index=True)
splits = splits.drop_duplicates(subset=["dataset", "split", "qrasar"], keep="last")

ad = pd.concat([pd.read_csv(f) for f in glob.glob(R + "/ad_*.csv")], ignore_index=True)
ad = ad.drop_duplicates(subset=["dataset", "lo"], keep="last")

util = pd.read_csv(R + "/utility_table.csv")
tr = pd.read_csv(R + "/transfer.csv") if os.path.exists(R + "/transfer.csv") else pd.DataFrame()

# dataset provenance: mine vs the report's stated counts
PROV = {
    "NA_influenza":     dict(label="NA 流感神经氨酸酶 pIC50", page="VC-AVI",
                             rep_n=1143, rep_scaf=318, rep_rho=0.681, rep_rho_qr=0.809),
    "PA_influenza":     dict(label="PA 核酸内切酶 pIC50", page="VC-AVI",
                             rep_n=252, rep_scaf=126, rep_rho=0.676, rep_rho_qr=0.734),
    "CC50_pan":         dict(label="泛细胞毒性 pCC50", page="VC-AVI / VC-TRIAGE",
                             rep_n=33641, rep_scaf=13678, rep_rho=0.535, rep_rho_qr=0.704),
    "MIC_H_influenzae": dict(label="流感嗜血杆菌 pMIC", page="VC-MIC2",
                             rep_n=3499, rep_scaf=1681, rep_rho=0.730, rep_rho_qr=None),
    "MIC_E_cloacae":    dict(label="阴沟肠杆菌 pMIC", page="VC-MIC2",
                             rep_n=2440, rep_scaf=822, rep_rho=0.785, rep_rho_qr=None),
    "hERG":             dict(label="hERG 抑制 pIC50", page="VC-ADMET",
                             rep_n=9716, rep_scaf=5106, rep_rho=0.606, rep_rho_qr=None),
}

for k, v in PROV.items():
    f = HERE + f"/curated/{k}.csv"
    if os.path.exists(f):
        d = pd.read_csv(f)
        v["my_n"] = int(len(d)); v["my_scaf"] = int(d.scaffold.nunique())

import narrative
import external as ext_mod

CHAINS = [
    dict(name="选择性指数 SI = CC50 / EC50",
         parts="CC50 模型 0.63 ⊕ NA 活性模型 1.10", sigma=1.27, fold="约 18 倍"),
    dict(name="抗菌窗口 = pMIC − pCC50（阴沟肠杆菌）",
         parts="CC50 模型 0.63 ⊕ MIC 模型 0.90", sigma=1.09, fold="约 12 倍"),
    dict(name="抗菌窗口 = pMIC − pCC50（流感嗜血杆菌）",
         parts="CC50 模型 0.63 ⊕ MIC 模型 0.74", sigma=0.97, fold="约 9 倍"),
    dict(name="EC50 → ODE 峰下降 %",
         parts="NA 活性模型 1.10 ⊕ 生化 IC50→细胞 EC50 系统偏移 ≈1.0", sigma=1.49, fold="约 31 倍"),
]

out = dict(
    splits=json.loads(splits.to_json(orient="records")),
    ad=json.loads(ad.to_json(orient="records")),
    utility=json.loads(util.to_json(orient="records")),
    transfer=json.loads(tr.to_json(orient="records")) if len(tr) else [],
    provenance=PROV,
    verdicts=narrative.VERDICTS,
    notes=narrative.NOTES,
    limits=narrative.LIMITS,
    chains=CHAINS,
    external=ext_mod.EXTERNAL,
)
with open(R + "/report_data.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)

# inject into the page template
tpl = open(HERE + "/audit.html", encoding="utf-8").read()
blob = json.dumps(out, ensure_ascii=False, separators=(",", ":"))
assert "/*__DATA__*/" in tpl, "template marker missing"
open(HERE + "/audit_final.html", "w", encoding="utf-8").write(
    tpl.replace("/*__DATA__*/ {}", blob))
print("page written:", HERE + "/audit_final.html",
      f"({len(blob)/1024:.0f} KB data)")
print("datasets:", sorted(splits.dataset.unique()))
print("splits rows:", len(splits), "ad rows:", len(ad), "transfer rows:", len(tr))
print(splits.pivot_table(index="dataset", columns=["split", "qrasar"],
                         values="spearman").round(3).to_string())
