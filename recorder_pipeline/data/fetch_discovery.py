"""Reverse pipeline: indication -> cells -> donor pseudobulk -> ratio nomination.

Instead of validating a given marker list, start from the tissue and let the
data nominate singles AND pairs, scored by the r2_pairing criterion.

SLE: 99 healthy + 162 SLE donors, one dataset (Perez et al.), classical
monocytes and B cells. Candidate panel = HPA "Secreted to blood" (measurability
gate applied BEFORE looking at effect sizes) plus the report's own complement
candidates for comparison.
"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, cellxgene_census

OUT = "/home/user/ClaudeCode/recorder_pipeline/data"
DS = "218acb0f-9f2f-4f76-b90b-15a4b7c7f629"

sec = pd.read_csv("secreted_blood.tsv", sep="\t")
PANEL = sorted(set(sec.Gene.dropna().tolist()) | {
    # the report's own SLE-relevant candidates, for head-to-head comparison
    "C3", "C4A", "C4B", "CR1", "C1QA", "C1QB", "CFB", "CFH", "SERPING1",
    # interferon axis (known SLE biology, acts as a positive control)
    "IFI44L", "ISG15", "IFI6", "MX1", "OAS1", "SIGLEC1", "USP18",
})
print("panel size:", len(PANEL))

CELLTYPES = ["classical monocyte", "B cell"]
with cellxgene_census.open_soma(census_version="2025-11-08") as census:
    a = cellxgene_census.get_anndata(
        census=census, organism="Homo sapiens", measurement_name="RNA",
        X_name="normalized",
        obs_value_filter=(f"tissue_general=='blood' and is_primary_data==True and "
                          f"dataset_id=='{DS}' and cell_type in {CELLTYPES}"),
        var_value_filter=f"feature_name in {sorted(PANEL)}",
        obs_column_names=["disease", "cell_type", "donor_id", "sex"])

X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
genes = list(a.var["feature_name"])
df = pd.DataFrame(X, columns=genes)
for c in ["disease", "cell_type", "donor_id"]:
    df[c] = a.obs[c].values
print("cells:", len(df), "genes recovered:", len(genes))

for ct in CELLTYPES:
    sub = df[df.cell_type == ct]
    if len(sub) < 5000:
        continue
    # donor pseudobulk: mean normalized expression per donor, then log
    pb = sub.groupby(["donor_id", "disease"], observed=True)[genes].mean().reset_index()
    pb = pb[pb.groupby("donor_id")["donor_id"].transform("size") == 1]  # 1 label/donor
    n_ok = pb.groupby("disease").size()
    print(f"\n{ct}: donors {dict(n_ok)}")
    L = np.log(pb[genes].values + 1e-6)
    lab = (pb.disease.str.contains("lupus")).values
    keep = (np.isfinite(L).all(0)) & (L.std(0) > 1e-9)
    L, gg = L[:, keep], [g for g, k in zip(genes, keep) if k]
    out = pd.DataFrame({"gene": gg,
                        "mean_ctrl": L[~lab].mean(0), "mean_dis": L[lab].mean(0),
                        "sd_ctrl": L[~lab].std(0, ddof=1), "sd_dis": L[lab].std(0, ddof=1)})
    out["delta"] = out.mean_dis - out.mean_ctrl
    out["sd_pool"] = np.sqrt((out.sd_ctrl**2 + out.sd_dis**2) / 2)
    out["d"] = out.delta / out.sd_pool
    out["cell_type"] = ct
    out.to_csv(f"{OUT}/discover_singles_{ct.replace(' ','_')}.tsv", sep="\t", index=False)
    np.save(f"/tmp/L_{ct.replace(' ','_')}.npy", L)
    np.save(f"/tmp/lab_{ct.replace(' ','_')}.npy", lab)
    pd.Series(gg).to_csv(f"/tmp/gg_{ct.replace(' ','_')}.csv", index=False, header=False)
    print(out.reindex(out.d.abs().sort_values(ascending=False).index).head(8)[
        ["gene", "delta", "d"]].to_string(index=False))
print("\nsaved")
