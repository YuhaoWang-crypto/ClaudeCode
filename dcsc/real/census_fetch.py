"""
Fetch a real CELLxGENE Census subset for DCSC tests.

8 blood datasets that each contain BOTH healthy donors and COVID-19 or lupus
(SLE) donors, spanning 5 assays (10x 3' v2/v3, 10x 5' v1/v2, Seq-Well), so
that disease effects can be tested inside a dataset and technical (dataset /
assay) effects across datasets. Up to 25 donors per disease per dataset,
400 random cells per donor, restricted to cell types that map onto the
2-level identity tree below.

    python3 -m dcsc.real.census_fetch      # -> data_cache/census_pbmc.h5ad
"""
import os

import numpy as np
import pandas as pd

CENSUS_VERSION = "2025-01-30"
DATASETS = {
    "218acb0f-9f2f-4f76-b90b-15a4b7c7f629": "SLE_1.2M (10x 3' v2)",
    "ebc2e1ff-c8f9-466a-acf4-9d291afaf8b3": "COMBAT (10x 5' v1)",
    "c7775e88-49bf-4ba2-a03b-93f00447c958": "Stephenson (10x 3')",
    "01ad3cd7-3929-4654-84c0-6db05bd5fd59": "IFN-autoAb (10x 5')",
    "5e717147-0f75-4de1-8bd2-6fda01b8d75f": "Schulte-Schrepping (10x 3' v3)",
    "de2c780c-1747-40bd-9ccf-9588ec186cee": "Lee COVID/flu (10x 3' v3)",
    "59b69042-47c2-47fd-ad03-d21beb99818f": "Arunachalam (10x 3' v3)",
    "456e8b9b-f872-488b-871d-94534090a865": "Wilk (Seq-Well)",
}
DISEASES = ["normal", "COVID-19", "systemic lupus erythematosus"]

# Census cell_type -> (lineage, leaf). Ambiguous labels (generic "T cell",
# gamma-delta, MAIT, NKT, progenitors, platelets...) are dropped.
TREE = {
    "Myeloid": ["CD14_Mono", "CD16_Mono", "cDC", "pDC"],
    "T": ["CD4_T", "CD8_T"],
    "NK": ["NK"],
    "B": ["B", "Plasma"],
}


def map_cell_type(ct):
    c = ct.lower()
    if "plasmacytoid" in c:
        return "pDC"
    if "dendritic" in c:
        return "cDC"
    if "monocyte" in c:
        if "non-classical" in c or "cd16-positive monocyte" in c or "cd14-low" in c:
            return "CD16_Mono"
        if c in ("classical monocyte", "cd14-positive monocyte", "monocyte",
                 "cd14-positive, cd16-negative classical monocyte"):
            return "CD14_Mono"
        return None
    if "natural killer" in c and "t cell" not in c:
        return "NK"
    if "gamma-delta" in c or "mucosal invariant" in c or "nk t" in c:
        return None
    if "cd4-positive" in c or c == "regulatory t cell":
        return "CD4_T"
    if "cd8-positive" in c:
        return "CD8_T"
    if "plasma" in c:
        return "Plasma"
    if "b cell" in c and "precursor" not in c and "pro-b" not in c and "pre-b" not in c:
        return "B"
    return None


LEAF2LIN = {leaf: lin for lin, leaves in TREE.items() for leaf in leaves}


def main(out="data_cache/census_pbmc.h5ad", per_donor=400, max_donors=25, seed=0):
    import cellxgene_census as cc
    rng = np.random.default_rng(seed)
    obs = pd.read_parquet("data_cache/census_blood_obs.parquet")
    obs = obs[obs.dataset_id.isin(DATASETS) & obs.disease.isin(DISEASES)].copy()
    obs["leaf"] = obs.cell_type.astype(str).map(map_cell_type)
    obs = obs[obs.leaf.notna()]
    pick = []
    for (ds, dis), g in obs.groupby(["dataset_id", "disease"], observed=True):
        donors = g.donor_id.unique()
        donors = rng.choice(donors, min(max_donors, len(donors)), replace=False)
        for d in donors:
            ids = g.soma_joinid[g.donor_id == d].values
            pick.append(rng.choice(ids, min(per_donor, len(ids)), replace=False))
    ids = np.sort(np.concatenate(pick))
    print(f"fetching {len(ids):,} cells")
    with cc.open_soma(census_version=CENSUS_VERSION) as census:
        ad = cc.get_anndata(census, "Homo sapiens", obs_coords=ids,
                            obs_column_names=["soma_joinid", "dataset_id", "assay", "cell_type",
                                              "disease", "donor_id", "suspension_type"])
    ad.obs["leaf"] = ad.obs.cell_type.astype(str).map(map_cell_type)
    ad.obs["lineage"] = ad.obs.leaf.map(LEAF2LIN)
    ad.obs["dataset"] = ad.obs.dataset_id.map(DATASETS)
    ad.var_names = ad.var.feature_name.astype(str).values
    ad.var_names_make_unique()
    keep = np.asarray((ad.X > 0).sum(0)).ravel() >= 50
    ad = ad[:, keep].copy()
    ad.write_h5ad(out, compression="gzip")
    print(ad)


if __name__ == "__main__":
    os.makedirs("data_cache", exist_ok=True)
    main()
