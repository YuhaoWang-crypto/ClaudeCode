"""
Modal app: Geneformer in-silico perturbation on the IPF atlas (Step 1-3).

Staged and verified one function at a time:
  inspect_census  -> discover IPF lung datasets + cell_type labels (pick data)
  fetch_ipf_data  -> download chosen dataset (raw counts) to a Volume as h5ad
  tokenize        -> Geneformer rank-value encoding (+ cell_state labels)
  perturb         -> InSilicoPerturber goal_state_shift -> ranked genes

The heavy CELLxGENE stack lives in its own image so it can't conflict with the
pinned geneformer/transformers env used for tokenize/perturb.
"""
import modal

app = modal.App("gf-ipf-perturb")

# --- data image: CELLxGENE Census + scanpy (no geneformer) -------------------
data_image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install("cellxgene-census", "scanpy", "anndata", "numpy", "pandas")
)

# --- geneformer image: identical to the verified smoke image (cache-shared) --
gf_image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("git", "git-lfs")
    .pip_install(
        "torch", "transformers==4.46.3", "datasets", "numpy", "pandas", "scipy",
        "scikit-learn", "loompy", "anndata", "scanpy", "pyarrow",
        "tdigest", "peft", "optuna", "ray",
    )
    .run_commands(
        "GIT_LFS_SKIP_SMUDGE=1 git clone "
        "https://huggingface.co/ctheodoris/Geneformer.git /opt/Geneformer",
        "cd /opt/Geneformer && git lfs pull "
        '--include="geneformer/*.pkl"',
        "pip install /opt/Geneformer",
    )
)

model_vol = modal.Volume.from_name("gf-models", create_if_missing=True)
data_vol = modal.Volume.from_name("gf-ipf-data", create_if_missing=True)
MODEL_DIR, DATA_DIR = "/models", "/data"
CKPT = "Geneformer-V2-104M"


@app.function(image=data_image, volumes={DATA_DIR: data_vol}, timeout=1800)
def inspect_census():
    """Discover IPF lung datasets and their cell_type labels in Census."""
    import cellxgene_census

    with cellxgene_census.open_soma(census_version="stable") as census:
        obs = cellxgene_census.get_obs(
            census, "homo_sapiens",
            value_filter="tissue_general == 'lung'",
            column_names=["dataset_id", "cell_type", "disease"],
        )
    print(f"lung cells in Census: {len(obs):,}")
    # observed=True so we only see labels actually present (not empty categories)
    dz = obs["disease"].value_counts()
    dz = dz[dz > 0]
    print("\n== disease labels present in lung ==")
    print(dz.head(30).to_string())
    # find the fibrosis-related label(s)
    fib = [d for d in dz.index if "fibros" in str(d).lower()]
    print("\n== fibrosis-related disease labels ==", fib)
    if fib:
        sub = obs[obs["disease"].isin(fib)]
        print(f"\n== fibrosis cells: {len(sub):,} ==")
        print("top datasets:")
        print(sub["dataset_id"].value_counts().head(6).to_string())
        kw = ["basal", "baso", "fibro", "aberr", "epitheli", "alveolar", "AT2", "AT1"]
        hits = sorted({c for c in sub["cell_type"].unique()
                       if any(k.lower() in str(c).lower() for k in kw)})
        print("\n== relevant cell_type labels in fibrosis cells ==")
        for h in hits:
            n = int((sub["cell_type"] == h).sum())
            print(f"   {h}: {n:,}")
    return {"disease_labels": list(map(str, dz.head(30).index)), "fibrosis": fib}


# Fibrosis lung atlas dataset in Census (from inspect_census).
PF_DATASET = "9f222629-9e39-47d0-b83f-e08d610c7479"
# Epithelial + fibroblast lineage = the whitespace compartments (basaloid proxy
# = respiratory basal cell) plus a healthy epithelial/fibroblast reference.
LINEAGE = [
    "respiratory basal cell",
    "pulmonary alveolar type 2 cell",
    "pulmonary alveolar type 1 cell",
    "myofibroblast cell",
    "alveolar adventitial fibroblast",
    "alveolar type 1 fibroblast cell",
    "bronchus fibroblast of lung",
    "epithelial cell of lower respiratory tract",
]
IPF_H5AD = f"{DATA_DIR}/ipf_lung_raw.h5ad"


@app.function(image=data_image, volumes={DATA_DIR: data_vol}, timeout=2400)
def fetch_ipf_data(max_per_state: int = 4000):
    """Fetch raw-count fibrosis + normal cells (epithelial/fibroblast lineage)."""
    import cellxgene_census
    import numpy as np

    lineage = "[" + ", ".join(f"'{c}'" for c in LINEAGE) + "]"
    vf = (
        f"dataset_id == '{PF_DATASET}' and "
        f"disease in ['pulmonary fibrosis', 'normal'] and "
        f"cell_type in {lineage}"
    )
    with cellxgene_census.open_soma(census_version="stable") as census:
        adata = cellxgene_census.get_anndata(
            census, organism="Homo sapiens", X_name="raw",
            obs_value_filter=vf,
            column_names={"obs": ["disease", "cell_type", "dataset_id"],
                          "var": ["feature_id", "feature_name"]},
        )
    print(f"fetched {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print("disease:\n", adata.obs["disease"].value_counts().to_string())

    # If this dataset lacks normals, pull a normal reference from the same lineage.
    if (adata.obs["disease"] == "normal").sum() < 200:
        print("[info] few normals in dataset; pulling normal reference broadly")
        vf_norm = (
            f"tissue_general == 'lung' and disease == 'normal' and "
            f"cell_type in {lineage}"
        )
        with cellxgene_census.open_soma(census_version="stable") as census:
            norm = cellxgene_census.get_anndata(
                census, organism="Homo sapiens", X_name="raw",
                obs_value_filter=vf_norm,
                column_names={"obs": ["disease", "cell_type", "dataset_id"],
                              "var": ["feature_id", "feature_name"]},
            )
        import anndata as ad
        # cap normals before concat to keep it bounded
        if norm.n_obs > max_per_state:
            idx = np.linspace(0, norm.n_obs - 1, max_per_state).astype(int)
            norm = norm[idx].copy()
        adata = ad.concat([adata[adata.obs["disease"] == "pulmonary fibrosis"], norm])

    # Balanced subsample per disease state.
    keep = []
    for dz in ["pulmonary fibrosis", "normal"]:
        idx = np.where(adata.obs["disease"].to_numpy() == dz)[0]
        if len(idx) > max_per_state:
            idx = np.random.RandomState(0).choice(idx, max_per_state, replace=False)
        keep.append(idx)
    keep = np.sort(np.concatenate(keep))
    adata = adata[keep].copy()

    # Geneformer needs Ensembl IDs + raw counts + n_counts.
    adata.var["ensembl_id"] = adata.var["feature_id"].values
    adata.obs["n_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()
    adata.obs["disease"] = adata.obs["disease"].astype(str)

    adata.write_h5ad(IPF_H5AD)
    data_vol.commit()
    print(f"\n[done] wrote {IPF_H5AD}")
    print("final disease balance:\n", adata.obs["disease"].value_counts().to_string())
    print("final cell_type:\n", adata.obs["cell_type"].value_counts().to_string())
    return {"n_obs": int(adata.n_obs), "n_vars": int(adata.n_vars)}


DATASET = f"{DATA_DIR}/ipf.dataset"
HVG_JSON = f"{DATA_DIR}/hvg_ensembl.json"
PERT_OUT = f"{DATA_DIR}/perturb"
STATES = {"state_key": "disease", "start_state": "pulmonary fibrosis",
          "goal_state": "normal", "alt_states": []}


def _gc104m_dicts():
    """Absolute paths to the gc104M dictionaries in the installed package."""
    import os
    import geneformer
    d = os.path.dirname(geneformer.__file__)
    return (f"{d}/gene_median_dictionary_gc104M.pkl",
            f"{d}/token_dictionary_gc104M.pkl",
            f"{d}/ensembl_mapping_dict_gc104M.pkl")


@app.function(image=gf_image, volumes={DATA_DIR: data_vol}, timeout=1800)
def tokenize(n_hvg: int = 256):
    """Geneformer rank-value tokenization + unbiased HVG selection."""
    import json
    import os
    import scanpy as sc
    from geneformer import TranscriptomeTokenizer

    med, tok, mapf = _gc104m_dicts()
    adata = sc.read_h5ad(IPF_H5AD)

    # Unbiased HVG selection (data-driven) to bound the perturbation set.
    a = adata.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=n_hvg)
    hvg = adata.var.loc[a.var["highly_variable"].to_numpy(), "ensembl_id"].tolist()
    json.dump(hvg, open(HVG_JSON, "w"))
    print(f"selected {len(hvg)} HVGs")

    tok_in = f"{DATA_DIR}/tok_in"
    os.makedirs(tok_in, exist_ok=True)
    adata.write_h5ad(f"{tok_in}/ipf.h5ad")

    tk = TranscriptomeTokenizer(
        custom_attr_name_dict={"disease": "disease", "cell_type": "cell_type"},
        nproc=4, model_input_size=4096, special_token=True,
        gene_median_file=med, token_dictionary_file=tok, gene_mapping_file=mapf,
    )
    tk.tokenize_data(tok_in, DATA_DIR, "ipf", file_format="h5ad")
    data_vol.commit()
    from datasets import load_from_disk
    ds = load_from_disk(DATASET)
    print(f"[done] dataset {DATASET}: {ds}")
    return {"n_cells": ds.num_rows, "n_hvg": len(hvg), "cols": ds.column_names}


@app.local_entrypoint()
def main():
    print(tokenize.remote())
