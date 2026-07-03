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


@app.local_entrypoint()
def main():
    print(inspect_census.remote())
