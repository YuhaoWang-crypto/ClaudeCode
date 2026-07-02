"""
Step 2 — Tokenize the atlas into Geneformer's rank-value encoding.

Geneformer represents each cell as its genes ranked by expression (normalized by
each gene's non-zero median across a large corpus). This step:
  1. ensures Ensembl IDs + a raw count layer (Geneformer requirements),
  2. writes a .loom / .h5ad that geneformer.TranscriptomeTokenizer consumes,
  3. emits a tokenized HF dataset at config.tokenize.out_dataset.

No GPU needed. Mirrors the preprocessing in bionemo-recipes/recipes/geneformer.
"""
from __future__ import annotations
import argparse
import os

import yaml


def load_config(path: str) -> dict:
    with open(path) as fh:
        return yaml.safe_load(fh)


def prepare_for_geneformer(adata, cfg: dict):
    """Add the obs/var fields Geneformer's tokenizer expects."""
    import numpy as np

    # Geneformer keys cells on raw counts in a layer, and genes on Ensembl IDs.
    if "ensembl_id" not in adata.var.columns:
        # CELLxGENE Census already stores Ensembl IDs in var index; adapt as needed.
        adata.var["ensembl_id"] = adata.var.get("feature_id", adata.var_names)
    adata.obs["n_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()

    # Tag each cell with which pathological / healthy state it belongs to, so
    # Step 3 can compute the goal-shift from pathological -> healthy.
    state_col = _assign_states(adata, cfg)
    adata.obs["cell_state"] = state_col
    return adata


def _assign_states(adata, cfg: dict):
    """Coarse state assignment by marker-score argmax over configured states.

    Real runs should use a proper annotation (Leiden + marker scoring or a
    reference-mapped label). This keeps the dependency surface minimal.
    """
    import numpy as np
    import scanpy as sc

    labels = np.array(["other"] * adata.n_obs, dtype=object)
    best = np.full(adata.n_obs, -np.inf)
    for state in cfg["target_cell_states"]:
        genes = [g for g in state["markers"] if g in adata.var_names]
        if not genes:
            continue
        sc.tl.score_genes(adata, genes, score_name="_s")
        s = adata.obs["_s"].to_numpy()
        take = s > best
        labels[take] = state["name"]
        best[take] = s[take]
    del adata.obs["_s"]
    return labels


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.yaml")
    args = ap.parse_args()
    cfg = load_config(args.config)

    import anndata as ad

    adata = ad.read_h5ad(cfg["data"]["processed_h5ad"])
    adata = prepare_for_geneformer(adata, cfg)

    tok_h5ad = cfg["data"]["processed_h5ad"].replace(".h5ad", ".tokin.h5ad")
    adata.write_h5ad(tok_h5ad)

    # ---- Geneformer tokenization -------------------------------------------
    from geneformer import TranscriptomeTokenizer

    out_ds = cfg["tokenize"]["out_dataset"]
    os.makedirs(os.path.dirname(out_ds), exist_ok=True)
    tk = TranscriptomeTokenizer(
        custom_attr_name_dict={"cell_state": "cell_state", "disease": "disease"},
        model_input_size=cfg["tokenize"]["model_input_size"],
    )
    tk.tokenize_data(
        data_directory=os.path.dirname(tok_h5ad),
        output_directory=os.path.dirname(out_ds),
        output_prefix=os.path.basename(out_ds).replace(".dataset", ""),
        file_format="h5ad",
    )
    print(f"[done] tokenized dataset -> {out_ds}")


if __name__ == "__main__":
    main()
