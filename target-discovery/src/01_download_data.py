"""
Step 1 — Acquire the IPF single-cell lung atlas.

Pulls a harmonized IPF vs. healthy lung dataset. Default path uses CELLxGENE
Census (no manual GEO wrangling); a GEO fallback is stubbed for Habermann 2020
(GSE135893) and Adams 2020 (GSE136831).

Output: an AnnData .h5ad with raw counts + a `disease` obs column, written to
config.data.processed_h5ad. Downstream steps assume raw counts in .X or .raw.
"""
from __future__ import annotations
import argparse
import os

import yaml


def load_config(path: str) -> dict:
    with open(path) as fh:
        return yaml.safe_load(fh)


def download_via_census(cfg: dict) -> "anndata.AnnData":
    import cellxgene_census

    d = cfg["data"]
    print(f"[census] opening Census {d['census_version']}")
    with cellxgene_census.open_soma(census_version=d["census_version"]) as census:
        # Pull IPF + healthy lung cells in one harmonized query.
        obs_filter = (
            f"tissue_general == 'lung' and "
            f"disease in ['{d['disease_label']}', '{d['healthy_label']}']"
        )
        adata = cellxgene_census.get_anndata(
            census,
            organism=d["organism"],
            obs_value_filter=obs_filter,
            column_names={"obs": ["cell_type", "disease", "donor_id", "assay"]},
        )
    print(f"[census] fetched {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    return adata


def download_via_geo(cfg: dict):
    """Fallback: download raw supplementary matrices from GEO.

    Left as a stub — implement per-accession parsing (10x mtx or processed
    h5ad) for cfg.data.geo_accessions. CELLxGENE is strongly preferred.
    """
    raise NotImplementedError(
        "GEO path not implemented; use source: cellxgene_census in config, "
        f"or add parsing for {cfg['data']['geo_accessions']}."
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.yaml")
    args = ap.parse_args()
    cfg = load_config(args.config)

    out = cfg["data"]["processed_h5ad"]
    os.makedirs(os.path.dirname(out), exist_ok=True)

    if cfg["data"]["source"] == "cellxgene_census":
        adata = download_via_census(cfg)
    else:
        adata = download_via_geo(cfg)

    # Basic sanity: keep protein-coding-ish genes, drop empty cells.
    import scanpy as sc

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=10)
    adata.write_h5ad(out)
    print(f"[done] wrote {out}  ({adata.n_obs:,} cells, {adata.n_vars:,} genes)")


if __name__ == "__main__":
    main()
