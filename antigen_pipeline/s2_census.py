"""Step 2 -- tumour single-cell expression from CZ CELLxGENE Census.

Discovers every whole-cell lung adenocarcinoma atlas in a pinned Census
release, samples each one to a fixed ceiling with a fixed seed, assigns every
cell to one of four tumour-microenvironment compartments, and reduces the
cohort to per-gene, per-compartment, per-dataset mean expression and detection
rate.  Cross-dataset consensus is computed on top of the per-dataset table so
that a target seen in one atlas only can never look reproducible.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from . import config as C


def assign_compartment(cell_type: str) -> str | None:
    ct = str(cell_type).lower()
    for pattern, compartment in C.COMPARTMENT_RULES:
        if pattern in ct:
            return compartment
    return None


def _compartment_stats(X, cell_idx: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Mean expression and detection rate over a set of rows of X."""
    sub = X[cell_idx]
    n = sub.shape[0]
    if n == 0:
        return np.zeros(X.shape[1]), np.zeros(X.shape[1])
    if hasattr(sub, "toarray"):
        mean = np.asarray(sub.mean(axis=0)).ravel()
        detected = np.asarray((sub > 0).sum(axis=0)).ravel() / n
    else:
        mean = sub.mean(axis=0)
        detected = (sub > 0).mean(axis=0)
    return mean, detected


def fetch(surfaceome: pd.DataFrame | None = None) -> pd.DataFrame:
    """Pull the cohort out of Census and write the per-dataset long table."""
    import cellxgene_census as cc

    if surfaceome is None:
        surfaceome = pd.read_csv(C.RESULTS / "s1_surfaceome.csv")
    wanted_genes = set(surfaceome.gene)

    print(f"[S2] CZ CELLxGENE Census {C.CENSUS_VERSION}: discovering "
          f"'{C.DISEASE_LABEL}' atlases")
    value_filter = (
        f"disease == '{C.DISEASE_LABEL}' and is_primary_data == True "
        f"and suspension_type == '{C.SUSPENSION_TYPE}'"
    )

    per_dataset_rows = []
    cohort_rows = []

    with cc.open_soma(census_version=C.CENSUS_VERSION) as census:
        human = census["census_data"]["homo_sapiens"]
        obs = human.obs.read(
            value_filter=value_filter,
            column_names=["soma_joinid", "dataset_id", "assay", "cell_type",
                          "tissue_general", "donor_id"],
        ).concat().to_pandas()
        print(f"  {len(obs):,} whole-cell {C.DISEASE_LABEL} cells "
              f"in {obs.dataset_id.nunique()} atlases")

        datasets_meta = census["census_info"]["datasets"].read().concat().to_pandas()
        meta_by_id = datasets_meta.set_index("dataset_id")

        obs["compartment"] = obs.cell_type.map(assign_compartment)

        rng = np.random.default_rng(C.SAMPLING_SEED)
        for dataset_id, block in obs.groupby("dataset_id", observed=True):
            if len(block) == 0:
                continue
            if len(block) > C.MAX_CELLS_PER_DATASET:
                take = rng.choice(len(block), C.MAX_CELLS_PER_DATASET, replace=False)
                block = block.iloc[np.sort(take)]
            usable = block[block.compartment.notna()]
            title = str(meta_by_id.loc[dataset_id, "dataset_title"]) if dataset_id in meta_by_id.index else ""
            collection = str(meta_by_id.loc[dataset_id, "collection_name"]) if dataset_id in meta_by_id.index else ""
            print(f"  {dataset_id[:8]}  sampled {len(block):,} -> "
                  f"{len(usable):,} in-compartment  | {title[:60]}")

            adata = cc.get_anndata(
                census=census,
                organism="Homo sapiens",
                measurement_name="RNA",
                X_name=C.CENSUS_LAYER,
                obs_coords=usable.soma_joinid.to_numpy(),
                column_names={"obs": ["soma_joinid", "cell_type"],
                              "var": ["feature_id", "feature_name"]},
            )
            var_names = adata.var["feature_name"].astype(str).to_numpy()
            keep = np.where(np.isin(var_names, list(wanted_genes)))[0]
            X = adata.X[:, keep]
            genes_here = var_names[keep]
            comp = usable.set_index("soma_joinid").loc[
                adata.obs.soma_joinid.to_numpy(), "compartment"
            ].to_numpy()

            counts = {}
            for compartment in C.COMPARTMENTS:
                idx = np.where(comp == compartment)[0]
                counts[compartment] = int(len(idx))
                mean, det = _compartment_stats(X, idx)
                per_dataset_rows.append(pd.DataFrame({
                    "gene": genes_here,
                    "dataset_id": dataset_id,
                    "compartment": compartment,
                    "mean_expression": mean,
                    "detection_rate": det,
                    "n_cells": len(idx),
                }))

            cohort_rows.append({
                "dataset_id": dataset_id,
                "collection": collection,
                "title": title,
                "assays": sorted(set(block.assay.astype(str))),
                "donors": int(block.donor_id.nunique()),
                "cells_available": int((obs.dataset_id == dataset_id).sum()),
                "cells_sampled": int(len(block)),
                "cells_in_compartments": int(len(usable)),
                "genes_measured": int(len(genes_here)),
                **{f"cells_{k}": v for k, v in counts.items()},
            })
            del adata, X

    long = pd.concat(per_dataset_rows, ignore_index=True)
    long.to_csv(C.RESULTS / "s2_expression_by_dataset.csv.gz", index=False,
                compression="gzip")

    cohort = pd.DataFrame(cohort_rows)
    cohort.to_csv(C.RESULTS / "s2_cohort.csv", index=False)
    (C.RESULTS / "s2_cohort.json").write_text(json.dumps(cohort_rows, indent=2, default=str))
    return long


def build_consensus(long: pd.DataFrame | None = None) -> pd.DataFrame:
    """Reduce the per-dataset long table to a cross-atlas consensus.

    Census stores library-size-normalised values, i.e. per-cell fractions of
    about 1e-5 for a typical gene.  Everything here works in counts per 10,000
    (CP10K) instead, so that the pseudocount in the specificity ratio is small
    relative to real expression rather than swamping it.
    """
    if long is None:
        long = pd.read_csv(C.RESULTS / "s2_expression_by_dataset.csv.gz")
    cohort = pd.read_csv(C.RESULTS / "s2_cohort.csv")

    long = long.copy()
    long["mean_expression"] = long.mean_expression * C.EXPRESSION_SCALE

    wide = long.pivot_table(index=["gene", "dataset_id"], columns="compartment",
                            values=["mean_expression", "detection_rate"])
    wide.columns = [f"{a}_{b}" for a, b in wide.columns]
    wide = wide.reset_index()

    other = ["mean_expression_caf", "mean_expression_immune", "mean_expression_endothelial"]
    wide["tme_max"] = wide[other].max(axis=1)
    # log2 fold change of malignant/epithelial over the loudest other compartment
    wide["log2fc_epi_vs_tme"] = np.log2(
        (wide["mean_expression_epithelial"] + C.EXPRESSION_PSEUDOCOUNT)
        / (wide["tme_max"] + C.EXPRESSION_PSEUDOCOUNT)
    )
    wide.to_csv(C.RESULTS / "s2_per_dataset_wide.csv.gz", index=False, compression="gzip")

    n_datasets = wide.dataset_id.nunique()
    agg = wide.groupby("gene").agg(
        epi_mean=("mean_expression_epithelial", "mean"),
        epi_detection=("detection_rate_epithelial", "mean"),
        caf_mean=("mean_expression_caf", "mean"),
        immune_mean=("mean_expression_immune", "mean"),
        endothelial_mean=("mean_expression_endothelial", "mean"),
        tme_max_mean=("tme_max", "mean"),
        log2fc_mean=("log2fc_epi_vs_tme", "mean"),
        log2fc_sd=("log2fc_epi_vs_tme", "std"),
        n_datasets=("dataset_id", "nunique"),
    ).reset_index()
    # reproducibility: in how many atlases is the target both enriched in the
    # epithelial compartment and actually detected there?
    repro = wide.assign(
        ok=lambda d: (d.log2fc_epi_vs_tme > C.CONSENSUS_MIN_LOG2FC)
        & (d.detection_rate_epithelial > C.CONSENSUS_MIN_DETECTION)
    ).groupby("gene")["ok"].sum().rename("n_datasets_reproducible")
    agg = agg.merge(repro, on="gene")
    agg["consensus_fraction"] = agg.n_datasets_reproducible / n_datasets
    agg.to_csv(C.RESULTS / "s2_consensus_expression.csv", index=False)

    total_available = int(cohort.cells_available.sum())
    total_sampled = int(cohort.cells_sampled.sum())
    total_used = int(cohort.cells_in_compartments.sum())
    print(f"  cohort: {total_available:,} available -> {total_sampled:,} sampled -> "
          f"{total_used:,} assigned to a compartment")

    C.write_provenance("s2", {
        "census_version": C.CENSUS_VERSION,
        "disease_label": C.DISEASE_LABEL,
        "suspension_type": C.SUSPENSION_TYPE,
        "max_cells_per_dataset": C.MAX_CELLS_PER_DATASET,
        "sampling_seed": C.SAMPLING_SEED,
        "layer": C.CENSUS_LAYER,
        "expression_units": f"counts per {int(C.EXPRESSION_SCALE):,} (CP10K)",
        "expression_pseudocount": C.EXPRESSION_PSEUDOCOUNT,
        "n_datasets": int(n_datasets),
        "cells_available": total_available,
        "cells_sampled": total_sampled,
        "cells_analysed": total_used,
        "cells_per_compartment": {
            c: int(cohort[f"cells_{c}"].sum()) for c in C.COMPARTMENTS
        },
        "genes_with_expression": int(agg.gene.nunique()),
    })
    return agg


def run(surfaceome: pd.DataFrame | None = None, use_cache: bool = True) -> pd.DataFrame:
    """Fetch (unless the long table is already on disk) and aggregate."""
    cached = C.RESULTS / "s2_expression_by_dataset.csv.gz"
    if use_cache and cached.exists():
        print(f"[S2] reusing cached cohort matrix {cached.name}")
        long = pd.read_csv(cached)
    else:
        long = fetch(surfaceome)
    return build_consensus(long)


if __name__ == "__main__":
    import sys
    run(use_cache="--refetch" not in sys.argv)
