"""Shared config: paths, held-out design, and the label harmonisation used for scoring."""
import os

DATA = os.environ.get("SCB_DATA", "data")          # downloaded inputs
WORK = os.environ.get("SCB_WORK", "work")          # intermediate h5ad / models
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")

CXG = "https://datasets.cellxgene.cziscience.com/"
DATASETS = {  # CELLxGENE dataset-version assets
    "ref_hca_bm": CXG + "4ce47ed6-f41f-4c35-8120-75bbb3f55004.h5ad",   # CellTypist-v2 / HCA bone marrow (Cell 2023)
    "q_triana_wta": CXG + "a58b36e0-fafe-40d6-a4ca-a659cfe9a1f8.h5ad",  # Triana 2021, BD Rhapsody WTA, donor BM3
    "q_ts_bm": CXG + "4aedb645-7487-4ab0-ad95-ee4db1b8eccd.h5ad",       # Tabula Sapiens bone marrow
}

# Study removed from the reference so it can serve as an unseen query.
HELD_OUT_STUDY = "Tabula Sapiens 2022"
# Lineages deleted from the reference before training -> must be flagged as "unknown" in queries.
HELD_OUT_LINEAGES = ["Plasma", "pDC"]

# CL cell_type label -> coarse lineage (shared by reference and both queries)
COARSE = {
    "HSPC": ["hematopoietic multipotent progenitor cell", "hematopoietic stem cell",
             "hematopoietic precursor cell", "common myeloid progenitor",
             "granulocyte monocyte progenitor cell", "early lymphoid progenitor",
             "lymphoid lineage restricted progenitor cell"],
    "Ery/Mk": ["megakaryocyte-erythroid progenitor cell", "erythroid lineage cell",
               "erythroid progenitor cell", "erythrocyte", "megakaryocyte progenitor cell",
               "megakaryocyte"],
    "Neutrophil lineage": ["promyelocyte", "early promyelocyte", "late promyelocyte", "myelocyte",
                           "neutrophil progenitor cell", "neutrophil", "granulocyte"],
    "Mono/Mac": ["classical monocyte", "non-classical monocyte", "intermediate monocyte",
                 "monocyte", "macrophage"],
    "cDC": ["conventional dendritic cell", "common dendritic progenitor", "myeloid dendritic cell"],
    "pDC": ["plasmacytoid dendritic cell", "plasmacytoid dendritic cell, human"],
    "Mast/Baso": ["mast cell", "basophil", "basophil mast progenitor cell"],
    "B progenitor": ["pro-B cell", "small pre-B-II cell", "precursor B cell", "immature B cell",
                     "fraction A pre-pro B cell"],
    "B mature": ["naive B cell", "memory B cell", "unswitched memory B cell",
                 "class switched memory B cell", "B cell"],
    "Plasma": ["plasma cell", "plasmablast"],
    "NK": ["CD16-positive, CD56-dim natural killer cell, human",
           "CD16-negative, CD56-bright natural killer cell, human", "natural killer cell",
           "pre-natural killer cell"],
    "T": ["central memory CD4-positive, alpha-beta T cell", "effector memory CD4-positive, alpha-beta T cell",
          "effector memory CD8-positive, alpha-beta T cell, terminally differentiated",
          "effector memory CD8-positive, alpha-beta T cell", "central memory CD8-positive, alpha-beta T cell",
          "gamma-delta T cell", "mature NK T cell", "mucosal-associated invariant T cell",
          "regulatory T cell", "T follicular helper cell", "thymocyte",
          "CD4-positive, alpha-beta memory T cell", "CD8-positive, alpha-beta memory T cell",
          "naive thymus-derived CD8-positive, alpha-beta T cell",
          "naive thymus-derived CD4-positive, alpha-beta T cell",
          "CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell", "T cell"],
    "Stroma": ["fibroblast", "endothelial cell"],   # absent from the reference by nature
}
FINE2COARSE = {f: c for c, fs in COARSE.items() for f in fs}

# Ordinal developmental potency used as ground truth for the CytoTRACE 2 test
# 2 = multipotent stem/progenitor, 1 = lineage-committed progenitor/precursor, 0 = mature
POTENCY = {
    2: ["hematopoietic multipotent progenitor cell", "hematopoietic stem cell",
        "hematopoietic precursor cell"],
    1: ["common myeloid progenitor", "granulocyte monocyte progenitor cell", "early lymphoid progenitor",
        "lymphoid lineage restricted progenitor cell", "megakaryocyte-erythroid progenitor cell",
        "erythroid progenitor cell", "megakaryocyte progenitor cell", "common dendritic progenitor",
        "basophil mast progenitor cell", "promyelocyte", "early promyelocyte", "late promyelocyte",
        "myelocyte", "neutrophil progenitor cell", "pro-B cell", "small pre-B-II cell",
        "precursor B cell", "fraction A pre-pro B cell", "pre-natural killer cell",
        "erythroid lineage cell"],
}
FINE2POTENCY = {f: p for p, fs in POTENCY.items() for f in fs}
