"""Step 5 -- normal-tissue safety baseline (therapeutic index).

Two independent Human Protein Atlas releases are used: consensus RNA (nTPM
across normal tissues) and normal-tissue immunohistochemistry (protein level
per tissue and cell type).  Each gene gets an organ-weighted risk from both,
and the safety coefficient is the *conservative minimum* of the two: a target
that looks clean in RNA but stains strongly in heart muscle is treated as
dangerous.

Genes with no HPA record at all get a flagged neutral default rather than a
free pass or a silent drop.
"""

from __future__ import annotations

import urllib.request
import zipfile

import numpy as np
import pandas as pd

from . import config as C

GOOD_RELIABILITY = {"Enhanced", "Supported", "Approved"}


def _download(url, path):
    if path.exists():
        return path
    print(f"  downloading {url}")
    urllib.request.urlretrieve(url, path)
    return path


def _read_zip_tsv(path) -> pd.DataFrame:
    with zipfile.ZipFile(path) as zf:
        name = zf.namelist()[0]
        with zf.open(name) as fh:
            return pd.read_csv(fh, sep="\t")


def organ_weight(tissue: str) -> float:
    t = str(tissue).strip().lower()
    w = C.VITAL_ORGAN_WEIGHTS.get(t)
    if w is None:
        w = next((v for k, v in C.VITAL_ORGAN_WEIGHTS.items() if k in t),
                 C.DEFAULT_ORGAN_WEIGHT)
    return min(w, C.MAX_ORGAN_WEIGHT)


def run(candidates: pd.DataFrame | None = None) -> pd.DataFrame:
    if candidates is None:
        candidates = pd.read_csv(C.RESULTS / "s4_annotated_candidates.csv")
    print("[S5] normal-tissue safety baseline from HPA")

    _download(C.HPA_RNA_URL, C.HPA_RNA_FILE)
    _download(C.HPA_IHC_URL, C.HPA_IHC_FILE)

    rna = _read_zip_tsv(C.HPA_RNA_FILE)
    ihc = _read_zip_tsv(C.HPA_IHC_FILE)
    print(f"  HPA consensus RNA: {rna['Gene name'].nunique():,} genes x "
          f"{rna.Tissue.nunique()} tissues")
    print(f"  HPA normal IHC:    {ihc['Gene name'].nunique():,} genes x "
          f"{ihc.Tissue.nunique()} tissues")

    # ---- RNA arm ---------------------------------------------------------
    rna = rna.rename(columns={"Gene name": "gene", "Tissue": "tissue", "nTPM": "ntpm"})
    rna["weight"] = rna.tissue.map(organ_weight)
    rna["risk"] = rna.weight * (
        np.log1p(rna.ntpm.clip(lower=0)) / np.log1p(C.RNA_SATURATION_NTPM)
    ).clip(0, 1)
    rna_g = rna.groupby("gene").agg(
        rna_risk=("risk", "max"),
        rna_max_ntpm=("ntpm", "max"),
    ).reset_index()
    worst_rna = rna.loc[rna.groupby("gene").risk.idxmax(), ["gene", "tissue"]]
    rna_g = rna_g.merge(worst_rna.rename(columns={"tissue": "rna_worst_tissue"}), on="gene")
    rna_g["safety_rna"] = (1.0 - rna_g.rna_risk).clip(0, 1)

    # ---- protein arm -----------------------------------------------------
    ihc = ihc.rename(columns={"Gene name": "gene", "Tissue": "tissue",
                              "Level": "level", "Reliability": "reliability"})
    ihc = ihc[ihc.reliability.isin(GOOD_RELIABILITY)]
    ihc["level_score"] = ihc.level.str.lower().map(C.IHC_LEVEL_SCORE)
    ihc = ihc.dropna(subset=["level_score"])
    # strongest-stained cell type defines the tissue's level
    per_tissue = ihc.groupby(["gene", "tissue"], as_index=False).level_score.max()
    per_tissue["weight"] = per_tissue.tissue.map(organ_weight)
    per_tissue["risk"] = per_tissue.weight * per_tissue.level_score
    prot_g = per_tissue.groupby("gene").agg(protein_risk=("risk", "max")).reset_index()
    worst_prot = per_tissue.loc[per_tissue.groupby("gene").risk.idxmax(),
                                ["gene", "tissue", "level_score"]]
    prot_g = prot_g.merge(worst_prot.rename(columns={"tissue": "protein_worst_tissue"}), on="gene")
    prot_g["safety_protein"] = (1.0 - prot_g.protein_risk).clip(0, 1)
    n_high_tissues = per_tissue[per_tissue.level_score >= 0.66].groupby("gene").size()
    prot_g["n_normal_tissues_medium_high"] = prot_g.gene.map(n_high_tissues).fillna(0).astype(int)

    df = candidates.merge(rna_g, on="gene", how="left").merge(prot_g, on="gene", how="left")

    have_any = df.safety_rna.notna() | df.safety_protein.notna()
    df["safety_source"] = np.where(
        df.safety_rna.notna() & df.safety_protein.notna(), "protein+RNA (min)",
        np.where(df.safety_protein.notna(), "protein only",
                 np.where(df.safety_rna.notna(), "RNA only", "default (no HPA record)")))
    df["safety_score"] = np.where(
        have_any,
        np.fmin(df.safety_protein.fillna(np.inf), df.safety_rna.fillna(np.inf)),
        C.SAFETY_MISSING_DEFAULT,
    )
    df["safety_is_default"] = ~have_any

    n_default = int(df.safety_is_default.sum())
    print(f"  safety scored for {len(df) - n_default:,} candidates; "
          f"{n_default} fell back to the flagged default {C.SAFETY_MISSING_DEFAULT}")
    print(f"  median safety coefficient {df.safety_score.median():.3f}")

    df.to_csv(C.RESULTS / "s5_safety_annotated.csv", index=False)

    val = df[df.gene.isin(C.VALIDATION_POSITIVES)][
        ["gene", "safety_score", "protein_worst_tissue", "rna_worst_tissue", "rna_max_ntpm"]]
    print("  validated antigens, safety coefficient: " +
          ", ".join(f"{r.gene}={r.safety_score:.2f}" for r in val.itertuples()))

    C.write_provenance("s5", {
        "hpa_rna_genes": int(rna.gene.nunique()),
        "hpa_rna_tissues": int(rna.tissue.nunique()),
        "hpa_ihc_genes": int(ihc.gene.nunique()),
        "hpa_ihc_tissues": int(ihc.tissue.nunique()),
        "candidates_scored": int(len(df) - n_default),
        "candidates_default_safety": n_default,
        "default_value": C.SAFETY_MISSING_DEFAULT,
        "median_safety": float(df.safety_score.median()),
        "validated_antigen_safety": {
            r.gene: round(float(r.safety_score), 3) for r in val.itertuples()
        },
        "validated_antigen_worst_normal_tissue": {
            r.gene: str(r.protein_worst_tissue) for r in val.itertuples()
        },
    })
    return df


if __name__ == "__main__":
    run()
