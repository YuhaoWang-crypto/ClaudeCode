"""Step 4 -- antibody druggability, localisation and essentiality annotation.

Open Targets supplies the antibody tractability bucket, an independent
subcellular localisation call and the list of drugs already aimed at the
target.  DepMap mean gene effect is attached as *annotation only*: it is never
used as a filter and never enters the score.

The reason is mechanical.  TROP2, c-MET and HER2 -- targets with approved
antibody drugs in this indication -- are all non-essential in DepMap.  Gating
on essentiality would drop them and enrich for housekeeping genes whose loss
kills every cell in the body, which is the opposite of what a surface-antigen
campaign wants.
"""

from __future__ import annotations

import json
import time

import numpy as np
import pandas as pd
import requests

from . import config as C

QUERY = """
query($ids:[String!]!){
  targets(ensemblIds:$ids){
    id
    approvedSymbol
    isEssential
    tractability{label modality value}
    subcellularLocations{location source}
    drugAndClinicalCandidates{count rows{drug{name maximumClinicalStage}}}
  }
}
"""

# Open Targets antibody tractability buckets, best first.
AB_BUCKET_SCORE = [
    ("Approved Drug", 1.00),
    ("Advanced Clinical", 0.90),
    ("Phase 1 Clinical", 0.80),
    ("UniProt loc high conf", 0.70),
    ("GO CC high conf", 0.65),
    ("UniProt loc med conf", 0.55),
    ("UniProt SigP or TMHMM", 0.45),
    ("GO CC med conf", 0.40),
    ("Human Protein Atlas loc", 0.35),
]
SURFACE_LOC_TOKENS = ("plasma membrane", "cell membrane", "cell surface",
                      "apical", "basolateral")


def _post(ids: list[str], retries: int = 4) -> list[dict]:
    for attempt in range(retries):
        try:
            r = requests.post(C.OPENTARGETS_GQL, json={"query": QUERY, "variables": {"ids": ids}},
                              timeout=120)
            r.raise_for_status()
            payload = r.json()
            if "errors" in payload and "data" not in payload:
                raise RuntimeError(payload["errors"])
            return [t for t in (payload.get("data", {}).get("targets") or []) if t]
        except Exception as exc:  # network hiccup -> exponential backoff
            if attempt == retries - 1:
                print(f"    open targets batch failed permanently: {exc}")
                return []
            time.sleep(2 ** (attempt + 1))
    return []


def fetch_open_targets(ensembl_ids: list[str], batch: int = 100) -> pd.DataFrame:
    rows = []
    ids = [i for i in ensembl_ids if isinstance(i, str) and i.startswith("ENSG")]

    # reuse anything already fetched in an earlier run; only ask for the rest
    cache_path = C.RESULTS / "s4_open_targets.csv"
    cached = pd.DataFrame()
    if cache_path.exists():
        cached = pd.read_csv(cache_path)
        have = set(cached.ensembl_gene)
        missing = [i for i in ids if i not in have]
        print(f"  open targets cache: {len(have & set(ids)):,} of {len(ids):,} targets on disk")
        ids = missing
        cached = cached[cached.ensembl_gene.isin(set(ensembl_ids))]

    for i in range(0, len(ids), batch):
        chunk = ids[i:i + batch]
        for t in _post(chunk):
            ab = {x["label"]: x["value"] for x in t["tractability"] if x["modality"] == "AB"}
            score, bucket = 0.0, "no antibody tractability bucket"
            for label, value in AB_BUCKET_SCORE:
                if ab.get(label):
                    score, bucket = value, label
                    break
            locs = [f"{x['source']}:{x['location']}" for x in (t["subcellularLocations"] or [])]
            surface_loc = any(any(tok in l.lower() for tok in SURFACE_LOC_TOKENS) for l in locs)
            drugs = t.get("drugAndClinicalCandidates") or {}
            names = [r["drug"]["name"] for r in (drugs.get("rows") or [])][:6]
            stages = [r["drug"].get("maximumClinicalStage") for r in (drugs.get("rows") or [])]
            stages = [s for s in stages if s]
            rows.append({
                "ensembl_gene": t["id"],
                "ot_symbol": t["approvedSymbol"],
                "ab_tractability_bucket": bucket,
                "ab_druggability": score,
                "ot_surface_localisation": surface_loc,
                "ot_locations": "; ".join(locs[:6]),
                "ot_is_essential": t.get("isEssential"),
                "ot_n_drugs": int(drugs.get("count") or 0),
                "ot_drugs": "; ".join(names),
                "ot_max_clinical_stage": max(stages) if stages else "",
            })
        print(f"  open targets: {min(i + batch, len(ids)):,}/{len(ids):,} newly fetched")
    fresh = pd.DataFrame(rows)
    if len(cached) and len(fresh):
        return pd.concat([cached, fresh], ignore_index=True).drop_duplicates("ensembl_gene")
    return fresh if len(fresh) else cached


def depmap_mean_effect() -> pd.DataFrame:
    """Mean Chronos gene effect across all DepMap cell lines. Annotation only."""
    if not C.DEPMAP_FILE.exists():
        print("  DepMap matrix absent -- essentiality annotation skipped")
        return pd.DataFrame(columns=["gene", "depmap_mean_gene_effect"])
    sums: dict[str, float] = {}
    counts: dict[str, int] = {}
    for chunk in pd.read_csv(C.DEPMAP_FILE, chunksize=200, index_col=0, low_memory=False):
        block_sum = chunk.sum(axis=0, numeric_only=True)
        block_n = chunk.notna().sum(axis=0)
        for k, v in block_sum.items():
            sums[k] = sums.get(k, 0.0) + float(v)
        for k, v in block_n.items():
            counts[k] = counts.get(k, 0) + int(v)
    rows = []
    for col, total in sums.items():
        n = counts.get(col, 0)
        if n == 0:
            continue
        gene = col.split(" (")[0]
        rows.append({"gene": gene, "depmap_mean_gene_effect": total / n,
                     "depmap_n_cell_lines": n})
    df = pd.DataFrame(rows).groupby("gene", as_index=False).mean(numeric_only=True)
    print(f"  DepMap: mean gene effect for {len(df):,} genes "
          f"across {int(df.depmap_n_cell_lines.max())} cell lines")
    return df


def run(candidates: pd.DataFrame | None = None) -> pd.DataFrame:
    if candidates is None:
        candidates = pd.read_csv(C.RESULTS / "s3_accessible_candidates.csv")
    print(f"[S4] druggability + essentiality annotation for {len(candidates)} candidates")

    ot = fetch_open_targets(sorted(set(candidates.ensembl_gene.dropna())))
    ot.to_csv(C.RESULTS / "s4_open_targets.csv", index=False)

    dm = depmap_mean_effect()
    dm.to_csv(C.RESULTS / "s4_depmap.csv", index=False)

    df = candidates.merge(ot, on="ensembl_gene", how="left")
    df = df.merge(dm[["gene", "depmap_mean_gene_effect"]], on="gene", how="left")
    df["ab_druggability"] = df.ab_druggability.fillna(0.30)  # unqueried -> conservative
    df["ab_tractability_bucket"] = df.ab_tractability_bucket.fillna("not in Open Targets")

    upgrade = (df.surface_confirmation.str.startswith("unconfirmed")) & (df.ot_surface_localisation == True)  # noqa: E712
    df.loc[upgrade, "surface_confirmation"] = "Open Targets confirmed"

    df.to_csv(C.RESULTS / "s4_annotated_candidates.csv", index=False)

    tiers = df.surface_confirmation.value_counts().to_dict()
    print("  surface confirmation tiers: " +
          ", ".join(f"{k}={v}" for k, v in tiers.items()))

    validated = [g for g in C.VALIDATION_POSITIVES if g in set(df.gene)]
    ess = df[df.gene.isin(validated)][["gene", "depmap_mean_gene_effect"]].dropna()
    print("  DepMap sanity check on validated antigens (annotation only): " +
          ", ".join(f"{r.gene}={r.depmap_mean_gene_effect:+.2f}" for r in ess.itertuples()))

    C.write_provenance("s4", {
        "candidates": int(len(df)),
        "open_targets_hits": int(ot.ensembl_gene.nunique()) if len(ot) else 0,
        "surface_confirmation_tiers": {str(k): int(v) for k, v in tiers.items()},
        "depmap_genes": int(len(dm)),
        "depmap_used_as_filter": False,
        "validated_antigen_gene_effect": {
            r.gene: round(float(r.depmap_mean_gene_effect), 3) for r in ess.itertuples()
        },
        "median_gene_effect_all_candidates": float(
            np.nanmedian(df.depmap_mean_gene_effect.to_numpy(dtype=float))
        ) if df.depmap_mean_gene_effect.notna().any() else None,
    })
    return df


if __name__ == "__main__":
    run()
