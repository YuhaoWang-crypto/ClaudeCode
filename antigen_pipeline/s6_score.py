"""Step 6 -- composite score and tiering.

    final score = tumour quality  x  safety coefficient  x  consensus multiplier

Tumour quality is a weighted sum of five bounded sub-scores.  Safety enters
*multiplicatively*, so a target that is loud in a vital organ is pushed down
however good its tumour profile looks; no additive weighting can buy that
back.  The consensus multiplier rewards targets that reproduce across
independent atlases.

Tiers are absolute thresholds, not percentiles, so the tier of a target does
not change when other targets are added or removed from the run.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import config as C


def _intensity_score(epi_mean: pd.Series) -> pd.Series:
    """log-scaled expression, anchored on the 99th percentile of the cohort."""
    x = np.log1p(epi_mean.clip(lower=0))
    anchor = np.log1p(epi_mean.clip(lower=0).quantile(0.99))
    if not np.isfinite(anchor) or anchor <= 0:
        return pd.Series(np.zeros(len(epi_mean)), index=epi_mean.index)
    return (x / anchor).clip(0, 1)


def build_subscores(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["score_specificity"] = (df.log2fc_mean / 3.0).clip(0, 1).fillna(0)
    df["score_intensity"] = _intensity_score(df.epi_mean.fillna(0))
    df["score_uniformity"] = df.epi_detection.fillna(0).clip(0, 1)
    df["score_accessibility"] = df.accessibility_score.fillna(0).clip(0, 1)
    df["score_druggability"] = df.ab_druggability.fillna(0.3).clip(0, 1)
    w = C.TUMOUR_QUALITY_WEIGHTS
    df["tumour_quality"] = (
        w["specificity"] * df.score_specificity
        + w["intensity"] * df.score_intensity
        + w["uniformity"] * df.score_uniformity
        + w["accessibility"] * df.score_accessibility
        + w["druggability"] * df.score_druggability
    )
    df["consensus_multiplier"] = 0.5 + 0.5 * df.consensus_fraction.fillna(0)
    df["multiplicative_score"] = df.tumour_quality * df.safety_score * df.consensus_multiplier
    df["final_score"] = composite_index(df.multiplicative_score)
    return df


def composite_index(product: pd.Series | np.ndarray):
    """Geometric mean of the three factors.

    This is the cube root of the product, so it is a strictly monotone
    transform: it cannot reorder two candidates, and every rank, recall and
    stability number is identical to the one the raw product gives. What it
    changes is readability -- a score of 0.55 means the target averages 0.55
    across tumour quality, safety and consensus, which is what the absolute
    tier thresholds are stated against.
    """
    return np.cbrt(np.asarray(product, dtype=float))


def assign_tier(score: float) -> str:
    if score >= C.TIER1_THRESHOLD:
        return "Tier 1"
    if score >= C.TIER2_THRESHOLD:
        return "Tier 2"
    return "Tier 3"


def run(annotated: pd.DataFrame | None = None,
        expression: pd.DataFrame | None = None) -> pd.DataFrame:
    if annotated is None:
        annotated = pd.read_csv(C.RESULTS / "s5_safety_annotated.csv")
    if expression is None:
        expression = pd.read_csv(C.RESULTS / "s2_consensus_expression.csv")
    print(f"[S6] scoring {len(annotated)} accessible candidates")

    merged = annotated.merge(expression, on="gene", how="inner")
    dropped = len(annotated) - len(merged)
    print(f"  {len(merged)} candidates have tumour expression data "
          f"({dropped} accessible genes are absent from the Census feature set)")

    # Antigen display requirement, applied before ranking. Without it a gene
    # that no tumour cell expresses wins on a safety coefficient earned purely
    # by being silent everywhere.
    eligible = merged.epi_detection.fillna(0) >= C.MIN_EPITHELIAL_DETECTION
    set_aside = merged[~eligible].copy()
    set_aside.sort_values("epi_detection", ascending=False).to_csv(
        C.RESULTS / "s6_below_display_threshold.csv", index=False)
    df = merged[eligible].copy()
    print(f"  {len(df)} candidates are displayed on >= "
          f"{C.MIN_EPITHELIAL_DETECTION:.0%} of epithelial/malignant cells "
          f"({len(set_aside)} set aside below that threshold)")

    df = build_subscores(df)
    df = df.sort_values("final_score", ascending=False).reset_index(drop=True)
    df["rank"] = np.arange(1, len(df) + 1)
    df["tier"] = df.final_score.map(assign_tier)

    cols = ["rank", "gene", "protein_name", "final_score", "tier",
            "multiplicative_score", "tumour_quality",
            "safety_score", "consensus_multiplier", "score_specificity",
            "score_intensity", "score_uniformity", "score_accessibility",
            "score_druggability", "epi_mean", "epi_detection",
            "caf_mean", "immune_mean", "endothelial_mean", "tme_max_mean",
            "log2fc_mean", "log2fc_sd", "n_datasets_reproducible", "n_datasets",
            "topology_class", "max_ecd_segment", "surface_confirmation",
            "ab_tractability_bucket", "ot_n_drugs", "ot_drugs",
            "depmap_mean_gene_effect", "safety_source", "safety_is_default",
            "protein_worst_tissue", "rna_worst_tissue", "rna_max_ntpm",
            "n_normal_tissues_medium_high", "almen_class", "uniprot", "ensembl_gene"]
    cols = [c for c in cols if c in df.columns]
    ranked = df[cols]
    ranked.to_csv(C.RESULTS / "s6_ranked_candidates.csv", index=False)

    counts = ranked.tier.value_counts().to_dict()
    print(f"  Tier 1 (>= {C.TIER1_THRESHOLD}): {counts.get('Tier 1', 0)} | "
          f"Tier 2 ({C.TIER2_THRESHOLD}-{C.TIER1_THRESHOLD}): {counts.get('Tier 2', 0)} | "
          f"Tier 3: {counts.get('Tier 3', 0)}")
    print("  top 12: " + ", ".join(ranked.gene.head(12)))

    C.write_provenance("s6", {
        "weights": C.TUMOUR_QUALITY_WEIGHTS,
        "tier1_threshold": C.TIER1_THRESHOLD,
        "tier2_threshold": C.TIER2_THRESHOLD,
        "candidates_scored": int(len(ranked)),
        "accessible_without_expression": int(dropped),
        "min_epithelial_detection": C.MIN_EPITHELIAL_DETECTION,
        "set_aside_below_display_threshold": int(len(set_aside)),
        "tier_counts": {str(k): int(v) for k, v in counts.items()},
        "top20": ranked.head(20)[["rank", "gene", "final_score", "tier"]].to_dict("records"),
        "score_formula": ("geometric mean of tumour_quality, safety_coefficient and "
                          "consensus_multiplier (cube root of their product; a monotone "
                          "transform, so ranks are those of the raw product)"),
    })
    return ranked


if __name__ == "__main__":
    run()
