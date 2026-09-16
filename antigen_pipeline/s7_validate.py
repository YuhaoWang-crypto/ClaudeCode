"""Step 7 -- validation and quality control.

Three checks, all reported whatever they say:

1. Recall of the pre-registered clinically validated LUAD antigens at the
   ranks fixed in step 1.
2. Negative controls.  A housekeeping pump or a pan-epithelial adhesion
   molecule ranking near the top means the safety or topology layer failed,
   so each control is given an explicit PASS/FAIL.
3. Rank stability under six alternative safety-aggregation and weighting
   schemes, to show whether the head of the ranking is an artefact of one
   particular parameterisation.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from . import config as C
from .s6_score import composite_index

ALTERNATIVES = {
    "A_safety_mean": "safety = mean(protein, RNA) instead of the conservative min",
    "B_safety_protein_only": "safety = IHC protein arm only",
    "C_safety_rna_only": "safety = consensus RNA arm only",
    "D_safety_sqrt": "safety enters as sqrt(safety): a softer toxicity penalty",
    "E_equal_weights": "all five tumour-quality sub-scores weighted 0.20",
    "F_specificity_heavy": "specificity 0.50, the other four 0.125 each",
}


def _alt_scores(df: pd.DataFrame) -> pd.DataFrame:
    """Re-score every candidate under each alternative, on the baseline scale."""
    out = pd.DataFrame(index=df.index)
    q = df.tumour_quality
    cons = df.consensus_multiplier
    arms = np.vstack([df.safety_protein.to_numpy(dtype=float),
                      df.safety_rna.to_numpy(dtype=float)])
    with np.errstate(invalid="ignore"):
        have = np.isfinite(arms).sum(axis=0)
        safety_mean = np.where(have > 0, np.nansum(np.nan_to_num(arms), axis=0)
                               / np.maximum(have, 1), C.SAFETY_MISSING_DEFAULT)
    out["A_safety_mean"] = q * safety_mean * cons
    out["B_safety_protein_only"] = q * df.safety_protein.fillna(C.SAFETY_MISSING_DEFAULT) * cons
    out["C_safety_rna_only"] = q * df.safety_rna.fillna(C.SAFETY_MISSING_DEFAULT) * cons
    out["D_safety_sqrt"] = q * np.sqrt(df.safety_score) * cons
    eq = (df.score_specificity + df.score_intensity + df.score_uniformity
          + df.score_accessibility + df.score_druggability) / 5.0
    out["E_equal_weights"] = eq * df.safety_score * cons
    heavy = (0.5 * df.score_specificity
             + 0.125 * (df.score_intensity + df.score_uniformity
                        + df.score_accessibility + df.score_druggability))
    out["F_specificity_heavy"] = heavy * df.safety_score * cons
    return out.apply(composite_index)


def run(ranked: pd.DataFrame | None = None,
        scored_full: pd.DataFrame | None = None) -> dict:
    if ranked is None:
        ranked = pd.read_csv(C.RESULTS / "s6_ranked_candidates.csv")
    if scored_full is None:
        merged = pd.read_csv(C.RESULTS / "s5_safety_annotated.csv").merge(
            pd.read_csv(C.RESULTS / "s2_consensus_expression.csv"), on="gene", how="inner")
        # the stability check must run on exactly the set that was ranked
        merged = merged[merged.gene.isin(set(ranked.gene))]
        from .s6_score import build_subscores
        scored_full = build_subscores(merged)

    registry = json.loads((C.RESULTS / "s1_validation_registry.json").read_text())
    n = len(ranked)
    print(f"[S7] validation against the set locked in step 1 ({n} ranked candidates)")

    rank_of = dict(zip(ranked.gene, ranked["rank"]))
    score_of = dict(zip(ranked.gene, ranked.final_score))
    tier_of = dict(zip(ranked.gene, ranked.tier))

    # ---- 1. recall of validated antigens ---------------------------------
    # detection rates for every candidate, so an unranked positive can be
    # explained rather than just left blank
    consensus = pd.read_csv(C.RESULTS / "s2_consensus_expression.csv")
    detection = dict(zip(consensus.gene, consensus.epi_detection))

    pos_rows = []
    for gene, meta in registry["positives"].items():
        row = ranked[ranked.gene == gene]
        det = detection.get(gene)
        # a reason code, so the report can render it in either language
        if len(row):
            code, status = "ranked", "ranked"
        elif det is not None and det < C.MIN_EPITHELIAL_DETECTION:
            code = "below_display_threshold"
            status = (f"not ranked: displayed on {det:.1%} of epithelial cells, below the "
                      f"{C.MIN_EPITHELIAL_DETECTION:.0%} antigen display threshold")
        elif det is None:
            code = "no_expression_data"
            status = "not ranked: no expression data in the Census feature set"
        else:
            code = "topology_gate"
            status = "not ranked: removed by the topology gate"
        pos_rows.append({
            "status_code": code,
            "status": status,
            "epi_detection": det,
            "gene": gene,
            "rationale": meta["rationale"],
            "in_search_space": meta["in_search_space"],
            "ranked": bool(len(row)),
            "rank": int(row["rank"].iloc[0]) if len(row) else None,
            "percentile": round(100.0 * (1 - (row["rank"].iloc[0] - 1) / n), 1) if len(row) else None,
            "final_score": float(row.final_score.iloc[0]) if len(row) else None,
            "tier": str(row.tier.iloc[0]) if len(row) else None,
            "tumour_quality": float(row.tumour_quality.iloc[0]) if len(row) else None,
            "safety_score": float(row.safety_score.iloc[0]) if len(row) else None,
            "worst_normal_tissue": str(row.protein_worst_tissue.iloc[0]) if len(row) else None,
        })
    pos = pd.DataFrame(pos_rows).sort_values("rank", na_position="last")
    pos.to_csv(C.RESULTS / "s7_validation_positives.csv", index=False)

    recall = {}
    for k in C.RECALL_K:
        hits = int(((pos["rank"].notna()) & (pos["rank"] <= k)).sum())
        recall[f"recall@{k}"] = f"{hits}/{len(pos)}"
        print(f"  recall@{k} = {hits}/{len(pos)}")
    ranks = pos["rank"].dropna()
    if len(ranks):
        print(f"  validated antigens rank {int(ranks.min())}-{int(ranks.max())}, "
              f"median {int(ranks.median())} of {n}")

    # ---- 2. negative controls --------------------------------------------
    neg_rows = []
    for gene, meta in registry["negative_controls"].items():
        r = rank_of.get(gene)
        if not meta["in_search_space"]:
            verdict, code = "PASS", "not_in_search_space"
            note = "never entered the search space (SURFY calls it non-surface)"
        elif r is None:
            verdict, code = "PASS", "removed_before_ranking"
            note = "removed before ranking (topology gate or no expression data)"
        elif r <= C.NEGATIVE_CONTROL_FAIL_RANK:
            verdict, code = "FAIL", "inside_fail_window"
            note = f"ranked {r} of {n}, inside the top {C.NEGATIVE_CONTROL_FAIL_RANK}"
        else:
            verdict, code = "PASS", "ranked_below_window"
            note = f"ranked {r} of {n}"
        neg_rows.append({"gene": gene, "rationale": meta["rationale"], "rank": r,
                         "final_score": score_of.get(gene), "tier": tier_of.get(gene),
                         "verdict": verdict, "reason_code": code, "note": note})
        print(f"  negative control {gene}: {verdict} -- {note}")
    neg = pd.DataFrame(neg_rows)
    neg.to_csv(C.RESULTS / "s7_negative_controls.csv", index=False)

    # ---- 3. rank stability -----------------------------------------------
    alt = _alt_scores(scored_full)
    base = scored_full.final_score
    baseline_order = scored_full.assign(s=base).sort_values("s", ascending=False)
    base_top20 = set(baseline_order.gene.head(20))
    stability_rows = []
    for name, description in ALTERNATIVES.items():
        s = alt[name]
        rho = spearmanr(base, s).statistic
        top20 = set(scored_full.assign(s=s).sort_values("s", ascending=False).gene.head(20))
        tier1 = int((s >= C.TIER1_THRESHOLD).sum())
        stability_rows.append({
            "scheme": name,
            "description": description,
            "spearman_vs_baseline": round(float(rho), 4),
            "top20_overlap": len(base_top20 & top20),
            "tier1_count": tier1,
        })
        print(f"  {name}: rho={rho:.3f}, top-20 overlap {len(base_top20 & top20)}/20")
    stab = pd.DataFrame(stability_rows)
    stab.to_csv(C.RESULTS / "s7_rank_stability.csv", index=False)

    # numbers the report needs to describe *why* the validated antigens land
    # where they land, computed here rather than asserted there
    top10 = ranked.head(10)
    pos_ranked = pos.dropna(subset=["rank"])
    mechanism = {
        "positives_ranked": int(len(pos_ranked)),
        "positives_total": int(len(pos)),
        "positive_best_rank": int(pos_ranked["rank"].min()) if len(pos_ranked) else None,
        "positive_worst_rank": int(pos_ranked["rank"].max()) if len(pos_ranked) else None,
        "positive_median_percentile_from_top": round(
            float(100 * pos_ranked["rank"].median() / n), 1) if len(pos_ranked) else None,
        "positive_median_tumour_quality": round(float(pos_ranked.tumour_quality.median()), 3)
        if len(pos_ranked) else None,
        "cohort_median_tumour_quality": round(float(ranked.tumour_quality.median()), 3),
        "positive_median_safety": round(float(pos_ranked.safety_score.median()), 3)
        if len(pos_ranked) else None,
        "cohort_median_safety": round(float(ranked.safety_score.median()), 3),
        "top10_median_safety": round(float(top10.safety_score.median()), 3),
        "top10_median_tumour_quality": round(float(top10.tumour_quality.median()), 3),
        "positives_not_ranked": {r.gene: r.status for _, r in pos.iterrows()
                                 if pd.isna(r["rank"])},
        "positives_not_ranked_codes": {r.gene: {"code": r.status_code,
                                                "epi_detection": (None if pd.isna(r.epi_detection)
                                                                  else float(r.epi_detection))}
                                       for _, r in pos.iterrows() if pd.isna(r["rank"])},
    }

    summary = {
        "n_ranked": int(n),
        "recall": recall,
        "mechanism": mechanism,
        "validated_antigen_ranks": {r.gene: (int(r["rank"]) if pd.notna(r["rank"]) else None)
                                    for _, r in pos.iterrows()},
        "negative_controls": {r.gene: r.verdict for _, r in neg.iterrows()},
        "negative_control_overall": "FAIL" if (neg.verdict == "FAIL").any() else "PASS",
        "rank_stability": stability_rows,
        "median_spearman": float(np.median([r["spearman_vs_baseline"] for r in stability_rows])),
        "min_top20_overlap": int(min(r["top20_overlap"] for r in stability_rows)),
    }
    C.write_provenance("s7", summary)
    print(f"  negative-control gate overall: {summary['negative_control_overall']}")
    return summary


if __name__ == "__main__":
    run()
