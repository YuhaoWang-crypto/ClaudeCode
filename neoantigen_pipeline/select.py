"""Step 6b: constrained selection of the <=34 neoantigens that go into the mRNA.

Ranking alone is not a selection. A payload of 34 slots is a budget, and the
public product description fixes its size but not how to spend it. The
constraints implemented here are the ones that a vaccine designer actually
argues about:

  * cap per gene           -- do not spend 6 slots on 6 epitopes of one gene
  * class-I allele spread  -- an escape variant that loses one HLA allele should
                              not silence the whole vaccine
  * clonal-first           -- subclonal neoantigens are present in only part of
                              the tumor (McGranahan 2016)
  * forced inclusions      -- a known driver mutation can be mandated
  * de-duplication         -- two variants whose best epitopes are the same
                              peptide get one slot

Greedy with a coverage bonus: at each step take the candidate that maximizes
score + lambda * (marginal allele coverage gained). This is the standard
submodular-coverage heuristic; it is deterministic and auditable.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence

import pandas as pd


def select_neoantigens(ranked: pd.DataFrame, rules, alleles: Sequence[str],
                       coverage_lambda: float = 0.15) -> pd.DataFrame:
    """-> selected rows in payload order, with `slot` and `why_selected`."""
    if ranked.empty:
        return ranked.assign(slot=[], why_selected=[])

    d = ranked.copy().sort_values("neo_score", ascending=False)
    alleles = list(alleles)
    gene_count: Dict[str, int] = {}
    allele_count: Dict[str, int] = {a: 0 for a in alleles}
    used_peptides = set()
    chosen: List[dict] = []

    forced = set(rules.force_include_genes or ())

    def _eligible(row) -> Optional[str]:
        g = row.get("gene")
        if gene_count.get(g, 0) >= rules.max_per_gene:
            return f"gene cap ({rules.max_per_gene}) reached for {g}"
        if row.get("mut_peptide") in used_peptides:
            return "duplicate epitope already in payload"
        if rules.per_allele_cap and allele_count.get(row.get("allele"), 0) >= rules.per_allele_cap:
            return f"allele cap reached for {row.get('allele')}"
        return None

    def _take(row, why):
        gene_count[row.get("gene")] = gene_count.get(row.get("gene"), 0) + 1
        allele_count[row.get("allele")] = allele_count.get(row.get("allele"), 0) + 1
        used_peptides.add(row.get("mut_peptide"))
        rec = row.to_dict()
        rec["why_selected"] = why
        chosen.append(rec)

    # 1. forced driver inclusions first
    for _, row in d.iterrows():
        if len(chosen) >= rules.max_neoantigens:
            break
        if row.get("gene") in forced and _eligible(row) is None:
            _take(row, f"forced inclusion (driver gene {row.get('gene')})")

    # 2. greedy with coverage bonus
    taken_ids = {c["var_id"] for c in chosen}
    pool = d[~d["var_id"].isin(taken_ids)].copy()
    while len(chosen) < rules.max_neoantigens and not pool.empty:
        covered = {a for a, n in allele_count.items() if n > 0}
        bonus = pool["allele"].map(lambda a: coverage_lambda if a not in covered else 0.0)
        clonal_bonus = 0.0
        if rules.prefer_clonal and "ccf" in pool:
            clonal_bonus = (pool["ccf"].fillna(0) >= 0.8).astype(float) * 0.05
        pool = pool.assign(_eff=pool["neo_score"] + bonus + clonal_bonus)
        pool = pool.sort_values("_eff", ascending=False)
        picked = None
        for i, row in pool.iterrows():
            reason = _eligible(row)
            if reason is None:
                picked = (i, row)
                break
        if picked is None:
            break
        i, row = picked
        why = "top-ranked"
        if row["allele"] not in covered:
            why = f"top-ranked + first epitope for {row['allele']}"
        _take(row, why)
        pool = pool.drop(index=i)

    out = pd.DataFrame(chosen)
    if out.empty:
        return out
    out = out.sort_values("neo_score", ascending=False).reset_index(drop=True)
    out.insert(0, "slot", range(1, len(out) + 1))
    return out


def coverage_report(selected: pd.DataFrame, alleles: Sequence[str],
                    rules=None) -> pd.DataFrame:
    rows = []
    for a in alleles:
        n = int((selected["allele"] == a).sum()) if not selected.empty else 0
        rows.append({"allele": a, "n_epitopes": n})
    df = pd.DataFrame(rows)
    return df.sort_values("n_epitopes", ascending=False).reset_index(drop=True)


def selection_qc(selected: pd.DataFrame, rules, alleles: Sequence[str]) -> Dict[str, object]:
    n = len(selected)
    cov = coverage_report(selected, alleles)
    covered = int((cov["n_epitopes"] > 0).sum())
    clonal = float((selected["ccf"].fillna(0) >= 0.8).mean()) if "ccf" in selected and n else float("nan")
    genes = selected["gene"].nunique() if n else 0
    return {
        "n_selected": n,
        "meets_min": n >= rules.min_neoantigens,
        "meets_max": n <= rules.max_neoantigens,
        "alleles_covered": covered,
        "meets_allele_spread": covered >= min(rules.min_alleles_covered, len(alleles)),
        "unique_genes": int(genes),
        "fraction_clonal": clonal,
        "median_rank": float(selected["mut_rank"].median()) if n else float("nan"),
        "n_strong_binders": int((selected["binder"] == "strong").sum()) if n else 0,
    }
