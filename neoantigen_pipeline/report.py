"""Markdown report generation, with explicit rigour labels.

Every number that comes out of this pipeline is one of three things and the
report says which:

  [computed]   produced by a real predictor / real data in this run
  [assumed]    a configuration choice standing in for data we do not have
               (tumor purity, an HLA type we could not obtain, a weight)
  [unverified] a prediction with no experimental validation behind it --
               which is *every* immunogenicity call this pipeline makes
"""

from __future__ import annotations

from typing import Dict, List, Optional

import pandas as pd


def _md_table(df: pd.DataFrame, max_rows: int = 40, floatfmt: str = "{:.3g}") -> str:
    if df is None or df.empty:
        return "_(empty)_\n"
    d = df.head(max_rows).copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda v: "" if pd.isna(v) else floatfmt.format(v))
    head = "| " + " | ".join(map(str, d.columns)) + " |"
    sep = "| " + " | ".join(["---"] * len(d.columns)) + " |"
    rows = ["| " + " | ".join("" if pd.isna(v) else str(v) for v in r) + " |"
            for r in d.itertuples(index=False)]
    extra = f"\n_({len(df)} rows total, showing {min(max_rows, len(df))})_\n" if len(df) > max_rows else ""
    return "\n".join([head, sep] + rows) + "\n" + extra


def build_report(res: Dict[str, object], cfg, assumptions: Optional[List[str]] = None,
                 title: str = "Personalized neoantigen selection") -> str:
    p = cfg.patient
    sel = res.get("selected", pd.DataFrame())
    qc = res.get("selection_qc", {})
    con = res.get("construct")
    L: List[str] = []

    L.append(f"# {title}\n")
    L.append(f"**Patient / sample:** `{p.patient_id}`  \n"
             f"**Class-I HLA:** {', '.join(p.hla_class1) or '_none_'}  \n"
             f"**Class-II HLA:** {', '.join(p.hla_class2) or '_none_'}  \n"
             f"**Tumor purity assumed:** {p.tumor_purity}\n")

    L.append("\n## 0. What is and is not real here\n")
    L.append("| label | meaning |\n| --- | --- |")
    L.append("| `[computed]` | produced in this run by a real predictor on real data |")
    L.append("| `[assumed]` | a configuration stand-in for data not available |")
    L.append("| `[unverified]` | a prediction with no experimental confirmation |\n")
    if assumptions:
        L.append("**Assumptions in force for this run:**\n")
        for a in assumptions:
            L.append(f"- `[assumed]` {a}")
        L.append("")

    L.append("\n## 1-3. Variant gating `[computed]`\n")
    L.append(_md_table(res.get("gate_waterfall")))

    skipped = res.get("peptides_skipped")
    if isinstance(skipped, pd.DataFrame) and not skipped.empty:
        counts = skipped["reason"].value_counts().reset_index()
        counts.columns = ["reason variants were not tiled", "n"]
        L.append("\n**Variants that could not be turned into peptides** "
                 "(reported, not silently dropped):\n")
        L.append(_md_table(counts))

    L.append("\n## 4-5. Peptides and HLA presentation `[computed]`\n")
    preds = res.get("predictions")
    if isinstance(preds, pd.DataFrame) and not preds.empty:
        by = preds.groupby("binder").size().reset_index(name="n")
        L.append(f"- peptides tiled: **{len(res.get('peptides', []))}**")
        L.append(f"- peptide x allele predictions: **{len(preds)}** "
                 f"(NetMHCpan-4.1 EL via the IEDB cloud API)")
        L.append(f"- IEDB positive-assay epitopes used as the TCR prior: "
                 f"**{res.get('n_iedb_reference_epitopes', 0)}**\n")
        L.append(_md_table(by))

    L.append("\n## 6. Ranking `[computed]` / `[unverified]`\n")
    L.append("Weights actually applied:\n")
    wd = cfg.weight_dict()
    L.append(_md_table(pd.DataFrame({"feature": list(wd), "weight": list(wd.values())})))
    ranked = res.get("ranked")
    if isinstance(ranked, pd.DataFrame) and not ranked.empty:
        cols = [c for c in ["rank", "gene", "protein_change", "allele", "mut_peptide",
                            "mut_rank", "wt_rank", "tpm", "ccf", "neo_score"] if c in ranked]
        L.append("\nTop candidates:\n")
        L.append(_md_table(ranked[cols], max_rows=20))

    L.append("\n## 7. Selected payload `[computed]` under `[assumed]` rules\n")
    L.append(f"- selected: **{qc.get('n_selected')}** "
             f"(rules: {cfg.selection.min_neoantigens}-{cfg.selection.max_neoantigens}, "
             f"max {cfg.selection.max_per_gene} per gene)")
    L.append(f"- class-I alleles covered: **{qc.get('alleles_covered')}** / {len(p.hla_class1)}")
    L.append(f"- unique genes: **{qc.get('unique_genes')}**, "
             f"strong binders: **{qc.get('n_strong_binders')}**, "
             f"median %rank: **{qc.get('median_rank')}**")
    fc = qc.get("fraction_clonal")
    if fc == fc:
        L.append(f"- clonal fraction (CCF>=0.8): **{fc:.0%}**")
    L.append("")
    if isinstance(sel, pd.DataFrame) and not sel.empty:
        cols = [c for c in ["slot", "gene", "protein_change", "allele", "mut_peptide",
                            "mut_rank", "tpm", "ccf", "neo_score", "why_selected"] if c in sel]
        L.append(_md_table(sel[cols], max_rows=40))
    L.append("\nPer-allele coverage:\n")
    L.append(_md_table(res.get("coverage")))

    if con:
        L.append("\n## 8. mRNA construct `[computed]`\n")
        jn = con.get("junction_scan")
        n_flag = int(jn["flagged"].sum()) if isinstance(jn, pd.DataFrame) and not jn.empty else 0
        L.append(f"- minigene length: {cfg.construct.epitope_length} aa, "
                 f"linker: `{cfg.construct.linker or 'none (direct fusion)'}`")
        L.append(f"- junction cost, input order -> optimized order: "
                 f"**{con.get('junction_cost_naive'):.3f} -> "
                 f"{con.get('junction_cost_optimized'):.3f}**")
        nb, ob = con.get("junction_binders_naive"), con.get("junction_binders_optimized")
        if nb is not None:
            L.append(f"- junction binders (<= {cfg.construct.junction_scan_rank}% rank) "
                     f"at the optimized lengths {con.get('junction_cost_lengths')}: "
                     f"**{nb} -> {ob}** by reordering alone")
        L.append(f"- final rescan over lengths {con.get('junction_scan_lengths')}: "
                 f"**{n_flag}** junction peptides bind at <= "
                 f"{cfg.construct.junction_scan_rank}% rank")
        if isinstance(jn, pd.DataFrame) and not jn.empty:
            by_len = (jn[jn["flagged"]].groupby("length").size()
                      .reset_index(name="flagged_peptides"))
            if not by_len.empty:
                L.append("")
                L.append("  Flagged by peptide length (lengths outside the "
                         "optimization objective are reported, not silently dropped):")
                L.append("")
                L.append(_md_table(by_len))
        q = con["qc"]
        L.append(f"- CDS: **{q['length_nt']} nt** ({q['length_aa']} aa), "
                 f"GC **{q['gc_percent']}%**, uridine **{q['uridine_percent']}%**")
        L.append(f"- synthesis QC: **{'PASS' if q['pass'] else 'FAIL'}**"
                 + ("" if q["pass"] else " - " + "; ".join(q["flags"])))
        L.append(f"- codon repairs applied: {len(con.get('codon_repairs', []))}\n")
        if isinstance(jn, pd.DataFrame) and not jn.empty:
            top = jn.sort_values("percentile_rank").head(10)
            L.append("Strongest junction peptides (these are *not* tumor-specific "
                     "and are the reason ordering matters):\n")
            L.append(_md_table(top))

    L.append("\n## 9. Limits `[unverified]`\n")
    L.append("- Every immunogenicity number here is a **prediction**. Nothing in this "
             "run has been tested against a patient's T cells.")
    L.append("- Presentation prediction is the strongest link in the chain; "
             "TCR recognition is the weakest. Published prospective benchmarks put "
             "the positive rate of top-ranked predicted neoantigens in the low tens "
             "of percent, not near 100%.")
    L.append("- Frameshift / neo-ORF and fusion neoantigens need transcript-level "
             "annotation; variants that needed it are listed above rather than dropped.")
    L.append("- The composite weights are an editable literature-grounded default, "
             "not a fitted model, unless you refit them with `benchmark.py`.")
    return "\n".join(L)
