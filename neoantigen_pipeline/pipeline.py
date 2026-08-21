"""End-to-end driver: variants in -> ranked neoantigens + mRNA construct out.

    cfg = PipelineConfig(patient=PatientConfig("PT-01", hla_class1=[...]))
    res = run_pipeline(variants_df, cfg)

`variants_df` is whatever `variants.from_cbioportal` / `variants.from_maf`
produced, already annotated with tpm and ccf.
"""

from __future__ import annotations

import os
from typing import Dict, List, Optional, Sequence

import pandas as pd

from . import features as F
from . import fetch, peptides, presentation, score, select, construct, selfindex
from . import variants as V
from .config import PipelineConfig


def run_pipeline(variant_table: pd.DataFrame, cfg: PipelineConfig,
                 proteome: Optional[Dict[str, str]] = None,
                 iedb_positive: Optional[Sequence[str]] = None,
                 build_construct: bool = True,
                 max_variants_to_predict: Optional[int] = None,
                 verbose: bool = True) -> Dict[str, object]:
    log = (lambda m: print(m, flush=True)) if verbose else (lambda m: None)
    out: Dict[str, object] = {}

    # -- 1-3. variant gates -------------------------------------------------
    v = variant_table.copy()
    if "var_id" not in v.columns:
        v["var_id"] = v["gene"].astype(str) + ":" + v["protein_change"].astype(str)
    v = V.apply_variant_gates(v, cfg.gates)
    out["variants"] = v
    out["gate_waterfall"] = V.gate_summary(v)
    keep = v[v["passes_variant_gates"]].copy()
    log(f"[1-3] {len(v)} somatic variants -> {len(keep)} pass expression/VAF/coding gates")

    if max_variants_to_predict and len(keep) > max_variants_to_predict:
        keep = keep.sort_values("tpm", ascending=False).head(max_variants_to_predict)
        log(f"      capped to the {len(keep)} highest-expressed for prediction")

    # -- 4. peptides --------------------------------------------------------
    proteome = proteome if proteome is not None else fetch.human_proteome()
    pep, skipped = peptides.generate_peptides(
        keep, proteome, lengths=cfg.patient.mhc1_lengths,
        class2_lengths=cfg.patient.mhc2_lengths)
    pep = peptides.annotate_anchor(pep)
    out["peptides"], out["peptides_skipped"] = pep, skipped
    log(f"[4]   {pep['var_id'].nunique()} variants tiled -> {len(pep)} peptides "
        f"({len(skipped)} variants skipped, see peptides_skipped)")

    # -- 5. presentation ----------------------------------------------------
    preds = presentation.predict_all(pep, cfg.patient.hla_class1,
                                     cfg.patient.hla_class2)
    out["predictions"] = preds
    joined = presentation.join_predictions(pep, preds)
    log(f"[5]   {len(preds)} peptide x allele predictions; "
        f"{int((preds['binder'] == 'strong').sum())} strong binders")

    c1 = joined[joined["mhc_class"] == "I"].copy()
    c2 = joined[joined["mhc_class"] == "II"].copy()

    # -- 6a. features + gates + score --------------------------------------
    self_kmers = selfindex.build_indexes(
        proteome, sorted(set(c1["length"].dropna().astype(int)))) if not c1.empty else {}
    if iedb_positive is None:
        try:
            rows = fetch.iedb_positive_epitopes(lengths=tuple(
                sorted(set(int(x) for x in c1["length"].dropna().unique()))) or (9,))
            iedb_positive = [r["linear_sequence"] for r in rows if r.get("linear_sequence")]
        except Exception as exc:            # noqa: BLE001 - optional feature
            log(f"      IEDB positive-epitope prior unavailable ({exc}); tcr_prior=0")
            iedb_positive = []
    out["n_iedb_reference_epitopes"] = len(iedb_positive)

    feats = F.compute_features(c1, keep, iedb_positive=iedb_positive,
                               self_kmers=self_kmers)
    feats = F.add_mhc2_support(feats, c2)
    gated = score.apply_peptide_gates(feats, cfg.gates)
    out["candidates"] = gated
    passing = gated[gated["passes"]].copy()
    log(f"[6a]  {len(gated)} candidate epitopes -> {len(passing)} pass peptide gates "
        f"({passing['var_id'].nunique()} distinct variants)")

    scored = score.composite_score(passing, cfg.weight_dict())
    best = score.best_per_variant(scored)
    out["ranked"] = score.rank_table(best)
    out["scored"] = scored

    # -- 6b. constrained selection -----------------------------------------
    sel = select.select_neoantigens(best, cfg.selection, cfg.patient.hla_class1)
    out["selected"] = sel
    out["coverage"] = select.coverage_report(sel, cfg.patient.hla_class1)
    out["selection_qc"] = select.selection_qc(sel, cfg.selection, cfg.patient.hla_class1)
    log(f"[6b]  selected {len(sel)} neoantigens across "
        f"{out['selection_qc']['alleles_covered']}/{len(cfg.patient.hla_class1)} class-I alleles")

    # -- 7. construct -------------------------------------------------------
    if build_construct and not sel.empty:
        log("[7]   building concatemer (junction-aware ordering)...")
        out["construct"] = construct.assemble(sel, proteome, cfg.patient.hla_class1,
                                              cfg.construct, verbose=verbose)
        q = out["construct"]["qc"]
        log(f"      CDS {q['length_nt']} nt, GC {q['gc_percent']}%, "
            f"QC {'PASS' if q['pass'] else 'FAIL: ' + '; '.join(q['flags'][:3])}")
    return out


def write_outputs(res: Dict[str, object], outdir: str) -> List[str]:
    os.makedirs(outdir, exist_ok=True)
    written = []
    for key in ("gate_waterfall", "ranked", "selected", "coverage",
                "peptides_skipped", "candidates"):
        obj = res.get(key)
        if isinstance(obj, pd.DataFrame) and not obj.empty:
            p = os.path.join(outdir, f"{key}.csv")
            obj.to_csv(p, index=False)
            written.append(p)
    con = res.get("construct")
    if con:
        con["minigenes"].to_csv(os.path.join(outdir, "minigenes.csv"), index=False)
        written.append(os.path.join(outdir, "minigenes.csv"))
        if isinstance(con.get("junction_scan"), pd.DataFrame) and not con["junction_scan"].empty:
            con["junction_scan"].to_csv(os.path.join(outdir, "junction_scan.csv"), index=False)
            written.append(os.path.join(outdir, "junction_scan.csv"))
        with open(os.path.join(outdir, "construct.fasta"), "w") as fh:
            fh.write(f">neoantigen_concatemer_protein len={len(con['protein'])}aa\n")
            for i in range(0, len(con["protein"]), 60):
                fh.write(con["protein"][i:i + 60] + "\n")
            fh.write(f">neoantigen_concatemer_cds len={len(con['cds'])}nt\n")
            for i in range(0, len(con["cds"]), 60):
                fh.write(con["cds"][i:i + 60] + "\n")
        written.append(os.path.join(outdir, "construct.fasta"))
    return written
