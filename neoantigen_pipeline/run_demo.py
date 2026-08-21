"""End-to-end public-data demo.

    python -m neoantigen_pipeline.run_demo --out demo_out
    python -m neoantigen_pipeline.run_demo --benchmark --out demo_out

Patient data: one TCGA-SKCM (cutaneous melanoma) tumor from the PanCanAtlas
cohort via the open cBioPortal API -- real somatic mutations with tumor read
counts, and the matched tumor RNA-seq abundance for every mutated gene.

HLA type is NOT available for TCGA samples through any open endpoint (the
PanImmune OptiType calls are controlled access), so the demo runs against a
declared, common class-I haplotype and says so in the report. Replace
`--hla` with the patient's real four-digit type for a real run; nothing else
in the pipeline changes.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import pandas as pd

from . import benchmark as B
from . import fetch, pipeline, report
from . import variants as V
from .config import (DEMO_HLA_CLASS1, DEMO_HLA_CLASS2, PatientConfig,
                     PipelineConfig)

STUDY = "skcm_tcga_pan_can_atlas_2018"
DEFAULT_SAMPLE = "TCGA-BF-AAP0-06"


def load_patient(study: str, sample: str, purity: float) -> pd.DataFrame:
    muts = fetch.cbio_mutations(study, sample)
    df = V.from_cbioportal(muts, sample)
    entrez = [int(e) for e in df["entrez"].dropna().unique()]
    expr = fetch.cbio_expression(study, sample, entrez)
    df = V.add_expression(df, expr)                 # RSEM normalized counts, TPM-like
    df = V.add_clonality(df, purity=purity)
    return df


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--study", default=STUDY)
    ap.add_argument("--sample", default=DEFAULT_SAMPLE)
    ap.add_argument("--hla", nargs="*", default=DEMO_HLA_CLASS1)
    ap.add_argument("--hla2", nargs="*", default=[])
    ap.add_argument("--purity", type=float, default=0.7)
    ap.add_argument("--lengths", nargs="*", type=int, default=[9, 10])
    ap.add_argument("--max-variants", type=int, default=120,
                    help="cap on variants sent to the predictor (highest-expressed first)")
    ap.add_argument("--max-neoantigens", type=int, default=34)
    ap.add_argument("--out", default="demo_out")
    ap.add_argument("--no-construct", action="store_true")
    ap.add_argument("--benchmark", action="store_true",
                    help="also run the IEDB-validated-neoantigen vs TCGA-decoy benchmark")
    ap.add_argument("--bench-ratio", type=int, default=5)
    ap.add_argument("--bench-samples", type=int, default=25)
    ap.add_argument("--bench-pool-ratio", type=int, default=250,
                    help="decoy pool size per positive for the presentation-controlled benchmark")
    ap.add_argument("--bench-min-allele-n", type=int, default=10,
                    help="alleles with fewer validated positives are dropped from benchmark B")
    ap.add_argument("--tesla", action="store_true",
                    help="also benchmark on the public TESLA mirror (real assay labels)")
    ap.add_argument("--no-figures", action="store_true")
    a = ap.parse_args(argv)

    os.makedirs(a.out, exist_ok=True)
    print(f"== loading {a.sample} from {a.study} (cBioPortal, open API)")
    variants = load_patient(a.study, a.sample, a.purity)
    print(f"   {len(variants)} somatic variants, "
          f"{int(variants['variant_class'].eq('Missense_Mutation').sum())} missense")

    print("== loading human reference proteome (UniProt, reviewed)")
    proteome = fetch.human_proteome()
    print(f"   {len(proteome)} canonical protein sequences")

    cfg = PipelineConfig(patient=PatientConfig(
        patient_id=a.sample, hla_class1=list(a.hla), hla_class2=list(a.hla2),
        tumor_purity=a.purity, mhc1_lengths=tuple(a.lengths),
        mhc2_lengths=(15,) if a.hla2 else ()))
    cfg.selection.max_neoantigens = a.max_neoantigens

    res = pipeline.run_pipeline(variants, cfg, proteome=proteome,
                                build_construct=not a.no_construct,
                                max_variants_to_predict=a.max_variants)
    written = pipeline.write_outputs(
        res, a.out, cfg=cfg,
        proteome_path=os.path.join(os.path.dirname(fetch.__file__), 'data',
                                   'cache', 'uniprot_human_reviewed.fasta.gz'))

    assumptions = [
        "HLA class-I type is a declared common haplotype, not this patient's real "
        "type: TCGA HLA calls are controlled access. Supply --hla for a real run.",
        f"Tumor purity fixed at {a.purity} (no ABSOLUTE purity call is exposed for "
        "this sample), so CCF is an estimate, not a measurement.",
        "Copy number is treated as 2 everywhere; a real run should use the "
        "segment-level CN at each locus.",
        "Gene-level RSEM values stand in for transcript-level TPM from the "
        "patient's own tumor RNA-seq.",
        "Composite weights are presentation-dominant defaults, chosen by hand "
        "after two benchmarks found the other peptide-intrinsic features at or "
        "below baseline. They are not a model fitted on outcome data, and the "
        "evidence behind them is 35 positives across 6 patients.",
    ]
    md = report.build_report(res, cfg, assumptions=assumptions,
                             title=f"Neoantigen selection demo - {a.sample}")

    bench_scored = None
    bench_a_scored = None
    if a.benchmark:
        print("== benchmark: mining IEDB for validated, mutation-shaped epitopes")
        pos = B.build_positive_set(proteome, lengths=tuple(a.lengths))
        print(f"   {len(pos)} validated positives across {pos['allele'].nunique()} alleles"
              if not pos.empty else "   no positives found")
        if not pos.empty:
            ref = [r["linear_sequence"] for r in
                   fetch.iedb_positive_epitopes(lengths=tuple(a.lengths))
                   if r.get("linear_sequence")]

            # A. allele/length-matched decoys only
            neg = B.build_decoy_set(pos, proteome, study=a.study,
                                    n_samples=a.bench_samples, ratio=a.bench_ratio,
                                    exclude_samples=[a.sample])
            bench_a = pd.concat([pos, neg], ignore_index=True)
            print(f"   [A] allele-matched: {int((bench_a['label']==1).sum())} positives "
                  f"+ {int((bench_a['label']==0).sum())} decoys")
            scored_a = B.score_benchmark(bench_a, ref, cfg.weight_dict())
            table_a = B.evaluate(scored_a)
            print(table_a.to_string(index=False))

            # B. presentation-controlled: restrict to the alleles that actually
            #    have enough validated positives, draw decoys from a much larger
            #    pool, and balance the two classes within each binding stratum.
            counts = pos["allele"].value_counts()
            keep_alleles = list(counts[counts >= a.bench_min_allele_n].index)
            pos_b = pos[pos["allele"].isin(keep_alleles)].copy()
            print(f"   [B] presentation-controlled on {len(keep_alleles)} alleles "
                  f"({', '.join(keep_alleles)}): {len(pos_b)} of {len(pos)} positives")
            pool = B.build_decoy_pool(pos_b, proteome, study=a.study,
                                      n_samples=a.bench_samples,
                                      pool_ratio=a.bench_pool_ratio,
                                      exclude_samples=[a.sample])
            print(f"   [B] decoy pool: {len(pool)} peptides, predicting ranks")
            pos_r = B.add_ranks(pos_b, verbose=False)
            pool_r = B.add_ranks(pool, verbose=False)
            bench_b = B.balance_by_stratum(pd.concat([pos_r, pool_r], ignore_index=True),
                                           ratio=a.bench_ratio)
            print(f"   [B] stratum-balanced: {int((bench_b['label']==1).sum())} positives "
                  f"+ {int((bench_b['label']==0).sum())} decoys")
            scored_b = B.score_benchmark(bench_b, ref, cfg.weight_dict())
            table_b = B.evaluate(scored_b)
            print(table_b.to_string(index=False))
            strat = B.stratified_evaluate(
                scored_b, ["score_netmhcpan_only", "score_composite_no_tcr",
                           "score_composite", "feat_agretopicity",
                           "feat_hydrophobicity", "feat_dissimilarity"])
            print(strat.to_string(index=False))

            bench_scored = scored_b
            bench_a_scored = scored_a
            for name, obj in (("benchmark_scored_allele_matched", scored_a),
                              ("benchmark_metrics_allele_matched", table_a),
                              ("benchmark_scored_presentation_controlled", scored_b),
                              ("benchmark_metrics_presentation_controlled", table_b),
                              ("benchmark_metrics_by_stratum", strat)):
                p = os.path.join(a.out, f"{name}.csv")
                obj.to_csv(p, index=False)
                written.append(p)

            def _rank_summary(df):
                p = df[df["label"] == 1]["mut_rank"].median()
                n = df[df["label"] == 0]["mut_rank"].median()
                return f"median %rank: positives {p:.2f}, decoys {n:.2f}"

            md += (
                "\n\n## 10. Benchmark: validated neoantigens vs decoys `[computed]`\n\n"
                f"- positives: **{int((pos['label']==1).sum())}** IEDB epitopes with a "
                "positive human T-cell assay, human source antigen, exactly one "
                "substitution from a reference-proteome k-mer (that difference is "
                "mutation-shaped and hands back the wild-type counterpart, so "
                "agretopicity is exact)\n"
                "- decoys: mutant peptides from real TCGA-SKCM missense mutations in "
                "**other** patients\n"
                "- decoys are **unlabelled, not verified negative** -> every AUC is a "
                "*lower bound* on separability\n"
                "- expression and clonality are excluded: they do not exist for an "
                "IEDB epitope, and inventing them would inflate the result\n\n"
                "### A. decoys matched on allele and length only\n\n"
                f"_{_rank_summary(scored_a)}_\n\n"
                + report._md_table(table_a) +
                "\n**Read this one as a warning, not a win.** IEDB epitopes were largely "
                "*discovered* because they bind well, so positives and random decoys "
                "differ mostly in predicted binding. A near-perfect AUC here measures "
                "the binding predictor, not the selection layer.\n\n"
                "### B. presentation-controlled\n\n"
                f"Restricted to the alleles with at least {a.bench_min_allele_n} "
                f"validated positives ({', '.join(keep_alleles)}), with decoys drawn "
                f"from a pool of {len(pool)} peptides and then subsampled so that each "
                "predicted-%rank stratum holds the same class mixture. "
                f"_{_rank_summary(scored_b)}_\n\n"
                + report._md_table(table_b) +
                "\nBinding strength is now (approximately) equalized between the "
                "classes, so an AUC above 0.5 here cannot be explained by presentation "
                "alone. This is the question the selection layer exists to answer -- "
                "among peptides that all reach the cell surface, which does a T cell "
                "see? -- and it is a much harder problem than benchmark A makes it "
                "look.\n\n"
                "**Within each binding stratum** (NaN = too few of one class in that "
                "stratum to report a number worth trusting):\n\n"
                + report._md_table(strat))

    if a.tesla:
        from . import tesla as T
        print("== TESLA mirror: real T-cell assay labels, real negatives")
        td = T.score(T.load(proteome=proteome), cfg.weight_dict())
        tt = T.evaluate(td)
        td.to_csv(os.path.join(a.out, "tesla_scored.csv"), index=False)
        tt.to_csv(os.path.join(a.out, "tesla_metrics.csv"), index=False)
        pp = T.per_patient_table(td, "score_netmhcpan_only", 20)
        pp.to_csv(os.path.join(a.out, "tesla_per_patient.csv"), index=False)
        written += [os.path.join(a.out, f) for f in
                    ("tesla_scored.csv", "tesla_metrics.csv", "tesla_per_patient.csv")]
        print(tt.to_string(index=False))
        md += (
            "\n\n## 11. TESLA mirror `[computed]`\n\n"
            f"522 peptide-HLA pairs, **{int(td['label'].sum())} experimentally immunogenic**, "
            f"{td['patient'].nunique()} patients. Unlike section 10, a label of 0 here means "
            "*assayed and negative*, not *never tested* - so these are real metrics, not "
            "lower bounds. Base rate "
            f"**{tt.attrs['baseline_AP']:.3f}**; average precision (PR-AUC) is the metric "
            "that moves at this class imbalance, and top-N recovery per patient is what a "
            "fixed vaccine budget actually cares about.\n\n"
            "The published per-model columns shipped with the mirror are scored on the "
            "identical rows, so the comparison is like for like.\n\n"
            + report._md_table(tt, max_rows=20) +
            "\nPer-patient recovery in a 20-slot budget, ranked by presentation:\n\n"
            + report._md_table(pp))

    if not a.no_figures:
        try:
            from . import figures
            p1 = figures.summary_figure(res, cfg, os.path.join(a.out, "summary.png"),
                                        bench=bench_scored)
            written.append(p1)
            p2 = figures.junction_figure(res, os.path.join(a.out, "junctions.png"))
            if p2:
                written.append(p2)
            md += ("\n\n## Figures\n\n"
                   "![summary](summary.png)\n\n![junctions](junctions.png)\n")
            if bench_a_scored is not None and bench_scored is not None:
                p3 = figures.benchmark_figure(bench_a_scored, bench_scored,
                                              os.path.join(a.out, "benchmark.png"))
                written.append(p3)
                md += "\n![benchmark](benchmark.png)\n"
        except Exception as exc:            # noqa: BLE001 - figures are optional
            print(f"   (figures skipped: {exc})")

    path = os.path.join(a.out, "REPORT.md")
    with open(path, "w") as fh:
        fh.write(md)
    written.append(path)
    with open(os.path.join(a.out, "selection_qc.json"), "w") as fh:
        json.dump(res.get("selection_qc", {}), fh, indent=2, default=str)

    print("\n== wrote:")
    for w in written:
        print("   ", w)
    return 0


if __name__ == "__main__":
    sys.exit(main())
