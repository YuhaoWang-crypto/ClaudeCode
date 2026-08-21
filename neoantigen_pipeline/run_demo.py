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
    ap.add_argument("--bench-samples", type=int, default=10)
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
    written = pipeline.write_outputs(res, a.out)

    assumptions = [
        "HLA class-I type is a declared common haplotype, not this patient's real "
        "type: TCGA HLA calls are controlled access. Supply --hla for a real run.",
        f"Tumor purity fixed at {a.purity} (no ABSOLUTE purity call is exposed for "
        "this sample), so CCF is an estimate, not a measurement.",
        "Copy number is treated as 2 everywhere; a real run should use the "
        "segment-level CN at each locus.",
        "Gene-level RSEM values stand in for transcript-level TPM from the "
        "patient's own tumor RNA-seq.",
        "Composite weights are literature-grounded defaults, not a model fitted "
        "on outcome data.",
    ]
    md = report.build_report(res, cfg, assumptions=assumptions,
                             title=f"Neoantigen selection demo - {a.sample}")

    if a.benchmark:
        print("== benchmark: mining IEDB for validated, mutation-shaped epitopes")
        pos = B.build_positive_set(proteome, lengths=tuple(a.lengths))
        print(f"   {len(pos)} validated positives across {pos['allele'].nunique()} alleles"
              if not pos.empty else "   no positives found")
        if not pos.empty:
            neg = B.build_decoy_set(pos, proteome, study=a.study,
                                    n_samples=a.bench_samples, ratio=a.bench_ratio,
                                    exclude_samples=[a.sample])
            bench = pd.concat([pos, neg], ignore_index=True)
            print(f"   {int((bench['label'] == 1).sum())} positives + "
                  f"{int((bench['label'] == 0).sum())} decoys")
            ref = [r["linear_sequence"] for r in
                   fetch.iedb_positive_epitopes(lengths=tuple(a.lengths))
                   if r.get("linear_sequence")]
            scored = B.score_benchmark(bench, ref, cfg.weight_dict())
            table = B.evaluate(scored)
            scored.to_csv(os.path.join(a.out, "benchmark_scored.csv"), index=False)
            table.to_csv(os.path.join(a.out, "benchmark_metrics.csv"), index=False)
            written += [os.path.join(a.out, "benchmark_scored.csv"),
                        os.path.join(a.out, "benchmark_metrics.csv")]
            print(table.to_string(index=False))
            md += ("\n\n## 10. Benchmark: validated neoantigens vs matched decoys "
                   "`[computed]`\n\n"
                   f"- positives: **{int((bench['label']==1).sum())}** IEDB epitopes with a "
                   "positive human T-cell assay, human source antigen, exactly one "
                   "substitution from a reference-proteome k-mer\n"
                   f"- decoys: **{int((bench['label']==0).sum())}** mutant peptides from real "
                   "TCGA-SKCM missense mutations in other patients, matched on allele "
                   "and length\n"
                   "- decoys are **unlabelled, not verified negative**, so every AUC below "
                   "is a *lower bound* on separability\n"
                   "- expression and clonality are excluded: they do not exist for an "
                   "IEDB epitope, and inventing them would inflate the result\n\n"
                   + report._md_table(table))

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
