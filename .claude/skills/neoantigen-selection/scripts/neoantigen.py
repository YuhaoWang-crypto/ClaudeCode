#!/usr/bin/env python3
"""Command-line entry point for the neoantigen-selection skill.

Works from anywhere:

    python scripts/neoantigen.py selftest
    python scripts/neoantigen.py demo --out demo_out --benchmark --tesla
    python scripts/neoantigen.py run --maf patient.maf --patient PT-014 \
        --hla HLA-A*02:01 HLA-A*24:02 HLA-B*07:02 HLA-B*44:02 \
        --expression tpm.csv --purity 0.62 --out PT-014_out

The `neoantigen_pipeline` package is looked for next to this file (the
distributed skill vendors it there) and then in the repository this skill was
developed in, so the same script runs in both layouts.
"""

from __future__ import annotations

import argparse
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
for _candidate in (_HERE,                                     # scripts/neoantigen_pipeline
                   os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))):
    if os.path.isdir(os.path.join(_candidate, "neoantigen_pipeline")):
        if _candidate not in sys.path:
            sys.path.insert(0, _candidate)
        break
else:                                                          # noqa: PLW0120
    sys.exit("cannot find the neoantigen_pipeline package next to this script "
             "or in the parent repository -- is the skill unpacked completely?")


def cmd_selftest(_args) -> int:
    from neoantigen_pipeline import selftest
    return selftest.main()


def cmd_demo(args) -> int:
    from neoantigen_pipeline import run_demo
    argv = ["--out", args.out]
    if args.benchmark:
        argv.append("--benchmark")
    if args.tesla:
        argv.append("--tesla")
    if args.sample:
        argv += ["--sample", args.sample]
    if args.max_neoantigens:
        argv += ["--max-neoantigens", str(args.max_neoantigens)]
    if args.no_figures:
        argv.append("--no-figures")
    return run_demo.main(argv)


def cmd_tesla(args) -> int:
    """Benchmark the shipped scoring on the public TESLA mirror alone."""
    import pandas as pd
    from neoantigen_pipeline import fetch, tesla
    from neoantigen_pipeline.config import PatientConfig, PipelineConfig

    cfg = PipelineConfig(patient=PatientConfig("tesla-benchmark"))
    d = tesla.score(tesla.load(proteome=fetch.human_proteome()), cfg.weight_dict())
    table = tesla.evaluate(d)
    os.makedirs(args.out, exist_ok=True)
    d.to_csv(os.path.join(args.out, "tesla_scored.csv"), index=False)
    table.to_csv(os.path.join(args.out, "tesla_metrics.csv"), index=False)
    tesla.per_patient_table(d, "score_netmhcpan_only", 20).to_csv(
        os.path.join(args.out, "tesla_per_patient.csv"), index=False)
    print(f"random baseline AP = {table.attrs['baseline_AP']}")
    print(table.to_string(index=False))
    print(f"\nwrote {args.out}/tesla_metrics.csv")
    return 0


def cmd_run(args) -> int:
    """Real patient run: a somatic MAF plus an HLA type."""
    import pandas as pd
    from neoantigen_pipeline import fetch, pipeline, report
    from neoantigen_pipeline import variants as V
    from neoantigen_pipeline.config import PatientConfig, PipelineConfig

    var = V.from_maf(args.maf, sample_id=args.patient)
    if args.expression:
        expr = pd.read_csv(args.expression)
        cols = {c.lower(): c for c in expr.columns}
        key = cols.get("entrez") or cols.get("entrez_gene_id") or list(expr.columns)[0]
        val = cols.get("tpm") or list(expr.columns)[1]
        var = V.add_expression(var, dict(zip(expr[key].astype(int), expr[val].astype(float))))
    else:
        print("WARNING: no --expression given. The expression gate cannot run; "
              "every candidate will be treated as expressed, which is exactly the "
              "failure mode this pipeline exists to prevent.", file=sys.stderr)
        var["tpm"] = float("nan")
    var = V.add_clonality(var, purity=args.purity)

    cfg = PipelineConfig(patient=PatientConfig(
        patient_id=args.patient, hla_class1=list(args.hla),
        hla_class2=list(args.hla2 or []), tumor_purity=args.purity))
    cfg.selection.max_neoantigens = args.max_neoantigens
    if args.force_genes:
        cfg.selection.force_include_genes = tuple(args.force_genes)
    if args.no_expression_gate:
        cfg.gates.min_tpm = 0.0

    res = pipeline.run_pipeline(var, cfg, max_variants_to_predict=args.max_variants)
    written = pipeline.write_outputs(res, args.out, cfg=cfg, input_path=args.maf)
    md = report.build_report(res, cfg, title=f"Neoantigen selection - {args.patient}")
    path = os.path.join(args.out, "REPORT.md")
    with open(path, "w") as fh:
        fh.write(md)
    print("\nwrote:")
    for w in written + [path]:
        print("   ", w)
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="neoantigen", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)

    s = sub.add_parser("selftest", help="43 offline checks, no network")
    s.set_defaults(func=cmd_selftest)

    d = sub.add_parser("demo", help="public-data end-to-end demo (TCGA-SKCM)")
    d.add_argument("--out", default="demo_out")
    d.add_argument("--sample", default=None)
    d.add_argument("--benchmark", action="store_true")
    d.add_argument("--tesla", action="store_true")
    d.add_argument("--max-neoantigens", type=int, default=None)
    d.add_argument("--no-figures", action="store_true")
    d.set_defaults(func=cmd_demo)

    t = sub.add_parser("tesla", help="benchmark the scoring on the public TESLA mirror")
    t.add_argument("--out", default="tesla_out")
    t.set_defaults(func=cmd_tesla)

    r = sub.add_parser("run", help="run on a real patient MAF")
    r.add_argument("--maf", required=True, help="somatic MAF from your tumor/normal caller")
    r.add_argument("--patient", required=True)
    r.add_argument("--hla", nargs="+", required=True,
                   help="four-digit class-I type from the NORMAL sample")
    r.add_argument("--hla2", nargs="*", default=None)
    r.add_argument("--expression", default=None,
                   help="CSV with entrez,tpm from the patient's tumor RNA-seq")
    r.add_argument("--purity", type=float, default=0.7)
    r.add_argument("--max-neoantigens", type=int, default=34)
    r.add_argument("--max-variants", type=int, default=None)
    r.add_argument("--force-genes", nargs="*", default=None)
    r.add_argument("--no-expression-gate", action="store_true",
                   help="explicitly disable the expression gate (not recommended)")
    r.add_argument("--out", default="patient_out")
    r.set_defaults(func=cmd_run)

    a = p.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
