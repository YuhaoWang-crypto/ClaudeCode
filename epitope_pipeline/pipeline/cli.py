"""Command-line entry point for the multi-model epitope pipeline.

Example
-------
    python -m pipeline.cli \
        --fasta examples/example.fasta \
        --models netmhcpan,netmhciipan,netctlpan,bepipred \
        --alleles-i HLA-A02:01,HLA-A01:01,HLA-B07:02 \
        --alleles-ii DRB1_0101,DRB1_0401 \
        --out results/

Tools that are not installed are skipped (not fatal), so you can dry-run the
whole orchestration before obtaining the DTU academic licences.
"""

from __future__ import annotations

import argparse
import os
import sys

from . import aggregate
from .models import REGISTRY, get_models
from .runner import run_models


def _csv(s: str):
    return [x.strip() for x in s.split(",") if x.strip()]


def _lens(s):
    return tuple(int(x) for x in _csv(s)) if s else None


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="epitope-pipeline",
        description="Run one protein through many DTU epitope/immunogenicity "
                    "models in parallel and merge the results.")
    p.add_argument("--fasta", required=True, help="input protein FASTA")
    p.add_argument("--models", default="netmhcpan,netmhciipan,netctlpan,bepipred",
                   help="comma list; known: " + ",".join(REGISTRY))
    p.add_argument("--alleles-i", default="HLA-A02:01,HLA-A01:01,HLA-B07:02",
                   help="MHC-I alleles (NetMHCpan/NetCTLpan)")
    p.add_argument("--alleles-ii", default="DRB1_0101,DRB1_0401",
                   help="MHC-II alleles (NetMHCIIpan)")
    p.add_argument("--lengths-i", default="8,9,10,11", help="MHC-I peptide lengths")
    p.add_argument("--lengths-ii", default="15", help="MHC-II peptide lengths")
    p.add_argument("--out", default="results", help="output directory")
    p.add_argument("--workers", type=int, default=8, help="parallel workers")
    p.add_argument("--timeout", type=int, default=3600, help="per-model timeout (s)")
    p.add_argument("--top", type=int, default=25, help="top-N consensus epitopes")
    p.add_argument("--binary", action="append", default=[],
                   metavar="name=path",
                   help="override a tool's executable, e.g. netMHCpan=/opt/nmp/netMHCpan")
    p.add_argument("--list-models", action="store_true",
                   help="print the model registry and exit")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    if args.list_models:
        print(f"{'key':<14}{'name':<18}{'category':<14}{'binary':<18}homepage")
        for k, m in REGISTRY.items():
            print(f"{k:<14}{m.name:<18}{m.category:<14}{m.binary:<18}{m.homepage}")
        return 0

    if not os.path.isfile(args.fasta):
        print(f"error: FASTA not found: {args.fasta}", file=sys.stderr)
        return 2

    try:
        specs = get_models(_csv(args.models))
    except KeyError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    overrides = {}
    for item in args.binary:
        if "=" in item:
            name, path = item.split("=", 1)
            overrides[name.strip()] = path.strip()

    os.makedirs(args.out, exist_ok=True)
    workdir = os.path.join(args.out, "work")

    results = run_models(
        specs, args.fasta,
        alleles_i=_csv(args.alleles_i), alleles_ii=_csv(args.alleles_ii),
        lengths_i=_lens(args.lengths_i), lengths_ii=_lens(args.lengths_ii),
        workdir=workdir, binary_overrides=overrides,
        max_workers=args.workers, timeout=args.timeout)

    long_csv = aggregate.write_long_table(results, os.path.join(args.out, "all_predictions.csv"))
    cons = aggregate.consensus(results, args.top)
    cons_csv = aggregate.write_consensus(results, os.path.join(args.out, "consensus_epitopes.csv"), args.top)
    summ = aggregate.write_summary_json(results, os.path.join(args.out, "run_summary.json"))

    aggregate.print_report(results, cons)
    print("\nOutputs:")
    print(f"  - {long_csv}")
    print(f"  - {cons_csv}")
    print(f"  - {summ}")

    ran = sum(1 for r in results if r.status == "ok")
    return 0 if ran or not specs else 0


if __name__ == "__main__":
    raise SystemExit(main())
