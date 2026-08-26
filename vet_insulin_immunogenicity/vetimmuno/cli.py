"""Command line entry point: ``python -m vetimmuno run --species dog``."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="vetimmuno",
        description="MHC-II immunogenicity workflow for veterinary insulins "
                    "(canine DLA / feline FLA).",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="run the workflow for one species")
    run.add_argument("--species", required=True,
                     help="dog | cat | path to a species YAML config")
    run.add_argument("--outdir", type=Path, default=None)
    run.add_argument("--backend", choices=("surrogate", "netmhciipan"), default="surrogate",
                     help="surrogate = untrained illustrative scorer (default); "
                          "netmhciipan = licensed NetMHCIIpan-4.3 binary via $NETMHCIIPAN")
    run.add_argument("--background", type=int, default=20000,
                     help="background peptides per molecule for %%rank calibration")
    run.add_argument("--offline", action="store_true",
                     help="fail rather than fetch anything not already cached")

    val = sub.add_parser("validate", help="run only the validation harness")
    val.add_argument("--species", required=True)
    val.add_argument("--offline", action="store_true")

    args = parser.parse_args(argv)
    from . import workflow

    if args.command == "run":
        result = workflow.run(args.species, outdir=args.outdir,
                              backend_name=args.backend, n_background=args.background,
                              offline=args.offline)
        vr = result["validation"]
        print(f"\nwrote {result['outdir']}")
        print(f"validation: {vr.passed}/{len(vr.checks)} passed")
        for check in vr.failed:
            print(f"  FAIL {check.id} {check.name}: {check.observed} "
                  f"(expected {check.expected})")
        return 1 if vr.failed else 0

    if args.command == "validate":
        result = workflow.run(args.species, offline=args.offline)
        vr = result["validation"]
        for check in vr.checks:
            print(f"{check.status:5s} {check.id:9s} {check.name}")
            print(f"                observed: {check.observed}")
        return 1 if vr.failed else 0

    return 2


if __name__ == "__main__":
    sys.exit(main())
