"""Local CLI (no Modal): run the same pipeline on a local GPU/CPU.

    python -m mdscreen.cli protein --pdb 3poz.pdb --ns 5
    python -m mdscreen.cli ligand  --ligand "CCO" --ns 5
    python -m mdscreen.cli complex --pdb 3poz.pdb --ligand ligand.sdf --ns 2
    python -m mdscreen.cli screen  --config examples/screen_config.yaml

This is handy for smoke-testing on a workstation; the production path is
``modal run modal_app.py`` which uses your Modal GPUs.
"""
from __future__ import annotations

import argparse
import sys

from .config import SimConfig, ScreenConfig
from .utils import log


def _add_sim_args(p):
    p.add_argument("--ns", type=float, default=None, help="production length (ns)")
    p.add_argument("--platform", default=None, help="CUDA | OpenCL | CPU")
    p.add_argument("--temperature", type=float, default=None)
    p.add_argument("--outdir", default="results")


def _cfg(args) -> SimConfig:
    cfg = SimConfig()
    if getattr(args, "ns", None) is not None:
        cfg.production_ns = args.ns
    if getattr(args, "platform", None):
        cfg.platform = args.platform
    if getattr(args, "temperature", None):
        cfg.temperature_k = args.temperature
    return cfg


def main(argv=None):
    parser = argparse.ArgumentParser("mdscreen")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("protein", help="protein-only MD")
    p.add_argument("--pdb", required=True)
    _add_sim_args(p)

    p = sub.add_parser("ligand", help="small-molecule-only MD")
    p.add_argument("--ligand", required=True, help="SMILES or SDF/MOL2 path")
    p.add_argument("--no-solvate", action="store_true")
    _add_sim_args(p)

    p = sub.add_parser("complex", help="protein-ligand MD + MM/GBSA")
    p.add_argument("--pdb", required=True)
    p.add_argument("--ligand", required=True)
    p.add_argument("--no-binding", action="store_true")
    _add_sim_args(p)

    p = sub.add_parser("screen", help="activity + allostery screen")
    p.add_argument("--config", required=True, help="YAML ScreenConfig")
    p.add_argument("--outdir", default="results")

    args = parser.parse_args(argv)
    from . import pipeline

    if args.cmd == "protein":
        out = pipeline.run_protein_md(args.pdb, _cfg(args), args.outdir)
    elif args.cmd == "ligand":
        out = pipeline.run_ligand_md(args.ligand, _cfg(args), args.outdir,
                                     solvate=not args.no_solvate)
    elif args.cmd == "complex":
        out = pipeline.run_complex_md(args.pdb, args.ligand, _cfg(args),
                                      args.outdir, compute_binding=not args.no_binding)
    elif args.cmd == "screen":
        import yaml
        with open(args.config) as fh:
            sc = ScreenConfig.from_dict(yaml.safe_load(fh))
        out = pipeline.run_screen(sc, args.outdir)
    else:  # pragma: no cover
        parser.error("unknown command")
        return 2

    log.info("Done. Results in %s", args.outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
