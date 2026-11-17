"""CLI entry point.

Dry run (no GPU, mock backends):
    python -m abproto.run --spec abproto/examples/il11_spec.json --backend mock

Real run (needs Boltz installed, and a fine-tuned Multiflow checkpoint):
    python -m abproto.run --spec my_spec.json --backend boltz-multiflow \
        --multiflow-ckpt /path/to/ckpt.pt
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from .beam_search import BeamConfig
from .pipeline import DesignPipeline, PipelineConfig
from .types import DesignSpec


def load_spec(path: Path) -> DesignSpec:
    data = json.loads(path.read_text())
    coords = data.get("antigen_coords")
    return DesignSpec(
        target_name=data["target_name"],
        antigen_seq=data["antigen_seq"],
        antigen_coords=np.asarray(coords, dtype=float) if coords else None,
        epitope_resids=tuple(data.get("epitope_resids", [])),
        heavy_framework=data.get("heavy_framework", ""),
        light_framework=data.get("light_framework", ""),
        cdr_lengths=data.get("cdr_lengths", {}),
    )


def build_backends(name: str, args):
    if name == "mock":
        from .backends.mock import MockFolder, MockGenerator
        return MockGenerator(), MockFolder()
    if name == "boltz-multiflow":
        from .backends.boltz import BoltzFolder
        from .backends.multiflow import MultiflowGenerator
        return MultiflowGenerator(checkpoint=args.multiflow_ckpt), BoltzFolder()
    raise SystemExit(f"unknown backend '{name}'")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Pythia-style antibody design prototype")
    p.add_argument("--spec", required=True, type=Path)
    p.add_argument("--backend", default="mock",
                   choices=["mock", "boltz-multiflow"])
    p.add_argument("--multiflow-ckpt", default=None)
    p.add_argument("--beam-width", type=int, default=16)
    p.add_argument("--branching", type=int, default=4)
    p.add_argument("--steps", type=int, default=3)
    p.add_argument("--init-pool", type=int, default=64)
    p.add_argument("--ctmc-variants", type=int, default=3)
    p.add_argument("--max-candidates", type=int, default=48)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out", type=Path, default=None, help="write results JSON here")
    args = p.parse_args(argv)

    spec = load_spec(args.spec)
    generator, folder = build_backends(args.backend, args)
    cfg = PipelineConfig(
        beam=BeamConfig(beam_width=args.beam_width, branching=args.branching,
                        steps=args.steps, init_pool=args.init_pool, seed=args.seed),
        ctmc_variants=args.ctmc_variants, max_candidates=args.max_candidates, seed=args.seed,
    )
    pipe = DesignPipeline(generator, folder, cfg)
    results = pipe.run(spec)

    print(f"\n=== {spec.target_name}: {len(results)} candidates selected ===")
    print("Beam Search trace (mean/best confidence per step):")
    for row in pipe.trace.get("beam", []):
        print(f"  step {row['step']}: best={row['best_score']}  mean={row['mean_score']}")
    print(f"generated={pipe.trace.get('generated')}  "
          f"passed_filters={pipe.trace.get('passed_filters')}  "
          f"clusters={pipe.trace.get('n_clusters')}  selected={pipe.trace.get('selected')}")
    print("\nTop 5 candidates:")
    for s in results[:5]:
        tag = "*" if s.is_representative else " "
        print(f" {tag} {s.design.design_id}  score={s.score:.3f}  "
              f"recovery={s.recovery_rmsd:.2f}A  cluster={s.cluster_id}  "
              f"prov={s.design.provenance}")

    payload = {
        "target": spec.target_name,
        "trace": pipe.trace,
        "candidates": [s.summary() for s in results],
    }
    if args.out:
        args.out.write_text(json.dumps(payload, indent=2))
        print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
