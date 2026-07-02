"""Command-line entry point for the design loop.

Runs the full mine → generate → score → select pipeline. By default it uses the
offline mock backends so it works with no dependencies or API keys; pass
``--backend real`` to use AlphaGenome + Evo 2 (requires the packages, a GPU for
Evo 2, and ``ALPHAGENOME_API_KEY``).

Examples::

    python -m alphagenome_evo2.cli --positive K562 --negative HepG2 --rounds 10
    python -m alphagenome_evo2.cli --seeds-fasta seeds.fa --json result.json
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from typing import List, Optional

from .config import CellContext, PipelineConfig  # re-exported for convenience
from .generator import MockGenerator
from .loop import DesignLoop
from .oracle import MockOracle
from .types import RoundResult, Sequence


def _read_fasta(path: str) -> List[Sequence]:
    seqs: List[Sequence] = []
    name: Optional[str] = None
    buf: List[str] = []

    def flush() -> None:
        if name is not None and buf:
            seqs.append(Sequence(id=name, seq="".join(buf)))

    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                name = line[1:].split()[0] or f"seq{len(seqs)}"
                buf = []
            else:
                buf.append(line)
    flush()
    return seqs


def _default_seeds(n: int, length: int, seed: int) -> List[Sequence]:
    """Generate random seed sequences (used when no FASTA is supplied)."""
    gen = MockGenerator(seed=seed, default_length=length)
    return gen.generate([], n)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="alphagenome-evo2", description=__doc__)
    p.add_argument("--positive", default="K562", help="cell line where the element should be ACTIVE")
    p.add_argument("--positive-ontology", default=None, help="ontology term for the positive context")
    p.add_argument("--negative", action="append", default=[], help="cell line(s) that must stay SILENT (repeatable)")
    p.add_argument("--rounds", type=int, default=8)
    p.add_argument("--population", type=int, default=32)
    p.add_argument("--children", type=int, default=64)
    p.add_argument("--temperature", type=float, default=0.8)
    p.add_argument("--mutation-rate", type=float, default=0.03)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--seeds-fasta", default=None, help="FASTA of seed sequences (else random seeds)")
    p.add_argument("--n-seeds", type=int, default=8, help="number of random seeds when no FASTA is given")
    p.add_argument("--seed-length", type=int, default=200)
    p.add_argument("--backend", choices=["mock", "real"], default="mock")
    p.add_argument("--json", default=None, help="write the full result as JSON to this path")
    p.add_argument("--top", type=int, default=5, help="how many final designs to print")
    return p


def _build_backends(args):
    if args.backend == "mock":
        return MockOracle(), MockGenerator(seed=args.seed, default_length=args.seed_length)
    # Real backends: imported lazily so mock runs need no heavy deps.
    from .oracle import AlphaGenomeOracle
    from .generator import Evo2Generator

    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    return AlphaGenomeOracle(api_key=api_key), Evo2Generator()


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)

    config = PipelineConfig(
        positive_context=CellContext(args.positive, args.positive_ontology),
        negative_contexts=[CellContext(n) for n in args.negative],
        rounds=args.rounds,
        population_size=args.population,
        children_per_round=args.children,
        temperature=args.temperature,
        mutation_rate=args.mutation_rate,
        seed=args.seed,
    )

    oracle, generator = _build_backends(args)

    if args.seeds_fasta:
        seeds = _read_fasta(args.seeds_fasta)
        if not seeds:
            print(f"no sequences found in {args.seeds_fasta}", file=sys.stderr)
            return 2
    else:
        seeds = _default_seeds(args.n_seeds, args.seed_length, args.seed)

    def report(rr: RoundResult) -> None:
        print(
            f"round {rr.round_index:>2}: best={rr.best_fitness:+.4f} "
            f"mean_elite={rr.mean_elite_fitness:+.4f} "
            f"generated={rr.n_generated} best_id={rr.best_candidate_id}",
            file=sys.stderr,
        )

    loop = DesignLoop(oracle, generator, config, on_round=report)
    result = loop.run(seeds)

    print(f"\nconverged={result.converged}  best_fitness={result.best_fitness:+.4f}")
    print(f"\nTop {args.top} designs:")
    for i, cand in enumerate(result.population[: args.top]):
        gap = cand.fitness_terms.get("specificity_gap", 0.0)
        print(f"  {i + 1}. {cand.id}  fitness={cand.fitness:+.4f}  "
              f"specificity_gap={gap:+.3f}  len={len(cand.sequence)}")
        print(f"     {cand.sequence.seq}")

    if args.json:
        _dump_json(result, config, args.json)
        print(f"\nwrote {args.json}")
    return 0


def _dump_json(result, config, path: str) -> None:
    payload = {
        "positive_context": config.positive_context.name,
        "negative_contexts": [c.name for c in config.negative_contexts],
        "converged": result.converged,
        "best": {
            "id": result.best.id,
            "sequence": result.best.sequence.seq,
            "fitness": result.best.fitness,
            "fitness_terms": dict(result.best.fitness_terms),
        },
        "history": [vars(r) for r in result.history],
        "population": [
            {
                "id": c.id,
                "sequence": c.sequence.seq,
                "fitness": c.fitness,
                "fitness_terms": dict(c.fitness_terms),
            }
            for c in result.population
        ],
    }
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2)


if __name__ == "__main__":
    raise SystemExit(main())
