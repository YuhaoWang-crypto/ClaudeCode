"""End-to-end demo of the AlphaGenome ↔ Evo 2 loop on mock backends.

Runs the full pipeline offline (no GPUs / API keys) and prints the fitness
trajectory plus the best cell-type-specific design found.

    python examples/run_demo.py
"""

from __future__ import annotations

from alphagenome_evo2 import (
    CellContext,
    DesignLoop,
    MockGenerator,
    MockOracle,
    PipelineConfig,
    Sequence,
    mine_elements,
)
from alphagenome_evo2.mining import MiningConfig


def main() -> None:
    # Design an element ACTIVE in a cancer line (K562), SILENT in a liver line.
    config = PipelineConfig(
        positive_context=CellContext("K562", "EFO:0002067"),
        negative_contexts=[CellContext("HepG2", "EFO:0001187")],
        rounds=15,
        population_size=24,
        children_per_round=48,
        seed=0,
    )

    oracle = MockOracle()
    generator = MockGenerator(seed=0, default_length=200)

    # Step 1 — mine candidate elements from a (synthetic) genome slice.
    genome = Sequence("chrDemo", "ACGTTGACTCAGATAAGCACGTGGGCGGCCAAT" * 90)
    seeds = mine_elements(
        [genome], oracle, config, MiningConfig(window=200, stride=100, top_n=12)
    )
    print(f"mined {len(seeds)} seed elements; "
          f"best seed fitness = {seeds[0].fitness:+.4f}\n")

    # Steps 2-4 — generate / score / select in a closed loop.
    def report(rr):
        bar = "#" * int(max(0, rr.best_fitness) * 20)
        print(f"round {rr.round_index:>2}  best={rr.best_fitness:+.4f}  {bar}")

    loop = DesignLoop(oracle, generator, config, on_round=report)
    result = loop.run(seeds)

    print(f"\nconverged={result.converged}  final best fitness={result.best_fitness:+.4f}")
    best = result.best
    gap = best.fitness_terms.get("specificity_gap", 0.0)
    print(f"best design: {best.id}  (K562-vs-HepG2 specificity gap = {gap:+.3f})")
    print(best.sequence.seq)


if __name__ == "__main__":
    main()
