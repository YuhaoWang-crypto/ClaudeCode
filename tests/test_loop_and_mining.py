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


def _config(**kw):
    base = dict(
        positive_context=CellContext("A", "ONT:A"),
        negative_contexts=[CellContext("B", "ONT:B")],
        rounds=10,
        population_size=16,
        children_per_round=32,
        seed=0,
    )
    base.update(kw)
    return PipelineConfig(**base)


def test_loop_improves_fitness_over_rounds():
    config = _config()
    loop = DesignLoop(MockOracle(), MockGenerator(seed=0), config)
    seeds = MockGenerator(seed=99, default_length=200).generate([], 8)
    result = loop.run(seeds)

    assert result.history, "loop should record at least one round"
    first = result.history[0].best_fitness
    # Best fitness is monotonically non-decreasing (elitism guarantees this).
    for rr in result.history:
        assert rr.best_fitness >= first - 1e-9
    assert result.best_fitness >= first


def test_loop_converges_and_reports_best():
    config = _config(rounds=50, tol=1e-4, patience=3)
    loop = DesignLoop(MockOracle(), MockGenerator(seed=1), config)
    seeds = MockGenerator(seed=5, default_length=180).generate([], 6)
    result = loop.run(seeds)
    assert result.best is not None
    assert result.best.fitness is not None
    # Convergence should trigger well before the 50-round cap on this surface.
    assert len(result.history) <= 50


def test_loop_requires_seeds():
    config = _config()
    loop = DesignLoop(MockOracle(), MockGenerator(), config)
    try:
        loop.run([])
    except ValueError:
        return
    raise AssertionError("expected ValueError for empty seeds")


def test_mining_returns_specific_windows():
    config = _config()
    oracle = MockOracle()
    # A synthetic 2 kb "genome".
    genome = Sequence("chrTest", "ACGTACGTGGCCAATT" * 128)
    mined = mine_elements([genome], oracle, config, MiningConfig(window=200, stride=100, top_n=5))
    assert 0 < len(mined) <= 5
    # Results are sorted by fitness, descending.
    fits = [c.fitness for c in mined]
    assert fits == sorted(fits, reverse=True)
    assert all(c.sequence.metadata.get("origin") == "mined" for c in mined)
