"""AlphaGenome ↔ Evo 2 integration: closed-loop regulatory element design.

A backend-agnostic active-learning framework that pairs Evo 2 (generation) with
AlphaGenome (oracle) to mine and de-novo design cell-type-specific promoters and
enhancers.

Quick start (offline, mock backends)::

    from alphagenome_evo2 import (
        CellContext, PipelineConfig, DesignLoop, MockOracle, MockGenerator, Sequence,
    )

    config = PipelineConfig(
        positive_context=CellContext("K562", "EFO:0002067"),
        negative_contexts=[CellContext("HepG2", "EFO:0001187")],
    )
    loop = DesignLoop(MockOracle(), MockGenerator(seed=0), config)
    result = loop.run([Sequence("seed", "ACGT" * 50)])
    print(result.best_fitness)

To use real models, swap ``MockOracle`` → ``AlphaGenomeOracle`` and
``MockGenerator`` → ``Evo2Generator``; everything else is unchanged.
"""

from .config import FitnessWeights, PipelineConfig
from .generator import Generator, MockGenerator
from .loop import DesignLoop, LoopResult
from .mining import MiningConfig, mine_elements
from .oracle import MockOracle, Oracle
from .scoring import apply_scores, score_candidate
from .types import (
    Candidate,
    CellContext,
    OraclePrediction,
    RoundResult,
    Sequence,
)

__version__ = "0.1.0"

__all__ = [
    "__version__",
    # config
    "PipelineConfig",
    "FitnessWeights",
    "MiningConfig",
    # backends
    "Oracle",
    "MockOracle",
    "Generator",
    "MockGenerator",
    # core
    "DesignLoop",
    "LoopResult",
    "mine_elements",
    "score_candidate",
    "apply_scores",
    # types
    "Sequence",
    "Candidate",
    "CellContext",
    "OraclePrediction",
    "RoundResult",
]
