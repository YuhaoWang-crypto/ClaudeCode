# alphagenome-evo2

A closed-loop framework for **mining and de-novo designing cell-type-specific
regulatory elements** (promoters / enhancers) by pairing two genomic foundation
models:

- **[AlphaGenome](https://deepmind.google.com/science/alphagenome)** (Google DeepMind) as the **oracle** — a sequence-to-function model that predicts, at single-base resolution, chromatin accessibility, TF binding, gene expression and RNA splicing in a **cell-line-specific** way.
- **[Evo 2](https://github.com/ArcInstitute/evo2)** (Arc Institute) as the **generator** — a 40B-parameter genomic language model that generates novel DNA de novo and via in-context directed evolution.

The two are wired into an **active-learning loop**: Evo 2 *designs*, AlphaGenome
*judges*, and the fittest designs are fed back to Evo 2 — iterating toward
elements that are **active in a target cell line and silent everywhere else**.

```
  ┌─────────────────────────────────────────────────────────────┐
  │  (1) MINE            (2) GENERATE          (3) SCORE          │
  │  AlphaGenome scans   Evo 2 proposes        AlphaGenome rates  │
  │  the genome for  ──▶  synthetic variants ──▶ each design in   │
  │  active-in-A /       from the elite pool   every cell context │
  │  silent-in-B tiles                                │           │
  │        ▲                                          ▼           │
  │        └────────────  (4) SELECT elites, iterate ─┘           │
  └─────────────────────────────────────────────────────────────┘
```

## Why it works this way

The design goal — an element that drives strong expression in **cell line A**
(e.g. a cancer line) while staying silent in **cell line B** (e.g. normal tissue)
— is a *specificity* objective. AlphaGenome's cell-line-specific predictions make
it a natural fitness oracle for exactly this, and Evo 2's generative capacity lets
you explore synthetic sequence space far beyond natural enhancers. The loop closes
the gap between "predict" and "design".

## Install

```bash
pip install -e .            # core loop (pure Python, no heavy deps)
pip install -e ".[real]"    # + AlphaGenome and Evo 2 adapters
pip install -e ".[dev]"     # + pytest
```

The **core loop and both mock backends have zero dependencies**, so you can develop,
test, and demo entirely offline. The real models are optional extras.

## Quick start (offline, mock backends)

```python
from alphagenome_evo2 import (
    CellContext, PipelineConfig, DesignLoop, MockOracle, MockGenerator, Sequence,
)

config = PipelineConfig(
    positive_context=CellContext("K562", "EFO:0002067"),   # active here
    negative_contexts=[CellContext("HepG2", "EFO:0001187")],  # silent here
    rounds=15,
)

loop = DesignLoop(MockOracle(), MockGenerator(seed=0), config)
result = loop.run(MockGenerator(seed=99).generate([], 8))   # 8 random seeds

print(result.best_fitness, result.best.sequence.seq)
```

Run the full worked demo (mining → loop → best design, with a fitness trajectory):

```bash
python examples/run_demo.py
```

### CLI

```bash
# Mock backend (default) — runs anywhere:
alphagenome-evo2 --positive K562 --negative HepG2 --negative GM12878 \
                 --rounds 10 --json result.json

# Seed from your own candidates and use the real models:
ALPHAGENOME_API_KEY=... alphagenome-evo2 --backend real \
                 --seeds-fasta seeds.fa --positive K562
```

## Using the real models

Swap the backends; **nothing else changes**:

```python
from alphagenome_evo2.oracle import AlphaGenomeOracle      # pip install alphagenome
from alphagenome_evo2.generator import Evo2Generator        # pip install evo2 (needs a GPU)

oracle = AlphaGenomeOracle(api_key=os.environ["ALPHAGENOME_API_KEY"])
generator = Evo2Generator("evo2_7b")
loop = DesignLoop(oracle, generator, config)
```

- `AlphaGenomeOracle` maps AlphaGenome's RNA-seq / ATAC / ChIP-TF / splice-site
  tracks onto the normalised `expression / accessibility / tf_binding /
  splice_anomaly` read-outs the fitness function consumes. It accepts an injected
  `client` for testing.
- `Evo2Generator` uses each elite as an in-context prompt for Evo 2 to extend, and
  sanitises the model's output back to ACGT. It accepts an injected `model` (or a
  hosted NVIDIA BioNeMo/NIM endpoint) for testing.

Both adapters import their heavy dependency lazily, so importing the package never
requires the model to be installed.

## The fitness function

`scoring.py` turns oracle read-outs into the scalar the loop maximises:

```
fitness =  w_expr   · expression(positive)
        +  w_acc    · accessibility(positive)
        +  w_tf     · tf_binding(positive)
        -  w_off    · max_over_negatives( expression(negative) )   # worst-case leakage
        -  w_splice · max_over_contexts( splice_anomaly )          # cryptic splicing
```

Weights live in `FitnessWeights` and are fully tunable via `PipelineConfig`. Using
the **max** over negative contexts (not the mean) enforces silence in the single
worst-offending line — which is what specificity actually demands.

## Package layout

| Module | Role |
|---|---|
| `types.py` | `Sequence`, `Candidate`, `CellContext`, `OraclePrediction` — the shared currency |
| `config.py` | `PipelineConfig`, `FitnessWeights` |
| `oracle/` | `Oracle` interface, `MockOracle` (motif-grammar), `AlphaGenomeOracle` (real) |
| `generator/` | `Generator` interface, `MockGenerator` (directed evolution), `Evo2Generator` (real) |
| `mining.py` | genome-wide sliding-window discovery of specific elements (step 1) |
| `scoring.py` | cell-type-specificity fitness |
| `selection.py` | top-k / tournament / elitist next-population |
| `loop.py` | `DesignLoop` — the generate→score→select active-learning cycle (steps 2–4) |
| `cli.py` | command-line runner |

## What the mock backends are (and are not)

`MockOracle` and `MockGenerator` are **deterministic, dependency-free stand-ins**,
not biological models. The mock oracle scores sequences against per-cell-context
motif grammars (so specificity is a real, optimisable signal); the mock generator
does point mutation, crossover, and motif transplantation. They exist so the entire
pipeline is runnable and testable without GPUs or API keys, and so the demo can show
fitness climbing over rounds. For real biology, use the real adapters — the loop
code is identical.

## Tests

```bash
pytest          # 19 tests, ~0.1s, no external dependencies
```

Coverage includes sequence validation, oracle determinism, specificity scoring,
generation invariants, loop improvement/convergence, mining, and that the real
adapters construct and run against injected fakes without the heavy packages
installed.

## License

MIT.
