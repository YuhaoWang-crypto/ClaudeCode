"""
Proto-language design pipeline for stimulus-responsive, cell-type-selective
promoters. Runs on a GPU host (your Modal environment), NOT in the assembler.

Pipeline per stimulus x target-cell:
  1. Build seed cassettes from elements.py (multimerised RE + minimal promoter).
  2. Evo2  -> sequence likelihood / "naturalness" (screen out un-synthesisable,
     repeat-prone, or cryptic-splice-site sequences).
  3. AlphaGenome (or Enformer/Borzoi) -> predicted expression track for the
     cassette in the TARGET cell type, scored twice:
        (a) baseline (unstimulated proxy)
        (b) stimulated proxy = inject the pathway TF motif context / use the
            cell+treatment track if available
     Inducibility = stimulated / baseline; specificity = target_cell / off_cells.
  4. Optimise: mutate spacers, copy number, and flanks; keep Pareto-best on
     (inducibility, target expression, low off-target, high Evo2 likelihood).

IMPORTANT / honesty note
------------------------
The proto-language *core* imports below follow the published docs
(proto.evodesign.org/docs/language/installation). The exact names of the
Evo2 / AlphaGenome scoring entrypoints in `proto_tools` vary by version, so the
two `score_*` functions are written as thin ADAPTERS with a documented contract.
Wire them to whatever your installed proto-tools / Modal endpoints expose, or to
the native evo2 / alphagenome Python APIs. Do not assume the function names here
are canonical without checking `python -c "import proto_tools; help(proto_tools)"`.

Install (per docs):
    pip install git+https://github.com/evo-design/proto-language.git
    export HF_TOKEN=...        # gated model access (Evo2, AlphaGenome)
    export PROTO_HOME=...      # optional persistent cache
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Callable

# --- proto-language core (per installation docs) ------------------------------
try:
    from proto_language.core import Segment, Construct
    from proto_language.generator import RandomNucleotideGenerator
    _PROTO = True
except Exception:  # allows dry-run / assembly on a non-GPU box
    _PROTO = False
    Segment = Construct = RandomNucleotideGenerator = None  # type: ignore

from elements import ELEMENTS, MINIMAL_PROMOTERS, CELL_CONTEXTS
from build_constructs import build, clean, multimerise


# ------------------------------------------------------------------ scoring ---
# ADAPTER CONTRACT (implement against your Modal / proto-tools endpoints):
#   score_evo2(seq) -> float
#       higher = more "natural"/likely under Evo2 (e.g. mean per-base log-lik).
#   score_alphagenome(seq, cell, treatment) -> float
#       predicted RNA/CAGE expression for `seq` in `cell`, under `treatment`
#       (treatment=None -> baseline). Use the cell/assay track matching your line.

def score_evo2(seq: str) -> float:
    raise NotImplementedError(
        "Wire to Evo2: modal run evo2_score.py  OR  from evo2 import Evo2; "
        "return mean log-likelihood of `seq`.")


def score_alphagenome(seq: str, cell: str, treatment: str | None = None) -> float:
    raise NotImplementedError(
        "Wire to AlphaGenome/Enformer: return predicted expression for `seq` in "
        "`cell` (RNA-seq/CAGE head). treatment=None -> baseline track; a stimulus "
        "string -> matching treated track if your model exposes one.")


# ------------------------------------------------------------------ types -----
@dataclass
class Candidate:
    stimulus: str
    cell: str
    seq: str
    evo2: float = 0.0
    baseline: float = 0.0
    stimulated: float = 0.0
    off_target: float = 0.0
    inducibility: float = 0.0     # stimulated / baseline
    specificity: float = 0.0      # target / off-target

    def score(self) -> float:
        # Composite objective; tune weights to taste.
        return (2.0 * self.inducibility + 1.0 * self.specificity
                + 0.5 * self.stimulated + 0.25 * self.evo2)


# --------------------------------------------------------------- evaluation ---
def evaluate(seq: str, stimulus: str, target_cell: str,
             off_cells: list[str],
             evo2: Callable = score_evo2,
             ag: Callable = score_alphagenome) -> Candidate:
    base = ag(seq, target_cell, treatment=None)
    stim = ag(seq, target_cell, treatment=ELEMENTS[stimulus]["stimulus"])
    off = max((ag(seq, c, treatment=ELEMENTS[stimulus]["stimulus"])
               for c in off_cells), default=0.0)
    c = Candidate(stimulus, target_cell, seq,
                  evo2=evo2(seq), baseline=base, stimulated=stim, off_target=off)
    c.inducibility = stim / base if base > 1e-9 else stim
    c.specificity = stim / off if off > 1e-9 else stim
    return c


# ---------------------------------------------------------------- optimise ----
def seed_variants(stimulus: str, copies_range=(2, 3, 4, 5),
                  promoters=("minCMV", "E1b_TATA", "YB_TATA")):
    """Combinatorial seeds over copy number and minimal promoter."""
    for n in copies_range:
        for p in promoters:
            yield build(stimulus, n, p)["seq"]


def optimise(stimulus: str, target_cell: str, off_cells: list[str],
             rounds: int = 3, keep: int = 8,
             evo2: Callable = score_evo2, ag: Callable = score_alphagenome):
    """
    Greedy round-based optimiser. Round 0 seeds from the module grid; later
    rounds mutate spacers/flanks of the survivors. If proto-language is
    installed, RandomNucleotideGenerator drives spacer mutagenesis; otherwise a
    minimal built-in mutator is used so the loop still runs.
    """
    pool = [evaluate(s, stimulus, target_cell, off_cells, evo2, ag)
            for s in seed_variants(stimulus)]
    pool.sort(key=Candidate.score, reverse=True)
    pool = pool[:keep]

    for _ in range(rounds):
        children = []
        for parent in pool:
            for mut in _mutate(parent.seq, n=4):
                children.append(
                    evaluate(mut, stimulus, target_cell, off_cells, evo2, ag))
        pool = sorted(pool + children, key=Candidate.score, reverse=True)[:keep]
    return pool


def _mutate(seq: str, n: int = 4) -> list[str]:
    """Spacer/flank mutagenesis. Uses proto-language generator if available."""
    seq = clean(seq)
    out = []
    if _PROTO:
        gen = RandomNucleotideGenerator()
        for _ in range(n):
            # Replace a short internal window with generated neutral sequence.
            i = len(seq) // 3
            window = str(gen.generate(length=6))  # per proto-language API
            out.append(seq[:i] + window + seq[i + 6:])
    else:
        bases = "ACGT"
        for k in range(n):
            i = (len(seq) // 3) + k
            j = i % len(seq)
            out.append(seq[:j] + bases[k % 4] + seq[j + 1:])
    return out


# --------------------------------------------------------------------- main ---
def run_all(rounds: int = 2):
    """Design one optimised promoter per (stimulus, recommended target cell)."""
    results = {}
    for stim in ELEMENTS:
        targets = [c for c, v in CELL_CONTEXTS.items() if stim in v["good_for"]]
        for cell in targets:
            off = [c for c in CELL_CONTEXTS if c != cell]
            best = optimise(stim, cell, off, rounds=rounds)
            if best:
                results[(stim, cell)] = best[0]
                top = best[0]
                print(f"{stim:<18} {cell:<8} induce={top.inducibility:6.2f} "
                      f"spec={top.specificity:6.2f} evo2={top.evo2:6.3f}")
    return results


if __name__ == "__main__":
    if not _PROTO:
        print("proto-language not importable here -> dry run only.\n"
              "This script is meant to run in your Modal GPU env with Evo2 + "
              "AlphaGenome. Implement score_evo2 / score_alphagenome first.")
    else:
        run_all()
