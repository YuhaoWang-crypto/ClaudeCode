"""
Proto-language design pipeline for stimulus-responsive, cell-type-selective
promoters. Written against the REAL proto_language API (v0.1.0), verified from
the installed package source (Segment/Construct/Constraint/Generator/Optimizer/
Program and the alphagenome_interval_track_constraint signature).

Runs on a GPU host with the proto_tools execution layer installed (your Modal
env). proto_tools is GitHub-only (evo-design/proto-tools) and provides the
Evo2 + AlphaGenome/Enformer/Borzoi wrappers; it is NOT on PyPI, so this file
imports lazily and degrades to a printable design plan when it is absent (as in
a CPU-only / sandbox box).

Architecture per promoter (ordered Segments -> Construct):

    [5' insulator]  fixed
    [ stimulus RE module ] xN     fixed (from elements.ELEMENTS[...].seed)
    [ spacer ]                    DESIGNABLE (length=k) -> generator fills it
    [ lineage RE module ] xM      fixed, ONLY for stimulus-AND-cell-type gates
    [ minimal promoter ]          fixed (elements.MINIMAL_PROMOTERS[...])
    [ Kozak/TSS stub + 3' handle] fixed

The DESIGNABLE spacer(s) are what Evo2 / the optimizer search over; the fixed RE
modules pin the biology. Scoring uses AlphaGenome interval-track signal with
`ontology_terms` = target cell and `contrastive_ontology_terms` = off-target
cells -> directly optimises cell-type SELECTIVITY (the AND gate), on top of the
stimulus responsiveness that the RE modules provide.

Install (Modal / GPU):
    pip install proto-language           # language layer (this box has it)
    pip install "git+https://github.com/evo-design/proto-tools.git"   # execution
    export HF_TOKEN=...                   # gated Evo2 / AlphaGenome weights
    python design_pipeline.py --stimulus interferon_typeII --target THP1
"""

from __future__ import annotations
import argparse

from elements import (ELEMENTS, LINEAGE_ELEMENTS, MINIMAL_PROMOTERS,
                      CELL_CONTEXTS)
from build_constructs import FLANK5, FLANK3, KOZAK_STUB, clean, multimerise

# proto_language imports its execution layer (proto_tools) at package import
# time, so guard the whole thing.
try:
    from proto_language.core import Segment, Construct, Constraint
    from proto_language import (
        Program,
        Evo2Generator, Evo2GeneratorConfig,
        RandomNucleotideGenerator, RandomNucleotideGeneratorConfig,
        alphagenome_interval_track_constraint,
    )
    from proto_language.optimizer import MCMCOptimizer, MCMCOptimizerConfig
    _READY = True
except Exception as exc:  # proto_tools missing / CPU-only box
    _READY = False
    _IMPORT_ERR = exc


# Deterministic neutral filler (balanced GC, 47-mer prime unit to avoid periodic
# motif artifacts) used to pad the AlphaGenome context window around the cassette.
_NEUTRAL_UNIT = "ACTGACAGTCGATCAGTACGATCGTACAGTCAGATCGTACGATCAGTA"  # 48 bp


def neutral_flank(n: int) -> str:
    if n <= 0:
        return ""
    reps = (n // len(_NEUTRAL_UNIT)) + 1
    return (_NEUTRAL_UNIT * reps)[:n]


# --------------------------------------------------------------- assembly -----
def make_segments(stimulus: str, target_cell: str, *,
                  stim_copies: int = 3, lineage: str | None = None,
                  lineage_copies: int = 2, min_promoter: str = "E1b_TATA",
                  spacer_len: int = 8, context_length: int = 16384):
    """Ordered list of Segments: fixed RE modules + one DESIGNABLE spacer.

    Returns (segments, designable_segment). `lineage` -> add a cell-lineage RE
    block to build a stimulus-AND-cell-type gate.
    """
    e = ELEMENTS[stimulus]
    stim_block = multimerise(e["seed"], e["spacer"], stim_copies)
    lin_block = (multimerise(LINEAGE_ELEMENTS[lineage]["seed"],
                             LINEAGE_ELEMENTS[lineage]["spacer"], lineage_copies)
                 if lineage else "")
    minp = clean(MINIMAL_PROMOTERS[min_promoter])
    tail = clean(KOZAK_STUB) + clean(FLANK3)
    f5 = clean(FLANK5)

    segs = [Segment(sequence=f5, sequence_type="dna", label="insulator5"),
            Segment(sequence=stim_block, sequence_type="dna", label=f"{stimulus}_x{stim_copies}")]
    spacer = Segment(length=spacer_len, sequence_type="dna", label="designable_spacer")
    segs.append(spacer)
    if lineage:
        segs.append(Segment(sequence=lin_block, sequence_type="dna",
                            label=f"lineage_{lineage}_x{lineage_copies}"))
    segs.append(Segment(sequence=minp, sequence_type="dna", label=f"minprom_{min_promoter}"))
    segs.append(Segment(sequence=tail, sequence_type="dna", label="tss_stub_3handle"))

    # AlphaGenome scores a fixed-size window (context_length, e.g. 16384 bp): the
    # designable spacer must be embedded with left+right context summing to
    # context_length - spacer_len. The functional cassette sits immediately beside
    # the spacer; the rest is neutral filler.
    upstream = f5 + stim_block                       # real cassette left of spacer
    downstream = lin_block + minp + tail             # real cassette right of spacer
    need = max(context_length - spacer_len, len(upstream) + len(downstream))
    left_len = need // 2
    right_len = need - left_len
    left_ctx = neutral_flank(max(left_len - len(upstream), 0)) + upstream
    right_ctx = downstream + neutral_flank(max(right_len - len(downstream), 0))
    meta = {
        "prompt_ctx": upstream,                      # Evo2 prompt = fixed upstream
        "left_ctx": left_ctx[-left_len:] if left_len else upstream,
        "right_ctx": right_ctx[:right_len] if right_len else downstream,
        "total_len": len(f5) + len(stim_block) + spacer_len + len(lin_block)
                     + len(minp) + len(tail),
        "spacer_len": spacer_len,
        "context_length": context_length,
    }
    return segs, spacer, meta


# --------------------------------------------------------------- scoring ------
def cell_type_constraint(spacer, target_cell: str, off_cells: list[str],
                         left_ctx: str, right_ctx: str, interval_len: int):
    """AlphaGenome interval-track constraint on the DESIGNABLE spacer, embedded in
    the fixed cassette (left/right context). Rewards RNA_SEQ signal in the TARGET
    cell ontology while penalising OFF-target ontologies (contrastive margin) ->
    stimulus-responsive AND cell-selective output.
    """
    target_term = CELL_CONTEXTS[target_cell]["ontology"]
    off_terms = [CELL_CONTEXTS[c]["ontology"] for c in off_cells]
    return Constraint(
        inputs=[spacer],                       # a Segment (has num_proposals)
        function=alphagenome_interval_track_constraint,
        function_config={
            "intervals": [(0, max(interval_len, 1))],   # within the spacer segment
            "left_context": left_ctx,
            "right_context": right_ctx,
            "requested_output": "RNA_SEQ",
            "ontology_terms": [target_term],
            "contrastive_ontology_terms": off_terms or None,
            "direction": "maximize",
            "maximize_inflection_value": 5.0,
            "margin_inflection_value": 0.0,   # reward any positive target-vs-off margin
            "margin_sigmoid_scale": 1.0,
            "context_length": 16384,          # pad to AlphaGenome input window
        },
    )


# --------------------------------------------------------------- optimise -----
def design(stimulus: str, target_cell: str, *, lineage: str | None = None,
           num_steps: int = 20, num_results: int = 4, use_evo2: bool = True,
           evo2_model: str = "evo2_1b_base", build_only: bool = False,
           contrastive: bool = True, spacer_len: int = 8):
    # contrastive=False -> maximize target-cell expression only (single ontology);
    # True -> also penalize off-target cells (the cell-type-specificity objective).
    # spacer_len = size of the designable proximal region (the AlphaGenome-scored,
    # Evo2-designed segment). Larger -> more leverage over cell-type specificity.
    off_cells = [c for c in CELL_CONTEXTS if c != target_cell] if contrastive else []
    segs, spacer, meta = make_segments(stimulus, target_cell, lineage=lineage,
                                       spacer_len=spacer_len)
    construct = Construct(segs, label=f"{stimulus}__{target_cell}"
                          + (f"__AND_{lineage}" if lineage else ""))

    # Generator drives the designable spacer. Evo2 = naturalness-aware fill;
    # RandomNucleotide = cheap baseline. Both are assigned to the spacer segment.
    if use_evo2:
        prompt = (meta["prompt_ctx"] or "ACGT")[:96]
        gen = Evo2Generator(Evo2GeneratorConfig(
            model_checkpoint=evo2_model,
            prompts=[prompt], temperature=0.7, top_k=4, prepend_prompt=False))
    else:
        gen = RandomNucleotideGenerator(RandomNucleotideGeneratorConfig())
    gen.assign(spacer)

    constraint = cell_type_constraint(spacer, target_cell, off_cells,
                                      meta["left_ctx"], meta["right_ctx"],
                                      meta["spacer_len"])

    optimizer = MCMCOptimizer(
        constructs=[construct],
        generators=[gen],
        constraints=[constraint],
        config=MCMCOptimizerConfig(num_results=num_results,
                                   proposals_per_result=16, num_steps=num_steps),
    )
    # build_only stops before the GPU-bound run so the wiring can be validated
    # on a cheap CPU container.
    if build_only:
        return {"ok": True, "total_len": meta["total_len"], "prompt": prompt if use_evo2 else None,
                "n_segments": len(segs), "spacer_len": meta["spacer_len"],
                "construct_label": construct.label}
    Program(optimizers=[optimizer], num_results=num_results).run()
    return optimizer


# --------------------------------------------------------------------- CLI ----
def _plan_only(stimulus, target_cell, lineage):
    """Print the concrete design plan when proto_tools is unavailable."""
    segs, spacer = ([], None)
    e = ELEMENTS[stimulus]
    off = [CELL_CONTEXTS[c]["ontology"] for c in CELL_CONTEXTS if c != target_cell]
    print(f"proto_tools (execution layer) not importable here:\n  {_IMPORT_ERR}\n")
    print("=== DESIGN PLAN (run this on your Modal GPU box) ===")
    print(f"stimulus         : {stimulus}  [{e['element']} / {e['tf']}]")
    print(f"target cell      : {target_cell}  ontology={CELL_CONTEXTS[target_cell]['ontology']}")
    if lineage:
        l = LINEAGE_ELEMENTS[lineage]
        print(f"AND lineage gate : {lineage}  [{l['tf']}]  ontology={l['ontology']}")
    print(f"generator        : Evo2Generator over a designable spacer segment")
    print(f"scorer           : alphagenome_interval_track_constraint")
    print(f"  ontology_terms            = ['{CELL_CONTEXTS[target_cell]['ontology']}']")
    print(f"  contrastive_ontology_terms= {off}")
    print(f"  direction='maximize', requested_output='RNA_SEQ'")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stimulus", default="interferon_typeII", choices=list(ELEMENTS))
    ap.add_argument("--target", default="THP1", choices=list(CELL_CONTEXTS))
    ap.add_argument("--lineage", default=None, choices=[None, *LINEAGE_ELEMENTS])
    ap.add_argument("--no-evo2", action="store_true")
    args = ap.parse_args()

    if not _READY:
        _plan_only(args.stimulus, args.target, args.lineage)
        return
    opt = design(args.stimulus, args.target, lineage=args.lineage,
                 use_evo2=not args.no_evo2)
    print("Top designs:")
    for i, seg in enumerate(getattr(opt, "segments", [])):
        print(f"  [{i}] {getattr(seg, 'sequence', '')}")


if __name__ == "__main__":
    main()
