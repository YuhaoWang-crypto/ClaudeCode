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


# --------------------------------------------------------------- assembly -----
def make_segments(stimulus: str, target_cell: str, *,
                  stim_copies: int = 3, lineage: str | None = None,
                  lineage_copies: int = 2, min_promoter: str = "E1b_TATA",
                  spacer_len: int = 8):
    """Ordered list of Segments: fixed RE modules + one DESIGNABLE spacer.

    Returns (segments, designable_segment). `lineage` -> add a cell-lineage RE
    block to build a stimulus-AND-cell-type gate.
    """
    e = ELEMENTS[stimulus]
    segs = []
    segs.append(Segment(sequence=clean(FLANK5), sequence_type="dna", label="insulator5"))
    segs.append(Segment(sequence=multimerise(e["seed"], e["spacer"], stim_copies),
                        sequence_type="dna", label=f"{stimulus}_x{stim_copies}"))
    # Designable spacer: this is what the generator/optimizer searches.
    spacer = Segment(length=spacer_len, sequence_type="dna", label="designable_spacer")
    segs.append(spacer)
    if lineage:
        l = LINEAGE_ELEMENTS[lineage]
        segs.append(Segment(sequence=multimerise(l["seed"], l["spacer"], lineage_copies),
                            sequence_type="dna", label=f"lineage_{lineage}_x{lineage_copies}"))
    segs.append(Segment(sequence=clean(MINIMAL_PROMOTERS[min_promoter]),
                        sequence_type="dna", label=f"minprom_{min_promoter}"))
    segs.append(Segment(sequence=clean(KOZAK_STUB) + clean(FLANK3),
                        sequence_type="dna", label="tss_stub_3handle"))
    return segs, spacer


# --------------------------------------------------------------- scoring ------
def cell_type_constraint(construct, target_cell: str, off_cells: list[str]):
    """AlphaGenome interval-track constraint that rewards RNA_SEQ signal in the
    TARGET cell ontology while penalising OFF-target ontologies (the contrastive
    margin). Maximising this yields stimulus-responsive AND cell-selective output.
    """
    target_term = CELL_CONTEXTS[target_cell]["ontology"]
    off_terms = [CELL_CONTEXTS[c]["ontology"] for c in off_cells]
    total_len = sum(len(s.sequence or "") or (s.length or 0)
                    for s in construct.segments) if hasattr(construct, "segments") else 0
    return Constraint(
        inputs=[construct],
        function=alphagenome_interval_track_constraint,
        function_config={
            # Score the whole assembled cassette window.
            "intervals": [(0, max(total_len, 1))],
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
           evo2_model: str = "evo2_1b_base"):
    off_cells = [c for c in CELL_CONTEXTS if c != target_cell]
    segs, spacer = make_segments(stimulus, target_cell, lineage=lineage)
    construct = Construct(segs, label=f"{stimulus}__{target_cell}"
                          + (f"__AND_{lineage}" if lineage else ""))

    # Generator drives the designable spacer. Evo2 = naturalness-aware fill;
    # RandomNucleotide = cheap baseline. Both are assigned to the spacer segment.
    if use_evo2:
        # Prompt Evo2 with the fixed upstream context (guaranteed non-empty).
        fixed_ctx = "".join((s.sequence or "") for s in segs if getattr(s, "sequence", None))
        prompt = (fixed_ctx or "ACGT")[:96]
        gen = Evo2Generator(Evo2GeneratorConfig(
            model_checkpoint=evo2_model,
            prompts=[prompt], temperature=0.7, top_k=4, prepend_prompt=False))
    else:
        gen = RandomNucleotideGenerator(RandomNucleotideGeneratorConfig())
    gen.assign(spacer)

    constraint = cell_type_constraint(construct, target_cell, off_cells)

    optimizer = MCMCOptimizer(
        constructs=[construct],
        generators=[gen],
        constraints=[constraint],
        config=MCMCOptimizerConfig(num_results=num_results,
                                   proposals_per_result=16, num_steps=num_steps),
    )
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
