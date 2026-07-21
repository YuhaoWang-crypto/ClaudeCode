"""
materials_discovery -- a GNoME-style crystalline-materials discovery pipeline.

Reproduces the *workflow* of Merchant et al., "Scaling deep learning for materials
discovery", Nature 2023 (doi:10.1038/s41586-023-06735-9) -- the GNoME engine:

    candidate generation (structural substitution + compositional)
      -> fast ML stability prescreen (energy / energy-above-hull)
      -> DFT verification of the top candidates (active-learning label)
      -> convex-hull decision (keep E_hull ~ 0)  -> retrain -> repeat.

Honesty split (same discipline as the sibling allosteric-biosensor skill):
  - the convex-hull math and the candidate generators are EXACT and unit-tested (offline).
  - the *energies* come from a pluggable scorer:
      * UMA (Meta FAIR fairchem) for the ML prescreen  [real, needs weights/GPU]
      * Quantum ESPRESSO on Modal for DFT verification  [real, needs Modal]
      * a deterministic surrogate for the OFFLINE demo    [⚠️ illustrative only]
"""
from .composition import Composition
from .hull import PhaseDiagram
