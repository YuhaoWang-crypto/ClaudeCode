"""abproto — minimal runnable prototype of a Pythia-1-style antibody design pipeline.

Pipeline (mirrors the technical plan):

    input spec
      -> (1) all-atom joint generation          [Generator backend]
      -> (2) confidence-guided Beam Search       [test-time compute scaling]
      -> (3) CTMC-style sequence redesign        [liability-constrained variants]
      -> (4) self-consistency refolding          [pose recovery via Folder backend]
      -> (5) physics / sequence filtering + 3A clustering
      -> <= N wet-lab candidates

The *orchestration logic* (beam search, Kabsch RMSD, self-consistency, liability
filters, agglomerative clustering) is real and runs on CPU with only numpy.

The two heavy learned components are hidden behind backend protocols:
  - Generator  : Multiflow / La-Proteina  (real adapter in backends/multiflow.py)
  - Folder     : Boltz-2 / AF3            (real adapter in backends/boltz.py)

A deterministic MockBackend lets the whole pipeline run end-to-end with no GPU,
no model weights and no network, so the scaffold is testable today.
"""

from .types import DesignSpec, Design, Pose, Confidence, ScoredDesign
from .pipeline import DesignPipeline, PipelineConfig

__all__ = [
    "DesignSpec",
    "Design",
    "Pose",
    "Confidence",
    "ScoredDesign",
    "DesignPipeline",
    "PipelineConfig",
]

__version__ = "0.1.0"
