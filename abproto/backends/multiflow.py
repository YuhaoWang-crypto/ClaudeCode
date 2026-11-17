"""Real Generator adapter around Multiflow (MIT) joint sequence+structure flow.

Multiflow (Campbell & Yim et al., ICML 2024, https://github.com/jasonkyuyim/multiflow)
co-generates backbone frames and sequence via a continuous structure flow plus a
discrete-sequence CTMC. This adapter documents how to drive it; the actual
sampling call depends on the checkpoint you train/fine-tune for antibody CDR
co-design against a conditioning antigen.

Because stock Multiflow is monomer co-design, turning it into a conditional
antibody generator means FINE-TUNING on SAbDab with the antigen + germline
framework as context (route B in the plan). Until that checkpoint exists, this
class raises a clear NotImplementedError pointing at the integration steps; use
MockGenerator for dry runs.

Integration steps (route B):
  1. git clone https://github.com/jasonkyuyim/multiflow && install its env
  2. fine-tune on SAbDab Ab-Ag complexes, framework+antigen as fixed context,
     CDR residues as the flow/CTMC-generated region
  3. implement ``_sample`` below to call the fine-tuned model and return
     per-chain CA coords + sequences + a confidence placeholder (real confidence
     comes from the Folder stage).
"""

from __future__ import annotations

from typing import Optional

from ..types import Design, DesignSpec


class MultiflowGenerator:
    def __init__(self, checkpoint: Optional[str] = None, device: str = "cuda"):
        self.checkpoint = checkpoint
        self.device = device
        self._model = None

    def _load(self):
        if self._model is not None:
            return
        try:
            import torch  # noqa: F401
        except ImportError as e:  # pragma: no cover
            raise RuntimeError("PyTorch required for MultiflowGenerator") from e
        raise NotImplementedError(
            "Wire up the fine-tuned Multiflow checkpoint here. See module docstring. "
            "For a dry run use abproto.backends.mock.MockGenerator."
        )

    def _sample(self, spec: DesignSpec, n: int, parents, seed) -> list[Design]:  # pragma: no cover
        self._load()
        # model.sample(...) -> frames + sequence logits -> argmax/sample
        # build Design(backbone={H,L,A}, heavy_seq, light_seq, cdr_spans, ...)
        raise NotImplementedError

    def generate(self, spec: DesignSpec, n: int, seed: Optional[int] = None) -> list[Design]:
        return self._sample(spec, n, None, seed)

    def refine(self, parents: list[Design], spec: DesignSpec, n_per_parent: int,
               seed: Optional[int] = None) -> list[Design]:
        # Beam-search refinement = partial-noise resampling around each parent
        # (keep framework + accepted CDRs, re-noise H3). Implement via the
        # model's inpainting / partial-diffusion entry point.
        return self._sample(spec, len(parents) * n_per_parent, parents, seed)
