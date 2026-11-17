"""Backend protocols. Anything implementing these plugs into the pipeline."""

from __future__ import annotations

from typing import Optional, Protocol, runtime_checkable

from ..types import Design, DesignSpec, Pose


@runtime_checkable
class Generator(Protocol):
    """Produces candidate antibody-antigen complexes conditioned on a spec.

    Real implementation: Multiflow / La-Proteina joint sequence+structure flow.
    """

    def generate(self, spec: DesignSpec, n: int, seed: Optional[int] = None) -> list[Design]:
        """Return ``n`` fresh candidate complexes."""
        ...

    def refine(self, parents: list[Design], spec: DesignSpec, n_per_parent: int,
               seed: Optional[int] = None) -> list[Design]:
        """Return children sampled in the neighbourhood of ``parents``.

        Used by Beam Search for test-time compute scaling: each surviving
        candidate seeds ``n_per_parent`` perturbed offspring. A trivial
        implementation may ignore parents and call ``generate``.
        """
        ...


@runtime_checkable
class Folder(Protocol):
    """Folds / refolds a sequence pair into a complex with confidence heads.

    Real implementation: Boltz-2 or AlphaFold3.
    """

    def fold(self, design: Design, spec: DesignSpec, seed: Optional[int] = None) -> Pose:
        """Unconditionally refold the design's sequences against the antigen.

        Returns the predicted complex geometry plus pLDDT/PAE/PDE confidence.
        This is the self-consistency check: does the sequence fold back to the
        intended binding pose?
        """
        ...
