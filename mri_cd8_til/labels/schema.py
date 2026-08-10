"""The label space, and which data source can actually populate it.

The reference study splits a three-class problem into two binary ones:

======================  ==========================================  ================
Model                   Contrast                                     Needs spatial?
======================  ==========================================  ================
``RM-whole``            inflamed  vs  (desert + excluded)            no
``RM-peri``             desert    vs  excluded                       **yes**
======================  ==========================================  ================

That table is the whole reason the multi-centre plan has two tiers.  A
transcriptomic immune score measures *how much* lymphocyte signal a tumour
carries, so it can stand in for ``RM-whole``.  It cannot stand in for
``RM-peri`` at all: immune-desert and immune-excluded tumours both have an
un-infiltrated core, and differ only in where the lymphocytes sit.  Bulk
expression of a core biopsy has no way to express that difference.

Reporting an ``RM-peri`` AUC against a transcriptomic label would therefore be
measuring something other than what the model claims to measure, which is the
same class of error this project exists to avoid.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class Immunophenotype(str, Enum):
    """CD8+ T-cell spatial immunophenotype as defined in the reference study."""

    DESERT = "immune-desert"
    EXCLUDED = "immune-excluded"
    INFLAMED = "inflamed"

    @property
    def is_inflamed(self) -> bool:
        return self is Immunophenotype.INFLAMED


class LabelTask(str, Enum):
    """The two binary tasks the reference study builds."""

    WHOLE = "inflamed_vs_rest"
    PERI = "desert_vs_excluded"

    @property
    def requires_spatial_label(self) -> bool:
        return self is LabelTask.PERI


def to_binary(phenotype: Immunophenotype, task: LabelTask) -> int | None:
    """Map a three-class phenotype onto one binary task.

    Returns ``None`` when the case is not part of that task's contrast (an
    inflamed case carries no information for ``desert_vs_excluded``), so callers
    must drop it rather than silently coding it as a negative.
    """
    if task is LabelTask.WHOLE:
        return int(phenotype.is_inflamed)
    if phenotype is Immunophenotype.INFLAMED:
        return None
    return int(phenotype is Immunophenotype.EXCLUDED)


@dataclass(frozen=True)
class LabelProvenance:
    """Where a label came from, carried alongside it through the pipeline.

    Mixing label sources without tracking them is how a "multi-centre" study
    quietly becomes a study of the label-generation process.
    """

    source: str
    modality: str
    spatially_resolved: bool
    cd8_specific: bool
    supports: tuple[LabelTask, ...]
    notes: str = ""


CD8_IHC = LabelProvenance(
    source="CD8 immunohistochemistry on whole surgical sections",
    modality="pathology",
    spatially_resolved=True,
    cd8_specific=True,
    supports=(LabelTask.WHOLE, LabelTask.PERI),
    notes="The reference standard used by the study being reproduced. Essentially "
    "never available in public repositories paired with MRI.",
)

HE_TIL_MAP = LabelProvenance(
    source="Deep-learning TIL maps over H&E whole-slide images (Saltz et al. 2018)",
    modality="pathology",
    spatially_resolved=True,
    cd8_specific=False,
    supports=(LabelTask.WHOLE, LabelTask.PERI),
    notes="Counts all lymphocytes, not CD8+ specifically. Public for TCGA-BRCA, "
    "which makes it the only spatially-resolved label paired with public MRI.",
)

TRANSCRIPTOMIC = LabelProvenance(
    source="Bulk expression immune signatures (CD8A/CD8B, cytolytic, TIS)",
    modality="transcriptome",
    spatially_resolved=False,
    cd8_specific=True,
    supports=(LabelTask.WHOLE,),
    notes="Available for 717 I-SPY2 patients with paired MRI. Supports the "
    "inflamed-vs-rest contrast only -- see the module docstring.",
)
