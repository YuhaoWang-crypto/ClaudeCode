"""Region-of-interest construction for whole-tumour and peritumoural models.

The reference study defines its "tumour periphery" as the **inner 2 mm rim** of
the segmented lesion -- a band *inside* the tumour boundary, not a ring of
surrounding breast tissue.  That distinction matters more than it looks:

* an **inner** rim contains only tumour voxels, so it is robust to segmentation
  leakage into fat and to breast-boundary effects, but it cannot see the
  peritumoural stroma where excluded lymphocytes actually sit;
* an **outer** ring does see that stroma, which is the biological rationale for
  the immune-excluded phenotype, but it is highly sensitive to how the lesion
  was segmented and will happily encode "how much fat is next to this tumour".

Both are implemented here, and any study comparing them should report both,
because a peritumoural model that beats a whole-tumour model on an outer ring
may simply be reading the fat/glandular background.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

import numpy as np

try:  # pragma: no cover - exercised implicitly wherever SimpleITK is installed
    import SimpleITK as sitk
except ImportError:  # pragma: no cover
    sitk = None  # type: ignore[assignment]


class RegionKind(str, Enum):
    WHOLE = "whole_tumour"
    INNER_RIM = "inner_rim"
    OUTER_RING = "outer_ring"
    CORE = "eroded_core"


@dataclass(frozen=True)
class RegionSpec:
    """How to build one ROI from a tumour mask."""

    kind: RegionKind
    width_mm: float = 2.0
    #: Only meaningful for OUTER_RING: clip the ring to a body/breast mask so it
    #: never extends into air outside the patient.
    clip_to_body: bool = True

    @property
    def name(self) -> str:
        if self.kind is RegionKind.WHOLE:
            return "whole_tumour"
        return f"{self.kind.value}_{self.width_mm:g}mm"


#: The reference study's two regions, plus the outer rings a multi-centre
#: replication should report alongside them.
DEFAULT_REGIONS: tuple[RegionSpec, ...] = (
    RegionSpec(RegionKind.WHOLE),
    RegionSpec(RegionKind.INNER_RIM, width_mm=2.0),
    RegionSpec(RegionKind.OUTER_RING, width_mm=3.0),
    RegionSpec(RegionKind.OUTER_RING, width_mm=5.0),
)


def build_region(
    mask: "sitk.Image",
    spec: RegionSpec,
    body_mask: "sitk.Image | None" = None,
) -> "sitk.Image":
    """Derive one ROI from a binary tumour mask, honouring voxel spacing.

    Morphology radii are computed per axis from the image spacing, so a 2 mm rim
    is 2 mm in all three directions even for the strongly anisotropic
    0.6 x 0.6 x 2.0 mm acquisition the reference study uses -- where an
    isotropic-in-voxels erosion would eat 3.3x more tissue in-plane than
    through-plane.
    """
    _require_sitk()
    mask = sitk.Cast(mask > 0, sitk.sitkUInt8)
    radius = _radius_in_voxels(mask, spec.width_mm)

    if spec.kind is RegionKind.WHOLE:
        return mask
    if spec.kind is RegionKind.INNER_RIM:
        eroded = sitk.BinaryErode(mask, radius)
        return sitk.Cast(mask - eroded > 0, sitk.sitkUInt8)
    if spec.kind is RegionKind.CORE:
        return sitk.BinaryErode(mask, radius)
    if spec.kind is RegionKind.OUTER_RING:
        dilated = sitk.BinaryDilate(mask, radius)
        ring = sitk.Cast(dilated - mask > 0, sitk.sitkUInt8)
        if spec.clip_to_body and body_mask is not None:
            ring = sitk.Cast(ring * sitk.Cast(body_mask > 0, sitk.sitkUInt8) > 0, sitk.sitkUInt8)
        return ring
    raise ValueError(f"unhandled region kind {spec.kind}")


def region_is_viable(region: "sitk.Image", min_voxels: int = 64) -> bool:
    """Reject ROIs too small for texture features to be defined.

    Texture matrices computed on a handful of voxels are numerically unstable
    and dominated by discretisation artefacts; small lesions eroded by a 2 mm
    rim routinely fall below this. Silently keeping them is a quiet source of
    the feature noise that later gets selected on.
    """
    _require_sitk()
    return int(np.count_nonzero(sitk.GetArrayFromImage(region))) >= min_voxels


def _radius_in_voxels(image: "sitk.Image", width_mm: float) -> list[int]:
    spacing = image.GetSpacing()
    return [max(1, int(round(width_mm / s))) for s in spacing]


def body_mask_from_image(image: "sitk.Image", percentile: float = 10.0) -> "sitk.Image":
    """Crude body/breast mask so outer rings never sample air.

    A percentile threshold plus largest-connected-component is enough for
    breast MRI, where the patient occupies a large, bright, connected region
    against a dark background.
    """
    _require_sitk()
    array = sitk.GetArrayFromImage(image).astype(float)
    threshold = float(np.percentile(array, percentile))
    binary = sitk.Cast(image > threshold, sitk.sitkUInt8)
    binary = sitk.BinaryMorphologicalClosing(binary, [3, 3, 1])
    components = sitk.ConnectedComponent(binary)
    return sitk.Cast(sitk.RelabelComponent(components, sortByObjectSize=True) == 1, sitk.sitkUInt8)


def _require_sitk() -> None:
    if sitk is None:  # pragma: no cover
        raise ImportError("SimpleITK is required for ROI construction: pip install SimpleITK")
