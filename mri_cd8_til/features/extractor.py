"""IBSI-oriented radiomic feature extraction reproducing the reference design.

Produces, per region per DCE phase, the reference study's ``1 + T + 12*T``
layout, where ``T`` is the size of the texture set.

**On the texture-set size.**  The reference study reports ``T = 64``, from
GLCM + GLRLM + GLSZM + NGTDM.  PyRadiomics' IBSI-aligned implementation of those
four classes yields ``T = 61`` (24 + 16 + 16 + 5).  The gap is three GLCM
features: the classic 27-feature Haralick set includes entries PyRadiomics
deliberately removes as mathematically redundant -- Sum Average, for instance,
is exactly ``2 x Joint Average`` for a symmetric GLCM and carries no independent
information.

So the honest statement is that ``T = 61`` here versus ``T = 64`` there, that the
difference is three redundant features, and that it has no bearing on any
conclusion.  Anyone reporting a reproduction should say which number they used
rather than quietly matching the published one.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from .wavelet_bank import PAPER_SUBBANDS, WaveletFilter, PAPER_FILTER_BANK, expected_feature_count

try:  # pragma: no cover
    import SimpleITK as sitk
    from radiomics import featureextractor
except ImportError:  # pragma: no cover
    sitk = None  # type: ignore[assignment]
    featureextractor = None  # type: ignore[assignment]

logging.getLogger("radiomics").setLevel(logging.ERROR)

TEXTURE_CLASSES: tuple[str, ...] = ("glcm", "glrlm", "glszm", "ngtdm")

#: Measured against PyRadiomics 3.0.1: glcm=24, glrlm=16, glszm=16, ngtdm=5.
IBSI_TEXTURE_COUNT = 61
#: As reported by the reference study (classic 27-feature Haralick GLCM).
PAPER_TEXTURE_COUNT = 64


@dataclass
class ExtractionConfig:
    """Extraction parameters.

    ``bin_width`` rather than ``bin_count`` is not a detail.  A fixed bin *count*
    makes the bin width depend on each lesion's own intensity range, so the same
    tissue discretised on two scanners with different gain lands in different
    bins -- re-introducing the scanner dependence that preprocessing removed.
    IBSI recommends fixed bin width for arbitrary-unit modalities such as MR,
    and multi-centre work needs it.
    """

    bin_width: float = 5.0
    #: Voxel array shift keeps discretisation stable when normalisation has
    #: produced negative intensities.
    voxel_array_shift: int = 300
    texture_classes: tuple[str, ...] = TEXTURE_CLASSES
    include_first_order: bool = False
    include_shape: bool = False
    include_volume_feature: bool = True
    filter_bank: tuple[WaveletFilter, ...] = PAPER_FILTER_BANK
    keep_subbands: tuple[str, ...] = PAPER_SUBBANDS
    force_2d: bool = False
    label: int = 1
    #: Extra settings passed straight through to PyRadiomics.
    extra_settings: dict[str, Any] = field(default_factory=dict)

    def base_settings(self) -> dict[str, Any]:
        settings: dict[str, Any] = {
            "binWidth": self.bin_width,
            "voxelArrayShift": self.voxel_array_shift,
            "label": self.label,
            "force2D": self.force_2d,
            # Resampling is handled upstream in imaging.preprocess so that the
            # mask and every DCE phase land on one shared grid.
            "resampledPixelSpacing": None,
            "normalize": False,
        }
        settings.update(self.extra_settings)
        return settings


@dataclass
class FeatureVector:
    """Extracted features for one case / region / phase."""

    case_id: str
    region: str
    phase: int
    values: dict[str, float]

    @property
    def n_features(self) -> int:
        return len(self.values)

    def prefixed(self) -> dict[str, float]:
        """Feature names namespaced by region and phase, ready for a wide table."""
        return {f"{self.region}|p{self.phase}|{k}": v for k, v in self.values.items()}


class RadiomicFeatureExtractor:
    """Extracts the ``1 + T + 12*T`` feature block for one region and phase."""

    def __init__(self, config: ExtractionConfig | None = None) -> None:
        self.config = config or ExtractionConfig()
        _require_radiomics()

    # ------------------------------------------------------------------ public

    def expected_n_features(self, texture_count: int = IBSI_TEXTURE_COUNT) -> int:
        return expected_feature_count(
            texture_count,
            full_3d=len(self.config.keep_subbands) > len(PAPER_SUBBANDS),
            include_volume=self.config.include_volume_feature,
        )

    def extract(
        self,
        image: "sitk.Image",
        mask: "sitk.Image",
        case_id: str,
        region: str,
        phase: int,
    ) -> FeatureVector:
        values: dict[str, float] = {}

        if self.config.include_volume_feature:
            values["shape_VoxelVolume"] = self._voxel_volume(mask)

        values.update(self._extract_original(image, mask))
        for wavelet_filter in self.config.filter_bank:
            values.update(self._extract_wavelet(image, mask, wavelet_filter))

        return FeatureVector(case_id=case_id, region=region, phase=phase, values=values)

    def extract_all_phases(
        self,
        phase_images: list["sitk.Image"],
        mask: "sitk.Image",
        case_id: str,
        region: str,
    ) -> dict[str, float]:
        """Flat, region- and phase-prefixed feature dict across all DCE phases."""
        out: dict[str, float] = {}
        for phase, image in enumerate(phase_images):
            out.update(self.extract(image, mask, case_id, region, phase).prefixed())
        return out

    # ----------------------------------------------------------------- private

    def _base_extractor(self, extra: dict[str, Any] | None = None) -> Any:
        settings = self.config.base_settings()
        if extra:
            settings.update(extra)
        extractor = featureextractor.RadiomicsFeatureExtractor(**settings)
        extractor.disableAllFeatures()
        for cls in self.config.texture_classes:
            extractor.enableFeatureClassByName(cls)
        if self.config.include_first_order:
            extractor.enableFeatureClassByName("firstorder")
        if self.config.include_shape:
            extractor.enableFeatureClassByName("shape")
        return extractor

    def _extract_original(self, image: "sitk.Image", mask: "sitk.Image") -> dict[str, float]:
        extractor = self._base_extractor()
        extractor.disableAllImageTypes()
        extractor.enableImageTypeByName("Original")
        return _to_floats(extractor.execute(image, mask), prefix="original_")

    def _extract_wavelet(
        self, image: "sitk.Image", mask: "sitk.Image", wavelet_filter: WaveletFilter
    ) -> dict[str, float]:
        extractor = self._base_extractor(
            {"wavelet": wavelet_filter.wavelet, "level": wavelet_filter.level, "start_level": 0}
        )
        extractor.disableAllImageTypes()
        extractor.enableImageTypeByName("Wavelet")
        raw = _to_floats(extractor.execute(image, mask), prefix="wavelet")

        keep = set(self.config.keep_subbands)
        out: dict[str, float] = {}
        for name, value in raw.items():
            subband = _subband_of(name)
            if subband in keep:
                out[f"wavelet-{wavelet_filter.name}-{subband}_{_feature_tail(name)}"] = value
        return out

    @staticmethod
    def _voxel_volume(mask: "sitk.Image") -> float:
        array = sitk.GetArrayFromImage(mask)
        spacing = mask.GetSpacing()
        return float(np.count_nonzero(array) * spacing[0] * spacing[1] * spacing[2])


def _to_floats(result: dict, prefix: str) -> dict[str, float]:
    out: dict[str, float] = {}
    for key, value in result.items():
        if not key.startswith(prefix):
            continue
        try:
            out[key] = float(value)
        except (TypeError, ValueError):
            continue
    return out


def _subband_of(feature_name: str) -> str:
    """``wavelet-LLH_glcm_Contrast`` -> ``LLH``."""
    head = feature_name.split("_", 1)[0]
    parts = head.split("-")
    return parts[1] if len(parts) > 1 else ""


def _feature_tail(feature_name: str) -> str:
    """``wavelet-LLH_glcm_Contrast`` -> ``glcm_Contrast``."""
    return feature_name.split("_", 1)[1] if "_" in feature_name else feature_name


def _require_radiomics() -> None:
    if featureextractor is None or sitk is None:  # pragma: no cover
        raise ImportError(
            "Feature extraction needs SimpleITK and pyradiomics:\n"
            "  pip install SimpleITK versioneer\n"
            "  pip install --no-build-isolation pyradiomics"
        )
