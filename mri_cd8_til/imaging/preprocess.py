"""DCE-MRI preprocessing for a *multi-centre* radiomics study.

For a single-centre study with one scanner and one protocol, preprocessing is
mostly cosmetic.  For a multi-centre study it is where most of the achievable
generalisation is won or lost, because MR intensities carry no physical units:
the same tissue in the same patient produces different numbers on a GE SIGNA
and a Siemens Skyra.  Left alone, a radiomic model trained on pooled data will
learn the scanner.

The steps here, in order, and what each is defending against:

1. **N4 bias-field correction** -- coil sensitivity, which varies by vendor and
   coil geometry and is a large low-frequency intensity gradient.
2. **Resampling to a common voxel grid** -- texture features are defined on a
   voxel lattice, so a model trained at 0.6 x 0.6 x 2.0 mm and applied at
   0.8 x 0.8 x 1.0 mm is being fed different features under the same names.
3. **Intensity normalisation** -- puts arbitrary MR units on a common scale.
4. **Fixed-bin-width discretisation** (applied at extraction time) -- with a
   fixed *count* of bins, the bin width depends on the intensity range of each
   individual lesion, which re-introduces exactly the scanner dependence steps
   1-3 removed.

None of this makes cross-scanner features identical; it makes the residual
difference small enough for ComBat to model.  Steps 1-4 plus ComBat plus
leave-one-centre-out validation is the minimum defensible protocol.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

try:  # pragma: no cover
    import SimpleITK as sitk
except ImportError:  # pragma: no cover
    sitk = None  # type: ignore[assignment]


@dataclass
class PreprocessConfig:
    """Preprocessing parameters, recorded so a run is reproducible."""

    #: The reference study's acquisition. Resampling every cohort onto this grid
    #: keeps features comparable with the published ones.
    target_spacing_mm: tuple[float, float, float] = (0.6, 0.6, 2.0)
    apply_n4: bool = True
    n4_shrink_factor: int = 4
    n4_iterations: tuple[int, ...] = (50, 50, 50, 50)
    #: "zscore" (whole-image), "roi_zscore", "percentile" or "none".
    normalization: str = "percentile"
    percentile_range: tuple[float, float] = (1.0, 99.0)
    interpolator: str = "linear"
    #: Mask interpolation must stay nearest-neighbour or labels get blended.
    mask_interpolator: str = "nearest"

    def as_dict(self) -> dict:
        return {
            "target_spacing_mm": list(self.target_spacing_mm),
            "apply_n4": self.apply_n4,
            "n4_shrink_factor": self.n4_shrink_factor,
            "n4_iterations": list(self.n4_iterations),
            "normalization": self.normalization,
            "percentile_range": list(self.percentile_range),
            "interpolator": self.interpolator,
        }


@dataclass
class DCESeries:
    """One patient's dynamic acquisition: a pre-contrast volume plus phases."""

    pre_contrast: "sitk.Image"
    post_contrast: list["sitk.Image"] = field(default_factory=list)
    phase_interval_s: float = 74.0

    @property
    def n_phases(self) -> int:
        return len(self.post_contrast)

    def subtraction(self, phase: int = 1) -> "sitk.Image":
        """Post-minus-pre subtraction image.

        The reference study segments on the subtraction of the **second**
        post-contrast phase (``phase=1`` zero-indexed), which is where invasive
        breast cancer typically peaks.
        """
        _require_sitk()
        if not 0 <= phase < self.n_phases:
            raise IndexError(f"phase {phase} out of range (n_phases={self.n_phases})")
        pre = sitk.Cast(self.pre_contrast, sitk.sitkFloat32)
        post = sitk.Cast(self.post_contrast[phase], sitk.sitkFloat32)
        return post - pre

    def wash_in_slope(self) -> "sitk.Image":
        """Voxelwise initial enhancement rate, a protocol-robust kinetic map.

        Ratios of phases cancel much of the scanner-dependent gain that raw
        intensities carry, so kinetic maps travel across centres considerably
        better than the raw phase images do.
        """
        _require_sitk()
        if self.n_phases < 1:
            raise ValueError("need at least one post-contrast phase")
        pre = sitk.Cast(self.pre_contrast, sitk.sitkFloat32)
        first = sitk.Cast(self.post_contrast[0], sitk.sitkFloat32)
        denominator = pre + sitk.Cast(sitk.Abs(pre) < 1e-6, sitk.sitkFloat32)
        return (first - pre) / denominator / self.phase_interval_s


def preprocess_volume(
    image: "sitk.Image",
    config: PreprocessConfig,
    mask: "sitk.Image | None" = None,
) -> "sitk.Image":
    """Run the full preprocessing chain on one volume."""
    _require_sitk()
    out = sitk.Cast(image, sitk.sitkFloat32)
    if config.apply_n4:
        out = correct_bias_field(out, config)
    out = resample_to_spacing(out, config.target_spacing_mm, config.interpolator)
    out = normalize_intensity(out, config, mask)
    return out


def correct_bias_field(image: "sitk.Image", config: PreprocessConfig) -> "sitk.Image":
    """N4ITK bias-field correction, fitted on a shrunk image for speed."""
    _require_sitk()
    image = sitk.Cast(image, sitk.sitkFloat32)
    head_mask = sitk.OtsuThreshold(image, 0, 1, 200)
    shrunk = sitk.Shrink(image, [config.n4_shrink_factor] * image.GetDimension())
    shrunk_mask = sitk.Shrink(head_mask, [config.n4_shrink_factor] * image.GetDimension())

    corrector = sitk.N4BiasFieldCorrectionImageFilter()
    corrector.SetMaximumNumberOfIterations(list(config.n4_iterations))
    corrector.Execute(shrunk, shrunk_mask)

    # Apply the full-resolution log bias field rather than the shrunk output, so
    # no resolution is lost to the speed-up.
    log_bias = corrector.GetLogBiasFieldAsImage(image)
    return image / sitk.Exp(log_bias)


def resample_to_spacing(
    image: "sitk.Image",
    spacing: tuple[float, float, float],
    interpolator: str = "linear",
) -> "sitk.Image":
    """Resample onto a fixed voxel grid, preserving physical extent."""
    _require_sitk()
    original_spacing = image.GetSpacing()
    original_size = image.GetSize()
    new_size = [
        max(1, int(round(osz * osp / nsp)))
        for osz, osp, nsp in zip(original_size, original_spacing, spacing)
    ]
    resampler = sitk.ResampleImageFilter()
    resampler.SetOutputSpacing(spacing)
    resampler.SetSize(new_size)
    resampler.SetOutputDirection(image.GetDirection())
    resampler.SetOutputOrigin(image.GetOrigin())
    resampler.SetTransform(sitk.Transform())
    resampler.SetDefaultPixelValue(0)
    resampler.SetInterpolator(_INTERPOLATORS[interpolator])
    return resampler.Execute(image)


def normalize_intensity(
    image: "sitk.Image",
    config: PreprocessConfig,
    mask: "sitk.Image | None" = None,
) -> "sitk.Image":
    """Put arbitrary MR units onto a common scale."""
    _require_sitk()
    if config.normalization == "none":
        return image

    array = sitk.GetArrayFromImage(image).astype(np.float64)
    if config.normalization == "roi_zscore":
        if mask is None:
            raise ValueError("roi_zscore normalization requires a mask")
        selection = array[sitk.GetArrayFromImage(mask) > 0]
    else:
        selection = array[array > 0]
    if selection.size == 0:
        return image

    if config.normalization in ("zscore", "roi_zscore"):
        centre, scale = selection.mean(), selection.std()
    elif config.normalization == "percentile":
        low, high = np.percentile(selection, config.percentile_range)
        centre, scale = low, (high - low)
    else:
        raise ValueError(f"unknown normalization {config.normalization!r}")

    scale = scale if scale > 1e-9 else 1.0
    out = sitk.GetImageFromArray((array - centre) / scale)
    out.CopyInformation(image)
    return sitk.Cast(out, sitk.sitkFloat32)


def resample_mask_like(mask: "sitk.Image", reference: "sitk.Image") -> "sitk.Image":
    """Put a mask onto a reference image's grid without blending labels."""
    _require_sitk()
    return sitk.Resample(
        sitk.Cast(mask > 0, sitk.sitkUInt8),
        reference,
        sitk.Transform(),
        sitk.sitkNearestNeighbor,
        0,
        sitk.sitkUInt8,
    )


_INTERPOLATORS = {
    "linear": getattr(sitk, "sitkLinear", None),
    "nearest": getattr(sitk, "sitkNearestNeighbor", None),
    "bspline": getattr(sitk, "sitkBSpline", None),
}


def _require_sitk() -> None:
    if sitk is None:  # pragma: no cover
        raise ImportError("SimpleITK is required for preprocessing: pip install SimpleITK")
