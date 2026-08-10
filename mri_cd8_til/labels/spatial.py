"""Turn a patch-level TIL map into the three-class spatial immunophenotype.

Input is the tabular form the Saltz TIL maps ship in: one row per image patch,
with a patch centroid in slide coordinates and a lymphocyte-infiltration
probability.  Optionally a per-patch tumour probability comes along too, which
is what lets us separate *tumour core* from *invasive margin* -- the distinction
the reference study makes by eye on CD8 IHC slides.

The geometry mirrors the pathology protocol:

* **core**    -- tumour patches at least ``margin_width_um`` inside the tumour
  boundary;
* **margin**  -- a band of ``margin_width_um`` either side of the boundary,
  standing in for the "invasion front";
* **stroma**  -- non-tumour tissue beyond the margin band.

The phenotype rule then follows the published definitions directly: prominent
infiltrate reaching the core is *inflamed*; infiltrate that stops at the margin
is *immune-excluded*; infiltrate absent from the margin as well is
*immune-desert*.

One deliberate design choice: the density cut-points are **fitted**, not
hard-coded.  A fitted cut-point is a parameter learned from data, so
:class:`PhenotypeThresholds` must be fit on training cases only.  Fitting it on
the pooled cohort and then evaluating on a held-out split would leak label
information across the split -- the exact failure mode this project measures in
``mri_cd8_til.experiments.optimism``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .schema import HE_TIL_MAP, Immunophenotype, LabelProvenance

#: Saltz TIL maps are computed on 50x50 micron patches at 20x.
DEFAULT_PATCH_UM = 50.0
#: Half-width of the invasive-margin band, in microns. 500 um is the value most
#: commonly used for "invasive margin" in the digital-pathology TIL literature.
DEFAULT_MARGIN_UM = 500.0


@dataclass(frozen=True)
class SpatialTILMetrics:
    """Per-case summary of a TIL map, before any thresholding."""

    case_id: str
    core_density: float
    margin_density: float
    stroma_density: float
    n_core_patches: int
    n_margin_patches: int
    n_tumour_patches: int

    @property
    def core_to_margin_ratio(self) -> float:
        """How far the infiltrate penetrates. Low ratio + high margin = excluded."""
        if self.margin_density <= 0:
            return float("inf") if self.core_density > 0 else 0.0
        return self.core_density / self.margin_density

    @property
    def is_evaluable(self) -> bool:
        """Guard against slides too small or too poorly segmented to phenotype."""
        return self.n_tumour_patches >= 100 and self.n_margin_patches >= 20


@dataclass(frozen=True)
class PhenotypeThresholds:
    """Fitted cut-points separating the three phenotypes.

    ``core_cut`` decides inflamed vs not; ``margin_cut`` then separates desert
    from excluded among the non-inflamed.
    """

    core_cut: float
    margin_cut: float

    @classmethod
    def fit(
        cls,
        metrics: list[SpatialTILMetrics],
        core_quantile: float = 0.60,
        margin_quantile: float = 0.50,
    ) -> "PhenotypeThresholds":
        """Fit cut-points on **training cases only**.

        Quantile defaults reproduce roughly the class balance reported in the
        reference study (a minority inflamed group, with the remainder split
        between desert and excluded); they are parameters, not facts, and any
        study using them should report a sensitivity analysis over them.
        """
        evaluable = [m for m in metrics if m.is_evaluable]
        if not evaluable:
            raise ValueError("no evaluable cases to fit thresholds on")
        core = np.asarray([m.core_density for m in evaluable], dtype=float)
        margin = np.asarray([m.margin_density for m in evaluable], dtype=float)
        return cls(
            core_cut=float(np.quantile(core, core_quantile)),
            margin_cut=float(np.quantile(margin, margin_quantile)),
        )

    def classify(self, metrics: SpatialTILMetrics) -> Immunophenotype:
        if metrics.core_density >= self.core_cut:
            return Immunophenotype.INFLAMED
        if metrics.margin_density >= self.margin_cut:
            return Immunophenotype.EXCLUDED
        return Immunophenotype.DESERT


def compute_spatial_metrics(
    case_id: str,
    x_um: np.ndarray,
    y_um: np.ndarray,
    til_prob: np.ndarray,
    tumour_prob: np.ndarray | None = None,
    til_threshold: float = 0.5,
    tumour_threshold: float = 0.5,
    margin_width_um: float = DEFAULT_MARGIN_UM,
    patch_um: float = DEFAULT_PATCH_UM,
) -> SpatialTILMetrics:
    """Compute core/margin/stroma TIL densities for one slide.

    Parameters
    ----------
    x_um, y_um:
        Patch centroid coordinates in microns.
    til_prob:
        Per-patch lymphocyte-infiltration probability.
    tumour_prob:
        Per-patch tumour probability.  When absent, the tumour region is taken
        to be the whole tissue and only a global density is meaningful -- the
        returned metrics will report ``n_margin_patches == 0`` and therefore
        be flagged non-evaluable, which is the honest outcome.
    """
    x_um = np.asarray(x_um, dtype=float)
    y_um = np.asarray(y_um, dtype=float)
    til = np.asarray(til_prob, dtype=float) >= til_threshold

    if tumour_prob is None:
        return SpatialTILMetrics(
            case_id=case_id,
            core_density=float(til.mean()) if til.size else 0.0,
            margin_density=0.0,
            stroma_density=0.0,
            n_core_patches=int(til.size),
            n_margin_patches=0,
            n_tumour_patches=int(til.size),
        )

    tumour = np.asarray(tumour_prob, dtype=float) >= tumour_threshold
    distance = _signed_distance_to_boundary(x_um, y_um, tumour, patch_um)

    core_mask = tumour & (distance >= margin_width_um)
    margin_mask = np.abs(distance) < margin_width_um
    stroma_mask = (~tumour) & (distance <= -margin_width_um)

    return SpatialTILMetrics(
        case_id=case_id,
        core_density=_density(til, core_mask),
        margin_density=_density(til, margin_mask),
        stroma_density=_density(til, stroma_mask),
        n_core_patches=int(core_mask.sum()),
        n_margin_patches=int(margin_mask.sum()),
        n_tumour_patches=int(tumour.sum()),
    )


def _density(til: np.ndarray, region: np.ndarray) -> float:
    n = int(region.sum())
    return float(til[region].sum() / n) if n else 0.0


def _signed_distance_to_boundary(
    x_um: np.ndarray, y_um: np.ndarray, tumour: np.ndarray, patch_um: float
) -> np.ndarray:
    """Distance from each patch to the nearest tumour-boundary patch.

    Positive inside the tumour, negative outside.  Patches are rasterised onto
    their native grid and the exact Euclidean distance transform is taken, so
    the result is in microns and independent of slide size.

    Falls back to a brute-force nearest-boundary search when SciPy is absent.
    """
    ix = np.rint((x_um - x_um.min()) / patch_um).astype(int)
    iy = np.rint((y_um - y_um.min()) / patch_um).astype(int)
    grid = np.zeros((iy.max() + 1, ix.max() + 1), dtype=bool)
    grid[iy, ix] = tumour

    try:
        from scipy.ndimage import distance_transform_edt

        inside = distance_transform_edt(grid) * patch_um
        outside = distance_transform_edt(~grid) * patch_um
        signed_grid = np.where(grid, inside, -outside)
    except ImportError:  # pragma: no cover - SciPy is a hard dependency in practice
        signed_grid = _brute_force_signed_distance(grid) * patch_um

    return signed_grid[iy, ix]


def _brute_force_signed_distance(grid: np.ndarray) -> np.ndarray:
    """O(n_boundary * n_cells) fallback; only used without SciPy."""
    from itertools import product

    boundary = []
    rows, cols = grid.shape
    for r, c in product(range(rows), range(cols)):
        if not grid[r, c]:
            continue
        for dr, dc in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            nr, nc = r + dr, c + dc
            if not (0 <= nr < rows and 0 <= nc < cols) or not grid[nr, nc]:
                boundary.append((r, c))
                break
    if not boundary:
        return np.where(grid, np.inf, -np.inf)
    bpts = np.asarray(boundary, dtype=float)
    rr, cc = np.meshgrid(np.arange(rows), np.arange(cols), indexing="ij")
    pts = np.stack([rr.ravel(), cc.ravel()], axis=1).astype(float)
    d = np.sqrt(((pts[:, None, :] - bpts[None, :, :]) ** 2).sum(-1)).min(axis=1)
    d = d.reshape(rows, cols)
    return np.where(grid, d, -d)


def provenance() -> LabelProvenance:
    """Provenance record to attach to labels produced by this module."""
    return HE_TIL_MAP
