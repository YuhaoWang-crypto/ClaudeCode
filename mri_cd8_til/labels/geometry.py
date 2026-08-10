"""TIL spatial-architecture descriptors that need no tumour mask.

The desert / excluded / inflamed trichotomy is defined against a tumour
boundary, and the published TIL maps do not carry one.  Rather than give up on
the spatial half of the problem, this module measures the *geometry* the
trichotomy is a description of:

============  ==========================================================
phenotype     what its geometry looks like
============  ==========================================================
desert        little infiltrate anywhere
excluded      infiltrate present but concentrated in bands or rims, not
              reaching the interior of the tissue it borders
inflamed      infiltrate distributed through the tissue
============  ==========================================================

Read that table again and notice that only the *middle* column mentions the
tumour.  "Concentrated in bands and not reaching the interior" is a statement
about the arrangement of one point pattern inside a tissue region, and it is
measurable directly.

**The central descriptor is penetration.**  For every tissue patch, take the
distance to the nearest TIL-positive patch.  Inflamed tissue has a short-tailed
distribution -- nowhere is far from lymphocytes.  Excluded tissue has a long
tail: the banded infiltrate leaves large interiors untouched.  Desert tissue has
a long tail too, but with almost no TIL patches to measure from, which is why
penetration has to be read alongside amount rather than instead of it.

Everything here is **scale-free**: distances are divided by the tissue's
equivalent radius and counts by tissue area, so a 400x450 slide and a 900x1100
slide are comparable.  Absolute pixel measures would otherwise encode specimen
size, which is a property of the grossing, not of the immune response.

A caveat that belongs in any paper using these: they *approximate* the published
phenotype, they do not reproduce it.  A band of lymphocytes at the tumour margin
and a band at the edge of a fat lobule look the same to a descriptor with no
tumour boundary.  What the transcriptome check in ``scripts/phenotype_til.py``
establishes is whether that approximation retains biological meaning.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

#: Radii, in patches (50 microns each), at which multi-scale statistics are
#: evaluated. 2 patches = 100 um is roughly a small lymphoid aggregate; 32
#: patches = 1.6 mm is the scale of a tumour nest plus its margin.
SCALE_RADII: tuple[int, ...] = (2, 4, 8, 16, 32)


@dataclass(frozen=True)
class TILGeometry:
    """Scale-free descriptors of one slide's TIL point pattern."""

    case_id: str

    # --- amount ----------------------------------------------------------
    til_fraction: float
    til_mean_prob: float
    n_tissue_patches: int

    # --- clustering ------------------------------------------------------
    morans_i: float
    largest_cluster_share: float
    n_clusters_per_area: float
    gini_cluster_size: float

    # --- penetration: the excluded-vs-inflamed axis ----------------------
    penetration_median: float
    penetration_p90: float
    penetration_gini: float
    unreached_fraction: float

    # --- shape of the infiltrate ----------------------------------------
    interface_density: float
    band_index: float
    cluster_elongation: float

    # --- multi-scale -----------------------------------------------------
    ripley_l: tuple[float, ...]
    lacunarity: tuple[float, ...]

    def as_dict(self) -> dict:
        record = asdict(self)
        for key in ("ripley_l", "lacunarity"):
            for radius, value in zip(SCALE_RADII, record.pop(key)):
                record[f"{key}_r{radius}"] = float(value)
        return record

    @property
    def is_evaluable(self) -> bool:
        return self.n_tissue_patches >= 1000


def describe(
    case_id: str,
    grid: np.ndarray,
    prob_cut: int = 128,
    radii: tuple[int, ...] = SCALE_RADII,
) -> TILGeometry:
    """Compute every descriptor for one TIL map.

    ``grid`` is the fractional TIL map: uint8 0-255 where 0 doubles as
    "no tissue", so tissue and TIL calls come from two thresholds on one array.
    """
    grid = np.asarray(grid)
    tissue = grid > 0
    til = tissue & (grid >= prob_cut)
    n_tissue = int(tissue.sum())

    if n_tissue < 1:
        raise ValueError(f"{case_id}: empty tissue mask")

    # Equivalent radius of a disc with the tissue's area, the natural length
    # unit for normalising distances on an irregular specimen outline.
    equivalent_radius = np.sqrt(n_tissue / np.pi)

    labels, sizes = _components(til)
    penetration = _penetration(tissue, til)

    return TILGeometry(
        case_id=case_id,
        til_fraction=float(til.sum() / n_tissue),
        til_mean_prob=float(grid[tissue].mean() / 255.0),
        n_tissue_patches=n_tissue,
        morans_i=_morans_i(grid, tissue),
        largest_cluster_share=float(sizes.max() / sizes.sum()) if sizes.size else 0.0,
        n_clusters_per_area=float(sizes.size / n_tissue * 1000.0),
        gini_cluster_size=_gini(sizes),
        penetration_median=float(np.median(penetration) / equivalent_radius)
        if penetration.size
        else 1.0,
        penetration_p90=float(np.quantile(penetration, 0.9) / equivalent_radius)
        if penetration.size
        else 1.0,
        penetration_gini=_gini(penetration),
        unreached_fraction=_unreached_fraction(penetration, equivalent_radius),
        interface_density=_interface_density(til, tissue),
        band_index=_band_index(labels, sizes),
        cluster_elongation=_cluster_elongation(labels, sizes),
        ripley_l=_ripley_l(til, tissue, radii),
        lacunarity=_lacunarity(til, tissue, radii),
    )


# ------------------------------------------------------------------ internals


def _components(mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    from scipy.ndimage import label

    labels, n = label(mask, structure=np.ones((3, 3)))
    if n == 0:
        return labels, np.array([], dtype=int)
    return labels, np.bincount(labels.ravel())[1:]


def _morans_i(field: np.ndarray, mask: np.ndarray) -> float:
    """Spatial autocorrelation of the TIL probability over tissue patches only.

    Computed over ``mask`` rather than the full grid: over the whole grid the
    background zeros dominate and the statistic reports that tissue is one
    contiguous blob, returning ~0.6 for every slide regardless of arrangement.
    """
    x = field.astype(float)
    mean = x[mask].mean()
    deviation = np.where(mask, x - mean, 0.0)
    denominator = float((deviation[mask] ** 2).sum())
    if denominator <= 0 or mask.sum() < 4:
        return 0.0

    numerator = 0.0
    n_pairs = 0
    for axis in (0, 1):
        head, tail = slice(1, None), slice(None, -1)
        if axis == 0:
            a, b, valid = deviation[head, :], deviation[tail, :], mask[head, :] & mask[tail, :]
        else:
            a, b, valid = deviation[:, head], deviation[:, tail], mask[:, head] & mask[:, tail]
        numerator += float((a[valid] * b[valid]).sum())
        n_pairs += int(valid.sum())

    if n_pairs == 0:
        return 0.0
    return float((int(mask.sum()) / n_pairs) * (numerator / denominator))


def _penetration(tissue: np.ndarray, til: np.ndarray) -> np.ndarray:
    """Distance from each tissue patch to the nearest TIL-positive patch.

    This is the excluded-vs-inflamed axis. Returned in patches; the caller
    normalises by the tissue's equivalent radius.
    """
    from scipy.ndimage import distance_transform_edt

    if not til.any():
        # No infiltrate at all: everything is maximally far. Use the tissue's
        # own extent rather than infinity so the statistic stays finite.
        return np.full(int(tissue.sum()), float(np.sqrt(tissue.sum() / np.pi)))
    distance = distance_transform_edt(~til)
    return distance[tissue]


def _unreached_fraction(penetration: np.ndarray, equivalent_radius: float) -> float:
    """Share of tissue more than 10% of the equivalent radius from any TIL.

    A single scalar for "how much of this specimen the infiltrate never gets
    near", which is the operational content of "excluded".
    """
    if penetration.size == 0:
        return 1.0
    return float((penetration > 0.10 * equivalent_radius).mean())


def _interface_density(til: np.ndarray, tissue: np.ndarray) -> float:
    """Boundary length of the TIL regions per unit TIL area.

    High for fragmented or banded infiltrate, low for compact blobs. Scale-free
    because it is a ratio of a perimeter to the area it encloses.
    """
    n_til = int(til.sum())
    if n_til == 0:
        return 0.0
    boundary = 0
    for axis in (0, 1):
        head, tail = slice(1, None), slice(None, -1)
        if axis == 0:
            a, b, valid = til[head, :], til[tail, :], tissue[head, :] & tissue[tail, :]
        else:
            a, b, valid = til[:, head], til[:, tail], tissue[:, head] & tissue[:, tail]
        boundary += int((a[valid] != b[valid]).sum())
    return float(boundary / n_til)


def _band_index(labels: np.ndarray, sizes: np.ndarray) -> float:
    """How band-like the infiltrate is, area-weighted over clusters.

    Uses the isoperimetric ratio ``P^2 / (4*pi*A)``: 1 for a disc, large for an
    elongated band. Immune exclusion presents as bands and rims, so this is the
    shape counterpart to the penetration statistic.
    """
    if sizes.size == 0:
        return 0.0
    from scipy.ndimage import find_objects

    total = 0.0
    weight = 0.0
    for index, box in enumerate(find_objects(labels), start=1):
        if box is None:
            continue
        area = float(sizes[index - 1])
        if area < 4:
            continue
        patch = labels[box] == index
        perimeter = _perimeter(patch)
        ratio = perimeter**2 / (4.0 * np.pi * area)
        total += ratio * area
        weight += area
    return float(total / weight) if weight else 0.0


def _perimeter(patch: np.ndarray) -> float:
    padded = np.pad(patch, 1)
    edges = 0
    for shift, axis in ((1, 0), (-1, 0), (1, 1), (-1, 1)):
        edges += int((padded & ~np.roll(padded, shift, axis=axis)).sum())
    return float(edges)


def _cluster_elongation(labels: np.ndarray, sizes: np.ndarray) -> float:
    """Area-weighted mean elongation of TIL clusters, from second moments.

    ``sqrt(lambda_max / lambda_min)`` of each cluster's coordinate covariance:
    1 for isotropic, larger for elongated. Complements the band index, which is
    boundary-based and therefore more sensitive to raggedness than to direction.
    """
    if sizes.size == 0:
        return 1.0
    from scipy.ndimage import find_objects

    total = 0.0
    weight = 0.0
    for index, box in enumerate(find_objects(labels), start=1):
        if box is None or sizes[index - 1] < 8:
            continue
        rows, cols = np.nonzero(labels[box] == index)
        coords = np.stack([rows, cols]).astype(float)
        covariance = np.cov(coords)
        if not np.all(np.isfinite(covariance)):
            continue
        eigenvalues = np.linalg.eigvalsh(covariance)
        eigenvalues = np.clip(eigenvalues, 1e-9, None)
        elongation = float(np.sqrt(eigenvalues.max() / eigenvalues.min()))
        area = float(sizes[index - 1])
        total += elongation * area
        weight += area
    return float(total / weight) if weight else 1.0


def _ripley_l(til: np.ndarray, tissue: np.ndarray, radii: tuple[int, ...]) -> tuple[float, ...]:
    """Besag's L(r) - r, on a lattice, restricted to tissue.

    Positive means clustering at scale ``r``, zero means the infiltrate is
    randomly placed within the tissue, negative means it is regularly spaced.
    Evaluating it at several radii separates tight aggregates from diffuse
    infiltration, which a single-scale statistic cannot do.
    """
    from scipy.ndimage import uniform_filter

    n_til = int(til.sum())
    n_tissue = int(tissue.sum())
    if n_til < 5 or n_tissue == 0:
        return tuple(0.0 for _ in radii)

    intensity = n_til / n_tissue
    til_float = til.astype(float)
    tissue_float = tissue.astype(float)

    out = []
    for radius in radii:
        window = 2 * radius + 1
        # Neighbour counts within a square window, with the tissue window as the
        # denominator so the specimen outline does not masquerade as structure.
        neighbours = uniform_filter(til_float, window, mode="constant") * window**2
        available = uniform_filter(tissue_float, window, mode="constant") * window**2
        valid = til & (available > 1)
        if not valid.any():
            out.append(0.0)
            continue
        observed = float(((neighbours[valid] - 1) / (available[valid] - 1)).mean())
        # K(r) expressed as a ratio to the expected count, then Besag's
        # variance-stabilising square root, on the same window geometry.
        k_ratio = observed / intensity if intensity > 0 else 0.0
        out.append(float(np.sqrt(max(k_ratio, 0.0)) - 1.0))
    return tuple(out)


def _lacunarity(til: np.ndarray, tissue: np.ndarray, radii: tuple[int, ...]) -> tuple[float, ...]:
    """Gliding-box lacunarity of the TIL field within tissue.

    Measures gappiness: how unevenly the infiltrate fills the space at each
    scale. Two slides with identical TIL fractions and identical Moran's I can
    still differ here if one has a few large voids and the other many small ones.
    """
    from scipy.ndimage import uniform_filter

    if not til.any() or not tissue.any():
        return tuple(0.0 for _ in radii)

    til_float = til.astype(float)
    tissue_float = tissue.astype(float)
    out = []
    for radius in radii:
        window = 2 * radius + 1
        counts = uniform_filter(til_float, window, mode="constant") * window**2
        available = uniform_filter(tissue_float, window, mode="constant") * window**2
        valid = tissue & (available > 0.5 * window**2)
        if valid.sum() < 10:
            out.append(0.0)
            continue
        values = counts[valid]
        mean = values.mean()
        out.append(float((values.var() + mean**2) / mean**2) if mean > 0 else 0.0)
    return tuple(out)


def _gini(values: np.ndarray) -> float:
    """Gini coefficient: 0 when all values are equal, →1 when concentrated."""
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return 0.0
    values = np.sort(np.clip(values, 0, None))
    total = values.sum()
    if total <= 0:
        return 0.0
    index = np.arange(1, values.size + 1)
    return float((2.0 * (index * values).sum()) / (values.size * total) - (values.size + 1) / values.size)


#: Descriptor names in the order :meth:`TILGeometry.as_dict` emits them, minus
#: the identifier -- convenient for building a feature matrix.
def descriptor_names(radii: tuple[int, ...] = SCALE_RADII) -> list[str]:
    base = [
        "til_fraction",
        "til_mean_prob",
        "n_tissue_patches",
        "morans_i",
        "largest_cluster_share",
        "n_clusters_per_area",
        "gini_cluster_size",
        "penetration_median",
        "penetration_p90",
        "penetration_gini",
        "unreached_fraction",
        "interface_density",
        "band_index",
        "cluster_elongation",
    ]
    return base + [f"ripley_l_r{r}" for r in radii] + [f"lacunarity_r{r}" for r in radii]
