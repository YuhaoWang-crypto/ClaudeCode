"""Client for the NCI Imaging Data Commons (IDC), used to fetch TIL maps.

The Saltz TIL maps are published on TCIA as a 73 GB bundle, which is an awkward
way to obtain the ~136 slides that matter here.  IDC re-publishes the same
analysis results as individual DICOM ``SEG`` objects in a public S3 bucket, with
a locally-queryable index (``idc-index``), so the maps for one patient list can
be pulled in seconds without touching the bundle.

Three map variants exist per slide; they are *not* interchangeable:

``Stony Brook Inception-V4 Binary TIL Map``
    Thresholded 0/1 per 50 micron patch. What the 2018 paper's figures show.
``Stony Brook Inception-V4 Fractional TIL Map``
    The underlying probability, quantised. Preferred here -- thresholding is a
    modelling decision that belongs downstream, not in the label loader.
``Stony Brook CNN-generated TIL Map``
    The earlier CNN. Present for only some slides.

**What these maps do not contain: a tumour region.** They mark lymphocyte
infiltration across the whole slide with no tumour/stroma boundary, so tumour
core and invasive margin cannot be separated from them alone.  That is the
single reason the tier-A cohort supports an overall-TIL target today and not the
desert-vs-excluded contrast -- see ``docs/TEST_SET.md``.
"""

from __future__ import annotations

import io
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .http import PublicDataSession

S3_BUCKET_URL = "https://idc-open-data.s3.amazonaws.com"

BINARY_TIL = "Stony Brook Inception-V4 Binary TIL Map"
FRACTIONAL_TIL = "Stony Brook Inception-V4 Fractional TIL Map"
CNN_TIL = "Stony Brook CNN-generated TIL Map"

_KEY_RE = re.compile(r"<Key>([^<]+)</Key>")


@dataclass(frozen=True)
class TILMap:
    """One slide's TIL map, as a patch grid."""

    patient_id: str
    series_uid: str
    description: str
    #: ``(rows, cols)`` array; binary maps are 0/1, fractional are quantised.
    grid: np.ndarray
    patch_um: float = 50.0

    @property
    def n_patches(self) -> int:
        return int(self.grid.size)

    @property
    def til_fraction(self) -> float:
        """Fraction of tissue patches classified as lymphocyte-infiltrated.

        Note this is over the *whole slide grid*, background included, so it is
        comparable across slides only if their tissue fraction is similar. Use
        :meth:`til_fraction_over_tissue` when that matters.
        """
        return float((self.grid > 0).mean())

    def til_fraction_over_tissue(self, tissue: np.ndarray | None = None) -> float:
        if tissue is None:
            return self.til_fraction
        n = int(tissue.sum())
        return float((self.grid > 0)[tissue].sum() / n) if n else 0.0

    def coordinates_um(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Patch centroids in microns plus the TIL value, ready for the
        spatial-metrics code in :mod:`mri_cd8_til.labels.spatial`."""
        rows, cols = np.indices(self.grid.shape)
        return (
            cols.ravel().astype(float) * self.patch_um,
            rows.ravel().astype(float) * self.patch_um,
            self.grid.ravel().astype(float),
        )


class IDCClient:
    """Fetches TIL-map DICOM objects from IDC's public S3 bucket."""

    def __init__(self, session: PublicDataSession | None = None) -> None:
        self.session = session or PublicDataSession()

    # ------------------------------------------------------------------ index

    @staticmethod
    def til_series(patient_ids: set[str] | None = None, description: str = FRACTIONAL_TIL):
        """Look up TIL-map series in the local ``idc-index`` snapshot.

        Returns a pandas DataFrame; raises ImportError with an actionable
        message if ``idc-index`` is absent, since there is no lightweight
        fallback (the alternative is enumerating 3,138 S3 objects).
        """
        try:
            from idc_index import IDCClient as _Index
        except ImportError as exc:  # pragma: no cover
            raise ImportError(
                "TIL-map lookup needs the idc-index package: pip install idc-index"
            ) from exc

        index = _Index().index
        brca = index[
            index.collection_id.astype(str).str.lower().str.replace("-", "_", regex=False)
            == "tcga_brca"
        ]
        series = brca[brca.SeriesDescription.astype(str) == description]
        if patient_ids is not None:
            series = series[series.PatientID.isin(patient_ids)]
        return series

    # --------------------------------------------------------------- download

    def fetch_map(self, crdc_series_uuid: str, patient_id: str, description: str) -> TILMap:
        """Download and decode one TIL-map series.

        Takes IDC's ``crdc_series_uuid``, **not** the SeriesInstanceUID: objects
        live under a per-series UUID prefix in the bucket, and listing by
        SeriesInstanceUID returns nothing at all.
        """
        import pydicom

        keys = self._object_keys(crdc_series_uuid)
        if not keys:
            raise FileNotFoundError(f"no S3 objects under prefix {crdc_series_uuid}")
        payload = self.session.get(f"{S3_BUCKET_URL}/{keys[0]}").content
        dataset = pydicom.dcmread(io.BytesIO(payload))
        grid = np.asarray(dataset.pixel_array)
        if grid.ndim == 3:  # multi-frame: collapse to the first frame
            grid = grid[0]
        return TILMap(
            patient_id=patient_id,
            series_uid=crdc_series_uuid,
            description=str(dataset.get("SeriesDescription", description)),
            grid=grid,
        )

    def _object_keys(self, prefix: str) -> list[str]:
        """S3 keys under an IDC prefix, via a plain list-objects call."""
        listing = self.session.get_text(
            f"{S3_BUCKET_URL}/?list-type=2&prefix={prefix}&max-keys=64"
        )
        return _KEY_RE.findall(listing)

    def fetch_by_crdc_prefix(self, prefix: str) -> bytes:
        keys = self._object_keys(prefix)
        if not keys:
            raise FileNotFoundError(f"no S3 objects under prefix {prefix}")
        return self.session.get(f"{S3_BUCKET_URL}/{keys[0]}").content


def cache_path(root: Path | str, patient_id: str) -> Path:
    return Path(root) / f"{patient_id}.npz"


def save_map(til_map: TILMap, root: Path | str) -> Path:
    path = cache_path(root, til_map.patient_id)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        grid=til_map.grid,
        patient_id=til_map.patient_id,
        series_uid=til_map.series_uid,
        description=til_map.description,
    )
    return path


def load_map(root: Path | str, patient_id: str) -> TILMap:
    data = np.load(cache_path(root, patient_id), allow_pickle=False)
    return TILMap(
        patient_id=patient_id,
        series_uid=str(data["series_uid"]),
        description=str(data["description"]),
        grid=data["grid"],
    )
