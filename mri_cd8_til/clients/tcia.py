"""Client for the public TCIA / NBIA REST API.

Only the *open-access* v1 endpoints are used, so no API key is required.  These
are the endpoints that let us answer, without downloading a single voxel, how
many patients a collection has, which scanners produced the images, and which
series are DCE-MRI.

API reference: https://wiki.cancerimagingarchive.net/x/ZoATB
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

from .http import PublicDataSession

BASE_URL = "https://services.cancerimagingarchive.net/nbia-api/services/v1"


@dataclass(frozen=True)
class SeriesRecord:
    """One DICOM series as described by the NBIA ``getSeries`` endpoint."""

    patient_id: str
    study_uid: str
    series_uid: str
    modality: str
    manufacturer: str
    model: str
    description: str
    image_count: int

    @property
    def scanner_key(self) -> str:
        """Vendor+model, used as the batch key for ComBat harmonisation."""
        return f"{self.manufacturer}|{self.model}"


@dataclass
class CollectionSummary:
    """Aggregate view of a TCIA collection, computed from series metadata."""

    collection: str
    n_patients: int
    modalities: list[str]
    n_series: int
    n_mr_series: int
    n_patients_with_mr: int
    scanner_counts: dict[str, int] = field(default_factory=dict)

    @property
    def n_distinct_scanners(self) -> int:
        return len(self.scanner_counts)


class TCIAClient:
    """Thin, read-only wrapper over the NBIA public API."""

    def __init__(self, session: PublicDataSession | None = None, base_url: str = BASE_URL) -> None:
        self.session = session or PublicDataSession()
        self.base_url = base_url.rstrip("/")

    # ---------------------------------------------------------------- queries

    def collections(self) -> list[str]:
        data = self.session.get_json(f"{self.base_url}/getCollectionValues")
        return sorted(item["Collection"] for item in data)

    def patients(self, collection: str) -> list[str]:
        data = self.session.get_json(f"{self.base_url}/getPatient", params={"Collection": collection})
        return sorted(item["PatientId"] for item in data)

    def modalities(self, collection: str) -> list[str]:
        data = self.session.get_json(
            f"{self.base_url}/getModalityValues", params={"Collection": collection}
        )
        return sorted(item["Modality"] for item in data)

    def series(self, collection: str, modality: str | None = None) -> list[SeriesRecord]:
        params: dict[str, str] = {"Collection": collection}
        if modality:
            params["Modality"] = modality
        data = self.session.get_json(f"{self.base_url}/getSeries", params=params)
        return [_parse_series(row) for row in data]

    # ------------------------------------------------------------- aggregates

    def summarize(self, collection: str, top_scanners: int | None = None) -> CollectionSummary:
        """Patient/series counts and scanner mix for one collection.

        ``top_scanners`` truncates the scanner histogram; ``None`` keeps all of
        it (which is what the harmonisation planning needs).
        """
        patients = self.patients(collection)
        modalities = self.modalities(collection)
        all_series = self.series(collection)
        mr_series = [s for s in all_series if s.modality == "MR"]
        counts = Counter(s.scanner_key for s in mr_series)
        scanner_counts = dict(counts.most_common(top_scanners) if top_scanners else counts)
        return CollectionSummary(
            collection=collection,
            n_patients=len(patients),
            modalities=modalities,
            n_series=len(all_series),
            n_mr_series=len(mr_series),
            n_patients_with_mr=len({s.patient_id for s in mr_series}),
            scanner_counts=scanner_counts,
        )

    def patients_with_modality(self, collection: str, modality: str = "MR") -> set[str]:
        return {s.patient_id for s in self.series(collection, modality=modality)}

    def download_series(self, series_uid: str, destination: Path | str) -> Path:
        """Download one series as DICOM files into ``destination``.

        The open-access ``getImage`` endpoint returns a zip of the whole series.
        Series sizes here run from a few MB (a single subtraction volume) to
        hundreds of MB (a full VIBRANT dynamic acquisition), so check
        ``SeriesRecord.image_count`` before pulling one.
        """
        import io
        import zipfile

        destination = Path(destination)
        destination.mkdir(parents=True, exist_ok=True)
        response = self.session.get(
            f"{self.base_url}/getImage", params={"SeriesInstanceUID": series_uid}
        )
        with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
            archive.extractall(destination)
        return destination

    def dce_series(
        self,
        collection: str,
        keywords: Iterable[str] = ("dce", "dyn", "vibrant", "ktrans", "post", "+c", "subtract"),
    ) -> list[SeriesRecord]:
        """Heuristic filter for dynamic contrast-enhanced series.

        TCIA series descriptions are free text and vary per site, so this is a
        recall-oriented pre-filter for a manifest that a human then curates --
        never treat it as a protocol-level guarantee.
        """
        wanted = tuple(k.lower() for k in keywords)
        return [
            s
            for s in self.series(collection, modality="MR")
            if any(k in s.description.lower() for k in wanted)
        ]


def _parse_series(row: dict[str, Any]) -> SeriesRecord:
    return SeriesRecord(
        patient_id=row.get("PatientID", ""),
        study_uid=row.get("StudyInstanceUID", ""),
        series_uid=row.get("SeriesInstanceUID", ""),
        modality=row.get("Modality", ""),
        manufacturer=row.get("Manufacturer") or "NA",
        model=row.get("ManufacturerModelName") or "NA",
        description=row.get("SeriesDescription") or "",
        image_count=int(row.get("ImageCount") or 0),
    )
