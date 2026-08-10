"""Minimal NCBI GEO client for linking imaging cohorts to expression cohorts.

The only thing we need from GEO at inventory time is the *sample-to-patient*
mapping, which the lightweight ``targ=gsm&form=text&view=brief`` view exposes
without downloading a multi-hundred-megabyte series matrix.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from .http import PublicDataSession

ACC_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"

_TITLE_RE = re.compile(r"^!Sample_title = (.+)$", re.M)
_ACCESSION_RE = re.compile(r"^!Sample_geo_accession = (\S+)$", re.M)
_CHARACTERISTIC_RE = re.compile(r"^!Sample_characteristics_ch1 = ([^:]+): (.*)$", re.M)


@dataclass
class GEOSeries:
    """Sample-level metadata for one GEO series."""

    accession: str
    titles: list[str] = field(default_factory=list)
    sample_accessions: list[str] = field(default_factory=list)
    characteristics: dict[str, list[str]] = field(default_factory=dict)

    @property
    def n_samples(self) -> int:
        return len(self.titles)

    def patient_ids(self, characteristic: str = "patient id") -> set[str]:
        """Patient identifiers taken from a ``!Sample_characteristics_ch1`` key.

        Falls back to digit runs inside the sample titles when the series does
        not carry an explicit patient-id characteristic.
        """
        values = self.characteristics.get(characteristic)
        if values:
            return {v.strip() for v in values if v.strip()}
        return {m for title in self.titles for m in re.findall(r"\d{3,}", title)}


class GEOClient:
    def __init__(self, session: PublicDataSession | None = None) -> None:
        self.session = session or PublicDataSession()

    def series(self, accession: str) -> GEOSeries:
        text = self.session.get_text(
            ACC_URL, params={"acc": accession, "targ": "gsm", "form": "text", "view": "brief"}
        )
        characteristics: dict[str, list[str]] = {}
        for key, value in _CHARACTERISTIC_RE.findall(text):
            characteristics.setdefault(key.strip().lower(), []).append(value.strip())
        return GEOSeries(
            accession=accession,
            titles=_TITLE_RE.findall(text),
            sample_accessions=_ACCESSION_RE.findall(text),
            characteristics=characteristics,
        )
