"""Client for the NCI Genomic Data Commons (GDC) public API.

Used to establish which TCGA-BRCA cases carry the *pathology* and
*transcriptome* halves of the label, so the imaging cohort can be intersected
against them without downloading anything.

API reference: https://docs.gdc.cancer.gov/API/Users_Guide/Getting_Started/
"""

from __future__ import annotations

import json
from dataclasses import dataclass

from .http import PublicDataSession

BASE_URL = "https://api.gdc.cancer.gov"
PAGE_SIZE = 1000


@dataclass(frozen=True)
class CaseAvailability:
    """Which data modalities exist for a set of TCGA cases."""

    project_id: str
    all_cases: frozenset[str]
    rna_seq: frozenset[str]
    diagnostic_slides: frozenset[str]
    tissue_slides: frozenset[str]

    @property
    def any_slide(self) -> frozenset[str]:
        return self.diagnostic_slides | self.tissue_slides

    def summary(self) -> dict[str, int]:
        return {
            "cases": len(self.all_cases),
            "rna_seq": len(self.rna_seq),
            "diagnostic_slides": len(self.diagnostic_slides),
            "tissue_slides": len(self.tissue_slides),
            "any_slide": len(self.any_slide),
        }


class GDCClient:
    def __init__(self, session: PublicDataSession | None = None, base_url: str = BASE_URL) -> None:
        self.session = session or PublicDataSession()
        self.base_url = base_url.rstrip("/")

    def _cases(self, project_id: str, extra_filters: list[dict]) -> set[str]:
        """Return submitter_ids (e.g. ``TCGA-AO-A03M``) matching the filters."""
        filters = {
            "op": "and",
            "content": [
                {
                    "op": "in",
                    "content": {"field": "cases.project.project_id", "value": [project_id]},
                }
            ]
            + extra_filters,
        }
        found: set[str] = set()
        offset = 0
        while True:
            payload = self.session.get_json(
                f"{self.base_url}/cases",
                params={
                    "filters": json.dumps(filters),
                    "fields": "submitter_id",
                    "format": "JSON",
                    "size": str(PAGE_SIZE),
                    "from": str(offset),
                },
            )
            hits = payload["data"]["hits"]
            found.update(hit["submitter_id"] for hit in hits)
            offset += len(hits)
            if not hits or offset >= payload["data"]["pagination"]["total"]:
                break
        return found

    def availability(self, project_id: str = "TCGA-BRCA") -> CaseAvailability:
        """Cases with RNA-seq / diagnostic (FFPE) WSI / tissue (frozen) WSI."""
        return CaseAvailability(
            project_id=project_id,
            all_cases=frozenset(self._cases(project_id, [])),
            rna_seq=frozenset(
                self._cases(
                    project_id,
                    [
                        {
                            "op": "in",
                            "content": {
                                "field": "files.data_type",
                                "value": ["Gene Expression Quantification"],
                            },
                        }
                    ],
                )
            ),
            diagnostic_slides=frozenset(
                self._cases(
                    project_id,
                    [
                        {
                            "op": "in",
                            "content": {
                                "field": "files.experimental_strategy",
                                "value": ["Diagnostic Slide"],
                            },
                        }
                    ],
                )
            ),
            tissue_slides=frozenset(
                self._cases(
                    project_id,
                    [
                        {
                            "op": "in",
                            "content": {
                                "field": "files.experimental_strategy",
                                "value": ["Tissue Slide"],
                            },
                        }
                    ],
                )
            ),
        )
