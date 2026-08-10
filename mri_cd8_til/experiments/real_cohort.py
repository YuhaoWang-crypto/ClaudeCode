"""Assemble the real 91-patient cohort: MRI radiomics x pathology TIL labels.

This closes the imaging-side gap that ``docs/TEST_SET.md`` recorded as open.
The missing piece was tumour delineation on TCGA-BRCA MRI -- no public
segmentations exist in the imaging collection itself, and MAMA-MIA covers
ISPY1/ISPY2/NACT/DUKE but not TCGA-BRCA.  They do exist, in a place the
collection page does not point at: the TCGA Breast Phenotype Research Group's
analysis result (TCIA, CC BY 3.0), which published radiologist-supervised lesion
segmentations for 91 cases together with the radiomic features computed from
them on the UChicago MRI Quantitative Radiomics workstation.

Using their published features rather than re-extracting has one real advantage
and one real cost, and both should be stated:

* advantage -- the segmentations are radiologist-supervised rather than a
  threshold surrogate, and the features are a fixed, documented set produced by
  one workstation, so nothing here depends on our own extraction choices;
* cost -- it is 36 features, not the reference study's 833 per phase, and it has
  no wavelet bank and no peritumoural region. So this cohort can test *whether
  DCE-MRI carries TIL signal at all*; it cannot reproduce the reference study's
  exact feature space.

The 36 features fall into four families, which is worth keeping track of because
they degrade very differently across scanners.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

TCIA_ATTACHMENTS = "https://wiki.cancerimagingarchive.net/download/attachments/19039112"
FEATURES_FILE = (
    "TCGA%20Run%202014_91cases_features_UChicago%20V2010%20MRI%20Workstation.xls?api=v2"
)
SEGMENTATIONS_FILE = "TCGA_Segmented_Lesions_UofC.zip?api=v2"

_PATIENT_RE = re.compile(r"(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})")

#: The four feature families in the UChicago set. Kinetic features come from the
#: DCE time course and are the most protocol-dependent; size and morphology are
#: the most protocol-robust but also the least specific to the microenvironment.
FEATURE_FAMILIES: dict[str, tuple[str, ...]] = {
    "kinetic": ("K1", "K2", "K3", "K4", "K5", "K6", "K7"),
    "enhancement_variance": ("E1", "E2", "E3", "E4"),
    "texture": tuple(f"T{i}" for i in range(1, 15)),
    "morphology": ("G1", "G2", "G3", "M1", "M2", "M3", "S1", "S2", "S3", "S4", "S5"),
}


@dataclass
class RealCohort:
    """The joined cohort, ready for modelling."""

    patient_ids: tuple[str, ...]
    feature_names: tuple[str, ...]
    X: np.ndarray = field(repr=False)
    #: Continuous label: TIL-positive share of tissue patches.
    til_fraction: np.ndarray = field(repr=False)
    #: Spatial arrangement label: Moran's I of the TIL field over tissue.
    til_moran: np.ndarray = field(repr=False)
    #: Scanner key (vendor|model) for the patient's MR series, as a batch label.
    scanner: np.ndarray = field(repr=False)

    @property
    def n(self) -> int:
        return len(self.patient_ids)

    def family_of(self, feature: str) -> str:
        code = _family_code(feature)
        for family, codes in FEATURE_FAMILIES.items():
            if code in codes:
                return family
        return "other"

    def family_indices(self, family: str) -> np.ndarray:
        return np.array(
            [i for i, name in enumerate(self.feature_names) if self.family_of(name) == family]
        )

    def binary_label(self, quantile: float = 2 / 3) -> np.ndarray:
        """TIL-high vs TIL-low -- the ``inflamed_vs_rest`` analogue.

        The cut is a *fitted* parameter. It is computed here on the whole cohort
        only to report class balance; every model fit re-derives it inside the
        training fold (see ``scripts/run_real_experiment.py``), because a cut
        chosen with the test labels in view is the same leak this project
        measures elsewhere.
        """
        return (self.til_fraction >= np.quantile(self.til_fraction, quantile)).astype(int)


def _family_code(feature: str) -> str:
    match = re.search(r"\(([KETGMS]\d{1,2})\)\s*$", feature)
    return match.group(1) if match else ""


def load_features(path: Path | str) -> tuple[list[str], list[str], np.ndarray]:
    """Read the UChicago feature workbook into (patients, names, matrix)."""
    import pandas as pd

    frame = pd.read_excel(path)
    patients = frame["Lesion Name"].astype(str).str.extract(_PATIENT_RE.pattern)[0]
    numeric = frame.drop(columns=["Lesion Name"])
    return list(patients), list(numeric.columns), numeric.to_numpy(dtype=float)


def load_til_labels(path: Path | str) -> dict[str, dict]:
    payload = json.loads(Path(path).read_text())
    return {record["patient_id"]: record for record in payload["records"]}


def build(
    features_path: Path | str,
    til_labels_path: Path | str,
    scanner_map: dict[str, str] | None = None,
) -> RealCohort:
    """Join features to labels, keeping only patients present in both."""
    patients, names, matrix = load_features(features_path)
    labels = load_til_labels(til_labels_path)

    keep = [i for i, pid in enumerate(patients) if pid in labels]
    if not keep:
        raise ValueError("no patient appears in both the feature table and the label set")

    kept_ids = [patients[i] for i in keep]
    X = matrix[keep]

    # Constant or all-NaN columns carry no information and break standardisation.
    finite = np.isfinite(X).all(axis=0)
    varying = np.nanstd(X, axis=0) > 0
    usable = finite & varying
    if not usable.all():
        dropped = [n for n, keep_it in zip(names, usable) if not keep_it]
        print(f"  dropping {len(dropped)} non-varying/non-finite feature(s): {dropped}")

    return RealCohort(
        patient_ids=tuple(kept_ids),
        feature_names=tuple(n for n, keep_it in zip(names, usable) if keep_it),
        X=X[:, usable],
        til_fraction=np.array([labels[p]["til_fraction_tissue"] for p in kept_ids]),
        til_moran=np.array([labels[p]["til_moran_i"] for p in kept_ids]),
        scanner=np.array([(scanner_map or {}).get(p, "unknown") for p in kept_ids]),
    )


def fetch_scanner_map(patient_ids: set[str]) -> dict[str, str]:
    """Dominant MR scanner per patient, from TCIA series metadata.

    Used as the batch key for leave-one-scanner-out. The TCGA-BRCA imaging came
    from several tissue source sites on several vendors, which is what makes a
    cross-scanner estimate possible at all in a cohort this size.
    """
    from collections import Counter

    from ..clients.tcia import TCIAClient

    counts: dict[str, Counter] = {}
    for series in TCIAClient().series("TCGA-BRCA", modality="MR"):
        if series.patient_id in patient_ids:
            counts.setdefault(series.patient_id, Counter())[series.scanner_key] += 1
    return {pid: counter.most_common(1)[0][0] for pid, counter in counts.items()}
