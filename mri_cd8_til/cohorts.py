"""Registry of public cohorts that can supply breast MRI + immune labels.

Every count in :data:`REGISTRY` was measured against the live TCIA, GDC and GEO
APIs (see ``scripts/refresh_inventory.py``, which re-measures them and reports
drift).  They are recorded here so the design can be read without network
access -- not as a substitute for re-verification.

The registry is organised by *label tier*, because the binding constraint on
this problem is not how many breast MRIs exist -- roughly 3,000 are public --
but how many of them have a usable immune label attached:

``A``  spatially resolved TIL label + imaging on the same patient.
``B``  quantitative but non-spatial (transcriptomic) immune label + imaging.
``C``  imaging only: useful for segmentation, harmonisation and domain shift,
       useless as a supervised immune-phenotype cohort.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum


class LabelTier(str, Enum):
    SPATIAL_TIL = "A"
    TRANSCRIPTOMIC = "B"
    IMAGING_ONLY = "C"


@dataclass(frozen=True)
class Cohort:
    """One public cohort, with everything needed to plan a multi-centre study."""

    key: str
    name: str
    tier: LabelTier
    n_patients: int
    n_with_mri: int
    #: Patients that also carry the label modality named in ``label_source``.
    n_paired: int
    label_source: str
    imaging: str
    centers: str
    scanner_vendors: tuple[str, ...]
    access: str
    identifier_link: str
    caveats: tuple[str, ...] = field(default_factory=tuple)

    @property
    def usable_for_supervision(self) -> bool:
        return self.tier is not LabelTier.IMAGING_ONLY and self.n_paired > 0


#: Verified 2026-08-10 against services.cancerimagingarchive.net, api.gdc.cancer.gov
#: and eutils.ncbi.nlm.nih.gov.  ``scripts/refresh_inventory.py`` re-checks these.
REGISTRY: dict[str, Cohort] = {
    "tcga_brca": Cohort(
        key="tcga_brca",
        name="TCGA-BRCA (TCIA imaging + GDC pathology/transcriptome)",
        tier=LabelTier.SPATIAL_TIL,
        n_patients=139,
        n_with_mri=137,
        n_paired=136,
        label_source=(
            "Saltz TIL maps (TIL-WSI-TCGA, DOI 10.7937/K9/TCIA.2018.Y75F9W1) over the "
            "diagnostic FFPE WSI, plus RNA-seq for CD8A/CD8B and deconvolution"
        ),
        imaging="Pre-operative breast DCE-MRI (1.5T and 3T, mixed protocols)",
        centers="MSKCC, Mayo Clinic, UPMC, Roswell Park and others (TCGA tissue source sites)",
        scanner_vendors=("GE SIGNA EXCITE", "GE SIGNA HDx", "SIEMENS Sonata", "Philips Achieva"),
        access="Fully open access; no data-use agreement required",
        identifier_link="TCIA PatientID == GDC submitter_id (e.g. TCGA-AO-A03M)",
        caveats=(
            "Small: 136 fully paired cases is a validation cohort, not a training cohort.",
            "Acquisition protocols are heterogeneous by design -- that is the point, "
            "but it means per-voxel intensity features need harmonisation.",
            "H&E-derived TIL is not CD8 IHC; it counts all lymphocytes.",
        ),
    ),
    "ispy2": Cohort(
        key="ispy2",
        name="I-SPY2 (TCIA ISPY2 imaging + GEO GSE194040 transcriptome)",
        tier=LabelTier.TRANSCRIPTOMIC,
        n_patients=719,
        n_with_mri=719,
        n_paired=717,
        label_source=(
            "GSE194040 pre-treatment Agilent 44K expression (n=987) -> CD8A/CD8B, "
            "immune signatures; plus pCR and treatment arm (incl. pembrolizumab, n=69)"
        ),
        imaging="Serial breast DCE-MRI across the neoadjuvant trial",
        centers="~20 US sites in the I-SPY2 TRIAL",
        scanner_vendors=("GE", "SIEMENS", "Philips"),
        access="TCIA open access + GEO open access",
        identifier_link="TCIA 'ISPY2-100899' == GEO !Sample_title 'ISPY2_100899'",
        caveats=(
            "Transcriptomic immune labels carry no spatial information at all: they "
            "cannot distinguish immune-excluded from inflamed, which is exactly the "
            "distinction the peritumoral model claims to make.",
            "Expression is from a pre-treatment core biopsy; the MRI images the whole "
            "tumour. Sampling mismatch is a real, non-trivial noise source.",
        ),
    ),
    "ispy1": Cohort(
        key="ispy1",
        name="I-SPY1 / ACRIN 6657 (TCIA ISPY1)",
        tier=LabelTier.TRANSCRIPTOMIC,
        n_patients=222,
        n_with_mri=222,
        n_paired=0,
        label_source="GSE22226 expression (n=150) exists but is not ID-linkable as published",
        imaging="Serial breast DCE-MRI, 4 timepoints",
        centers="9 ACRIN sites",
        scanner_vendors=("GE GENESIS_SIGNA", "SIEMENS Sonata", "SIEMENS MAGNETOM VISION", "Philips Intera"),
        access="TCIA open access",
        identifier_link=(
            "BROKEN: TCIA uses 'ISPY1_1001' (4-digit), GSE22226 uses 6-digit study IDs; "
            "measured overlap by naive numeric match is 0. Needs the ACRIN crosswalk."
        ),
        caveats=(
            "Counts as tier C in practice until the identifier crosswalk is obtained.",
        ),
    ),
    "acrin_6698": Cohort(
        key="acrin_6698",
        name="ACRIN 6698 (I-SPY2 DWI sub-study)",
        tier=LabelTier.IMAGING_ONLY,
        n_patients=385,
        n_with_mri=385,
        n_paired=0,
        label_source="none directly; would inherit I-SPY2 labels if IDs could be linked",
        imaging="Breast diffusion-weighted MRI (ADC)",
        centers="I-SPY2 sites",
        scanner_vendors=("GE", "SIEMENS", "Philips"),
        access="TCIA open access",
        identifier_link=(
            "BROKEN: 'ACRIN-6698-102212' vs 'ISPY2-100899' -- measured overlap of the "
            "numeric parts is 0, so the two TCIA collections are separately de-identified."
        ),
        caveats=("ADC would add real biological signal for TIL if the link were recoverable.",),
    ),
    "duke": Cohort(
        key="duke",
        name="Duke-Breast-Cancer-MRI",
        tier=LabelTier.IMAGING_ONLY,
        n_patients=922,
        n_with_mri=922,
        n_paired=0,
        label_source="Clinical/pathology tables (receptor status, stage); no TIL, no CD8",
        imaging="Pre-operative axial breast DCE-MRI, 4 post-contrast phases",
        centers="Single institution (Duke), 6+ distinct scanners",
        scanner_vendors=(
            "GE SIGNA HDx",
            "GE Signa HDxt",
            "GE Optima MR450w",
            "SIEMENS Avanto",
            "SIEMENS TrioTim",
            "SIEMENS Skyra",
        ),
        access="TCIA open access",
        identifier_link="n/a",
        caveats=(
            "Best available single source of scanner-induced batch effect within one "
            "institution: ideal for measuring how much of a radiomic signature is scanner.",
        ),
    ),
    "advanced_mri_breast": Cohort(
        key="advanced_mri_breast",
        name="Advanced-MRI-Breast-Lesions",
        tier=LabelTier.IMAGING_ONLY,
        n_patients=632,
        n_with_mri=632,
        n_paired=0,
        label_source="Lesion-level annotations; no immune label",
        imaging="Breast DCE-MRI with expert segmentations",
        centers="Single institution",
        scanner_vendors=("GE Signa HDxt",),
        access="TCIA open access",
        identifier_link="n/a",
        caveats=("Single scanner model -- useful as a homogeneous control arm.",),
    ),
    "nact_pilot": Cohort(
        key="nact_pilot",
        name="Breast-MRI-NACT-Pilot",
        tier=LabelTier.IMAGING_ONLY,
        n_patients=64,
        n_with_mri=64,
        n_paired=0,
        label_source="Response outcome only",
        imaging="Breast DCE-MRI across neoadjuvant chemotherapy",
        centers="Single institution (UCSF)",
        scanner_vendors=("GE GENESIS_SIGNA",),
        access="TCIA open access",
        identifier_link="n/a",
        caveats=(
            "Same size (n=64) as the neoadjuvant cohort in the paper being reproduced, "
            "which makes it a natural stand-in for that arm.",
        ),
    ),
    "qin_breast": Cohort(
        key="qin_breast",
        name="QIN-BREAST",
        tier=LabelTier.IMAGING_ONLY,
        n_patients=68,
        n_with_mri=41,
        n_paired=0,
        label_source="none",
        imaging="Breast DCE-MRI + PET/CT",
        centers="QIN network",
        scanner_vendors=("Philips Achieva", "GE Discovery STE"),
        access="TCIA open access",
        identifier_link="n/a",
        caveats=("Small; multimodal.",),
    ),
    "mama_mia": Cohort(
        key="mama_mia",
        name="MAMA-MIA (curated segmentations over ISPY1/ISPY2/NACT/DUKE)",
        tier=LabelTier.IMAGING_ONLY,
        n_patients=1506,
        n_with_mri=1506,
        n_paired=0,
        label_source="Expert tumour segmentations + 49 harmonised clinical variables",
        imaging="Pre-treatment breast DCE-MRI pooled from four TCIA collections",
        centers="Multi-centre by construction (the four source collections)",
        scanner_vendors=("GE", "SIEMENS", "Philips"),
        access="Open (Synapse/GitHub, CC-BY); see github.com/LidiaGarrucho/MAMA-MIA",
        identifier_link="Retains the source-collection patient IDs, so it joins to ISPY2/DUKE",
        caveats=(
            "Derived from the same patients as ispy1/ispy2/nact_pilot/duke -- do NOT add "
            "its 1506 to the other counts, it is a re-annotation, not new patients.",
            "Its value here is removing the segmentation bottleneck and supplying a "
            "pretrained nnU-Net, not supplying labels.",
        ),
    ),
}


def by_tier(tier: LabelTier) -> list[Cohort]:
    return [c for c in REGISTRY.values() if c.tier is tier]


def supervised_cohorts() -> list[Cohort]:
    """Cohorts that can actually train or validate an immune-phenotype model."""
    return [c for c in REGISTRY.values() if c.usable_for_supervision]


#: Collections whose patients are re-annotations or sub-studies of patients
#: already counted elsewhere in the registry.
DERIVED_OR_OVERLAPPING = {
    # MAMA-MIA re-annotates ispy1/ispy2/nact_pilot/duke; not new patients.
    "mama_mia": "re-annotation of four collections already counted",
    # ACRIN-6698 is the I-SPY2 DWI sub-study, so its patients are very likely
    # the same individuals as ispy2 entries under different de-identified IDs.
    # The IDs cannot be matched (measured overlap 0), so this cannot be proven
    # either way -- which is precisely why it is excluded from the lower bound.
    "acrin_6698": "I-SPY2 DWI sub-study; probably the same individuals as ispy2",
}


def total_public_mri(count_overlapping: bool = False) -> int:
    """Public breast-MRI patient records.

    Returns a **lower bound** by default, excluding collections that are
    re-annotations or sub-studies of patients already counted.  Pass
    ``count_overlapping=True`` for the upper bound -- the raw sum of records,
    which double-counts anyone appearing in two collections under different
    de-identified IDs.

    The two numbers differ by ACRIN-6698's 385 patients. Reporting the upper
    bound as if it were a distinct-patient count would be the same kind of
    unstated-assumption error this project exists to avoid.
    """
    return sum(
        cohort.n_with_mri
        for cohort in REGISTRY.values()
        # MAMA-MIA is always excluded: it adds no patients under any reading.
        if cohort.key != "mama_mia"
        and (count_overlapping or cohort.key not in DERIVED_OR_OVERLAPPING)
    )


#: The cohort in the paper this project reproduces, for side-by-side reporting.
REFERENCE_STUDY = {
    "citation": (
        "Frontiers in Immunology 2022;13:1080048 -- Radiomic models based on MRI predict "
        "the spatial distribution of CD8+ tumor-infiltrating lymphocytes in breast cancer"
    ),
    "n_total": 182,
    "n_train": 137,
    "n_validation": 45,
    "n_nact": 64,
    "centers": 1,
    "external_validation": False,
    "field_strength": "3T Philips Achieva",
    "voxel_mm": (0.6, 0.6, 2.0),
    "n_post_contrast_phases": 4,
    "phase_interval_s": 74,
    "features_per_region_per_phase": 833,  # 1 volume + 64 texture + 768 wavelet
    "wavelet_filters": 12,
    "reported_auc": {
        "whole_tumour_train": 0.973,
        "whole_tumour_validation": 0.985,
        "peritumoural_train": 0.993,
        "peritumoural_validation": 0.984,
    },
}
