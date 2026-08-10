"""The wavelet filter bank that turns 64 texture features into 768.

The reference study reports 833 features per region per DCE phase, decomposed as
``1 volume + 64 texture + 768 wavelet``, and describes the wavelet stage as
"Daubechies 2 (ratio = 1/2, 2/3, 3/2, or 2), Coiflets 1, and Symlets 4",
generating **12** transformed images with 64 texture features on each
(12 x 64 = 768).

Two things about that description need stating plainly rather than papering
over, because they change what a replication is actually replicating.

**1. The count.**  Six named filter configurations produce twelve images only if
each yields two.  The reading implemented here -- and the only one consistent
with 12 -- is one low-pass (approximation) and one high-pass (detail) image per
filter.  A full 3-D discrete wavelet transform instead yields eight sub-bands
per filter (LLL...HHH), which would give 48 images and 3,072 features.  Set
``full_3d=True`` for that variant; the default reproduces the published count.

**2. What 768 features from 182 patients means.**  With four phases and two
regions the reference design reaches 833 x 4 = 3,332 candidate features per
region against 182 patients: roughly **18 features per patient**.  Feature count
is not itself the problem -- penalised regression handles p >> n routinely.  The
problem is that a selection step run on all 182 patients gets 3,332 chances to
find a feature that happens to separate the 45 validation cases, and roughly
that many chances is enough to find one by luck alone.  That is what
``mri_cd8_til.experiments.optimism`` measures.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class WaveletFilter:
    """One entry in the filter bank."""

    name: str
    #: PyWavelets/PyRadiomics wavelet identifier.
    wavelet: str
    #: Decomposition level, standing in for the paper's unexplained "ratio".
    level: int = 1
    note: str = ""


#: The six filter configurations named in the reference study.  The four
#: Daubechies-2 entries are distinguished by the paper only as "ratio = 1/2, 2/3,
#: 3/2, or 2"; no toolkit exposes a wavelet parameter by that name, so they are
#: mapped here onto decomposition levels 1-4, which is the interpretation that
#: preserves both the filter identity and the count.  Flagged as an
#: interpretation, not a reconstruction of the original code.
PAPER_FILTER_BANK: tuple[WaveletFilter, ...] = (
    WaveletFilter("db2_r1", "db2", 1, "paper 'ratio = 1/2' -- interpreted as level 1"),
    WaveletFilter("db2_r2", "db2", 2, "paper 'ratio = 2/3' -- interpreted as level 2"),
    WaveletFilter("db2_r3", "db2", 3, "paper 'ratio = 3/2' -- interpreted as level 3"),
    WaveletFilter("db2_r4", "db2", 4, "paper 'ratio = 2'   -- interpreted as level 4"),
    WaveletFilter("coif1", "coif1", 1, ""),
    WaveletFilter("sym4", "sym4", 1, ""),
)

#: Sub-bands kept per filter under the 12-image reading.
PAPER_SUBBANDS: tuple[str, ...] = ("LLL", "HHH")

#: All sub-bands of a 3-D single-level DWT, used when ``full_3d=True``.
FULL_3D_SUBBANDS: tuple[str, ...] = (
    "LLL",
    "LLH",
    "LHL",
    "LHH",
    "HLL",
    "HLH",
    "HHL",
    "HHH",
)


def derived_image_names(full_3d: bool = False) -> list[str]:
    """Names of every wavelet-derived image the bank produces."""
    subbands = FULL_3D_SUBBANDS if full_3d else PAPER_SUBBANDS
    return [f"wavelet-{f.name}-{sb}" for f in PAPER_FILTER_BANK for sb in subbands]


def n_derived_images(full_3d: bool = False) -> int:
    return len(derived_image_names(full_3d))


def expected_feature_count(n_texture: int, full_3d: bool = False, include_volume: bool = True) -> int:
    """Total features per region per phase for a given texture-set size.

    With the paper's ``n_texture=64`` and the 12-image bank this returns 833,
    matching the published number exactly.
    """
    return int(include_volume) + n_texture + n_texture * n_derived_images(full_3d)


def pyradiomics_image_types(full_3d: bool = False) -> list[dict]:
    """Per-filter ``imageType`` settings blocks for PyRadiomics.

    PyRadiomics computes all eight 3-D sub-bands per invocation, so the caller
    filters the returned feature names down to :data:`PAPER_SUBBANDS` unless
    ``full_3d`` is set.  Each dict is one extractor configuration; they are run
    sequentially because ``wavelet`` is a scalar setting, not a list.
    """
    return [
        {
            "filter": f.name,
            "settings": {"wavelet": f.wavelet, "level": f.level, "start_level": 0},
            "keep_subbands": list(FULL_3D_SUBBANDS if full_3d else PAPER_SUBBANDS),
            "note": f.note,
        }
        for f in PAPER_FILTER_BANK
    ]
