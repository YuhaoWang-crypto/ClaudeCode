"""fftplsr -- FFT-PLSR combinatorial mutation modelling for enzyme engineering.

A cross-platform reimplementation of the FP4COM method (Hu et al., Nat Commun 16,
2025; doi:10.1038/s41467-025-61952-2), plus the baselines and held-out tests
needed to know when to trust it.

Typical use, starting from measured single-mutant activities:

    from fftplsr import datasets, design

    table, parent = datasets.trainset(1)
    report = design.design_round(
        parent=parent,
        measured=dict(zip(table["Variants"], table["Fitness"])),
        sites=datasets.ROUND12_SITES,
        n_rounds=1,
        pick=8,
    )
    print(report.picks)
"""

from . import baselines, datasets, design, encode, model, variants
from .encode import MutationEncoder, encode_sequences, spectrum
from .model import cv_score, fit_predict, nested_cv_r2, screen_indices
from .variants import apply_variant, enumerate_combinations, parse_variant, saturation_scan

__version__ = "0.1.0"

__all__ = [
    "baselines",
    "datasets",
    "design",
    "encode",
    "model",
    "variants",
    "MutationEncoder",
    "encode_sequences",
    "spectrum",
    "cv_score",
    "fit_predict",
    "nested_cv_r2",
    "screen_indices",
    "apply_variant",
    "enumerate_combinations",
    "parse_variant",
    "saturation_scan",
]
