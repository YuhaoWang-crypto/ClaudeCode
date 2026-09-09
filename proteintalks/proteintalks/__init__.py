"""Reproduction of ProteinTalks, a perturbation-proteomics virtual cell model.

Sun R., Qian L., Li Y. et al., bioRxiv 2025.02.07.637070;
Nature, doi:10.1038/s41586-026-11001-9.
"""

from .data import (
    PerturbationDataset,
    SyntheticPerturbationProteome,
    load_real,
    make_splits,
)
from .evaluate import bootstrap_ci, classification_metrics, trajectory_metrics
from .model import ProteinTalks
from .multitask import GradientConflictWeighter
from .train import predict, train_proteintalks

__version__ = "0.1.0"

__all__ = [
    "ProteinTalks",
    "GradientConflictWeighter",
    "PerturbationDataset",
    "SyntheticPerturbationProteome",
    "load_real",
    "make_splits",
    "train_proteintalks",
    "predict",
    "classification_metrics",
    "trajectory_metrics",
    "bootstrap_ci",
]
