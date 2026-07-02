"""mdscreen — run making-it-rain-style molecular dynamics on Modal GPUs to
screen small molecules for activity, allosteric effects, and enzyme
inhibition/activation.

Public modules:
  config     - SimConfig / ScreenConfig dataclasses
  prepare    - protein/ligand parameterisation & solvation (OpenMM/OpenFF)
  simulate   - OpenMM MD engine (minimise -> NVT -> NPT -> production)
  analyze    - RMSD/RMSF/Rg/contacts/H-bonds trajectory analysis
  binding    - end-state MM/GBSA binding free energy for ranking
  allostery  - DCCM + apo/holo comparison for allosteric detection
  pipeline   - high-level workflows and screening campaign
"""
from .config import SimConfig, ScreenConfig  # noqa: F401

__version__ = "0.1.0"
__all__ = ["SimConfig", "ScreenConfig"]
