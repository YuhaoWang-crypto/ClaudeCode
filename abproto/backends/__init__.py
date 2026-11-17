"""Pluggable learned backends.

Two roles, two protocols (see base.py):
  - Generator : produces candidate antibody-antigen complexes (Multiflow / La-Proteina)
  - Folder    : (re)folds a sequence into a complex + confidence (Boltz-2 / AF3)

Import the concrete backend you want:
    from abproto.backends.mock import MockGenerator, MockFolder
    from abproto.backends.boltz import BoltzFolder
    from abproto.backends.multiflow import MultiflowGenerator
"""

from .base import Generator, Folder

__all__ = ["Generator", "Folder"]
