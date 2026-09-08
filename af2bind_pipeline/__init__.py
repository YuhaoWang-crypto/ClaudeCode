"""AF2BIND: small-molecule binding-site prediction from the AlphaFold2 pair representation.

Headless reimplementation of the method described in:

    Gazizov A, Lian A, Goverde CA, Mou J, Ovchinnikov S, Polizzi NF.
    "AF2BIND: predicting small-molecule binding sites using the pair
    representation of AlphaFold2." Nature Methods (2026).
    doi:10.1038/s41592-026-03011-2   (preprint doi:10.1101/2023.10.15.562410)

Upstream reference implementation (MIT): https://github.com/sokrypton/af2bind

The pipeline is split so that the expensive part runs once:

    structure --(GPU, AlphaFold2 forward pass)--> pair features (L, 5120)
    pair features --(pure numpy logistic head)--> p(bind) per residue

`core` and everything downstream of it is numpy-only, so cached features can be
re-scored with different weights/seeds on any machine with no GPU and no jax.
"""

__all__ = ["core", "weights", "structure", "pockets", "validate"]

__version__ = "1.0.0"
