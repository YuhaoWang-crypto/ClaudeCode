r"""Epitope-level homology (optimisation #7).

Antigen-sequence clustering alone leaks when two different antigens share the
epitope fold. When antibody + antigen coordinates are available, we define the
epitope as the antigen residues contacting the antibody and cluster on that
epitope's residues, so two complexes binding the same site land on the same side
of the split even across antigens.

Coordinates are CA traces here; a production system uses heavy-atom contacts.
"""

from __future__ import annotations

import numpy as np


def epitope_residues(antibody_ca: np.ndarray, antigen_ca: np.ndarray,
                     contact: float = 10.0) -> np.ndarray:
    """Indices into the antigen of residues within ``contact`` A of the antibody."""
    antibody_ca = np.asarray(antibody_ca, float)
    antigen_ca = np.asarray(antigen_ca, float)
    if antibody_ca.size == 0 or antigen_ca.size == 0:
        return np.array([], dtype=int)
    d = np.linalg.norm(antigen_ca[:, None, :] - antibody_ca[None, :, :], axis=-1)
    return np.where(d.min(axis=1) <= contact)[0]


def epitope_signature(antigen_seq: str, epitope_idx: np.ndarray) -> str:
    """The epitope residues as a position-sorted string (for clustering)."""
    idx = sorted(int(i) for i in epitope_idx if 0 <= int(i) < len(antigen_seq))
    return "".join(antigen_seq[i] for i in idx)


def epitope_key(antibody_ca: np.ndarray, antigen_ca: np.ndarray, antigen_seq: str,
                contact: float = 10.0) -> str:
    """Convenience: epitope signature straight from coordinates + antigen sequence."""
    return epitope_signature(antigen_seq, epitope_residues(antibody_ca, antigen_ca, contact))
