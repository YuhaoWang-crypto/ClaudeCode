"""
design.py -- the deterministic chimera-design engine.  ✅ RIGOROUS

Implements, verbatim, the recipe from the paper's Methods section
"Design of chimeric proteins":

    "a[n] arbitrary residue located within a loop region was chosen as the
     permutation site to generate the circularly permuted binder. This amino
     acid was removed to establish novel N and C termini. The native N and C
     termini of the binder were joined by a flexible glycine-serine linker of
     sufficient length. This permuted binder was inserted at position 253 of
     beta-lactamase TEM-1, with a single glycine residue used as a linker
     between TEM-1 and the binder."

Everything here is pure string manipulation -> 100% reproducible.  No
randomness, no external services.  Given the same inputs it always returns
the same sequence.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


# ----------------------------------------------------------------------------
# 1. Circular permutation of the receptor                              ✅
# ----------------------------------------------------------------------------
def circular_permute(seq: str, site: int, gs_linker: str) -> str:
    """Return the circularly permuted sequence of ``seq``.

    Parameters
    ----------
    seq : str
        The receptor (ligand-binding domain) amino-acid sequence, 1-letter.
    site : int
        0-indexed position of the loop residue chosen as the permutation
        site.  This residue is *removed* (it becomes the break point that
        creates the new N and C termini), exactly as in the paper.
    gs_linker : str
        Flexible Gly/Ser linker that joins the *native* N- and C-termini of
        the receptor "of sufficient length".

    Returns
    -------
    str
        cpR = seq[site+1:]  +  gs_linker  +  seq[:site]

        i.e. the chain is reopened just after ``site``:
          * new N-terminus  = residue (site+1) of the original,
          * new C-terminus  = residue (site-1) of the original,
          * the original termini (seq[0] and seq[-1]) are now internal and
            are bridged by ``gs_linker``.
    """
    if not 0 <= site < len(seq):
        raise ValueError(f"permutation site {site} out of range 0..{len(seq)-1}")
    if not seq[:site] or not seq[site + 1:]:
        raise ValueError(
            "permutation site too close to a terminus; both flanking "
            "segments must be non-empty"
        )
    n_part = seq[site + 1:]      # becomes the new N-terminal segment
    c_part = seq[:site]          # becomes the new C-terminal segment (native N-term at its start)
    return n_part + gs_linker + c_part


# ----------------------------------------------------------------------------
# 2. Insertion of the permuted receptor into the reporter enzyme       ✅
# ----------------------------------------------------------------------------
@dataclass
class Chimera:
    """A fully specified single-component allosteric switch construct."""
    name: str
    sequence: str
    # provenance / design record (so every residue is traceable)
    reporter_name: str
    receptor_name: str
    insertion_index: int          # 0-indexed cut point in the reporter
    permutation_site: int         # 0-indexed removed residue in the receptor
    gs_linker: str
    flank_linker: str             # linker between reporter and receptor (paper: single "G")
    receptor_domain_span: tuple   # (start, end) 0-indexed, inclusive/exclusive of the cpReceptor within the chimera
    meta: dict = field(default_factory=dict)

    @property
    def length(self) -> int:
        return len(self.sequence)

    def fasta(self) -> str:
        return f">{self.name}\n{self.sequence}\n"


def build_chimera(
    *,
    name: str,
    reporter_name: str,
    reporter_seq: str,
    receptor_name: str,
    receptor_seq: str,
    insertion_index: int,
    permutation_site: int,
    gs_linker: str = "GGSGGSGGS",
    flank_linker: str = "G",
) -> Chimera:
    """Build one chimeric biosensor sequence.

    chimera = reporter[:q]  +  L  +  cpReceptor  +  L  +  reporter[q:]

    where q = insertion_index and L = flank_linker (paper default: a single
    glycine).  cpReceptor is produced by :func:`circular_permute`.

    Returns a :class:`Chimera` carrying the full design provenance.
    """
    if not 0 < insertion_index < len(reporter_seq):
        raise ValueError(
            f"insertion index {insertion_index} out of range 1..{len(reporter_seq)-1}"
        )
    cp = circular_permute(receptor_seq, permutation_site, gs_linker)

    head = reporter_seq[:insertion_index]
    tail = reporter_seq[insertion_index:]
    insert = flank_linker + cp + flank_linker

    seq = head + insert + tail
    dom_start = len(head) + len(flank_linker)
    dom_end = dom_start + len(cp)

    return Chimera(
        name=name,
        sequence=seq,
        reporter_name=reporter_name,
        receptor_name=receptor_name,
        insertion_index=insertion_index,
        permutation_site=permutation_site,
        gs_linker=gs_linker,
        flank_linker=flank_linker,
        receptor_domain_span=(dom_start, dom_end),
        meta={
            "reporter_len": len(reporter_seq),
            "receptor_len": len(receptor_seq),
            "cp_len": len(cp),
        },
    )


# ----------------------------------------------------------------------------
# 3. Round-trip check -- proves the construction is lossless           ✅
# ----------------------------------------------------------------------------
def verify_chimera(ch: Chimera, reporter_seq: str, receptor_seq: str) -> dict:
    """Structural-accounting checks on a built chimera.

    These are *exact* invariants (not heuristics): the chimera must contain
    the reporter split cleanly around the insert, and the circularly-permuted
    receptor must contain every original receptor residue exactly once
    (minus the single removed permutation-site residue).
    """
    q = ch.insertion_index
    head, tail = reporter_seq[:q], reporter_seq[q:]

    # reporter is preserved and contiguous on both sides of the insert
    reporter_ok = ch.sequence.startswith(head) and ch.sequence.endswith(tail)

    # the cp receptor sits in the recorded span
    ds, de = ch.receptor_domain_span
    cp = ch.sequence[ds:de]
    cp_expected = circular_permute(receptor_seq, ch.permutation_site, ch.gs_linker)
    cp_ok = cp == cp_expected

    # residue conservation: cp receptor = all receptor residues except the
    # removed permutation site, plus the GS linker
    removed = receptor_seq[ch.permutation_site]
    body = cp.replace(ch.gs_linker, "", 1)
    conserved = sorted(body) == sorted(
        receptor_seq[:ch.permutation_site] + receptor_seq[ch.permutation_site + 1:]
    )

    return {
        "reporter_preserved": reporter_ok,
        "cp_receptor_correct": cp_ok,
        "residues_conserved": conserved,
        "removed_residue": removed,
        "chimera_len": ch.length,
        "all_ok": reporter_ok and cp_ok and conserved,
    }
