"""Peptide tiling, 9-mer core enumeration and the self-tolerance filter.

The reasoning this module encodes
---------------------------------
A dog or cat is centrally tolerant to its *own* insulin. When a therapeutic
insulin is administered, the only MHC-II cores that can present a genuinely
novel signal to a naive T cell are the cores the recipient's own insulin does
not contain. Everything else has been seen in the thymus.

So the risk unit is: a 9-mer register inside the administered sequence that is
absent from the recipient's own insulin. This is set arithmetic on the
sequences -- no model, no training data, no species extrapolation. It is the
one part of the pipeline that is exactly as reliable for a cat as for a human.

Two deliberate scope limits:
  * Tiling is per chain. A two-chain insulin is two polypeptides held together
    by disulfides, so peptides spanning the B/A junction are not produced by
    antigen processing.
  * The tolerance reference is the recipient's insulin, not its whole proteome.
    A core that happens to recur elsewhere in the recipient proteome would also
    be tolerated; adding a proteome reference only ever *removes* risk calls,
    so this stays on the conservative side.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Sequence, Set

from .insulin import Insulin, diff

CORE_LENGTH = 9
PEPTIDE_LENGTH = 15


@dataclass(frozen=True)
class Core:
    """A 9-mer MHC-II binding core located on a specific chain."""

    chain: str
    start: int          # 1-based position of core residue P1 within the chain
    sequence: str

    @property
    def label(self) -> str:
        return f"{self.chain}{self.start}-{self.chain}{self.start + len(self.sequence) - 1}"

    def positions(self) -> List[int]:
        return list(range(self.start, self.start + len(self.sequence)))


@dataclass
class Peptide:
    """A tiled 15-mer, with the registers it can adopt."""

    chain: str
    start: int
    sequence: str
    cores: List[Core] = field(default_factory=list)


def tile(chain: str, seq: str, length: int = PEPTIDE_LENGTH) -> List[Peptide]:
    """Overlapping windows of ``length``; short chains yield the whole chain."""
    peptides: List[Peptide] = []
    if len(seq) <= length:
        windows = [(1, seq)]
    else:
        windows = [(i + 1, seq[i:i + length]) for i in range(len(seq) - length + 1)]
    for start, window in windows:
        pep = Peptide(chain, start, window)
        for off in range(len(window) - CORE_LENGTH + 1):
            pep.cores.append(Core(chain, start + off, window[off:off + CORE_LENGTH]))
        peptides.append(pep)
    return peptides


def all_cores(ins: Insulin) -> List[Core]:
    """Every 9-mer register in both chains, in order."""
    out: List[Core] = []
    for ch in ("A", "B"):
        seq = ins.chain(ch)
        for i in range(len(seq) - CORE_LENGTH + 1):
            out.append(Core(ch, i + 1, seq[i:i + CORE_LENGTH]))
    return out


def self_core_set(self_insulin: Insulin) -> Set[str]:
    """The 9-mers the recipient is centrally tolerant to."""
    return {c.sequence for c in all_cores(self_insulin)}


@dataclass
class NeoCore:
    """A core in the administered insulin that the recipient has never seen."""

    core: Core
    foreign_positions: List[int]         # chain positions differing from self
    foreign_residues: List[str]          # 'A8 T' style labels

    @property
    def n_foreign(self) -> int:
        return len(self.foreign_positions)


def neo_cores(drug: Insulin, self_insulin: Insulin) -> List[NeoCore]:
    """Cores present in ``drug`` but absent from the recipient's own insulin.

    A core is kept only if its exact 9-mer occurs nowhere in the recipient's
    insulin -- a shifted register that happens to reproduce a self 9-mer is
    tolerated and correctly dropped.
    """
    tolerated = self_core_set(self_insulin)
    differences = {(d.chain, d.position): d for d in diff(self_insulin, drug)}
    out: List[NeoCore] = []
    for core in all_cores(drug):
        if core.sequence in tolerated:
            continue
        pos = [p for p in core.positions() if (core.chain, p) in differences]
        labels = [differences[(core.chain, p)].label for p in pos]
        out.append(NeoCore(core, pos, labels))
    return out


def neo_peptides(drug: Insulin, self_insulin: Insulin,
                 length: int = PEPTIDE_LENGTH) -> List[Peptide]:
    """15-mers that contain at least one non-self core -- the NetMHCIIpan input set."""
    novel = {(nc.core.chain, nc.core.start) for nc in neo_cores(drug, self_insulin)}
    keep: List[Peptide] = []
    for ch in ("A", "B"):
        for pep in tile(ch, drug.chain(ch), length):
            if any((c.chain, c.start) in novel for c in pep.cores):
                keep.append(pep)
    return keep


def residue_coverage(drug: Insulin, cores: Sequence[NeoCore]) -> Dict[str, List[int]]:
    """How many non-self cores each residue of the drug participates in."""
    cov = {"A": [0] * len(drug.A), "B": [0] * len(drug.B)}
    for nc in cores:
        for p in nc.core.positions():
            if p <= len(cov[nc.core.chain]):
                cov[nc.core.chain][p - 1] += 1
    return cov


def foreign_burden(drug: Insulin, self_insulin: Insulin) -> Dict[str, object]:
    """Headline, model-free summary of how novel a product is to a recipient."""
    differences = diff(self_insulin, drug)
    cores = neo_cores(drug, self_insulin)
    total = len(all_cores(drug))
    return {
        "product": drug.name,
        "n_differences": len(differences),
        "differences": [d.label for d in differences],
        "n_cores_total": total,
        "n_neo_cores": len(cores),
        "frac_neo_cores": len(cores) / total if total else 0.0,
        "neo_core_labels": [c.core.label for c in cores],
    }
