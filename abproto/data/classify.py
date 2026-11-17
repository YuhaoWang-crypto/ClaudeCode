"""Classify polymer entities into heavy / light / antigen and assemble
antibody-antigen complexes.

Classification is description + length based. It is deliberately conservative:
a complex is kept only if it has at least one antibody chain AND one non-antibody
protein antigen. CDR extraction is left to ANARCI in production (see note); for
split/dedup we operate on full VH/VL sequences, which needs no numbering.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Optional

from .rcsb import Entry

_HEAVY = re.compile(r"\b(heavy chain|vhh|nanobody|heavy|hc\b)", re.I)
_LIGHT = re.compile(r"\b(light chain|lambda|kappa|light|lc\b)", re.I)
_ANTIBODY = re.compile(r"(fab|fv|scfv|antibody|immunoglobulin|nanobody|vhh|igg)", re.I)
_GPCR = re.compile(r"(receptor|gpcr|rhodopsin|chemokine receptor|adrenergic|opioid|"
                   r"muscarinic|adenosine receptor|angiotensin)", re.I)


@dataclass
class Chain:
    entity_id: str
    role: str          # "heavy" | "light" | "antigen"
    description: str
    sequence: str


@dataclass
class AbAgComplex:
    pdb_id: str
    deposit_date: str
    resolution: Optional[float]
    heavy: list[Chain] = field(default_factory=list)
    light: list[Chain] = field(default_factory=list)
    antigen: list[Chain] = field(default_factory=list)
    is_gpcr: bool = False

    @property
    def antibody_seq(self) -> str:
        return "".join(c.sequence for c in (self.heavy + self.light))

    @property
    def antigen_seq(self) -> str:
        # Longest antigen chain is the representative antigen.
        return max((c.sequence for c in self.antigen), key=len, default="")

    def summary(self) -> dict:
        return {
            "pdb_id": self.pdb_id,
            "deposit_date": self.deposit_date,
            "resolution": self.resolution,
            "n_heavy": len(self.heavy), "n_light": len(self.light),
            "antigen": (self.antigen[0].description if self.antigen else ""),
            "is_gpcr": self.is_gpcr,
            "antibody_len": len(self.antibody_seq),
            "antigen_len": len(self.antigen_seq),
        }


def _looks_antibody(desc: str, seq: str) -> bool:
    if _ANTIBODY.search(desc):
        return True
    # Fallback: Fab chains are ~210-235 aa; VHH ~120. Weak signal, needs the
    # description to also not look like a named receptor/enzyme.
    return False


def classify_entry(entry: Entry) -> Optional[AbAgComplex]:
    """Return an AbAgComplex if the entry is a valid Ab-Ag complex, else None."""
    cx = AbAgComplex(pdb_id=entry.pdb_id, deposit_date=entry.deposit_date,
                     resolution=entry.resolution)
    for ent in entry.entities:
        if ent.polymer_type and "Protein" not in ent.polymer_type:
            continue  # ignore nucleic acids, sugars
        desc, seq = ent.description, ent.sequence
        if not seq:
            continue
        if _looks_antibody(desc, seq):
            chain = Chain(ent.entity_id, "?", desc, seq)
            if _HEAVY.search(desc) and not _LIGHT.search(desc):
                chain.role = "heavy"; cx.heavy.append(chain)
            elif _LIGHT.search(desc) and not _HEAVY.search(desc):
                chain.role = "light"; cx.light.append(chain)
            else:
                # Ambiguous Fab/Fv entity: bucket by length heuristic.
                chain.role = "heavy" if len(seq) >= 215 else "light"
                (cx.heavy if chain.role == "heavy" else cx.light).append(chain)
        else:
            cx.antigen.append(Chain(ent.entity_id, "antigen", desc, seq))
            if _GPCR.search(desc):
                cx.is_gpcr = True

    has_ab = bool(cx.heavy or cx.light)
    has_ag = bool(cx.antigen)
    return cx if (has_ab and has_ag) else None


def to_complexes(entries: list[Entry]) -> tuple[list[AbAgComplex], list[str]]:
    """Return (valid complexes, skipped_pdb_ids)."""
    good, skipped = [], []
    for e in entries:
        cx = classify_entry(e)
        if cx is None:
            skipped.append(e.pdb_id)
        else:
            good.append(cx)
    return good, skipped
