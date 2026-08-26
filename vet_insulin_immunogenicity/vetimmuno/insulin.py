"""Insulin drug substances and the cross-species residue difference map.

The mature A/B chains of insulin are 21 and 30 residues and are strictly
colinear across mammals, so "alignment" here is positional comparison -- no
gap penalties, no ambiguity. Engineered analogues are expressed as edits on a
named parent so that the provenance of every residue stays explicit.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Sequence, Tuple

from . import data

_EDIT = re.compile(r"^(?P<chain>[AB])(?P<pos>\d+)(?P<ref>[A-Z*])>(?P<alt>[A-Z-]+)$")


@dataclass
class Insulin:
    """A two-chain insulin, mature A and B chains only."""

    name: str
    A: str
    B: str
    source: str = ""
    note: str = ""

    def chain(self, which: str) -> str:
        return self.A if which == "A" else self.B

    def residues(self) -> List[Tuple[str, int, str]]:
        return [("A", i, a) for i, a in enumerate(self.A, 1)] + [
            ("B", i, a) for i, a in enumerate(self.B, 1)
        ]


@dataclass
class Difference:
    chain: str
    position: int
    self_aa: str
    drug_aa: str

    @property
    def label(self) -> str:
        return f"{self.chain}{self.position}{self.self_aa}>{self.drug_aa}"


@dataclass
class ProductSpec:
    """Declarative product definition, as it appears in the species config."""

    name: str
    parent: str
    edits: Sequence[str] = field(default_factory=tuple)
    note: str = ""


def natural_insulin(species: str) -> Insulin:
    rec = data.insulin_chains(species)
    return Insulin(
        name=f"{species} insulin",
        A=str(rec["A"]),
        B=str(rec["B"]),
        source=f"UniProt {rec['accession']} ({rec['entry_name']})",
    )


def apply_edits(parent: Insulin, edits: Sequence[str], name: str, note: str = "") -> Insulin:
    """Apply ``A21N>G`` / ``B30T>-`` / ``B31*>RR`` style edits to a parent insulin.

    ``>-`` deletes the residue (des-B30 analogues); a ``*`` reference residue
    means "append here" and is used for C-terminal extensions such as the
    glargine B31-B32 Arg-Arg.
    """
    chains = {"A": list(parent.A), "B": list(parent.B)}
    for edit in edits:
        m = _EDIT.match(edit)
        if not m:
            raise ValueError(f"malformed edit {edit!r} (expected e.g. A21N>G)")
        chain, pos = m["chain"], int(m["pos"])
        ref, alt = m["ref"], m["alt"]
        seq = chains[chain]
        if ref == "*":
            if pos != len(seq) + 1:
                raise ValueError(f"{edit}: extension must start at position {len(seq) + 1}")
            seq.extend(alt)
            continue
        if not 1 <= pos <= len(seq):
            raise ValueError(f"{edit}: position outside chain {chain} (len {len(seq)})")
        if seq[pos - 1] != ref:
            raise ValueError(f"{edit}: parent has {seq[pos - 1]} at {chain}{pos}, not {ref}")
        seq[pos - 1] = "" if alt == "-" else alt
    return Insulin(
        name=name,
        A="".join(chains["A"]),
        B="".join(chains["B"]),
        source=f"{parent.source} + {', '.join(edits)}" if edits else parent.source,
        note=note,
    )


def build_product(spec: ProductSpec) -> Insulin:
    parent = natural_insulin(spec.parent)
    if not spec.edits:
        return Insulin(spec.name, parent.A, parent.B, parent.source, spec.note)
    return apply_edits(parent, spec.edits, spec.name, spec.note)


def diff(self_insulin: Insulin, drug: Insulin) -> List[Difference]:
    """Residues in ``drug`` that the recipient's immune system has never seen.

    Positions beyond the length of the recipient's own chain (analogue
    extensions) count as differences; a residue deleted in the drug cannot be
    presented and so is not a difference.
    """
    out: List[Difference] = []
    for ch in ("A", "B"):
        s, d = self_insulin.chain(ch), drug.chain(ch)
        for i in range(len(d)):
            self_aa = s[i] if i < len(s) else "-"
            if self_aa != d[i]:
                out.append(Difference(ch, i + 1, self_aa, d[i]))
    return out


def difference_matrix(species: Sequence[str]) -> Dict[Tuple[str, str], List[Difference]]:
    """All pairwise natural-insulin differences, keyed ``(recipient, donor)``."""
    ins = {sp: natural_insulin(sp) for sp in species}
    return {
        (rec, don): diff(ins[rec], ins[don])
        for rec in species
        for don in species
        if rec != don
    }


def foreign_positions(self_insulin: Insulin, drug: Insulin) -> Dict[str, set]:
    out: Dict[str, set] = {"A": set(), "B": set()}
    for d in diff(self_insulin, drug):
        out[d.chain].add(d.position)
    return out


def linear_sequence(ins: Insulin, linker: str = "") -> Tuple[str, List[Tuple[str, int]]]:
    """Concatenate B then A into one string, with a position index.

    NOT used by the default pipeline. A two-chain insulin drug substance is two
    separate polypeptides held together by disulfides, so peptides spanning the
    B/A junction do not exist after processing and tiling is done per chain.
    This helper exists for proinsulin-style analyses, where the chains really
    are contiguous; ``linker`` stands in for the C-peptide and its positions are
    indexed as ``("-", 0)``.
    """
    index: List[Tuple[str, int]] = []
    seq = []
    for i, aa in enumerate(ins.B, 1):
        seq.append(aa)
        index.append(("B", i))
    for aa in linker:
        seq.append(aa)
        index.append(("-", 0))
    for i, aa in enumerate(ins.A, 1):
        seq.append(aa)
        index.append(("A", i))
    return "".join(seq), index
