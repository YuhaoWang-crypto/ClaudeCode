"""
architectures.py -- extra biosensor architectures integrated from the
luminescent/complementation reporter literature (see
reference/reporter-architectures.md):

  #4  split-fragment complementation (NanoBiT): two constructs
      (LgBiT-fusion + SmBiT-fusion) whose luminescence reports LIGAND-INDUCED
      PROXIMITY -- the two-component / auxiliary-binding-domain design.
  #5  ligand-gated modules (TrpR / DtxR / MetJ) as orthogonal inputs, and
      multi-input INTRAMOLECULAR AND gates (the paper's logic-gate architecture:
      two orthogonal receptors inserted at two reporter loops).

Everything here is deterministic sequence construction (✅). Whether a given
fusion actually switches / computes the logic is ⚠️ until a wet-lab luminescence
assay -- the same honesty contract as the rest of the package.
"""

from __future__ import annotations
from dataclasses import dataclass, field

from .design import circular_permute


# ===========================================================================
# #4  Split-fragment complementation (NanoBiT)
# ===========================================================================
@dataclass
class SplitConstruct:
    name: str
    sequence: str
    fusion_partner: str          # binder name
    fragment: str                # "LgBiT" | "SmBiT"
    orientation: str             # "N" (fragment-linker-binder) | "C" (binder-linker-fragment)
    length: int


@dataclass
class SplitPair:
    analyte: str
    linker: str
    lgbit_construct: SplitConstruct
    smbit_construct: SplitConstruct
    note: str = ""


def _fuse(fragment_seq, binder_seq, linker, orientation):
    return (fragment_seq + linker + binder_seq) if orientation == "N" \
        else (binder_seq + linker + fragment_seq)


def build_split_complementation(
    *, analyte: str, binder_a_name: str, binder_a_seq: str,
    binder_b_name: str, binder_b_seq: str,
    lgbit_seq: str, smbit_seq: str, linker: str = "GGSGGSGG",
    orientation_a: str = "C", orientation_b: str = "C",
) -> SplitPair:
    """Two NanoBiT constructs. Binder A gets LgBiT, binder B gets SmBiT.

    When the analyte brings binder A and binder B into proximity (a sandwich on
    one target, or a ligand-induced dimer such as FKBP·rapamycin·FRB), LgBiT and
    SmBiT reconstitute an active luciferase -> luminescence. The fragments are
    engineered to associate only weakly on their own (KD ~190 uM), so signal is
    dictated by the binders' interaction, not the tags.
    """
    a = SplitConstruct(
        name=f"{binder_a_name}-LgBiT",
        sequence=_fuse(lgbit_seq, binder_a_seq, linker, orientation_a),
        fusion_partner=binder_a_name, fragment="LgBiT",
        orientation=orientation_a, length=0)
    a.length = len(a.sequence)
    b = SplitConstruct(
        name=f"{binder_b_name}-SmBiT",
        sequence=_fuse(smbit_seq, binder_b_seq, linker, orientation_b),
        fusion_partner=binder_b_name, fragment="SmBiT",
        orientation=orientation_b, length=0)
    b.length = len(b.sequence)
    return SplitPair(analyte=analyte, linker=linker, lgbit_construct=a, smbit_construct=b,
                     note="⚠️ construction is exact; ligand-induced complementation "
                          "(dynamic range) is a wet-lab luminescence readout.")


def verify_split(pair: SplitPair, binder_a_seq, binder_b_seq, lgbit_seq, smbit_seq) -> dict:
    a, b = pair.lgbit_construct, pair.smbit_construct
    return {
        "lgbit_has_binder": binder_a_seq in a.sequence,
        "lgbit_has_fragment": lgbit_seq in a.sequence,
        "smbit_has_binder": binder_b_seq in b.sequence,
        "smbit_has_fragment": smbit_seq in b.sequence,
        "all_ok": all([binder_a_seq in a.sequence, lgbit_seq in a.sequence,
                       binder_b_seq in b.sequence, smbit_seq in b.sequence]),
    }


# ===========================================================================
# #5  Multi-input INTRAMOLECULAR logic gate (paper's AND-gate architecture)
# ===========================================================================
@dataclass
class LogicGate:
    name: str
    reporter: str
    sequence: str
    length: int
    inputs: list                 # [(module_name, ligand, insertion_label), ...]
    gate: str                    # "AND" | "YES"
    meta: dict = field(default_factory=dict)


def build_multi_insertion(reporter_seq: str, inserts, flank: str = "G") -> str:
    """Insert several cp-modules into the reporter at given 0-idx cut points.

    inserts: list of (cp_seq, index0). Processed high->low index so earlier cut
    points are not shifted. Each insert is flanked by `flank` on both sides.
    """
    seq = reporter_seq
    for cp_seq, idx in sorted(inserts, key=lambda t: t[1], reverse=True):
        seq = seq[:idx] + flank + cp_seq + flank + seq[idx:]
    return seq


def build_and_gate(
    *, name: str, reporter_name: str, reporter_seq: str,
    modules,                     # [(module_name, ligand, module_seq, perm_site, insertion_index), ...]
    gs_linker: str = "GGSGGSGGS", flank: str = "G",
) -> LogicGate:
    """Two (or more) orthogonal ligand-gated modules, each circularly permuted and
    inserted at a distinct reporter loop -> an intramolecular AND gate: the
    reporter is ON only when ALL inputs are liganded (paper: loops 41 & 197).
    """
    inserts, inputs = [], []
    for mname, ligand, mseq, psite, idx in modules:
        cp = circular_permute(mseq, psite, gs_linker)
        inserts.append((cp, idx))
        inputs.append((mname, ligand, idx))
    seq = build_multi_insertion(reporter_seq, inserts, flank)
    return LogicGate(
        name=name, reporter=reporter_name, sequence=seq, length=len(seq),
        inputs=[(m, l) for m, l, _ in inputs], gate="AND" if len(modules) > 1 else "YES",
        meta={"n_inputs": len(modules), "insertion_indices": [i for _, _, i in inputs]},
    )


def truth_table(gate: LogicGate) -> list:
    """Design-intent truth table (⚠️ hypothesis): an AND gate is ON only when
    every input's ligand is present. Construction is ✅; the actual ligand-
    dependent activity must be measured."""
    import itertools
    ligs = [l for _, l in gate.inputs]
    rows = []
    for combo in itertools.product([0, 1], repeat=len(ligs)):
        on = all(combo) if gate.gate == "AND" else any(combo)
        rows.append({**{ligs[i]: combo[i] for i in range(len(ligs))}, "reporter_ON": int(on)})
    return rows
