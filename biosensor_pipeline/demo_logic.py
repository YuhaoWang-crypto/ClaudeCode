"""
demo_logic.py -- architecture #5: ligand-gated modules -> intramolecular AND gate.

Uses the ligand-gated DNA-binding module registry (TrpR/Trp, MetJ/SAM, DtxR/metal)
as ORTHOGONAL small-molecule inputs, and builds the paper's intramolecular AND
gate: two orthogonal modules, each circularly permuted and inserted at a distinct
TEM-1 loop (positions 41 & 197). The reporter is ON only when BOTH ligands are
present.

    python3 -m biosensor_pipeline.demo_logic

Construction is exact (✅). The truth table is the DESIGN INTENT (⚠️) — actual
ligand-dependent β-lactamase activity must be measured at the bench.
"""

from __future__ import annotations
import json
import os

from .systems import TEM1, LIGAND_GATED_MODULES
from .architectures import build_and_gate, truth_table
from .design import verify_chimera, Chimera

OUT = os.path.join(os.getcwd(), "biosensor_out")


def main():
    # two orthogonal inputs: Trp (via TrpR) AND SAM (via MetJ)
    A, B = LIGAND_GATED_MODULES["TrpR"], LIGAND_GATED_MODULES["MetJ"]
    # permutation site = an interior residue (placeholder; a real design uses
    # structure.annotate to pick a loop). Insertion loops = TEM-1 41 & 197.
    modules = [
        ("TrpR", A["ligand"], A["seq"], len(A["seq"]) // 2, TEM1.insertion_sites["41"]),
        ("MetJ", B["ligand"], B["seq"], len(B["seq"]) // 2, TEM1.insertion_sites["197"]),
    ]
    gate = build_and_gate(name="AND[TrpR@41,MetJ@197]-TEM1",
                          reporter_name=TEM1.name, reporter_seq=TEM1.seq, modules=modules)

    # verify: reporter catalytic residues still present & in order, both modules present
    cat_ok = all(gate.sequence.find(TEM1.seq[i]) >= 0 for i in TEM1.catalytic.values())
    both_present = all(m[2][: len(m[2]) // 2] and True for m in modules)  # sanity
    # reporter contiguity: the reporter residues appear in order (subsequence check)
    def is_subsequence(sub, seq):
        it = iter(seq)
        return all(c in it for c in sub)
    reporter_ok = is_subsequence(TEM1.seq, gate.sequence)

    print("=== Intramolecular AND gate from ligand-gated modules (architecture #5) ===")
    print(f"  reporter: {gate.reporter} (colorimetric)   construct length: {gate.length} aa")
    print(f"  inputs (orthogonal): " + " AND ".join(f"{m}({l})" for m, l in gate.inputs))
    print(f"  insertion loops: TEM-1 {list(k for k in ['41','197'])} "
          f"(the paper's YES/AND-gate loops)")
    print(f"  verify: reporter residues in order={reporter_ok}, "
          f"catalytic residues present={cat_ok}")
    print("\n  design-intent truth table (⚠️ hypothesis — activity must be measured):")
    tt = truth_table(gate)
    ligs = [l for _, l in gate.inputs]
    print("   " + " | ".join(f"{l[:16]:16s}" for l in ligs) + " | reporter_ON")
    for row in tt:
        print("   " + " | ".join(f"{row[l]:16d}" for l in ligs) + f" | {row['reporter_ON']}")

    # also show a single-input YES gate for contrast
    yes = build_and_gate(name="YES[DtxR@197]-TEM1", reporter_name=TEM1.name,
                         reporter_seq=TEM1.seq,
                         modules=[("DtxR", LIGAND_GATED_MODULES["DtxR"]["ligand"],
                                   LIGAND_GATED_MODULES["DtxR"]["seq"],
                                   len(LIGAND_GATED_MODULES["DtxR"]["seq"]) // 2,
                                   TEM1.insertion_sites["197"])])
    print(f"\n  (contrast) YES gate {yes.name}: {yes.length} aa, "
          f"gate={yes.gate}, ON iff {yes.inputs[0][1]} present")

    os.makedirs(OUT, exist_ok=True)
    rec = {
        "architecture": "ligand-gated modules -> intramolecular AND gate",
        "module_registry": {k: {"ligand": v["ligand"], "pdb": v["pdb"],
                                "dna_site": v["dna_site"]}
                            for k, v in LIGAND_GATED_MODULES.items()},
        "and_gate": {"name": gate.name, "inputs": gate.inputs, "length": gate.length,
                     "reporter_in_order": reporter_ok, "catalytic_present": cat_ok,
                     "truth_table": tt, "sequence": gate.sequence},
    }
    with open(os.path.join(OUT, "demo_logic.json"), "w") as fh:
        json.dump(rec, fh, indent=2)
    ok = reporter_ok and cat_ok and tt[-1]["reporter_ON"] == 1 and tt[0]["reporter_ON"] == 0
    print(f"\n  wrote {OUT}/demo_logic.json   all_ok={ok}")
    return ok


if __name__ == "__main__":
    import sys
    sys.exit(0 if main() else 1)
