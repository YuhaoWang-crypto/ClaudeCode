"""
demo_nanobit.py -- architecture #4: split-fragment complementation (NanoBiT).

Builds the two-component, ligand-induced-proximity sensor: binder A fused to
LgBiT and binder B fused to SmBiT. When the analyte brings A and B together,
LgBiT+SmBiT reconstitute an active luciferase -> luminescence.

Demo analyte: rapamycin (FKBP12 + FRB), the canonical NanoBiT positive control —
rapamycin bridges FKBP12·FRB, forcing complementation. This is the two-component
counterpart to the paper's single-component switches / auxiliary-binding designs.

    python3 -m biosensor_pipeline.demo_nanobit
"""

from __future__ import annotations
import json
import os

from .systems import NANOBIT, PROXIMITY_PAIRS
from .architectures import build_split_complementation, verify_split

OUT = os.path.join(os.getcwd(), "biosensor_out")


def main():
    analyte = "rapamycin"
    pair = PROXIMITY_PAIRS[analyte]
    (a_name, a_seq) = pair["binder_a"]
    (b_name, b_seq) = pair["binder_b"]
    lg, sm = NANOBIT["LgBiT"], NANOBIT["SmBiT"]

    sp = build_split_complementation(
        analyte=analyte, binder_a_name=a_name, binder_a_seq=a_seq,
        binder_b_name=b_name, binder_b_seq=b_seq, lgbit_seq=lg, smbit_seq=sm)
    chk = verify_split(sp, a_seq, b_seq, lg, sm)

    print("=== NanoBiT split-complementation sensor (architecture #4) ===")
    print(f"  analyte: {analyte}  |  {pair['note']}")
    print(f"  split site: NanoLuc {NANOBIT['split_site']}  ({NANOBIT['ref']})")
    print(f"  construct A: {sp.lgbit_construct.name}  ({sp.lgbit_construct.length} aa, "
          f"{a_name}+LgBiT[{len(lg)}aa])")
    print(f"  construct B: {sp.smbit_construct.name}  ({sp.smbit_construct.length} aa, "
          f"{b_name}+SmBiT[{len(sm)}aa peptide])")
    print(f"  linker: {sp.linker}")
    print(f"  verify: {chk}")
    print(f"  mechanism: no rapamycin -> A,B apart -> LgBiT/SmBiT weak (KD~190uM) -> DARK;")
    print(f"             + rapamycin -> FKBP12·rapamycin·FRB -> fragments reconstituted -> LIGHT")
    print(f"  {sp.note}")

    os.makedirs(OUT, exist_ok=True)
    rec = {
        "architecture": "split_complementation (NanoBiT)", "analyte": analyte,
        "split_site": NANOBIT["split_site"], "linker": sp.linker, "verify": chk,
        "construct_A": {"name": sp.lgbit_construct.name, "len": sp.lgbit_construct.length,
                        "sequence": sp.lgbit_construct.sequence},
        "construct_B": {"name": sp.smbit_construct.name, "len": sp.smbit_construct.length,
                        "sequence": sp.smbit_construct.sequence},
    }
    with open(os.path.join(OUT, "demo_nanobit.json"), "w") as fh:
        json.dump(rec, fh, indent=2)
    print(f"\n  wrote {OUT}/demo_nanobit.json   all_ok={chk['all_ok']}")
    return chk["all_ok"]


if __name__ == "__main__":
    import sys
    sys.exit(0 if main() else 1)
