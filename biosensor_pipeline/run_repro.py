"""
run_repro.py -- deterministic reproduction driver.

Runs the *offline, 100%-reproducible* half of the workflow for every system:
build the focused circular-permutation library, verify each construction, and
emit the sequences + Boltz payloads that the in-silico validation half will
submit.

    python -m biosensor_pipeline.run_repro            # both systems
    python -m biosensor_pipeline.run_repro dig        # one system

Outputs land in ./biosensor_out/.
"""

from __future__ import annotations
import json
import os
import sys

from .systems import SYSTEMS, get_system
from .screen import build_library, library_summary
from . import boltz_io

OUT = os.path.join(os.getcwd(), "biosensor_out")


def run_system(key: str) -> dict:
    sysm = get_system(key)
    lib = build_library(sysm)
    summ = library_summary(sysm, lib)

    os.makedirs(OUT, exist_ok=True)
    # write FASTA of the whole library
    with open(os.path.join(OUT, f"{key}_library.fasta"), "w") as fh:
        for c in lib:
            fh.write(c.fasta())

    # pick the "primary" construct = middle loop site (a neutral default; the
    # in-silico proxy re-ranks the library once Boltz metrics are in).
    primary = lib[len(lib) // 2]

    payloads = {
        "holo": boltz_io.holo_payload(primary.sequence, sysm.receptor.ligand_smiles),
        "apo": boltz_io.apo_payload(primary.sequence),
        "control": boltz_io.control_payload(sysm.receptor.seq, sysm.receptor.ligand_smiles),
    }
    with open(os.path.join(OUT, f"{key}_boltz_payloads.json"), "w") as fh:
        json.dump({"primary_construct": primary.name, "payloads": payloads}, fh, indent=2)

    summ["primary_construct"] = primary.name
    summ["primary_length"] = primary.length
    with open(os.path.join(OUT, f"{key}_summary.json"), "w") as fh:
        json.dump(summ, fh, indent=2)
    return summ


def main(argv):
    keys = argv[1:] if len(argv) > 1 else list(SYSTEMS)
    allsumm = {}
    for k in keys:
        s = run_system(k)
        allsumm[k] = s
        print(f"\n=== system '{k}' ({s['role']}) ===")
        print(f"  receptor={s['receptor']}  reporter={s['reporter']}  analyte={s['analyte']}")
        print(f"  insertion={s['insertion']}  n_variants={s['n_variants']}  "
              f"all_valid={s['all_constructions_valid']}")
        print(f"  primary construct: {s['primary_construct']} ({s['primary_length']} aa)")
        for v in s["variants"]:
            print(f"    - {v['name']}: site={v['permutation_site']} "
                  f"(removed {v['removed_residue']}), len={v['length']}, valid={v['valid']}")
    print(f"\nwrote outputs to {OUT}/")
    return allsumm


if __name__ == "__main__":
    main(sys.argv)
