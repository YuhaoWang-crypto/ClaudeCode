"""
tune_linker.py -- vary the receptor<->reporter linker of one chimera.

The paper tunes "the linker connecting the last alpha-helix of the receptor and
beta-strand 8 of beta-lactamase": stepwise glycine lengthening lowers dynamic
range but raises kcat, and serine rigidification reverts it. This module builds
that flank-linker series for a chosen system's primary construct so each can be
Boltz-validated.

    python3 -m biosensor_pipeline.tune_linker vitd

Emits sequences + a manifest to biosensor_out/<key>_linker_scan.json.
"""

from __future__ import annotations
import json
import os
import sys

from .systems import get_system
from .screen import build_library
from .design import build_chimera, verify_chimera

OUT = os.path.join(os.getcwd(), "biosensor_out")

# flank linker series: 1 Gly (baseline) -> longer flexible Gly -> Ser-rigidified
FLANK_SERIES = ["G", "GGG", "GGGGG", "GSGSG"]


def scan(key: str, flanks=FLANK_SERIES, gs_linker: str = "GGSGGSGGS") -> dict:
    sysm = get_system(key)
    rc, rp = sysm.receptor, sysm.reporter
    q = rp.insertion_sites[sysm.primary_insertion]
    # use the same permutation site as the registered primary (middle loop)
    site = rc.loop_sites[len(rc.loop_sites) // 2]

    variants = []
    for fl in flanks:
        ch = build_chimera(
            name=f"cp{rc.name}-{site}_{rp.name}-{sysm.primary_insertion}_L{len(fl)}{'S' if 'S' in fl else 'G'}",
            reporter_name=rp.name, reporter_seq=rp.seq,
            receptor_name=rc.name, receptor_seq=rc.seq,
            insertion_index=q, permutation_site=site,
            gs_linker=gs_linker, flank_linker=fl,
        )
        chk = verify_chimera(ch, rp.seq, rc.seq)
        variants.append({
            "name": ch.name, "flank_linker": fl, "flank_len": len(fl),
            "rigidified": "S" in fl, "length": ch.length,
            "valid": chk["all_ok"], "sequence": ch.sequence,
        })

    out = {"system": key, "analyte": rc.ligand_name, "permutation_site": site,
           "ligand_smiles": rc.ligand_smiles, "variants": variants}
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, f"{key}_linker_scan.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    return out


def main(argv):
    key = argv[1] if len(argv) > 1 else "vitd"
    out = scan(key)
    print(f"linker scan for '{key}' (analyte={out['analyte']}, site={out['permutation_site']})")
    for v in out["variants"]:
        print(f"  {v['name']:44s} flank={v['flank_linker']:6s} len={v['length']} "
              f"rigid={v['rigidified']} valid={v['valid']}")
    print(f"\nwrote {OUT}/{key}_linker_scan.json  (submit each 'sequence' as a Boltz holo)")


if __name__ == "__main__":
    main(sys.argv)
