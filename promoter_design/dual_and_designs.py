"""
Composite designs: (1) ISRE+GAS pan-interferon sensor, and (2) stimulus-AND-
cell-type gated promoters. Deterministic assembly (no GPU); writes FASTA that
design_pipeline.py then refines/validates with Evo2 + AlphaGenome.

Design logic
------------
pan-interferon sensor:
    [insulator] - (ISRE - GAS)xN - [minimal promoter] - [Kozak stub] - [3' handle]
    Fires on type I IFN (ISGF3 -> ISRE) OR type II IFN (STAT1 -> GAS).

stimulus-AND-cell-type gate (analog AND):
    [insulator] - (stimulus RE)xN - (lineage RE)xM - [minimal promoter] - ...
    Strong output needs BOTH: the stimulus signal (stimulus RE occupied) AND the
    right cell (lineage TF present -> lineage RE occupied). This is analog/graded
    AND from enhancer synergy. For DIGITAL AND use a two-component relay (e.g. a
    stimulus-driven split-transactivator whose halves are each cell-restricted);
    noted in README. Cell selectivity is then verified in silico by AlphaGenome
    contrastive scoring (target ontology vs off-target ontologies).

Usage:
    python dual_and_designs.py
"""

import re
from elements import (ELEMENTS, LINEAGE_ELEMENTS, MINIMAL_PROMOTERS, COMPOSITE,
                      CELL_CONTEXTS)
from build_constructs import FLANK5, FLANK3, KOZAK_STUB, clean, multimerise


def build_composite(name: str, copies: int, min_promoter: str) -> dict:
    spec = COMPOSITE[name]
    mods = [ELEMENTS[k] for k in spec["modules"]]
    # Interleave one monomer of each module per repeat unit.
    unit = "".join(clean(m["seed"]) + clean(m["spacer"]) for m in mods)
    enhancer = "".join([unit] * copies)
    core = clean(MINIMAL_PROMOTERS[min_promoter])
    seq = (clean(FLANK5) + enhancer + core + clean(KOZAK_STUB) + clean(FLANK3)).upper()
    return {"name": f"{name}_{copies}x_{min_promoter}",
            "kind": "composite", "desc": spec["description"],
            "enh": len(enhancer), "seq": seq}


def build_and_gate(stim_key: str, lineage_key: str, stim_copies: int,
                   lineage_copies: int, min_promoter: str) -> dict:
    s, l = ELEMENTS[stim_key], LINEAGE_ELEMENTS[lineage_key]
    stim_block = multimerise(s["seed"], s["spacer"], stim_copies)
    lin_block = multimerise(l["seed"], l["spacer"], lineage_copies)
    core = clean(MINIMAL_PROMOTERS[min_promoter])
    seq = (clean(FLANK5) + stim_block + "gaattc" + lin_block + core
           + clean(KOZAK_STUB) + clean(FLANK3)).upper()
    return {"name": f"AND_{stim_key}_x{lineage_key}_{min_promoter}",
            "kind": "and_gate",
            "desc": (f"{s['element']} ({s['tf']}) AND lineage {l['tf']} "
                     f"-> ontology {l['ontology']}"),
            "enh": len(stim_block) + len(lin_block), "seq": seq}


def write_fasta(records, path):
    with open(path, "w") as fh:
        for r in records:
            fh.write(f">{r['name']} | {r['kind']} | {r['desc']} | "
                     f"enhancer={r['enh']}bp | total={len(r['seq'])}bp\n")
            for i in range(0, len(r["seq"]), 70):
                fh.write(r["seq"][i:i + 70] + "\n")


def main():
    records = []
    # 1) pan-interferon sensor (a few copy numbers, low-baseline E1b core)
    for n in (2, 3, 4):
        records.append(build_composite("pan_interferon", n, "E1b_TATA"))

    # 2) worked AND gates
    #    (a) interferon AND myeloid   -> IFN-responsive ONLY in monocyte/macrophage
    records.append(build_and_gate("interferon_typeII", "myeloid", 3, 3, "E1b_TATA"))
    records.append(build_and_gate("interferon_typeI",  "myeloid", 3, 3, "E1b_TATA"))
    #    (b) hypoxia AND hepatocyte   -> hypoxia-responsive ONLY in liver cells
    records.append(build_and_gate("hypoxia", "hepatocyte", 4, 2, "minCMV"))
    #    (c) oxidative stress AND hepatocyte
    records.append(build_and_gate("oxidative_stress", "hepatocyte_hnf1", 3, 2, "minCMV"))
    #    (d) cAMP AND neuronal
    records.append(build_and_gate("camp", "neuronal", 4, 3, "E1b_TATA"))

    out = "designs/composite_and_designs.fasta"
    write_fasta(records, out)
    print(f"Wrote {len(records)} composite/AND designs to {out}\n")
    for r in records:
        print(f"  {r['name']:<42} {len(r['seq'])}bp  ({r['kind']})")


if __name__ == "__main__":
    main()
