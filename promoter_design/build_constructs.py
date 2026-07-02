"""
Deterministic assembler: response-element modules -> full synthetic promoters.

Runs WITHOUT a GPU or any model. It multimerises each stimulus module, appends
a minimal promoter, and pads with neutral flanks so the cassette can later be
dropped into an AlphaGenome/Enformer scoring window. Output is FASTA that you
can order/clone directly, and that the design pipeline (design_pipeline.py)
reads back in for Evo2 / AlphaGenome scoring.

Usage:
    python build_constructs.py                 # default: 3 copies, minCMV
    python build_constructs.py --copies 5 --min_promoter E1b_TATA
"""

import argparse
import re
from elements import ELEMENTS, MINIMAL_PROMOTERS, CELL_CONTEXTS

# Neutral, low-complexity-avoiding flank (GC ~50%, no obvious TFBS). Kept short.
FLANK5 = "gctagcAGGACCGGT"          # includes NheI + AgeI cloning handles
FLANK3 = "GCTAGCggtaccGAattc"       # KpnI + EcoRI handles
KOZAK_STUB = "gccaccATGg"           # optional TSS/Kozak stub after core promoter


def clean(seq: str) -> str:
    return re.sub(r"[^ACGTacgtN]", "", seq)


def multimerise(seed: str, spacer: str, copies: int) -> str:
    seed, spacer = clean(seed), clean(spacer)
    return spacer.join([seed] * copies)


def build(stim_key: str, copies: int, min_promoter: str) -> dict:
    e = ELEMENTS[stim_key]
    enhancer = multimerise(e["seed"], e["spacer"], copies)
    core = clean(MINIMAL_PROMOTERS[min_promoter])
    cassette = clean(FLANK5) + enhancer + core + clean(KOZAK_STUB) + clean(FLANK3)
    return {
        "key": stim_key,
        "name": f"{stim_key}_{copies}x_{min_promoter}",
        "element": e["element"],
        "tf": e["tf"],
        "stimulus": e["stimulus"],
        "enhancer_len": len(enhancer),
        "seq": cassette.upper(),
    }


def which_cells(stim_key: str):
    return [c for c, v in CELL_CONTEXTS.items() if stim_key in v["good_for"]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--copies", type=int, default=3)
    ap.add_argument("--min_promoter", default="minCMV",
                    choices=list(MINIMAL_PROMOTERS))
    ap.add_argument("--out", default="designs/stimulus_responsive_promoters.fasta")
    args = ap.parse_args()

    records = [build(k, args.copies, args.min_promoter) for k in ELEMENTS]

    with open(args.out, "w") as fh:
        for r in records:
            cells = ",".join(which_cells(r["key"])) or "generic"
            header = (f">{r['name']} | {r['element']} | TF={r['tf']} | "
                      f"stimulus={r['stimulus']} | enhancer={r['enhancer_len']}bp | "
                      f"total={len(r['seq'])}bp | recommended_cells={cells}")
            fh.write(header + "\n")
            for i in range(0, len(r["seq"]), 70):
                fh.write(r["seq"][i:i + 70] + "\n")

    print(f"Wrote {len(records)} cassettes to {args.out} "
          f"({args.copies}x, {args.min_promoter})")
    for r in records:
        print(f"  {r['name']:<34} {len(r['seq'])}bp  cells=[{','.join(which_cells(r['key'])) or 'generic'}]")


if __name__ == "__main__":
    main()
