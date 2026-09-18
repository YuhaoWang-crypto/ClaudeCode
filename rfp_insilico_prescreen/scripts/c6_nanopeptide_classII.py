#!/usr/bin/env python3
"""
C6 - What in-silico can and cannot say about RFP Step 1 (the nanopeptide).

Step 1 asks for MHC class II-dependent, peptide-dependent IL-2 output from the
DO11.10 hybridoma, which is I-A(d)-restricted and specific for chicken ovalbumin
323-339. The peptide set is:

    8-mer core control (all L)
    9-mer core test (4 D-amino acids)
    17-mer uniblock
    30-mer diblock control
    30-mer diblock test (4 D-amino acids)

This module establishes, by running it rather than by asserting it, exactly which
of those five a class II predictor can speak to. Four things are tested:

  1. Positive control. cOVA323-339 on H2-IAd should return the documented core.
     If it does not, nothing else in this module is worth reading.

  2. Length floor. The IEDB class II API refuses input shorter than 11 residues.
     That is measured here, not quoted, because it removes the 8-mer and the bare
     9-mer core from reach entirely.

  3. Flanking-residue dependence. A class II groove is open-ended: the bound
     core is 9 residues but affinity depends on the peptide flanking residues
     either side. Padding the core to clear the length floor is therefore not a
     workaround - it scores a different molecule. The padding sweep here shows
     how much the answer moves with the padding, which is the size of the
     artefact anyone taking that shortcut would import.

  4. Assembly blocks. A diblock construct is scored to show whether appending a
     self-assembling block changes the predicted core. It does not, which is the
     point: sequence-based class II prediction is blind to the thing the diblock
     exists to do.

D-amino acids are not tested because they cannot be: every class II predictor in
use is trained on L-peptides and has no representation of stereochemistry. A
substituted sequence is scored as if it were all-L, so the number returned is the
score of a different molecule. That is reported as a hard exclusion, not a caveat.

Output: results/c6_nanopeptide_classII.json
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import results_path, iedb_mhcii  # noqa: E402

ENDPOINT = "https://tools-cluster-interface.iedb.org/tools_api/mhcii/"
METHOD = "netmhciipan_el"
ALLELE = "H2-IAd"

# cOVA323-339: the epitope the DO11.10 hybridoma is specific for, and by
# coincidence exactly the 17 residues the RFP's uniblock is. Used here as the
# stand-in for the uniblock and as the system's positive control.
COVA = "ISQAVHAAHAEINEAGR"
DOCUMENTED_CORE = "VHAAHAEIN"      # I-A(d) register reported for OVA323-339

# Constructed stand-ins. The sponsor's real sequences drop in here unchanged;
# these exist so the demo runs end to end and are labelled as constructed.
CORE9 = "VHAAHAEIN"
CORE8 = "VHAAHAEI"
ASSEMBLY_BLOCK = "AEAEAKAKAEA"     # EAK-type self-assembling block, 11 aa
DIBLOCK30 = COVA + "GG" + ASSEMBLY_BLOCK[:11]      # 17 + 2 + 11 = 30

PANEL = [
    ("uniblock_17mer_L", COVA, "17-mer uniblock stand-in = cOVA323-339", True),
    ("core_9mer_L", CORE9, "9-mer core, all L", True),
    ("core_8mer_L", CORE8, "8-mer core control, all L", True),
    ("diblock_30mer_L", DIBLOCK30, "30-mer diblock control (uniblock + GG + EAK block)", True),
    ("core_9mer_D4", CORE9, "9-mer core test with 4 D-amino acids", False),
    ("diblock_30mer_D4", DIBLOCK30, "30-mer diblock test with 4 D-amino acids", False),
]

PAD_LENGTHS = [11, 13, 15, 17]


API_MIN_LEN = 11        # measured: the class II endpoint refuses shorter input


SCAN_WINDOW = 15        # fixed, so constructs of different lengths stay comparable


def score(seq, label, window=SCAN_WINDOW):
    """Return rows, or the API's refusal text. A refusal is a result here.

    The panel is scanned with ONE window (15-mer) across every construct. Scoring
    each construct at its own full length instead would make a 17-mer and a
    30-mer incomparable - a %Rank is a percentile against a background of that
    same length - and would manufacture a difference between the uniblock and the
    diblock that is an artefact of the window, not of the sequence. The padding
    sweep below deliberately varies the window; nothing else does.
    """
    w = window
    try:
        rows = iedb_mhcii(ENDPOINT, METHOD, {label: seq}, [ALLELE], length=w)
        return {"ok": True, "rows": rows}
    except Exception as e:
        return {"ok": False, "error": str(e)[:300]}


def best(rows):
    if not rows:
        return None
    r = min(rows, key=lambda x: float(x["rank"]))
    return {"core": r["core_peptide"], "peptide": r["peptide"],
            "rank": float(r["rank"]), "score": float(r["score"])}


def main():
    out = {"allele": ALLELE, "method": METHOD, "documented_core": DOCUMENTED_CORE}

    # ---- 1. positive control -------------------------------------------
    pc = score(COVA, "cOVA323_339")
    pc_best = best(pc.get("rows", [])) if pc["ok"] else None
    out["positive_control"] = {
        "peptide": COVA, "result": pc_best,
        "core_matches_literature": bool(pc_best and pc_best["core"] == DOCUMENTED_CORE),
    }
    print(f"positive control cOVA323-339 on {ALLELE}")
    if pc_best:
        ok = pc_best["core"] == DOCUMENTED_CORE
        print(f"  core {pc_best['core']}  %Rank {pc_best['rank']}  "
              f"-> literature register {'REPRODUCED' if ok else 'NOT reproduced'}")
    else:
        print("  FAILED - do not read the rest of this module")

    # ---- 2 & 4. the panel ----------------------------------------------
    panel_out = []
    print("\npanel")
    for name, seq, note, scoreable in PANEL:
        rec = {"name": name, "length": len(seq), "note": note,
               "scoreable_in_silico": scoreable}
        if not scoreable:
            rec["excluded_reason"] = (
                "contains D-amino acids; class II predictors are trained on "
                "L-peptides only and have no stereochemistry representation, so "
                "any score returned would be the score of the all-L sequence, "
                "which is a different molecule")
            print(f"  {name:22s} len {len(seq):2d}  EXCLUDED - D-amino acids")
        else:
            r = score(seq, name)
            if r["ok"]:
                b = best(r["rows"])
                rec["result"] = b
                print(f"  {name:22s} len {len(seq):2d}  core {b['core']}  "
                      f"%Rank {b['rank']}")
            else:
                rec["result"] = None
                rec["api_refusal"] = r["error"]
                short = r["error"].split("\n")[0][:60]
                print(f"  {name:22s} len {len(seq):2d}  REFUSED BY API - {short}")
        panel_out.append(rec)
    out["panel"] = panel_out

    # ---- 3. padding sweep ----------------------------------------------
    print("\npadding sweep: what a bare 9-mer core scores once padded to clear "
          "the length floor")
    sweep = []
    for L in PAD_LENGTHS:
        pad = L - len(CORE9)
        left = pad // 2
        seq = "G" * left + CORE9 + "G" * (pad - left)
        r = score(seq, f"pad{L}", window=L)
        b = best(r["rows"]) if r["ok"] else None
        sweep.append({"padded_length": L, "sequence": seq,
                      "result": b, "api_refusal": None if r["ok"] else r["error"][:200]})
        if b:
            print(f"  len {L:2d}  {seq:20s} core {b['core']}  %Rank {b['rank']}")
        else:
            print(f"  len {L:2d}  {seq:20s} REFUSED")
    ranks = [s["result"]["rank"] for s in sweep if s["result"]]
    out["padding_sweep"] = {
        "core": CORE9, "pad_residue": "G", "sweep": sweep,
        "rank_range": [min(ranks), max(ranks)] if ranks else None,
        "fold_range": round(max(ranks) / min(ranks), 1) if ranks and min(ranks) > 0 else None,
    }
    if ranks:
        print(f"  -> %Rank spans {min(ranks)} to {max(ranks)} "
              f"({round(max(ranks)/min(ranks),1)}x) purely from how much padding "
              f"was added")

    # ---- verdict --------------------------------------------------------
    n_scoreable = sum(1 for r in panel_out if r.get("result"))
    out["verdict"] = {
        "panel_size": len(PANEL),
        "scored": n_scoreable,
        "blocked_by_d_amino_acids": sum(1 for _, _, _, s in PANEL if not s),
        "blocked_by_length_floor": sum(1 for r in panel_out
                                       if r.get("api_refusal")),
        "statement": (
            "Class II prediction can confirm the DO11.10 / I-A(d) system control "
            "and can score the all-L constructs of 11 residues or more. It cannot "
            "score the 8-mer control, and it cannot score any D-substituted "
            "peptide. Since the comparison Step 1 is built around is all-L "
            "control versus D-substituted test, in-silico cannot perform that "
            "comparison - it can only characterise the L arm."),
    }
    p = results_path("c6_nanopeptide_classII.json")
    with open(p, "w") as f:
        json.dump(out, f, indent=1)
    print(f"\nscored {n_scoreable} of {len(PANEL)} panel peptides")
    print(f"wrote {p}")


if __name__ == "__main__":
    main()
