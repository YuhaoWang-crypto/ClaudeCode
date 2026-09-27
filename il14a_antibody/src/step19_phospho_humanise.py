"""Step 19 - Humanise the phospho arm and merge it into the deliverable.

Step 16's designs are made on the same murine SM3 scaffold as step 15's, so
they need the same treatment: the designed loops onto a human germline
acceptor that actually splits and carries exactly the two canonical
cysteines. Only the phosphate-visible condition is carried forward - the
context-off and no-phosphate runs are controls and are not candidates.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import RESULTS
from humanness import germline_9mer_coverage, load_germline_kmers
from step17_humanise_designs import (graft, liabilities, pick_acceptor,
                                     split_heavy, split_light)

TEMPLATE_CHAINS = ["H", "L"]


def main():
    KM, _ = load_germline_kmers()
    arms = json.loads((RESULTS / "step16_phospho.json").read_text())
    out, fasta = [], []
    for arm in arms:
        W = arm["window"]
        # the first stored record is the native input sequence, which
        # LigandMPNN echoes and which carries no id= field
        designs = [d for d in arm["designs"]["pS_context_ON"]
                   if "id=" in d["header"]]
        print("=" * 96)
        print(f"{W['id']}  TXLNA {W['start']}-{W['start']+len(W['window'])-1}"
              f"  {W['window']}  pS515 at position {W['idx']}")
        print(f"  pocket {', '.join(arm['chosen_rotamer']['near'])}")
        print(f"  context control: {arm['context_control']['verdict']}")
        print("=" * 96)
        seen, rank = set(), 0
        for d in designs:
            vh_m, vl_m = d["seqs"][0], d["seqs"][1]
            ph, pl = split_heavy(vh_m), split_light(vl_m)
            if ph is None or pl is None:
                continue
            if any("C" in ph[k] or "C" in pl[k]
                   for k in ("CDR1", "CDR2", "CDR3")):
                continue
            key = (ph["CDR1"], ph["CDR2"], ph["CDR3"],
                   pl["CDR1"], pl["CDR2"], pl["CDR3"])
            if key in seen:
                continue
            seen.add(key)
            gh = pick_acceptor(vh_m[:130], "IGHV", "H")
            gk = pick_acceptor(vl_m[:130], "IGKV", "L")
            gl = pick_acceptor(vl_m[:130], "IGLV", "L")
            cands = [c for c in (gk, gl) if c]
            if gh is None or not cands:
                continue
            gsel = max(cands, key=lambda c: c["identity"])
            lkind = "K" if gsel is gk else "L"
            vh_h, _ = graft(ph, gh["sequence"], "H")
            vl_h, _ = graft(pl, gsel["sequence"], lkind)
            if vh_h is None or vl_h is None:
                continue
            lia = {**{f"VH:{k}": v for k, v in liabilities(vh_h).items()},
                   **{f"VL:{k}": v for k, v in liabilities(vl_h).items()}}
            ch, _ = germline_9mer_coverage(vh_h, KM)
            cl, _ = germline_9mer_coverage(vl_h, KM)
            rank += 1
            print(f"  design {rank}  CDR-H3 {ph['CDR3']:14} "
                  f"CDR-L3 {pl['CDR3']:12} acceptor {gh['gene']}/"
                  f"{gsel['gene']}  9mer {ch:.3f}/{cl:.3f}  "
                  f"liabilities {lia or 'none'}")
            name = f"{W['id']}_d{rank}"
            fasta += [
                f">{name}_VH | phospho arm | acceptor {gh['gene']} | epitope "
                f"TXLNA {W['start']}-{W['start']+len(W['window'])-1} "
                f"{W['window']} with pSer at position {W['idx']} | "
                f"CDRH3 {ph['CDR3']} | 9mer_cov {ch:.3f}", vh_h,
                f">{name}_VL | acceptor {gsel['gene']} | CDRL3 {pl['CDR3']} | "
                f"9mer_cov {cl:.3f}", vl_h]
            out.append({"id": name, "arm": W["id"], "epitope": W["window"],
                        "txlna_window": [W["start"],
                                         W["start"] + len(W["window"]) - 1],
                        "phospho_position_in_window": W["idx"],
                        "pocket": arm["chosen_rotamer"]["near"],
                        "cdrs": {"H": {k: ph[k] for k in
                                       ("CDR1", "CDR2", "CDR3")},
                                 "L": {k: pl[k] for k in
                                       ("CDR1", "CDR2", "CDR3")}},
                        "acceptor": {"VH": gh["gene"], "VL": gsel["gene"]},
                        "humanised": {"VH": vh_h, "VL": vl_h},
                        "germline_9mer_coverage": {"VH": round(ch, 4),
                                                   "VL": round(cl, 4)},
                        "liabilities": lia})
            if rank >= 3:
                break
        print()

    (RESULTS / "step19_phospho_humanised.json").write_text(json.dumps(out, indent=2))

    # merge into one deliverable
    pan = (RESULTS / "IL14A_ctd_designed_humanised.fasta").read_text().rstrip()
    merged = pan + "\n#\n# --- phospho arm (LigandMPNN, phosphate as ligand context) ---\n" \
        + "\n".join(fasta) + "\n"
    (RESULTS / "IL14A_ctd_designed_humanised.fasta").write_text(merged)
    n = merged.count(">")
    print(f"wrote {len(out)} phospho designs; deliverable now holds "
          f"{n//2} humanised VH/VL pairs")


if __name__ == "__main__":
    main()
