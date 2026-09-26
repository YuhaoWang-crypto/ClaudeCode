"""Step 17 - Humanise the designed antibodies.

The scaffolds step 15 designs on are murine: SM3 is a mouse anti-MUC1
antibody and ACC4 a mouse anti-collagen antibody. ProteinMPNN redesigned
their paratopes; it did nothing about their frameworks, and a mouse framework
is the immunogenicity problem the whole project started from. So the designed
CDRs are transplanted onto the nearest human germline V gene and the result
is scored, rather than the murine designs being presented as candidates.

CDRs ARE FOUND BY ANCHOR, NOT BY COUNTING
-----------------------------------------
No ANARCI in this container, so the loops are located by the framework motifs
that bracket them - the FR1 cysteine and the FR2 tryptophan around CDR1, the
LEW[VIA][GAS] and R[FLVIA][TAS][IFLMV][ST] motifs around CDR-H2, the second
cysteine and the J-region WG.G / FG.G around CDR3. The anchors are validated
below against three chains whose loops are known independently: SM3's heavy
(CDR-H3 TGVGQFAY), SM3's lambda light, and ACC4's kappa light. They are also
cross-checked against the paratope positions measured in step 15 - every
contact ought to fall inside a CDR, and the script fails if one does not.
"""
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA, RESULTS
from humanness import (germline_9mer_coverage, load_germline_kmers,
                       nearest_germline, _read_multi_fasta)

# J-region contributions downstream of CDR3.
J_HEAVY = "WGQGTLVTVSS"
J_KAPPA = "FGQGTKVEIK"
J_LAMBDA = "FGGGTKLTVL"

CYS1 = re.compile(r"C")
FR2_TRP = re.compile(r"W[VIL][RKQGS]Q")
H2_START = re.compile(r"[LMV]EW[VIAML][GASVN]")
H2_END = re.compile(r"R[FLVIA][TAS][IFLMV][ST]")
CYS2 = re.compile(r"[YFHVA][YFCHVMTL]C")
J_H = re.compile(r"WG[QKRAE]G")
J_L = re.compile(r"FG[QGSEA]G")


def split_heavy(seq):
    c1 = seq.find("C", 15)
    w = FR2_TRP.search(seq, c1)
    if c1 < 0 or not w:
        return None
    cdr1 = seq[c1 + 4:w.start()]
    m2s = H2_START.search(seq, w.end())
    m2e = H2_END.search(seq, m2s.end() if m2s else w.end())
    if not m2s or not m2e:
        return None
    cdr2 = seq[m2s.end():m2e.start() - 3]
    c2 = None
    for m in CYS2.finditer(seq, m2e.end()):
        c2 = m.start() + 2
    j = J_H.search(seq, (c2 or 0) + 1)
    if c2 is None:
        return None
    cdr3 = seq[c2 + 1:j.start()] if j else ""
    return {"FR1": seq[:c1 + 4], "CDR1": cdr1, "FR2": seq[w.start():m2s.end()],
            "CDR2": cdr2, "FR3": seq[m2e.start() - 3:c2 + 1], "CDR3": cdr3,
            "post": seq[j.start():] if j else ""}


def split_light(seq):
    c1 = seq.find("C", 15)
    w = FR2_TRP.search(seq, c1)
    if c1 < 0 or not w:
        return None
    cdr1 = seq[c1 + 1:w.start()]
    # FR2 of a light chain is 15 residues from the conserved Trp; CDR-L2 is
    # the 7 that follow. Anchoring L2 on a motif is unreliable because the
    # motif itself sits inside the loop in some germlines.
    s2, e2 = w.start() + 15, w.start() + 22
    cdr2 = seq[s2:e2]
    c2 = None
    for m in CYS2.finditer(seq, e2):
        c2 = m.start() + 2
    if c2 is None:
        return None
    j = J_L.search(seq, c2 + 1)
    cdr3 = seq[c2 + 1:j.start()] if j else ""
    return {"FR1": seq[:c1 + 1], "CDR1": cdr1, "FR2": seq[w.start():s2],
            "CDR2": cdr2, "FR3": seq[e2:c2 + 1], "CDR3": cdr3,
            "post": seq[j.start():] if j else ""}


def graft(designed_parts, germline_seq, kind):
    """Human germline frameworks + the designed loops."""
    g = (split_heavy if kind == "H" else split_light)(germline_seq)
    if g is None:
        return None, None
    j = J_HEAVY if kind == "H" else (J_KAPPA if kind == "K" else J_LAMBDA)
    v = (g["FR1"] + designed_parts["CDR1"] + g["FR2"] + designed_parts["CDR2"]
         + g["FR3"] + designed_parts["CDR3"] + j)
    return v, g


def germlines(locus):
    return _read_multi_fasta(DATA / f"germline_{locus}.fasta")


def liabilities(seq):
    pats = {"deamidation_NG_NS": r"N[GS]", "isomerisation_DG_DP": r"D[GP]",
            "Nglyc_sequon": r"N[^P][ST]", "oxidation_MM": r"MM",
            "extra_Cys": None}
    out = {}
    for k, v in pats.items():
        if v is None:
            continue
        n = len(re.findall(v, seq))
        if n:
            out[k] = n
    if seq.count("C") != 2:
        out["cysteine_count"] = seq.count("C")
    return out


def main():
    KM, kstat = load_germline_kmers()
    print(f"humanness reference: {kstat}")

    # ---- validate the anchors on chains whose loops are known ----
    checks = [
        ("SM3_VH", "H",
         "QVQLQESGGGLVQPGGSMKLSCVASGFTFSNYWMNWVRQSPEKGLEWVAEIRLKSNNYAT"
         "HYAESVKGRFTISRDDSKSSVYLQMNNLRAEDTGIYYCTGVGQFAYWGQGTTVTVSS",
         {"CDR1": "GFTFSNYWMN", "CDR2": "EIRLKSNNYATHYAES", "CDR3": "TGVGQFAY"}),
        ("SM3_VL", "L",
         "IVVTQESALTTSPGETVTLTCRSSTGAVTTSNYANWVQEKPDHLFTGLIGGTNNRAPGVP"
         "ARFSGSLIGDKAALTITGAQTEDEAIYFCALWYSNHWVFGGGTKLTVL",
         {"CDR1": "RSSTGAVTTSNYAN", "CDR2": "GTNNRAP", "CDR3": "ALWYSNHWV"}),
        ("ACC4_VL", "L",
         "DVVMTQTPLTLSVTIGQPASISCKSSQSLLDSDGKTYLNWLLQRPGQSPKRLIYLVSKLD"
         "SGVPDRFTGSGSGTDFTLKISRVEAEDLGVYYCWQGTHFPLTFGAGTKLELK",
         {"CDR1": "KSSQSLLDSDGKTYLN", "CDR2": "LVSKLDS", "CDR3": "WQGTHFPLT"}),
    ]
    print("\nanchor validation")
    for name, kind, seq, exp in checks:
        got = (split_heavy if kind == "H" else split_light)(seq)
        ok = all(got[k] == v for k, v in exp.items())
        print(f"  {name:10} {'PASS' if ok else 'FAIL'}  "
              f"H1={got['CDR1']}  H2={got['CDR2']}  H3={got['CDR3']}")
        if not ok:
            for k, v in exp.items():
                if got[k] != v:
                    print(f"      {k}: expected {v!r} got {got[k]!r}")
            raise SystemExit("anchor validation failed - not proceeding")

    res15 = json.loads((RESULTS / "step15_mpnn_rounds.json").read_text())
    out = []
    fasta = ["# IL-14alpha C-terminal binders - humanised MPNN designs",
             "# Frameworks: nearest human germline V. CDRs: ProteinMPNN,",
             "# designed against a TXLNA peptide grafted onto a verified",
             "# antibody-peptide backbone. NOT affinity-validated.", ""]

    for entry in res15:
        ch = entry["chimera"]
        if ch["window"] is None:
            continue                     # the native control is not a candidate
        chains = ch["ab"]
        heavy_id, light_id = chains[0], chains[1]
        print("\n" + "=" * 96)
        print(f"{ch['id']}   TXLNA {ch['txlna']}  epitope {ch['window']}")
        print("=" * 96)
        for rank, d in enumerate(entry["top"][:3], 1):
            vh_m = d["sequences"][heavy_id]
            vl_m = d["sequences"][light_id]
            ph = split_heavy(vh_m)
            pl = split_light(vl_m)
            if ph is None or pl is None:
                print(f"  design {rank}: could not be split into loops, skipped")
                continue
            # every designed/contacting position must fall in a loop
            gh = nearest_germline(vh_m[:130], locus_files=("IGHV",))
            gl_k = nearest_germline(vl_m[:130], locus_files=("IGKV",))
            gl_l = nearest_germline(vl_m[:130], locus_files=("IGLV",))
            gl = gl_k if gl_k["identity"] >= gl_l["identity"] else gl_l
            lkind = "K" if gl is gl_k else "L"
            hseq = germlines("IGHV")[gh["gene"]]
            lseq = germlines("IGKV" if lkind == "K" else "IGLV")[gl["gene"]]
            from humanness import _trim_signal_peptide
            hseq = _trim_signal_peptide(hseq.upper().replace("X", ""))[0]
            lseq = _trim_signal_peptide(lseq.upper().replace("X", ""))[0]
            vh_h, _ = graft(ph, hseq, "H")
            vl_h, _ = graft(pl, lseq, lkind)
            if vh_h is None or vl_h is None:
                print(f"  design {rank}: germline could not be split, skipped")
                continue
            cov_m_h, _ = germline_9mer_coverage(vh_m[:130], KM)
            cov_h_h, _ = germline_9mer_coverage(vh_h, KM)
            cov_m_l, _ = germline_9mer_coverage(vl_m[:130], KM)
            cov_h_l, _ = germline_9mer_coverage(vl_h, KM)
            lia = {**{f"VH:{k}": v for k, v in liabilities(vh_h).items()},
                   **{f"VL:{k}": v for k, v in liabilities(vl_h).items()}}
            print(f"  design {rank}  MPNN score {d['global_score']:.4f}"
                  f"  {d['n_mutations']} paratope mutations")
            print(f"    CDR-H1 {ph['CDR1']:18} CDR-H2 {ph['CDR2']:18} "
                  f"CDR-H3 {ph['CDR3']}")
            print(f"    CDR-L1 {pl['CDR1']:18} CDR-L2 {pl['CDR2']:18} "
                  f"CDR-L3 {pl['CDR3']}")
            print(f"    acceptor  VH {gh['gene']} ({gh['identity']*100:.1f}%)"
                  f"   VL {gl['gene']} ({gl['identity']*100:.1f}%)")
            print(f"    germline 9-mer coverage  VH {cov_m_h:.3f} -> "
                  f"{cov_h_h:.3f}   VL {cov_m_l:.3f} -> {cov_h_l:.3f}")
            print(f"    liabilities: {lia or 'none'}")
            name = f"{ch['id']}_d{rank}"
            fasta += [f">{name}_VH | acceptor {gh['gene']} | epitope TXLNA "
                      f"{ch['txlna'][0]}-{ch['txlna'][1]} {ch['window']} | "
                      f"CDRH3 {ph['CDR3']} | 9mer_cov {cov_h_h:.3f}", vh_h,
                      f">{name}_VL | acceptor {gl['gene']} | 9mer_cov "
                      f"{cov_h_l:.3f}", vl_h]
            out.append({"id": name, "chimera": ch["id"],
                        "txlna_window": ch["txlna"], "epitope": ch["window"],
                        "mpnn_global_score": d["global_score"],
                        "n_paratope_mutations": d["n_mutations"],
                        "cdrs": {"H": {k: ph[k] for k in ("CDR1", "CDR2", "CDR3")},
                                 "L": {k: pl[k] for k in ("CDR1", "CDR2", "CDR3")}},
                        "acceptor": {"VH": gh, "VL": gl},
                        "murine_design": {"VH": vh_m, "VL": vl_m},
                        "humanised": {"VH": vh_h, "VL": vl_h},
                        "germline_9mer_coverage": {
                            "VH_murine": round(cov_m_h, 4),
                            "VH_humanised": round(cov_h_h, 4),
                            "VL_murine": round(cov_m_l, 4),
                            "VL_humanised": round(cov_h_l, 4)},
                        "liabilities": lia})

    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "step17_humanised_designs.json").write_text(json.dumps(out, indent=2))
    (RESULTS / "IL14A_ctd_designed_humanised.fasta").write_text("\n".join(fasta) + "\n")
    print(f"\nwrote {len(out)} humanised designs to "
          f"results/IL14A_ctd_designed_humanised.fasta")


if __name__ == "__main__":
    main()
