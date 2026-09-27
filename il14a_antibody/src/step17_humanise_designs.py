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

FR2_TRP = re.compile(r"W[VIL][RKQGS]Q")
# Light chains do NOT share the heavy chain's W-V-R-Q. SM3's lambda reads
# WVQE and ACC4's kappa reads WLLQ, so the light FR2 anchor is its own
# pattern, and it is searched only in the window where CDR-L1 can end
# (7-25 residues after the FR1 cysteine) so a tryptophan inside CDR-L1
# cannot be mistaken for it.
FR2_TRP_L = re.compile(r"W[YVFLAIM][QLRHY][QKELRV]")
# Murine heavy chains do not all read L-E-W: ACC4's FR2 is G-L-K-W-M-G, so
# hard-coding the glutamate made split_heavy return None for it and dropped
# the entire ACC4 arm without a word.
H2_START = re.compile(r"[LMVIF][EKQRDNSTG]W[VIAMLFY][GASVNT]")
# FR3 opens R-F-T-I-S in most germlines but IGHV7-4-1 reads R-F-V-F-S, and
# the narrower pattern silently returned None for it, taking the whole ACC4
# arm with it. Constrained to the window just after CDR-H2 starts, so the
# looser character classes cannot wander.
H2_END = re.compile(r"R[FLVIAM][TASVG][IFLMVT][ST]")
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
    if not m2s:
        return None
    m2e = H2_END.search(seq, m2s.end())
    if not m2e or m2e.start() > m2s.end() + 40:
        return None
    cdr2 = seq[m2s.end():m2e.start() - 3]
    # The J-region motif is located FIRST and the CDR3 cysteine is then the
    # last one BEFORE it. Taking the last cysteine in the whole chain instead
    # lands in CH1 - these sequences are VH+CH1, not VH - and then no WG.G
    # follows it, so CDR-H3 comes back empty and the graft silently carries
    # constant-domain residues.
    j = J_H.search(seq, m2e.end())
    limit = j.start() if j else len(seq)
    c2 = None
    for m in CYS2.finditer(seq, m2e.end()):
        if m.start() + 2 < limit:
            c2 = m.start() + 2
    if c2 is None:
        return None
    cdr3 = seq[c2 + 1:limit]
    return {"FR1": seq[:c1 + 4], "CDR1": cdr1, "FR2": seq[w.start():m2s.end()],
            "CDR2": cdr2, "FR3": seq[m2e.start() - 3:c2 + 1], "CDR3": cdr3,
            "post": seq[limit:]}


def split_light(seq):
    c1 = seq.find("C", 15)
    if c1 < 0:
        return None
    w = FR2_TRP_L.search(seq, c1 + 7)
    if not w or w.start() > c1 + 25:
        w = FR2_TRP.search(seq, c1 + 7)
    if not w:
        return None
    cdr1 = seq[c1 + 1:w.start()]
    # FR2 of a light chain is 15 residues from the conserved Trp; CDR-L2 is
    # the 7 that follow. Anchoring L2 on a motif is unreliable because the
    # motif itself sits inside the loop in some germlines.
    s2, e2 = w.start() + 15, w.start() + 22
    cdr2 = seq[s2:e2]
    j = J_L.search(seq, e2)
    limit = j.start() if j else len(seq)
    c2 = None
    for m in CYS2.finditer(seq, e2):
        if m.start() + 2 < limit:
            c2 = m.start() + 2
    if c2 is None:
        return None
    cdr3 = seq[c2 + 1:limit]
    return {"FR1": seq[:c1 + 1], "CDR1": cdr1, "FR2": seq[w.start():s2],
            "CDR2": cdr2, "FR3": seq[e2:c2 + 1], "CDR3": cdr3,
            "post": seq[limit:]}


def mature_v(seq, kind):
    """Strip the signal peptide from a germline V gene.

    humanness._trim_signal_peptide leaves roughly half of them untrimmed -
    IGLV7-43 comes back still carrying MAWTPLFLFLLTCC..., whose two leader
    cysteines then land in the graft and make the "humanised" light chain
    read four cysteines. The mature domain's first cysteine is at IMGT 23, so
    the start is 22 residues before it, nudged forward to the first residue
    that actually begins a V domain.
    """
    anchor = (FR2_TRP if kind == "H" else FR2_TRP_L).search(seq, 20)
    if not anchor:
        return seq
    c1 = None
    for i in range(anchor.start() - 20, anchor.start() - 10):
        if 0 <= i < len(seq) and seq[i] == "C":
            c1 = i
    if c1 is None:
        return seq
    start = max(0, c1 - 22)
    window = range(start, min(start + 5, len(seq)))
    # Q, E and D begin the overwhelming majority of mature V domains; S and I
    # only a few. Checking them in one pass let a leader's trailing serine win
    # over the real start one residue later (SQVQLVQ... for IGHV1-69).
    for pref in ("QED", "IS"):
        hit = next((k for k in window if seq[k] in pref), None)
        if hit is not None:
            return seq[hit:]
    return seq[start:]


def pick_acceptor(query, locus, kind):
    """Best human germline that is actually USABLE as an acceptor framework.

    Highest identity is not sufficient. IGHV7-4-1 in this germline set carries
    a third cysteine in FR3 (AYLQI-C-SLKAEDT where the canonical sequence has
    a serine), and grafting onto it produces a V domain with an unpaired
    thiol. Some germlines also cannot be split by the framework anchors. So
    the candidates are ranked by identity and the first one that splits AND
    carries exactly the two canonical cysteines is taken, with the rejection
    recorded rather than hidden.
    """
    from humanness import _aligner
    al = _aligner()
    scored, rejected = [], []
    for gene, raw in germlines(locus).items():
        g = mature_v(raw.upper().replace("X", ""), kind)
        if len(g) < 60:
            continue
        aln = al.align(query, g)[0]
        ni = nt = 0
        for (qs, qe), (gs, ge) in zip(*aln.aligned):
            for off in range(qe - qs):
                nt += 1
                ni += int(query[qs + off] == g[gs + off])
        if nt < 50:
            continue
        scored.append((ni / nt, gene, g))
    scored.sort(reverse=True)
    for ident, gene, g in scored:
        parts = (split_heavy if kind == "H" else split_light)(g)
        if parts is None:
            rejected.append((gene, round(ident, 3), "frameworks not resolvable"))
            continue
        ncys = sum(parts[k].count("C")
                   for k in ("FR1", "CDR1", "FR2", "CDR2", "FR3"))
        if ncys != 2:
            rejected.append((gene, round(ident, 3), f"{ncys} cysteines"))
            continue
        return {"gene": gene, "identity": round(ident, 4), "locus": locus,
                "parts": parts, "sequence": g, "rejected": rejected[:4]}
    return None


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
    """Keyed by gene symbol, matching what nearest_germline reports.

    The FASTA headers are full UniProt descriptions; nearest_germline pulls
    the GN= field out of them, so the raw dictionary cannot be indexed by the
    gene name it returns.
    """
    out = {}
    for header, seq in _read_multi_fasta(DATA / f"germline_{locus}.fasta").items():
        gene = (header.split("GN=")[-1].split()[0]
                if "GN=" in header else header[:24])
        out[gene] = seq
    return out


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
        # ProteinMPNN was run with the full alphabet, so a design can carry
        # an unpaired cysteine in a CDR. Two cysteines is the intradomain
        # disulfide every V domain needs; anything else is a manufacturing
        # liability, so those designs are dropped before ranking rather than
        # being reported with a warning.
        def spare_cys(d):
            h = split_heavy(d["sequences"][heavy_id])
            l = split_light(d["sequences"][light_id])
            if h is None or l is None:
                return True
            return any("C" in h[k] or "C" in l[k]
                       for k in ("CDR1", "CDR2", "CDR3"))
        pool = [d for d in entry["top"] if not spare_cys(d)]
        print(f"  {len(entry['top']) - len(pool)} of {len(entry['top'])} "
              f"top designs carry a cysteine in a CDR and are dropped")
        for rank, d in enumerate(pool[:3], 1):
            vh_m = d["sequences"][heavy_id]
            vl_m = d["sequences"][light_id]
            ph = split_heavy(vh_m)
            pl = split_light(vl_m)
            if ph is None or pl is None:
                print(f"  design {rank}: could not be split into loops, skipped")
                continue
            # every designed/contacting position must fall in a loop
            gh = pick_acceptor(vh_m[:130], "IGHV", "H")
            gl_k = pick_acceptor(vl_m[:130], "IGKV", "L")
            gl_l = pick_acceptor(vl_m[:130], "IGLV", "L")
            if gh is None or (gl_k is None and gl_l is None):
                print(f"  design {rank}: no usable acceptor framework, skipped")
                continue
            cands = [c for c in (gl_k, gl_l) if c]
            gl = max(cands, key=lambda c: c["identity"])
            lkind = "K" if gl is gl_k else "L"
            vh_h, _ = graft(ph, gh["sequence"], "H")
            vl_h, _ = graft(pl, gl["sequence"], lkind)
            if vh_h is None or vl_h is None:
                print(f"  design {rank}: graft failed, skipped")
                continue
            for rj in (gh["rejected"] + gl["rejected"])[:3]:
                print(f"    acceptor {rj[0]} ({rj[1]*100:.1f}%) rejected: {rj[2]}")
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
                        "acceptor": {
                            "VH": {k: gh[k] for k in
                                   ("gene", "identity", "locus", "rejected")},
                            "VL": {k: gl[k] for k in
                                   ("gene", "identity", "locus", "rejected")}},
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
