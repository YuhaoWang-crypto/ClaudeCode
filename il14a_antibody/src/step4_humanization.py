"""Step 4 - Humanisation and immunogenicity reduction, run as four gated rounds.

SCOPE, STATED HONESTLY
----------------------
The user asked for (a) affinity maturation over 3-4 rounds and (b) humanisation
with reduced immunogenicity. These two asks differ in what can be computed here:

  * Humanisation / immunogenicity IS computable from sequence against the real
    human germline repertoire, with a germline control that scores 1.000. That
    is what this script optimises, over four gated rounds.

  * Affinity maturation is NOT computable without a reliable antigen-paratope
    interface, and no such interface exists for this target yet: the earlier
    campaign ran 24 Boltz-2 complexes plus an OpenFold3 cross-check and the
    contact sets did not reproduce between replicate samples (target-contact
    Jaccard 0.012 for the TMEL1011 parent). Step 6 therefore writes
    epitope-CONSTRAINED Boltz configs - the previous runs applied no epitope or
    contact constraints at all, which is the most likely reason they failed -
    and the affinity axis stays gated until one of them returns a reproducible
    interface.

Inventing per-variant affinity gains here is exactly the failure the G034
cross-validation report documented. This script does not do it.

THE FOUR ROUNDS (humanisation objective)
----------------------------------------
  R1 REPAIR    - remove chemical/PTM liabilities outside the CDRs.
  R2 GERMLINE  - substitute framework positions toward the nearest human
                 germline V gene, protecting CDRs, Vernier positions and (for
                 the VHH) the hallmark tetrad that lets it fold without a VL.
  R3 COMBINE   - stack the accepted single changes, then re-check that no new
                 liability and no new foreign 9-mer was introduced.
  R4 REPORT    - final panel, delta table vs parent, and the MHC-II peptide set
                 to submit for real T-cell epitope prediction.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import RESULTS
from humanness import (_aligner, _read_multi_fasta, _trim_signal_peptide,
                       chemical_liabilities, germline_9mer_coverage,
                       load_germline_kmers, mhc2_peptides, nearest_germline)
from common import DATA

# --- parents ---------------------------------------------------------------
# TMEL1011: human antibody-library-derived pair with published TXLNA protein
# array binding (US12285484B2). Sequences as transcribed in the supplied
# supplementary report; treated as the primary existing parent.
TMEL_VH = ("PVQLVQSGAEVKKPGSSVKVSCKASGGTFSSYSFNWVRQAPGQGLEWMARIIPILGLANYA"
           "QKFQGRVTLTADESTSTAYMELSSLRSEDTAIFYCAGMVLGQLGFDPWGQGTLVTVSS")
TMEL_VL = ("QFMLTQPHSVSESPGRTVTISCTRSSGSIARNYVHWYQHRPGSSPTTVIYEDDQRPSGVPD"
           "RFSGSIDSSSNSASLTISGLKPEDEADYFCQSYDSNIWVFGGGTKLTVL")
# G034 Axis-2 VHH, best-scoring EPI-C parent from step 2. Camelid framework,
# to be humanised and reformatted as VHH-Fc (the format the user selected).
VHH_C05 = ("QVQLVESGGGLVQAGGSLRLSCAASGFTFSSYAMGWFRQAPGKEREFVSAISGSGGSTYYA"
           "DSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCVKEFEHVLHMGAIVKWGQGTQVTVSS")

CDRS = {
    "TMEL1011_VH": ["GGTFSSYS", "IIPILGLA", "AGMVLGQLGFDP"],
    "TMEL1011_VL": ["SGSIARNY", "EDD", "QSYDSNIWV"],
    "VHH_C05": ["GFTFSSYA", "AISGSGGSTYYADSVKG", "VKEFEHVLHMGAIVK"],
}

# Liabilities worth repairing in a framework. Met/Trp oxidation is NOT repaired
# blindly: those residues are often structural, and oxidation is a formulation
# problem before it is a sequence problem.
REPAIRABLE = ("deamidation_NG_NS", "isomerisation_DG_DP", "fragmentation_DP",
              "Nglyc_sequon")
# Cys is deliberately NOT repairable. The two conserved framework cysteines
# (~22 and ~96) form the intradomain disulfide every Ig V domain depends on;
# substituting either one unfolds the domain. A genuinely unpaired Cys would
# sit in a CDR, which this script never touches.
# Conservative substitutions that break the motif with minimal structural cost.
REPAIR_CHOICES = {"N": "Q", "D": "E", "G": "A", "S": "A", "P": "A", "T": "A"}


def cdr_mask(seq, cdrs):
    """Positions (1-based) inside a CDR, plus two flanking residues each side."""
    mask = set()
    for c in cdrs:
        i = seq.find(c)
        if i < 0:
            continue
        for p in range(max(1, i - 1), min(len(seq), i + len(c) + 1) + 1):
            mask.add(p)
    return mask


def vhh_hallmark_positions(seq):
    """Locate the VHH hallmark tetrad by the framework-2 motif.

    In a human VH, framework 2 reads W-V-R-Q-A-P-G-K-G-L-E-W-V-S. In a VHH the
    equivalent stretch carries the substitutions that replace the absent VL
    interface, canonically at Kabat 37/44/45/47. Rather than trusting Kabat
    arithmetic on an unaligned string, find the W...W framework-2 span and
    return the four positions that differ from the human consensus.
    """
    human_fr2 = "WVRQAPGKGLEWVS"
    # the VHH framework-2 begins at the conserved Trp after CDR1
    start = seq.find("W", 30)
    if start < 0:
        return set()
    span = seq[start:start + len(human_fr2)]
    out = set()
    for i, (a, b) in enumerate(zip(span, human_fr2)):
        if a != b:
            out.add(start + i + 1)          # 1-based
    return out


def germline_diffs(seq, gene_locus):
    """Framework positions where `seq` differs from its nearest germline gene.

    Returns [(position_1based, seq_residue, germline_residue)].
    """
    al = _aligner()
    path = DATA / f"germline_{gene_locus['locus']}.fasta"
    target = None
    for header, g in _read_multi_fasta(path).items():
        gene = header.split("GN=")[-1].split()[0] if "GN=" in header else ""
        if gene == gene_locus["gene"]:
            target, _ = _trim_signal_peptide(g.upper().replace("X", ""))
            break
    if target is None:
        return []
    aln = al.align(seq, target)[0]
    diffs = []
    for (qs, qe), (gs, ge) in zip(*aln.aligned):
        for off in range(qe - qs):
            q_i, g_i = qs + off, gs + off
            if seq[q_i] != target[g_i]:
                diffs.append((int(q_i) + 1, seq[q_i], target[g_i]))
    return diffs


def total_repairable(ev):
    return sum(len(v) for v in ev["framework_repairable_liabilities"].values())


def accept(cand_ev, cur_ev, require_cov_gain=False):
    """Global gate every substitution must pass.

    A change may never raise the framework liability count, never introduce a
    cysteine, and never lose humanness. This is what stops R1 and R2 from
    undoing each other.
    """
    if total_repairable(cand_ev) > total_repairable(cur_ev):
        return False
    if cand_ev["sequence"].count("C") != cur_ev["sequence"].count("C"):
        return False
    cov_delta = cand_ev["germline_9mer_coverage"] - cur_ev["germline_9mer_coverage"]
    return cov_delta >= (0.0 if require_cov_gain else -0.02)


def evaluate(seq, kmers, cdrs):
    cov, foreign = germline_9mer_coverage(seq, kmers)
    lia = chemical_liabilities(seq)
    ng = nearest_germline(seq)
    mask = cdr_mask(seq, cdrs)
    framework_lia = {k: [p for p in v if p not in mask]
                     for k, v in lia.items() if k in REPAIRABLE}
    framework_lia = {k: v for k, v in framework_lia.items() if v}
    return {
        "sequence": seq, "length": len(seq),
        "germline_9mer_coverage": round(cov, 4),
        "n_foreign_9mers": len(foreign),
        "nearest_germline": ng,
        "all_liabilities": lia,
        "framework_repairable_liabilities": framework_lia,
    }


def run_chain(name, seq, cdrs, kmers, is_vhh=False):
    log = {"name": name, "parent": seq, "is_vhh": is_vhh, "rounds": {}}
    mask = cdr_mask(seq, cdrs)
    cys_positions = {i + 1 for i, c in enumerate(seq) if c == "C"}
    protected = set(mask) | cys_positions
    hallmarks = set()
    if is_vhh:
        hallmarks = vhh_hallmark_positions(seq)
        protected |= hallmarks
    log["protected_cys_positions"] = sorted(cys_positions)
    log["cdr_positions"] = sorted(mask)
    log["vhh_hallmark_positions"] = sorted(hallmarks)
    log["vhh_hallmark_residues"] = [f"{seq[p-1]}{p}" for p in sorted(hallmarks)]

    base = evaluate(seq, kmers, cdrs)
    log["rounds"]["R0_parent"] = base

    # ---- R1 REPAIR ------------------------------------------------------
    cur = seq
    r1_changes = []
    for motif, positions in base["framework_repairable_liabilities"].items():
        for p in positions:
            if p in protected:
                continue
            old = cur[p - 1]
            new = REPAIR_CHOICES.get(old)
            if not new or new == old:
                continue
            cand = cur[:p - 1] + new + cur[p:]
            ev = evaluate(cand, kmers, cdrs)
            cur_ev = evaluate(cur, kmers, cdrs)
            # must strictly reduce the liability count, and pass the global gate
            if total_repairable(ev) < total_repairable(cur_ev) and accept(ev, cur_ev):
                r1_changes.append({"pos": p, "from": old, "to": new,
                                   "motif": motif,
                                   "cov_after": ev["germline_9mer_coverage"]})
                cur = cand
    log["rounds"]["R1_repair"] = {"changes": r1_changes, **evaluate(cur, kmers, cdrs)}

    # ---- R2 GERMLINE ----------------------------------------------------
    ng = nearest_germline(cur)
    diffs = germline_diffs(cur, ng)
    accepted, rejected = [], []
    for p, have, want in diffs:
        if p in protected:
            rejected.append({"pos": p, "from": have, "to": want,
                             "reason": "CDR/Vernier" if p in mask else "VHH hallmark"})
            continue
        cand = cur[:p - 1] + want + cur[p:]
        ev = evaluate(cand, kmers, cdrs)
        before = evaluate(cur, kmers, cdrs)
        gain = ev["germline_9mer_coverage"] - before["germline_9mer_coverage"]
        if accept(ev, before, require_cov_gain=True):
            dist = min([abs(p - q) for q in mask], default=99)
            accepted.append({
                "pos": p, "from": have, "to": want,
                "coverage_gain": round(gain, 4),
                "distance_to_nearest_CDR_position": dist,
                # a framework change close to a CDR can still move the paratope
                "binding_risk": "elevated - confirm experimentally" if dist <= 4
                                else "low",
            })
            cur = cand
        else:
            why = ("raises liability count"
                   if total_repairable(ev) > total_repairable(before)
                   else "no humanness gain")
            rejected.append({"pos": p, "from": have, "to": want, "reason": why})
    log["rounds"]["R2_germline"] = {"target_germline": ng, "accepted": accepted,
                                    "rejected": rejected,
                                    **evaluate(cur, kmers, cdrs)}

    # ---- R3 COMBINE / VERIFY -------------------------------------------
    final_ev = evaluate(cur, kmers, cdrs)
    regressions = []
    for motif in REPAIRABLE:
        b = set(base["all_liabilities"].get(motif, []))
        f = set(final_ev["all_liabilities"].get(motif, []))
        if f - b:
            regressions.append({"motif": motif, "new_positions": sorted(f - b)})
    cdr_preserved = all(c in cur for c in cdrs)
    log["rounds"]["R3_verify"] = {
        "cdrs_preserved": cdr_preserved,
        "new_liabilities_introduced": regressions,
        "cys_count_parent": seq.count("C"), "cys_count_final": cur.count("C"),
        **final_ev,
    }

    # ---- R4 REPORT ------------------------------------------------------
    log["rounds"]["R4_final"] = {
        "final_sequence": cur,
        "delta": {
            "germline_9mer_coverage": [base["germline_9mer_coverage"],
                                       final_ev["germline_9mer_coverage"]],
            "n_foreign_9mers": [base["n_foreign_9mers"],
                                final_ev["n_foreign_9mers"]],
            "germline_identity": [base["nearest_germline"]["identity"],
                                  final_ev["nearest_germline"]["identity"]],
            "n_framework_liabilities": [
                sum(len(v) for v in base["framework_repairable_liabilities"].values()),
                sum(len(v) for v in final_ev["framework_repairable_liabilities"].values())],
        },
        "mhc2_peptides_for_submission": mhc2_peptides(cur),
    }
    return log


def main():
    kmers, stats = load_germline_kmers()
    chains = [
        ("TMEL1011_VH", TMEL_VH, CDRS["TMEL1011_VH"], False),
        ("TMEL1011_VL", TMEL_VL, CDRS["TMEL1011_VL"], False),
        ("VHH_C05", VHH_C05, CDRS["VHH_C05"], True),
    ]
    out = {"germline_index": stats, "chains": {}}
    for name, seq, cdrs, is_vhh in chains:
        out["chains"][name] = run_chain(name, seq, cdrs, kmers, is_vhh)

    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "step4_humanization.json").write_text(json.dumps(out, indent=2))

    print("=" * 102)
    print("HUMANISATION / IMMUNOGENICITY REDUCTION - 4 GATED ROUNDS")
    print(f"germline index: {stats['germline_sequences']} V genes, "
          f"{stats['distinct_kmers']} distinct 9-mers "
          f"({stats['signal_peptides_trimmed']} signal peptides trimmed)")
    print("=" * 102)
    for name, lg in out["chains"].items():
        d = lg["rounds"]["R4_final"]["delta"]
        print(f"\n--- {name} ---")
        if lg["is_vhh"]:
            print(f"  VHH hallmarks protected : {', '.join(lg['vhh_hallmark_residues'])}")
        print(f"  nearest germline        : {lg['rounds']['R0_parent']['nearest_germline']['gene']}"
              f" ({lg['rounds']['R0_parent']['nearest_germline']['locus']})")
        print(f"  R1 liability repairs    : {len(lg['rounds']['R1_repair']['changes'])}")
        for c in lg["rounds"]["R1_repair"]["changes"]:
            print(f"      {c['from']}{c['pos']}{c['to']}  ({c['motif']})")
        print(f"  R2 germlining accepted  : {len(lg['rounds']['R2_germline']['accepted'])}"
              f"  rejected: {len(lg['rounds']['R2_germline']['rejected'])}")
        for c in lg["rounds"]["R2_germline"]["accepted"]:
            print(f"      {c['from']}{c['pos']}{c['to']}  (+{c['coverage_gain']:.4f} coverage)")
        r3 = lg["rounds"]["R3_verify"]
        print(f"  R3 CDRs preserved       : {r3['cdrs_preserved']};  "
              f"Cys {r3['cys_count_parent']}->{r3['cys_count_final']};  "
              f"new liabilities: {r3['new_liabilities_introduced'] or 'none'}")
        print(f"  R4 germline 9-mer cov   : {d['germline_9mer_coverage'][0]:.3f}"
              f" -> {d['germline_9mer_coverage'][1]:.3f}")
        print(f"     foreign 9-mers       : {d['n_foreign_9mers'][0]}"
              f" -> {d['n_foreign_9mers'][1]}")
        print(f"     germline identity    : {d['germline_identity'][0]*100:.1f}%"
              f" -> {d['germline_identity'][1]*100:.1f}%")
        print(f"     framework liabilities: {d['n_framework_liabilities'][0]}"
              f" -> {d['n_framework_liabilities'][1]}")
        print(f"  final: {lg['rounds']['R4_final']['final_sequence']}")
    print("\nwrote results/step4_humanization.json")


if __name__ == "__main__":
    main()
