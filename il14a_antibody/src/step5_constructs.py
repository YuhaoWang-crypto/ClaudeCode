"""Step 5 - Assemble the deliverable constructs and QC them.

Lead format is the humanised VHH-Fc the user selected: the humanised VHH from
step 4 fused to a human IgG1 hinge-CH2-CH3, giving a bivalent ~80 kDa molecule
with IgG-like pharmacokinetics. Humanised VHH-Fc / VHH-based biologics have
regulatory precedent (caplacizumab, ozoralizumab).

A humanised full-length IgG1/lambda from the TMEL1011 parent is also assembled,
because TMEL1011 is the one existing parent with BOTH a public paired sequence
and reported TXLNA binding.

CONSTANT REGION PROVENANCE - and one trap
-----------------------------------------
IGHG1 is UniProt P01857, IGLC2 is P0DOY2, IGKC is P01834. P01857 is the
MEMBRANE-BOUND heavy chain entry: after the CH3 domain it continues into the
transmembrane and cytoplasmic exons
(ELQLEESCAEAQDGELDGLWTTITIFITLFLLSVCYSATV...). A secreted antibody must stop at
the CH3 C-terminus. Every constant region here is therefore truncated at the
'...HYTQKSLSLSP' motif and completed with the canonical 'GK', and the assembled
sequences are asserted to contain none of the membrane tail.
"""
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA, RESULTS
from humanness import (_read_multi_fasta, chemical_liabilities,
                       germline_9mer_coverage, load_germline_kmers)

MEMBRANE_TAIL_START = "ELQLEESCAEAQ"
CH3_END_MOTIF = "SLSLSP"
HINGE_START = "EPKSCDKTHTCPPCP"
CH1_START = "ASTKGPSVFPLAP"

# Average residue masses (Da) for a quick MW estimate, plus water.
AA_MASS = {
    "A": 71.08, "R": 156.19, "N": 114.10, "D": 115.09, "C": 103.14, "Q": 128.13,
    "E": 129.12, "G": 57.05, "H": 137.14, "I": 113.16, "L": 113.16, "K": 128.17,
    "M": 131.19, "F": 147.18, "P": 97.12, "S": 87.08, "T": 101.10, "W": 186.21,
    "Y": 163.18, "V": 99.13,
}
PKA = {"D": 3.9, "E": 4.3, "C": 8.3, "Y": 10.1, "H": 6.0, "K": 10.5, "R": 12.5}


def load_constant(acc):
    seq = next(iter(_read_multi_fasta(DATA / f"C_{acc}.fasta").values())).upper()
    return seq


def truncate_secreted(seq):
    """Cut a heavy-chain constant region at the secreted CH3 C-terminus."""
    j = seq.find(CH3_END_MOTIF)
    if j < 0:
        raise ValueError("CH3 end motif not found")
    return seq[:j + len(CH3_END_MOTIF)] + "GK"


def mw_kda(seq):
    return round((sum(AA_MASS[c] for c in seq) + 18.02) / 1000.0, 2)


def isoelectric_point(seq):
    """Bisection on net charge using standard side-chain pKa values."""
    counts = {aa: seq.count(aa) for aa in PKA}

    def charge(ph):
        z = 1.0 / (1.0 + 10 ** (ph - 8.0))          # N-terminus
        z -= 1.0 / (1.0 + 10 ** (3.1 - ph))         # C-terminus
        for aa, n in counts.items():
            if not n:
                continue
            if aa in "KRH":
                z += n / (1.0 + 10 ** (ph - PKA[aa]))
            else:
                z -= n / (1.0 + 10 ** (PKA[aa] - ph))
        return z

    lo, hi = 2.0, 13.0
    for _ in range(80):
        mid = (lo + hi) / 2
        if charge(mid) > 0:
            lo = mid
        else:
            hi = mid
    return round((lo + hi) / 2, 2)


def qc(name, seq, kmers, variable_len=None):
    """QC an assembled chain.

    Two checks that look like defects but are not, and are handled explicitly:

    * Cys parity. An antibody chain normally carries an ODD number of
      cysteines, because the hinge cysteines and the CH1<->light-chain cysteine
      form INTERCHAIN disulfides with the partner chain. Parity per chain is
      therefore meaningless. What is meaningful is that each VARIABLE domain
      carries exactly the two conserved cysteines of its intradomain disulfide.

    * N-glycosylation sequon. Wild-type IgG1 Fc has one conserved sequon at
      Asn297 (in the 'QYNSTYR' motif) and it is required, not a liability. It
      is reported separately from any sequon in a variable domain, which would
      be a genuine finding.
    """
    cov = None
    var_cys = None
    if variable_len:
        var = seq[:variable_len]
        cov = round(germline_9mer_coverage(var, kmers)[0], 4)
        var_cys = var.count("C")
    lia = chemical_liabilities(seq)
    sequons = lia.get("Nglyc_sequon", [])
    n297 = seq.find("QYNSTYR")
    n297_pos = n297 + 3 if n297 >= 0 else None      # the Asn in QY-N-STYR
    variable_sequons = [p for p in sequons
                        if variable_len and p <= variable_len]
    return {
        "name": name, "length": len(seq), "mw_kDa_monomer": mw_kda(seq),
        "predicted_pI": isoelectric_point(seq),
        "cys_total": seq.count("C"),
        "cys_in_variable_domain": var_cys,
        "variable_domain_disulfide_intact": var_cys == 2 if var_cys is not None else None,
        "n_glyc_sequons_all": sequons,
        "conserved_Fc_N297_sequon_at": n297_pos,
        "n_glyc_sequons_in_variable_domain": variable_sequons,
        "deamidation_NG_NS": lia.get("deamidation_NG_NS", []),
        "isomerisation_DG_DP": lia.get("isomerisation_DG_DP", []),
        "variable_domain_germline_9mer_coverage": cov,
        "contains_membrane_tail": MEMBRANE_TAIL_START in seq,
        "sequence": seq,
    }


def main():
    kmers, _ = load_germline_kmers()
    hum = json.loads((RESULTS / "step4_humanization.json").read_text())["chains"]
    vhh = hum["VHH_C05"]["rounds"]["R4_final"]["final_sequence"]
    vh = hum["TMEL1011_VH"]["rounds"]["R4_final"]["final_sequence"]
    vl = hum["TMEL1011_VL"]["rounds"]["R4_final"]["final_sequence"]

    ighg1_full = load_constant("P01857")
    iglc2 = load_constant("P0DOY2")
    heavy_constant = truncate_secreted(ighg1_full[ighg1_full.find(CH1_START):])
    fc_hinge_ch2_ch3 = truncate_secreted(ighg1_full[ighg1_full.find(HINGE_START):])

    # Effector-silenced Fc. LALA (L234A/L235A, EU numbering) is anchored here by
    # the unambiguous 'APELLGG' hinge-CH2 junction motif rather than by position
    # arithmetic. P329G and N297A/aglycosyl are further options but their EU
    # positions must be confirmed against a numbered reference before ordering.
    assert fc_hinge_ch2_ch3.count("APELLGG") == 1
    fc_lala = fc_hinge_ch2_ch3.replace("APELLGG", "APEAAGG")

    constructs = {
        "IL14A-VHHFc-01_hIgG1wt": vhh + fc_hinge_ch2_ch3,
        "IL14A-VHHFc-01_hIgG1_LALA": vhh + fc_lala,
        "IL14A-TMEL1011hz_HeavyChain_IgG1": vh + heavy_constant,
        "IL14A-TMEL1011hz_LightChain_lambda": vl + iglc2,
    }
    var_len = {
        "IL14A-VHHFc-01_hIgG1wt": len(vhh),
        "IL14A-VHHFc-01_hIgG1_LALA": len(vhh),
        "IL14A-TMEL1011hz_HeavyChain_IgG1": len(vh),
        "IL14A-TMEL1011hz_LightChain_lambda": len(vl),
    }

    report = {
        "constant_region_provenance": {
            "IGHG1": "UniProt P01857 (membrane-bound entry; truncated at CH3 + GK)",
            "IGLC2": "UniProt P0DOY2",
            "membrane_tail_excluded": MEMBRANE_TAIL_START,
        },
        "variable_domains": {"humanised_VHH": vhh,
                             "humanised_TMEL1011_VH": vh,
                             "humanised_TMEL1011_VL_lambda": vl},
        "constructs": {},
    }
    print("=" * 100)
    print("ASSEMBLED CONSTRUCTS")
    print("=" * 100)
    print(f"{'construct':38} {'len':5} {'kDa':7} {'pI':6} {'Cys':5} {'V-SS':5} {'V-sequon':9} {'memTail':7}")
    print("-" * 100)
    all_ok = True
    for name, seq in constructs.items():
        rec = qc(name, seq, kmers, var_len[name])
        report["constructs"][name] = rec
        # meaningful gates: conserved variable-domain disulfide present, no
        # sequon introduced into a variable domain, no membrane tail.
        ok = (rec["variable_domain_disulfide_intact"]
              and not rec["n_glyc_sequons_in_variable_domain"]
              and not rec["contains_membrane_tail"])
        all_ok &= bool(ok)
        print(f"{name:38} {rec['length']:<5} {rec['mw_kDa_monomer']:<7.1f} "
              f"{rec['predicted_pI']:<6.2f} {rec['cys_total']:<5} "
              f"{str(rec['variable_domain_disulfide_intact']):5} "
              f"{len(rec['n_glyc_sequons_in_variable_domain']):<9} "
              f"{str(rec['contains_membrane_tail']):7}")
    print()
    print(f"VHH-Fc bivalent assembled mass  : "
          f"~{2 * report['constructs']['IL14A-VHHFc-01_hIgG1wt']['mw_kDa_monomer']:.0f} kDa (homodimer)")
    h = report["constructs"]["IL14A-TMEL1011hz_HeavyChain_IgG1"]["mw_kDa_monomer"]
    l = report["constructs"]["IL14A-TMEL1011hz_LightChain_lambda"]["mw_kDa_monomer"]
    print(f"full IgG1 assembled mass        : ~{2 * (h + l):.0f} kDa (H2L2)")
    fc = report["constructs"]["IL14A-VHHFc-01_hIgG1wt"]
    print(f"conserved Fc N297 sequon at    : residue {fc['conserved_Fc_N297_sequon_at']} "
          f"(required for wild-type IgG1, not a liability)")
    print(f"odd total Cys counts are expected: hinge and H-L cysteines form "
          f"INTERCHAIN disulfides")
    print(f"all constructs pass V-domain disulfide / V-sequon / membrane-tail "
          f"checks: {all_ok}")

    # --- FASTA deliverable ----------------------------------------------
    lines = ["# IL-14alpha (TXLNA) humanised blocking antibody candidates",
             "# Variable domains humanised in step4; constant regions from UniProt.",
             "# Affinity is NOT experimentally or structurally validated - see REPORT.md",
             ""]
    for name, seq in report["variable_domains"].items():
        lines.append(f">{name} | variable domain only")
        lines += [seq[i:i + 60] for i in range(0, len(seq), 60)]
    for name, rec in report["constructs"].items():
        lines.append(f">{name} | len={rec['length']} pI={rec['predicted_pI']} "
                     f"MW={rec['mw_kDa_monomer']}kDa")
        s = rec["sequence"]
        lines += [s[i:i + 60] for i in range(0, len(s), 60)]
    (RESULTS / "IL14A_humanised_candidates.fasta").write_text("\n".join(lines) + "\n")
    (RESULTS / "step5_constructs.json").write_text(json.dumps(report, indent=2))
    print("\nwrote results/step5_constructs.json and results/IL14A_humanised_candidates.fasta")


if __name__ == "__main__":
    main()
