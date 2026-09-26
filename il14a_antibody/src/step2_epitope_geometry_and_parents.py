"""Step 2 - Epitope 3D footprint + parent (existing antibody) triage.

Part A. Epitope geometry
------------------------
The lead epitope IL14A-EPI-N1 (aa 367-382) sits on a coiled-coil helix. Before
designing a paratope we measure the real 3D footprint from the AlphaFold model:
how far apart the specificity-determining residues (SDRs) and the functional
hotspots actually are. That distance sets how long a CDR3 has to be to touch
them all, instead of guessing a length range.

Part B. Parent triage
---------------------
The user asked to optimise an EXISTING antibody rather than go de novo. Of the
candidates in the G034 package, only one series is usable:

  * Axis-1 IgG4 series (5 leads) - all bind EPI-B (aa 58-98). EPI-B is absent
    from secreted IL-14alpha (starts at Met176) AND is 80% identical to a
    beta-taxilin 15-mer. RETIRED, not matured. Their CDR-H3s are also strongly
    cationic (net +2.5 to +5.5, up to 46% R+K), a known polyreactivity and
    fast-clearance liability.
  * Axis-2 VHH series - the EPI-C members (aa 358-384) overlap the new lead
    epitope aa 367-382 and are therefore genuine parents. They were written as
    camelid intrabodies, so they need humanising and reformatting, which is
    exactly what steps 4-5 do.

This script scores the EPI-C VHH parents on the criteria that matter for the
NEW epitope and picks the ones to carry into maturation.
"""
import json
import re
import sys
from itertools import combinations

from Bio.PDB import PDBParser

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from common import DATA, RESULTS, gravy, net_charge

# --- lead epitope, fixed by step 1b -----------------------------------------
EPITOPE = (367, 382)
EPITOPE_SEQ = "VESQRMCELMKQQETH"
# exposed TXLNA-unique positions inside the epitope (paralog discrimination)
SDR = [367, 370, 371, 372, 373, 375, 382]
# residues that also form the syntaxin-H3 interface (function / blocking)
HOTSPOT = [371, 372, 375]
# paralog residues at the SDR positions, from the step1b global alignment
PARALOG_AT_SDR = {
    367: {"txlnb": "A", "txlng": "T"},
    370: {"txlnb": "K", "txlng": "R"},
    371: {"txlnb": "L", "txlng": "H"},
    372: {"txlnb": "Q", "txlng": "K"},
    373: {"txlnb": "A", "txlng": "Y"},
    375: {"txlnb": "V", "txlng": "Q"},
    382: {"txlnb": "V", "txlng": "Q"},
}

# --- the EPI-C-directed VHH parents from G034_blocking_antibody_Axis2.json ---
VHH_FRAMEWORK_PRE = ("QVQLVESGGGLVQAGGSLRLSCAASGFTFSSYAMGWFRQAPGKEREFVSAIS"
                     "GSGGSTYYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYC")
VHH_FRAMEWORK_POST = "WGQGTQVTVSS"

PARENTS = {
    "G034-Ab-Ax2-C-01": {"cdr3": "YQEHAGDKQKFAKVYKEKI", "dg": -29.94},
    "G034-Ab-Ax2-C-02": {"cdr3": "AYTVGATAYIAVEG", "dg": -30.20},
    "G034-Ab-Ax2-C-03": {"cdr3": "IANFGQHRQNLYGNEAYLD", "dg": -28.50},
    "G034-Ab-Ax2-C-04": {"cdr3": "ACTANKVRGYHFTTSDEA", "dg": -27.60},
    "G034-Ab-Ax2-C-05": {"cdr3": "VKEFEHVLHMGAIVK", "dg": -32.71},
    "G034-Ab-Ax2-C-06": {"cdr3": "IQFDVGLLLATYWD", "dg": -32.60},
    "G034-Ab-Ax2-C-07": {"cdr3": "MCFFIAAGQSGINFK", "dg": -31.91},
    "G034-Ab-Ax2-C-08": {"cdr3": "RFAIWVDDHFTHHKN", "dg": -31.43},
    "G034-Ab-Ax2-C-09": {"cdr3": "QRRNELTHHMVESYFVDKHL", "dg": -30.53},
}

# hard developability filters (applied to CDR3 only; framework is fixed)
LIABILITY_PATTERNS = {
    "free_Cys": r"C",
    "deamidation": r"N[GS]",
    "isomerisation": r"D[GPST]",
    "Nglyc_sequon": r"N[^P][ST]",
    "oxidation_MW": r"[MW]",
    "hydrophobic_run": r"[ILVFWMY]{4,}",
}


def liability_report(seq):
    out = {}
    for name, pat in LIABILITY_PATTERNS.items():
        hits = [m.start() + 1 for m in re.finditer(pat, seq)]
        if hits:
            out[name] = hits
    return out


def epitope_geometry(pdb_path):
    """Measure the real 3D footprint of the epitope on the AlphaFold model."""
    structure = PDBParser(QUIET=True).get_structure("af", str(pdb_path))
    chain = structure[0]["A"]

    def anchor(resnum):
        """CB where present, else CA - the side-chain projection point."""
        res = chain[resnum]
        return res["CB"].get_coord() if "CB" in res else res["CA"].get_coord()

    pts = {p: anchor(p) for p in SDR}
    dists = {}
    for a, b in combinations(sorted(SDR), 2):
        d = float(((pts[a] - pts[b]) ** 2).sum() ** 0.5)
        dists[f"{a}-{b}"] = round(d, 2)
    span = max(dists.values())

    def cluster_span(positions):
        p = {q: anchor(q) for q in positions}
        return max(float(((p[a] - p[b]) ** 2).sum() ** 0.5)
                   for a, b in combinations(sorted(positions), 2))

    hot_span = cluster_span(HOTSPOT)
    # CORE footprint = the syntaxin hotspots plus the SDRs within the same one
    # or two helix turns. This is the sub-surface a single CDR3 can realistically
    # cover while delivering BOTH blocking (hotspots) and paralog discrimination.
    CORE = [370, 371, 372, 373, 375]
    core_span = cluster_span(CORE)

    def cdr3_length_for(span_a):
        """Residues of CDR3 needed to cover a surface span of `span_a` Angstrom.

        A CDR3 is a loop anchored at BOTH ends in the framework beta-sheet: it
        exits, traverses the epitope, and returns. So only about half its
        residues advance across the surface, at ~3.3 A rise per residue in an
        extended loop conformation. Four residues are spent on the take-off and
        re-entry turns and contribute no traversal.

            span = ((L - 4) / 2) * 3.3   =>   L = 2*span/3.3 + 4

        A naive span/3.3 (one-way) model is wrong here and returns lengths far
        below the 14-22 residues actually observed for VHH CDR3s.
        """
        return int(round(2.0 * span_a / 3.3)) + 4

    return {
        "sdr_pairwise_CB_distances_A": dists,
        "max_sdr_span_A": round(span, 2),
        "hotspot_span_A": round(hot_span, 2),
        "core_footprint_positions": CORE,
        "core_footprint_span_A": round(core_span, 2),
        "loop_model": "L = 2*span/3.3 + 4 (out-and-back loop, 3.3 A/residue)",
        "cdr3_len_for_core_footprint": cdr3_length_for(core_span),
        "cdr3_len_for_full_sdr_span": cdr3_length_for(span),
        "cdr3_len_for_hotspots_only": cdr3_length_for(hot_span),
    }


def score_parent(name, cdr3, geom):
    """Score an existing VHH parent against the NEW epitope requirements."""
    lia = liability_report(cdr3)
    q = net_charge(cdr3)
    g = gravy(cdr3)
    # epitope net charge over the contact face -> complementary paratope charge
    epi_q = net_charge(EPITOPE_SEQ)
    # ideal paratope charge is roughly the negative of the epitope charge, but
    # capped: strongly cationic CDRs are the Axis-1 failure mode.
    charge_fit = max(0.0, 1.0 - abs(q - (-epi_q)) / 4.0)
    # length fit against the measured footprint
    want = geom["cdr3_len_for_core_footprint"]
    want_max = geom["cdr3_len_for_full_sdr_span"]
    if want <= len(cdr3) <= want_max:
        length_fit = 1.0
    else:
        off = min(abs(len(cdr3) - want), abs(len(cdr3) - want_max))
        length_fit = max(0.0, 1.0 - off / 6.0)
    # hydrophobicity window: too hydrophobic aggregates, too polar loses affinity
    gravy_fit = max(0.0, 1.0 - abs(g - (-0.2)) / 2.0)
    # blocking hard filters
    blockers = [k for k in ("free_Cys", "Nglyc_sequon", "hydrophobic_run") if k in lia]
    return {
        "id": name, "cdr3": cdr3, "cdr3_len": len(cdr3),
        "legacy_dG_REU": PARENTS[name]["dg"],
        "net_charge": round(q, 2), "gravy": round(g, 2),
        "charge_fit": round(charge_fit, 3),
        "length_fit": round(length_fit, 3),
        "gravy_fit": round(gravy_fit, 3),
        "liabilities": lia,
        "hard_blockers": blockers,
        "parent_fitness": round(charge_fit * length_fit * gravy_fit, 4),
        "usable_as_parent": not blockers,
    }


def main():
    geom = epitope_geometry(DATA / "AF-P40222.pdb")
    scored = [score_parent(n, v["cdr3"], geom) for n, v in PARENTS.items()]
    scored.sort(key=lambda r: (-r["usable_as_parent"], -r["parent_fitness"]))
    selected = [r for r in scored if r["usable_as_parent"]][:4]

    out = {
        "lead_epitope": {
            "id": "IL14A-EPI-N1", "range": list(EPITOPE), "sequence": EPITOPE_SEQ,
            "sdr_positions": SDR, "syntaxin_hotspots": HOTSPOT,
            "paralog_residues_at_sdr": PARALOG_AT_SDR,
            "epitope_net_charge_pH74": round(net_charge(EPITOPE_SEQ), 2),
        },
        "geometry": geom,
        "retired_series": {
            "Axis-1_IgG4_EPI-B": "EPI-B absent from secreted IL-14a (Met176 start) "
                                 "and 80% identical to a beta-taxilin 15-mer; all 5 "
                                 "leads also strongly cationic (net +2.5..+5.5).",
        },
        "parent_triage": scored,
        "selected_parents": [r["id"] for r in selected],
    }
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "step2_parents.json").write_text(json.dumps(out, indent=2))

    print("=" * 100)
    print("A. MEASURED 3D FOOTPRINT OF IL14A-EPI-N1 (aa 367-382) on AF-P40222 v6")
    print("=" * 100)
    print(f"  epitope sequence           : {EPITOPE_SEQ}")
    print(f"  epitope net charge (pH7.4) : {net_charge(EPITOPE_SEQ):+.1f}")
    print(f"  SDR positions              : {SDR}")
    print(f"  syntaxin-H3 hotspots       : {HOTSPOT}")
    print(f"  max SDR CB-CB span         : {geom['max_sdr_span_A']} A")
    print(f"  hotspot-only span          : {geom['hotspot_span_A']} A")
    print(f"  core footprint {geom['core_footprint_positions']} span : {geom['core_footprint_span_A']} A")
    print(f"  loop model                 : {geom['loop_model']}")
    print(f"  => CDR3 length window      : "
          f"{geom['cdr3_len_for_core_footprint']}-"
          f"{geom['cdr3_len_for_full_sdr_span']} residues"
          f"  (hotspots alone would need {geom['cdr3_len_for_hotspots_only']})")
    print()
    print("=" * 100)
    print("B. TRIAGE OF EXISTING EPI-C VHH PARENTS AGAINST THE NEW EPITOPE")
    print("=" * 100)
    print(f"{'id':20} {'len':4} {'netQ':6} {'GRAVY':7} {'fit':7} {'use?':5} blockers / liabilities")
    print("-" * 100)
    for r in scored:
        lia = ",".join(r["hard_blockers"]) or ",".join(r["liabilities"].keys()) or "clean"
        print(f"{r['id']:20} {r['cdr3_len']:<4} {r['net_charge']:<6.1f} {r['gravy']:<7.2f} "
              f"{r['parent_fitness']:<7.3f} {'yes' if r['usable_as_parent'] else 'NO':5} {lia}")
    print()
    print(f"selected parents for maturation: {', '.join(r['id'] for r in selected)}")
    print("wrote results/step2_parents.json")


if __name__ == "__main__":
    main()
