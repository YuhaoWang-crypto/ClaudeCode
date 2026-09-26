"""Step 1 - Target (epitope) selection on the SECRETED IL-14alpha species.

Why this step exists
--------------------
The earlier G034 pipeline nominated four epitopes (EPI-A..D) and then built its
entire extracellular IgG4 lead series against EPI-B (aa 58-98). Two independent
checks invalidate that choice:

  (1) Accessibility. Secreted IL-14alpha is translated from an alternate start
      at Met176, so the extracellular species spans aa 176-546 only. Epitopes
      at aa 1-17 (EPI-A) and aa 58-98 (EPI-B) are simply NOT PRESENT on the
      antigen a systemic IgG can reach.
  (2) Paralog specificity. EPI-B contains SEELSRQLEDILSTY (aa 72-86), which is
      12/15 identical to beta-taxilin SEELNRQLEDIINTY. Re-derived max 15-mer
      identity for EPI-B is 0.80 against TXLNB - the WORST of the four windows,
      not the best (the legacy report recorded 0.33).

This script therefore re-scans the secreted region from scratch using real
structure and real paralog sequences, and emits a ranked epitope table.

Scoring (all terms in [0,1], transparent and reproducible):
    exposure    = mean relative SASA over the window (Shrake-Rupley on AF model)
    confidence  = mean pLDDT / 100 over the window
    specificity = 1 - max(contiguous 15-mer identity vs TXLNB, vs TXLNG)
    cleanliness = 1 - (chemical liability motifs per residue, capped)
    composite   = specificity^1.5 * exposure * confidence * cleanliness

specificity is raised to 1.5 because paralog cross-reactivity is the dominant
de-risking axis for this target family (3 highly similar coiled-coil paralogs,
all broadly expressed); a beautiful but cross-reactive epitope is worthless.
"""
import json
import re
import sys

from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from common import (DATA, LEGACY_EPITOPES, MAX_ASA, PATENT_ZONE, RESULTS,
                    SECRETED_END, SECRETED_START, load_taxilins)

THREE2ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}

# Chemical/PTM liability motifs. Relevant in an antigen because they make the
# epitope heterogeneous (assay drift, lot-to-lot ELISA variation).
LIABILITY_MOTIFS = {
    "deamidation_NG": r"N[GS]",
    "isomerisation_DG": r"D[GPST]",
    "oxidation_M": r"M",
    "oxidation_W": r"W",
    "free_Cys": r"C",
    "Nglyc_sequon": r"N[^P][ST]",
}


def per_residue_structure(pdb_path):
    """Return {resnum: (aa, rel_sasa, plddt)} from an AlphaFold PDB model."""
    structure = PDBParser(QUIET=True).get_structure("af", str(pdb_path))
    ShrakeRupley().compute(structure[0], level="R")
    out = {}
    for res in structure[0]["A"]:
        name = res.get_resname()
        if name not in THREE2ONE:
            continue
        aa = THREE2ONE[name]
        num = res.id[1]
        plddt = sum(a.get_bfactor() for a in res) / len(res)
        out[num] = (aa, min(res.sasa / MAX_ASA[aa], 1.0), plddt)
    return out


def max_kmer_identity(frag, other, k=15):
    """Worst-case contiguous identity: max over all k-mer pairs.

    This is the metric that matters for cross-reactivity - a single highly
    identical contiguous stretch is enough for a paratope to bind both.
    """
    if len(frag) < k:
        k = len(frag)
    best = 0.0
    best_pair = None
    for i in range(len(frag) - k + 1):
        w = frag[i:i + k]
        for j in range(len(other) - k + 1):
            o = other[j:j + k]
            m = sum(1 for a, b in zip(w, o) if a == b)
            if m > best:
                best = m
                best_pair = (i, j, w, o)
    return best / k, best_pair


def liabilities(frag):
    hits = []
    for name, pat in LIABILITY_MOTIFS.items():
        for m in re.finditer(pat, frag):
            hits.append(f"{name}@{m.start() + 1}")
    return hits


def overlaps(a, b, c, d):
    return not (b < c or d < a)


def score_window(seq, struct, beta, gamma, start, end):
    """start/end are 1-based inclusive residue numbers on TXLNA."""
    frag = seq[start - 1:end]
    res = [struct[i] for i in range(start, end + 1) if i in struct]
    if len(res) != len(frag):
        return None
    exposure = sum(r[1] for r in res) / len(res)
    confidence = sum(r[2] for r in res) / len(res) / 100.0
    id_b, pair_b = max_kmer_identity(frag, beta)
    id_g, pair_g = max_kmer_identity(frag, gamma)
    worst_id = max(id_b, id_g)
    specificity = 1.0 - worst_id
    lia = liabilities(frag)
    cleanliness = max(0.0, 1.0 - 0.5 * len(lia) / len(frag))
    composite = (specificity ** 1.5) * exposure * confidence * cleanliness
    return {
        "start": start, "end": end, "length": len(frag), "sequence": frag,
        "mean_rel_sasa": round(exposure, 4),
        "mean_plddt": round(confidence * 100, 2),
        "max_id_TXLNB": round(id_b, 4),
        "max_id_TXLNG": round(id_g, 4),
        "worst_paralog_id": round(worst_id, 4),
        "specificity": round(specificity, 4),
        "liabilities": lia,
        "cleanliness": round(cleanliness, 4),
        "composite": round(composite, 6),
        "in_secreted_species": SECRETED_START <= start and end <= SECRETED_END,
        "hits_patent_zone": overlaps(start, end, *PATENT_ZONE),
        "worst_TXLNB_window": None if not pair_b else
            {"txlna": pair_b[2], "txlnb": pair_b[3]},
        "worst_TXLNG_window": None if not pair_g else
            {"txlna": pair_g[2], "txlng": pair_g[3]},
    }


def main():
    alpha, beta, gamma = load_taxilins()
    struct = per_residue_structure(DATA / "AF-P40222.pdb")
    assert len(alpha) == 546, len(alpha)

    # --- A. audit the four legacy epitopes -------------------------------
    legacy = {}
    for name, (s, e) in LEGACY_EPITOPES.items():
        rec = score_window(alpha, struct, beta, gamma, s, e)
        legacy[name] = rec

    # --- B. de novo scan of the secreted species -------------------------
    candidates = []
    for length in (14, 16, 18, 20):
        for start in range(SECRETED_START, SECRETED_END - length + 2):
            end = start + length - 1
            rec = score_window(alpha, struct, beta, gamma, start, end)
            if rec is None or rec["hits_patent_zone"]:
                continue
            candidates.append(rec)
    candidates.sort(key=lambda r: -r["composite"])

    # Greedy non-redundant selection: no two picks may overlap.
    picked = []
    for c in candidates:
        if all(not overlaps(c["start"], c["end"], p["start"], p["end"])
               for p in picked):
            picked.append(c)
        if len(picked) == 8:
            break

    RESULTS.mkdir(exist_ok=True)
    out = {
        "target": "IL-14alpha (secreted TXLNA, alternate start Met176)",
        "antigen_span_for_systemic_IgG": [SECRETED_START, SECRETED_END],
        "patent_zone_excluded": list(PATENT_ZONE),
        "legacy_epitope_audit": legacy,
        "denovo_ranked_epitopes": picked,
        "n_windows_scanned": len(candidates),
    }
    (RESULTS / "step1_epitope_scan.json").write_text(json.dumps(out, indent=2))

    # --- report ----------------------------------------------------------
    print("=" * 96)
    print("A. AUDIT OF LEGACY G034 EPITOPES (recomputed from primary sources)")
    print("=" * 96)
    hdr = f"{'epitope':12} {'range':11} {'secreted?':10} {'relSASA':8} {'pLDDT':7} {'paralogID':10} {'composite':10}"
    print(hdr)
    print("-" * 96)
    for name, r in legacy.items():
        flag = "YES" if r["in_secreted_species"] else "*** NO ***"
        print(f"{name:12} {str((r['start'], r['end'])):11} {flag:10} "
              f"{r['mean_rel_sasa']:<8.3f} {r['mean_plddt']:<7.1f} "
              f"{r['worst_paralog_id']:<10.2f} {r['composite']:<10.4f}")
    print()
    print("Key corrections vs the legacy report:")
    b = legacy["G034-EPI-B"]
    print(f"  * EPI-B ({b['start']}-{b['end']}) is absent from secreted IL-14alpha "
          f"(starts at Met{SECRETED_START}).")
    print(f"  * EPI-B worst paralog identity is {b['worst_paralog_id']:.2f}, not 0.33:")
    print(f"      TXLNA {b['worst_TXLNB_window']['txlna']}")
    print(f"      TXLNB {b['worst_TXLNB_window']['txlnb']}")
    print()
    print("=" * 96)
    print("B. DE NOVO RANKED EPITOPES WITHIN THE SECRETED SPECIES (aa 176-546)")
    print("=" * 96)
    print(f"{'rank':5} {'range':11} {'relSASA':8} {'pLDDT':7} {'specif':7} {'composite':10} sequence")
    print("-" * 96)
    for i, r in enumerate(picked, 1):
        print(f"{i:<5} {str((r['start'], r['end'])):11} {r['mean_rel_sasa']:<8.3f} "
              f"{r['mean_plddt']:<7.1f} {r['specificity']:<7.2f} "
              f"{r['composite']:<10.4f} {r['sequence']}")
    print()
    print(f"scanned {len(candidates)} windows; wrote results/step1_epitope_scan.json")


if __name__ == "__main__":
    main()
