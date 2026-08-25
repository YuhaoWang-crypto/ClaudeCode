#!/usr/bin/env python3
"""
M0 - Fetch every sequence used in the demo from a public, citable source.

Nothing in this pipeline is typed from memory: each test article, benchmark
ligand and assay control is pulled live from RCSB PDB or UniProt and written
to data/sequences.fasta with its accession recorded in
data/sequences_metadata.tsv.

Replace the TEST-ARTICLE entry with your own de-identified FASTA to run the
pipeline on a proprietary ligand; everything downstream is source-agnostic.
"""
import json
import os
import sys
import time
import urllib.parse
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")


def get(url, tries=4):
    last = None
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read().decode("utf-8")
        except Exception as e:  # network flakiness in the sandbox
            last = e
            time.sleep(2 ** i)
    raise RuntimeError(f"failed to fetch {url}: {last}")


def pdb_entity(pdb_id, entity):
    d = json.loads(get(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity}"))
    seq = d["entity_poly"]["pdbx_seq_one_letter_code_can"].replace("\n", "").strip()
    desc = d["rcsb_polymer_entity"]["pdbx_description"]
    return seq, desc


def uniprot(acc):
    txt = get(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta")
    lines = txt.strip().split("\n")
    return "".join(lines[1:]), lines[0][1:]


def strip_tag(seq):
    """Remove common expression tags that are not part of the ligand itself."""
    for tag in ("ALEHHHHHH", "LEHHHHHH", "HHHHHH"):
        if seq.endswith(tag):
            seq = seq[: -len(tag)]
    if seq.startswith("MHHHHHHAM"):
        seq = seq[9:]
    return seq


def main():
    os.makedirs(DATA, exist_ok=True)
    records = []   # (id, role, seq, source, note)

    # ---------------------------------------------------------------- test
    seq, desc = pdb_entity("9DC3", 2)
    records.append(("AAVX_VHH", "test_article", seq, "PDB 9DC3 entity 2",
                    "AAVX affinity ligand (camelid VHH) in complex with AAV8; "
                    "public stand-in for a proprietary AAV affinity ligand"))

    # ----------------------------------------------------------- benchmarks
    seq, _ = pdb_entity("1Q2N", 1)
    records.append(("ProteinA_Z", "benchmark_ligand", seq, "PDB 1Q2N",
                    "Engineered Z domain of staphylococcal protein A - the "
                    "ligand unit of alkali-stable rProtein A resins; leachate "
                    "has decades of clinical exposure history"))

    spa, _ = uniprot("P38507")
    # Native domain B of protein A. Located by its own motif rather than by a
    # hard-coded offset so a UniProt sequence-version bump cannot silently
    # shift the slice.
    key = "ADNKFNKEQQNAFYEILHLPNLNEEQRNGFIQSLKDDPSQSANLLAEAKKLNDAQAPK"
    i = spa.find(key)
    if i < 0:
        sys.exit("protein A domain B motif not found in P38507 - check sequence version")
    records.append(("ProteinA_B_native", "benchmark_ligand", key, f"UniProt P38507 {i+1}-{i+len(key)}",
                    "Native (non-engineered) protein A domain B"))

    pl, _ = uniprot("Q51918")
    key = "EEVTIKANLIFADGSTQNAEFKGTFAKAVSDAYAYADALKKDNGEYTVDVADKGLTLNIKFAG"
    i = pl.find(key)
    if i < 0:
        sys.exit("protein L B1 motif not found in Q51918 - check sequence version")
    records.append(("ProteinL_B1", "benchmark_ligand", key, f"UniProt Q51918 {i+1}-{i+len(key)}",
                    "Ig kappa-light-chain-binding domain B1 of Finegoldia magna "
                    "protein L - ligand of CaptureSelect/Capto L type resins"))

    # ------------------------------------------------- VHH class comparators
    seq, _ = pdb_entity("7EOW", 2)
    seq = strip_tag(seq)
    if seq.startswith("M"):
        seq = seq[1:]
    records.append(("Caplacizumab_VHH", "clinical_anchor", seq, "PDB 7EOW entity 2",
                    "Caplacizumab - humanised bivalent anti-vWF VHH, EMA/FDA "
                    "approved; the only nanobody-class molecule in this set "
                    "with published clinical ADA rates"))

    seq, _ = pdb_entity("4KRL", 1)
    records.append(("VHH_7D12", "class_comparator", strip_tag(seq),
                    "PDB 4KRL entity 1",
                    "Non-humanised camelid anti-EGFR VHH 7D12 - VHH-scaffold "
                    "background control"))

    # ------------------------------------------------------------- controls
    ighv, _ = uniprot("P01764")
    mature = ighv[ighv.find("EVQLVESGGG"):]           # drop the signal peptide
    records.append(("HumanVH3_23_germline", "negative_control_self", mature,
                    "UniProt P01764 (IGHV3-23), mature V region",
                    "Human germline VH3-23 - the human framework a VHH is "
                    "closest to; tolerised floor for framework-derived hits"))

    alb, _ = uniprot("P02768")
    records.append(("HSA_D1", "negative_control_self", alb[24:24 + 130],
                    "UniProt P02768 residues 25-154 (mature albumin D1)",
                    "Human serum albumin domain I - self-protein floor"))

    tt, _ = uniprot("P04958")
    p2 = tt[809:859]     # 810-859, brackets the p2 epitope at 830-844
    p30 = tt[926:976]    # 927-976, brackets the p30 epitope at 947-967
    assert "QYIKANSKFIGITEL" in p2, "TT p2 epitope not where expected"
    assert "FNNFTVSFWLRVPKVSASHLE" in p30, "TT p30 epitope not where expected"
    records.append(("TT_p2_region", "positive_control", p2,
                    "UniProt P04958 residues 810-859",
                    "Tetanus toxin region containing universal T-helper "
                    "epitope p2 (830-844) - promiscuity positive control"))
    records.append(("TT_p30_region", "positive_control", p30,
                    "UniProt P04958 residues 927-976",
                    "Tetanus toxin region containing universal T-helper "
                    "epitope p30 (947-967) - promiscuity positive control"))

    fa = os.path.join(DATA, "sequences.fasta")
    with open(fa, "w") as f:
        for sid, role, seq, src, note in records:
            f.write(f">{sid}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")

    with open(os.path.join(DATA, "sequences_metadata.tsv"), "w") as f:
        f.write("id\trole\tlength\tsource\tnote\tsequence\n")
        for sid, role, seq, src, note in records:
            f.write(f"{sid}\t{role}\t{len(seq)}\t{src}\t{note}\t{seq}\n")

    print(f"wrote {len(records)} sequences to {fa}")
    for sid, role, seq, src, _ in records:
        print(f"  {sid:24s} {role:22s} {len(seq):4d} aa   {src}")


if __name__ == "__main__":
    main()
