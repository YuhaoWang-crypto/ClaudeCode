#!/usr/bin/env python
"""In-silico pre-screen for the peptide / HLA validation RFP.

What this is
------------
The RFP is a request for wet-lab work at a CRO. Most of it cannot be replaced
computationally. This script does the part that genuinely can be: predicting
peptide-MHC binding for the three alleles the RFP names, so that the peptides
sent for tetramer assembly and cell assays are prioritised rather than chosen
blind.

It is a *triage* step that runs in minutes and costs nothing. It does not
replace a single assay in the RFP.

What it runs
------------
1. **Allele support check.** Confirms the three alleles the RFP names are
   available in the prediction backend, against the live allele list.

2. **Positive controls, on the RFP's own system.** The predictor is only
   trustworthy here if it reproduces immunology that is already known for
   *these* alleles:
     - OVA 323-339 on mouse H2-IAd. This is the exact antigen/restriction of
       the DO11.10 hybridoma assay in RFP Step 1.
     - A documented HLA-B*57:01-restricted 9-mer, which should bind B*57:01 and
       not A*32:01.
   If these fail, nothing downstream is worth reading.

3. **Shuffle control.** Each positive control is compared against
   composition-matched shuffles of itself. A predictor that scores the real
   peptide no better than its own anagrams is keying on amino-acid composition,
   not on the binding motif.

4. **9-mer scan** of a source antigen across HLA-A*32:01 and HLA-B*57:01,
   producing a ranked shortlist of the kind the RFP's Step 2 needs.

5. **The stereochemistry blind spot.** The decisive limitation, demonstrated
   rather than asserted. See ``d_amino_acid_blind_spot``.

Backend: IEDB NetMHCpan / NetMHCIIpan cloud REST (no licence, no token).
"""
from __future__ import annotations

import argparse
import io
import json
import os
import random
import sys

import pandas as pd
import requests

MHCI_URL = "https://tools-cluster-interface.iedb.org/tools_api/mhci/"
MHCII_URL = "https://tools-cluster-interface.iedb.org/tools_api/mhcii/"

# The three alleles named in the RFP.
RFP_ALLELES_I = ["HLA-A*32:01", "HLA-B*57:01"]
RFP_ALLELE_II = "H2-IAd"

# IEDB / DTU %rank thresholds.
STRONG_I, WEAK_I = 0.5, 2.0
STRONG_II, WEAK_II = 2.0, 10.0

# The DO11.10 antigen. Note the numbering: the immunology literature numbers
# ovalbumin on the MATURE protein (no initiator Met), so the classic
# "OVA 323-339" is UniProt P01012 residues 324-340. Ordering peptides by
# position rather than by sequence gives a one-residue frameshift.
OVA_323_339 = "ISQAVHAAHAEINEAGR"

# A documented HLA-B*57:01-restricted 9-mer (HIV-1 Gag p24, "IW9").
# Used only as a positive control for the predictor, not as a study peptide.
B5701_POSITIVE = "ISPRTLNAW"


# --------------------------------------------------------------------------- #
def _post(url: str, data: dict, timeout: int = 300) -> pd.DataFrame:
    r = requests.post(url, data=data, timeout=timeout)
    r.raise_for_status()
    txt = r.text.strip()
    if not txt or txt.lower().startswith(("error", "<!doctype")):
        raise RuntimeError(f"backend returned no table: {txt[:200]}")
    return pd.read_csv(io.StringIO(txt), sep="\t")


def predict_mhci(peptides, alleles, method="netmhcpan_el") -> pd.DataFrame:
    """Score fixed-length peptides against MHC-I alleles. Returns %rank per pair."""
    rows = []
    for pep in peptides:
        df = _post(MHCI_URL, {
            "method": method,
            "sequence_text": pep,
            "allele": ",".join(alleles),
            "length": ",".join([str(len(pep))] * len(alleles)),
        })
        df["input_peptide"] = pep
        rows.append(df)
    out = pd.concat(rows, ignore_index=True)
    out = out.rename(columns={"percentile_rank": "rank"})
    out["call"] = out["rank"].apply(
        lambda r: "strong" if r <= STRONG_I else "weak" if r <= WEAK_I else "none")
    return out


def predict_mhcii(sequence: str, allele=RFP_ALLELE_II,
                  method="netmhciipan_el") -> pd.DataFrame:
    """Scan a sequence against an MHC-II allele; returns the best-ranked register."""
    df = _post(MHCII_URL, {
        "method": method, "sequence_text": sequence, "allele": allele}, timeout=600)
    df = df.rename(columns={"rank": "rank"})
    df["call"] = df["rank"].apply(
        lambda r: "strong" if r <= STRONG_II else "weak" if r <= WEAK_II else "none")
    return df.sort_values("rank")


def scan_9mers(protein: str, alleles=RFP_ALLELES_I) -> pd.DataFrame:
    """All 9-mers of a protein against the MHC-I alleles, in one request."""
    df = _post(MHCI_URL, {
        "method": "netmhcpan_el",
        "sequence_text": protein,
        "allele": ",".join(alleles),
        "length": ",".join(["9"] * len(alleles)),
    }, timeout=900)
    df = df.rename(columns={"percentile_rank": "rank"})
    df["call"] = df["rank"].apply(
        lambda r: "strong" if r <= STRONG_I else "weak" if r <= WEAK_I else "none")
    return df


# --------------------------------------------------------------------------- #
def d_amino_acid_blind_spot(peptide: str = B5701_POSITIVE,
                            d_positions=(2, 4, 6, 8)) -> dict:
    """Demonstrate that the predictor cannot see D-amino acid substitution.

    The RFP's central comparison is an all-L peptide against the same peptide
    carrying four D-amino acids. Every sequence-based MHC predictor in routine
    use — NetMHCpan, NetMHCIIpan, MHCflurry — consumes a 20-letter alphabet with
    no stereochemistry channel. A D-residue and its L-enantiomer are the same
    character, so the two molecules are the same input string.

    The prediction is therefore not merely uncertain for the D-variant. It is
    *identical to the L-variant by construction*, which means it carries no
    information about the one variable the study exists to measure, while
    looking exactly like a normal result.

    MHC class I and II both bind peptide backbones through a stereochemically
    specific groove, and D-substitution is used precisely because it changes
    backbone presentation and protease resistance. So the quantity the tool
    silently assumes is unchanged is the quantity under test.
    """
    l_form = peptide
    # A D-substituted peptide has the SAME one-letter sequence; the difference
    # is chirality at those positions, which the encoding cannot express.
    d_form_as_encoded = peptide
    return {
        "peptide": peptide,
        "d_positions": list(d_positions),
        "L_form_string": l_form,
        "D_form_string_as_the_model_sees_it": d_form_as_encoded,
        "strings_identical": l_form == d_form_as_encoded,
        "implication": ("Any sequence-based MHC predictor returns the same number "
                        "for both, so it cannot inform the L-vs-D comparison."),
        "also_unrepresentable": ["N-terminal acetylation", "C-terminal amidation",
                                 "TFA-to-acetate counter-ion form",
                                 "diblock/uniblock assembly state"],
    }


# --------------------------------------------------------------------------- #
def shuffle_control(peptide: str, n: int = 20, seed: int = 0) -> list[str]:
    """Composition-matched anagrams, to test that the motif and not the
    amino-acid composition is driving the score."""
    rng = random.Random(seed)
    out, seen = [], {peptide}
    while len(out) < n:
        s = list(peptide)
        rng.shuffle(s)
        s = "".join(s)
        if s not in seen:
            seen.add(s)
            out.append(s)
    return out


def check_allele_support() -> dict:
    """Confirm the RFP's alleles exist in the live backend allele lists."""
    support = {}
    ri = requests.post(MHCI_URL, data={"method": "netmhcpan_el",
                                       "species": "human"}, timeout=180).text
    for a in RFP_ALLELES_I:
        support[a] = a in ri
    rii = requests.post(MHCII_URL, data={"method": "netmhciipan_el"},
                        timeout=180).text
    support[RFP_ALLELE_II] = RFP_ALLELE_II in rii
    return support


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--antigen-fasta", default=None,
                    help="source protein to scan for 9-mers (default: ovalbumin P01012)")
    ap.add_argument("--peptides", nargs="*", default=None,
                    help="your own 9-mers; replaces the ovalbumin scan")
    ap.add_argument("--top", type=int, default=12)
    ap.add_argument("--out", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "results"))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    report = {}

    print("=" * 74)
    print("1. ALLELE SUPPORT  (the three alleles named in the RFP)")
    print("=" * 74)
    sup = check_allele_support()
    report["allele_support"] = sup
    for a, ok in sup.items():
        print(f"   {a:14s} {'supported' if ok else 'NOT SUPPORTED'}")
    if not all(sup.values()):
        print("\n   Unsupported allele: stop here, the rest is meaningless.")
        return 1

    print()
    print("=" * 74)
    print("2. POSITIVE CONTROLS  (does the predictor reproduce known immunology")
    print("   for THESE alleles, not just in general?)")
    print("=" * 74)

    ii = predict_mhcii(OVA_323_339)
    best = ii.iloc[0]
    print(f"\n   RFP Step 1 antigen: OVA 323-339 = {OVA_323_339}")
    print(f"   on {RFP_ALLELE_II} (the DO11.10 restriction element)")
    print(f"      best register {best['peptide']}  core {best['core_peptide']}")
    print(f"      %rank {best['rank']}  -> {best['call']} binder")
    report["control_ova_iad"] = {"rank": float(best["rank"]),
                                 "core": str(best["core_peptide"]),
                                 "call": str(best["call"])}

    pi = predict_mhci([B5701_POSITIVE], RFP_ALLELES_I)
    print(f"\n   Documented HLA-B*57:01-restricted 9-mer: {B5701_POSITIVE}")
    for _, r in pi.iterrows():
        print(f"      {r['allele']:14s} %rank {r['rank']:>6}  -> {r['call']}")
    report["control_b5701"] = pi[["allele", "rank", "call"]].to_dict("records")

    print()
    print("=" * 74)
    print("3. SHUFFLE CONTROL  (motif, or just composition?)")
    print("=" * 74)
    shuf = shuffle_control(B5701_POSITIVE, n=20)
    sh = predict_mhci(shuf, ["HLA-B*57:01"])
    real = float(pi[pi["allele"] == "HLA-B*57:01"]["rank"].iloc[0])
    better = int((sh["rank"] < real).sum())
    print(f"\n   real peptide  %rank {real}")
    print(f"   20 anagrams   %rank median {sh['rank'].median():.2f}  "
          f"best {sh['rank'].min():.2f}")
    print(f"   anagrams scoring better than the real peptide: {better}/20")
    report["shuffle_control"] = {"real_rank": real,
                                 "shuffle_median": float(sh["rank"].median()),
                                 "shuffle_best": float(sh["rank"].min()),
                                 "n_better": better, "n_shuffles": 20}

    print()
    print("=" * 74)
    print("4. 9-MER RANKING  (the shortlist RFP Step 2 would load onto tetramers)")
    print("=" * 74)
    if args.peptides:
        scan = predict_mhci(args.peptides, RFP_ALLELES_I)
        scan = scan.rename(columns={"input_peptide": "peptide"})
        src = "user-supplied peptides"
    else:
        if args.antigen_fasta:
            seq = "".join(l.strip() for l in open(args.antigen_fasta)
                          if not l.startswith(">"))
            src = os.path.basename(args.antigen_fasta)
        else:
            fa = requests.get("https://rest.uniprot.org/uniprotkb/P01012.fasta",
                              timeout=120).text
            seq = "".join(l.strip() for l in fa.splitlines() if not l.startswith(">"))
            src = "ovalbumin P01012 (example source antigen)"
        scan = scan_9mers(seq)
    print(f"\n   source: {src}")
    for allele in RFP_ALLELES_I:
        sub = scan[scan["allele"] == allele].nsmallest(args.top, "rank")
        n_strong = int((scan[scan["allele"] == allele]["rank"] <= STRONG_I).sum())
        print(f"\n   {allele}   strong binders: {n_strong}")
        print(f"      {'peptide':12s} {'%rank':>7s}  call")
        for _, r in sub.iterrows():
            print(f"      {r['peptide']:12s} {r['rank']:>7}  {r['call']}")
    scan.to_csv(os.path.join(args.out, "mhci_scan.csv"), index=False)

    print()
    print("=" * 74)
    print("5. THE STEREOCHEMISTRY BLIND SPOT  (why this cannot replace the assay)")
    print("=" * 74)
    bs = d_amino_acid_blind_spot()
    report["d_amino_acid_blind_spot"] = bs
    print(f"\n   L-form, as the model sees it : {bs['L_form_string']}")
    print(f"   D-form, as the model sees it : {bs['D_form_string_as_the_model_sees_it']}")
    print(f"   identical input strings      : {bs['strings_identical']}")
    print(f"\n   {bs['implication']}")
    print("\n   Also outside the encoding entirely:")
    for x in bs["also_unrepresentable"]:
        print(f"      - {x}")

    with open(os.path.join(args.out, "report.json"), "w") as f:
        json.dump(report, f, indent=2)
    print(f"\n[saved] {args.out}/report.json and mhci_scan.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
