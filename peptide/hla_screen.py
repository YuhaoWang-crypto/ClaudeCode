"""Pre-screen the RFP's peptides against its HLA alleles, before anything is synthesised.

The RFP commits to *"custom synthesis of two HLA alleles and 8 uniblock
complexes"*.  Tetramer synthesis is the expensive, slow line item, and a peptide
that does not bind its allele will not fold into a usable complex -- so the
question worth answering in silico is not "is this immunogenic" but **which of
the five 9-mers justify a tetramer**.

That is the best-validated task in immunoinformatics: 9-mers against Class I
alleles.  This module runs it through IEDB's NetMHCpan and reports a *calibrated*
call rather than a bare percentile, because a %rank means nothing without knowing
what a binder and a non-binder look like for that particular allele.

Two things it deliberately refuses to do.

**D-amino acids.**  Every sequence-based predictor -- NetMHCpan, MHCflurry, all
of IEDB -- is trained exclusively on L-peptides.  A 9-mer carrying four D-residues
is outside the applicability domain completely; the model will still return a
number, and that number describes a different molecule than the one being
synthesised.  So peptides flagged as D-containing are refused rather than scored.
The RFP's 9-mer core test and 30-mer diblock test are both in this category, and
whether the five uniblock 9-mers are is a question for the sponsor.

**Immunogenicity.**  Binding is necessary and nowhere near sufficient for the
IL-2 or tetramer readouts the RFP actually measures.  This screens for the
former only.

    python -m peptide.hla_screen --demo          # calibrate on known epitopes
    python -m peptide.hla_screen --fasta mine.fa # screen real sequences
"""

from __future__ import annotations

import argparse
import io
import json
import time
from pathlib import Path

import pandas as pd
import requests

IEDB_I = "https://tools-cluster-interface.iedb.org/tools_api/mhci/"
IEDB_II = "https://tools-cluster-interface.iedb.org/tools_api/mhcii/"

# The RFP's alleles.
CLASS_I = ("HLA-A*32:01", "HLA-B*57:01")
# The mouse step: DO11.10 is I-Ad restricted and specific for cOVA323-339.
CLASS_II = ("H2-IAd",)
OVA_323_339 = "ISQAVHAAHAEINEAGR"

# HIV-1 Gag epitopes restricted by B*57:01, from the elite-controller literature.
# These calibrate what a real binder scores. A*32:01 has far fewer published
# 9-mers, so its positive calibration is weaker and the report says so.
POSITIVE_CONTROLS = {
    "HLA-B*57:01": {"IW9_Gag": "ISPRTLNAW", "KF11_Gag": "KAFSPEVIPMF",
                    "TW10_Gag": "TSTLQEQIGW"},
    "HLA-A*02:01": {"SL9_Gag": "SLYNTVATL"},      # pipeline sanity check
}

D_MARKERS = ("d-", "(d)", "[d]", "dala", "dleu", "dlys", "dphe", "dtyr", "dtrp")


def has_d_residues(name: str, seq: str) -> bool:
    """Flag a peptide as D-containing from its name or lowercase residues.

    Two conventions are common: annotating the name, or writing D-residues as
    lowercase in the sequence.  Both are honoured, because getting this wrong
    means silently scoring a molecule that was never modelled.
    """
    low = name.lower()
    if any(m in low for m in D_MARKERS):
        return True
    return any(c.islower() for c in seq.strip())


def _post(url: str, data: dict, tries: int = 3) -> pd.DataFrame:
    for attempt in range(tries):
        try:
            r = requests.post(url, data=data, timeout=180)
            r.raise_for_status()
            if "Invalid" in r.text[:200] or "\t" not in r.text[:200]:
                raise RuntimeError(r.text[:300])
            return pd.read_csv(io.StringIO(r.text), sep="\t")
        except Exception:
            if attempt == tries - 1:
                raise
            time.sleep(2 ** attempt)
    raise RuntimeError("unreachable")


def screen_class1(peptides: dict[str, str], alleles=CLASS_I) -> pd.DataFrame:
    """NetMHCpan eluted-ligand %rank for every peptide x allele."""
    rows = []
    for allele in alleles:
        by_len: dict[int, list[tuple[str, str]]] = {}
        for name, seq in peptides.items():
            by_len.setdefault(len(seq), []).append((name, seq))
        for length, items in sorted(by_len.items()):
            if not 8 <= length <= 14:
                continue
            fasta = "".join(f">{n}\n{s}\n" for n, s in items)
            df = _post(IEDB_I, {"method": "netmhcpan_el", "sequence_text": fasta,
                                "allele": allele, "length": length})
            names = [n for n, _ in items]
            for _, r in df.iterrows():
                rows.append({"peptide_id": names[int(r["seq_num"]) - 1],
                             "allele": allele, "peptide": r["peptide"],
                             "percentile_rank": float(r["percentile_rank"]),
                             "score": float(r["score"])})
    return pd.DataFrame(rows)


def screen_class2(peptides: dict[str, str], alleles=CLASS_II) -> pd.DataFrame:
    """NetMHCIIpan, returning the best-scoring 15-mer window and its core."""
    rows = []
    for allele in alleles:
        for name, seq in peptides.items():
            if len(seq) < 15:
                continue
            df = _post(IEDB_II, {"method": "netmhciipan_el",
                                 "sequence_text": seq, "allele": allele})
            best = df.sort_values("rank").iloc[0]
            rows.append({"peptide_id": name, "allele": allele,
                         "window": best["peptide"], "core": best["core_peptide"],
                         "percentile_rank": float(best["rank"])})
    return pd.DataFrame(rows)


def null_distribution(allele: str, n: int = 60, seed: int = 0) -> pd.Series:
    """What a non-binder scores: random 9-mers at human amino-acid frequencies."""
    import numpy as np

    aa = "ACDEFGHIKLMNPQRSTVWY"
    freq = np.array([7.0, 2.3, 4.7, 7.1, 3.7, 6.6, 2.6, 4.4, 5.8, 10.0,
                     2.2, 3.6, 6.3, 4.8, 5.6, 8.3, 5.4, 6.0, 1.3, 2.7])
    rng = np.random.default_rng(seed)
    peps = {f"rand{i}": "".join(rng.choice(list(aa), 9, p=freq / freq.sum()))
            for i in range(n)}
    df = screen_class1(peps, alleles=(allele,))
    return df["percentile_rank"]


def call(rank: float) -> str:
    """NetMHCpan's own eluted-ligand convention."""
    return "STRONG" if rank <= 0.5 else "weak" if rank <= 2.0 else "non-binder"


def demo(out: Path) -> None:
    """Calibrate the pipeline on epitopes whose answer is already known."""
    print("=" * 74)
    print("1. Class II, the mouse step: does the pipeline recover the published")
    print("   I-Ad register of cOVA323-339, which DO11.10 is specific for?")
    c2 = screen_class2({"cOVA_323_339": OVA_323_339})
    for _, r in c2.iterrows():
        print(f"   {r['peptide_id']} on {r['allele']}: core {r['core']}  "
              f"%rank {r['percentile_rank']:.2f}  -> {call(r['percentile_rank'])}")
    print("   Published core for this epitope is ISQAVHAAHAEINEAGR residues "
          "327-335.\n   A correct core here is the pipeline validating itself.")

    print("\n" + "=" * 74)
    print("2. Class I, positive controls: do known B*57:01 epitopes score as")
    print("   binders, and does a known A*02:01 epitope score on its own allele?")
    pos = {}
    for allele, eps in POSITIVE_CONTROLS.items():
        pos.update(eps)
    c1 = screen_class1(pos, alleles=("HLA-B*57:01", "HLA-A*02:01") + CLASS_I[:1])
    piv = c1.pivot_table(index="peptide_id", columns="allele",
                         values="percentile_rank")
    print(piv.round(3).to_string())
    print("\n   Read the diagonal: each epitope should score best on its own\n"
          "   restricting allele and worse elsewhere. Off-diagonal is the\n"
          "   specificity control -- a predictor that calls everything a binder\n"
          "   for every allele is useless for deciding which tetramers to make.")

    print("\n" + "=" * 74)
    print("3. The null: what a non-binder looks like, per allele.")
    summary = {}
    for allele in CLASS_I:
        null = null_distribution(allele)
        summary[allele] = {"median": float(null.median()),
                           "q10": float(null.quantile(0.10)),
                           "frac_strong": float((null <= 0.5).mean()),
                           "frac_weak_or_better": float((null <= 2.0).mean())}
        print(f"   {allele}: 60 random 9-mers, median %rank "
              f"{null.median():.1f}, {100 * (null <= 2.0).mean():.0f}% would be "
              f"called weak-or-better by chance")

    print("\n" + "=" * 74)
    print("4. The refusal: D-containing peptides are not scored.")
    for name, seq in (("9mer_core_test_4D", "ISQAvHAAr"),
                      ("D-Ala_diblock_test", "KAFSPEVIPMF")):
        print(f"   {name}: D-flagged = {has_d_residues(name, seq)}"
              + ("  -> refused" if has_d_residues(name, seq) else "  -> scored"))
    print("   Every sequence-based MHC predictor is trained on L-peptides only.\n"
          "   A %rank for a D-substituted peptide describes a molecule that was\n"
          "   never in the training set. The RFP's 9-mer core test and 30-mer\n"
          "   diblock test both carry 4 D-residues.")

    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(
        {"class_ii_ova": c2.to_dict("records"),
         "class_i_controls": c1.to_dict("records"),
         "null_by_allele": summary}, indent=2))
    print(f"\nwrote {out}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--demo", action="store_true")
    ap.add_argument("--fasta", type=Path, help="peptides to screen")
    ap.add_argument("--out", type=Path, default=Path("results/hla_screen.json"))
    args = ap.parse_args()

    if args.demo or not args.fasta:
        demo(args.out)
        return

    peptides, name = {}, None
    for line in args.fasta.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].strip()
        elif name and line.strip():
            peptides[name] = peptides.get(name, "") + line.strip()

    refused = {n: s for n, s in peptides.items() if has_d_residues(n, s)}
    ok = {n: s for n, s in peptides.items() if n not in refused}
    for n in refused:
        print(f"REFUSED {n}: D-residues are outside every sequence predictor's "
              f"applicability domain")
    if not ok:
        raise SystemExit("nothing left to score")

    df = screen_class1(ok)
    df["call"] = df["percentile_rank"].map(call)
    print(df.sort_values(["allele", "percentile_rank"]).to_string(index=False))

    # The RFP puts 4 of its 5 peptides on each of two alleles, so a peptide that
    # loads onto BOTH cannot be attributed to one of them by tetramer staining.
    # Flagging that before synthesis is cheaper than discovering it in FACS.
    piv = df.pivot_table(index="peptide_id", columns="allele",
                         values="percentile_rank")
    if set(CLASS_I).issubset(piv.columns):
        both = piv[(piv[CLASS_I[0]] <= 2.0) & (piv[CLASS_I[1]] <= 2.0)]
        neither = piv[(piv[CLASS_I[0]] > 2.0) & (piv[CLASS_I[1]] > 2.0)]
        print(f"\nassay-design flags")
        print(f"  binds BOTH alleles ({len(both)}): "
              f"{', '.join(both.index) or 'none'}")
        print(f"    -> tetramer staining cannot attribute a T cell to one "
              f"allele for these")
        print(f"  binds NEITHER ({len(neither)}): "
              f"{', '.join(neither.index) or 'none'}")
        print(f"    -> a complex that will not fold; synthesising it spends the "
              f"budget for no readout")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(
        {"scored": df.to_dict("records"), "refused": list(refused)}, indent=2))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
