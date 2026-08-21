"""Offline self-test: `python -m neoantigen_pipeline.selftest`.

No network, no cached data. Checks the parts that are pure logic, because
those are the parts where a silent bug would corrupt a real patient's payload
without anything looking wrong: peptide tiling, wild-type pairing, the self
k-mer index, junction enumeration, codon optimization's translation invariant,
the ordering optimizer, and the AUC.
"""

from __future__ import annotations

import sys

import numpy as np
import pandas as pd

from . import benchmark as B
from . import construct as C
from . import features as F
from . import peptides as P
from . import score as S
from . import select as SEL
from . import selfindex
from .config import Gates, SelectionRules

FAILS = []


def check(name, cond, detail=""):
    status = "ok  " if cond else "FAIL"
    print(f"  [{status}] {name}" + (f"   {detail}" if detail and not cond else ""))
    if not cond:
        FAILS.append(name)


def test_peptides():
    print("peptides")
    proteome = {"GENE1": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"}
    assert proteome["GENE1"][9] == "R"          # the variant below must match
    var = pd.DataFrame([{"gene": "GENE1", "protein_change": "R10E", "ref_aa": "R",
                         "aa_pos": 10, "alt_aa": "E",
                         "variant_class": "Missense_Mutation"}])
    pep, skipped = P.generate_peptides(var, proteome, lengths=(9,), class2_lengths=())
    check("9 windows cover the mutation", len(pep) == 9, f"got {len(pep)}")
    if pep.empty:
        FAILS.append("peptide tiling produced nothing; later checks are vacuous")
        return
    check("every mutant peptide carries the mutation",
          all(r.mut_peptide[r.mut_offset] == "E" for r in pep.itertuples()))
    check("wild-type counterpart differs at exactly one position",
          all(sum(a != b for a, b in zip(r.mut_peptide, r.wt_peptide)) == 1
              for r in pep.itertuples()))
    check("wild-type peptide is a real substring of the reference",
          all(w in proteome["GENE1"] for w in pep["wt_peptide"]))

    bad = pd.DataFrame([{"gene": "GENE1", "protein_change": "A10E", "ref_aa": "A",
                         "aa_pos": 10, "alt_aa": "E",   # proteome has R at 10
                         "variant_class": "Missense_Mutation"}])
    _, sk = P.generate_peptides(bad, proteome, lengths=(9,), class2_lengths=())
    check("reference mismatch is reported, not silently applied",
          len(sk) == 1 and "mismatch" in sk.iloc[0]["reason"])

    fs = pd.DataFrame([{"gene": "GENE1", "protein_change": "Q10fs", "ref_aa": None,
                        "aa_pos": None, "alt_aa": None,
                        "variant_class": "Frame_Shift_Del"}])
    _, sk2 = P.generate_peptides(fs, proteome, lengths=(9,), class2_lengths=())
    check("frameshift is reported as needing transcript annotation",
          len(sk2) == 1 and "frameshift" in sk2.iloc[0]["reason"])

    anch = P.annotate_anchor(pep)
    check("anchor annotation marks P2 and P-omega",
          set(anch.loc[anch["mut_at_anchor"], "mut_offset"]) <= {0, 1, 8})


def test_selfindex():
    print("self k-mer index")
    proteome = {"A": "MKTAYIAKQRQISFVKSHFSRQ", "B": "GILGFVFTLAAAAAAAA"}
    idx = selfindex.SelfKmerIndex(proteome, 9)
    check("self peptide is found", "GILGFVFTL" in idx)
    check("non-self peptide is not found", "WWWWWWWWW" not in idx)
    check("one-mismatch lookup finds the wild type",
          idx.one_mismatch("GILGFVFTW") == ["GILGFVFTL"])
    check("a peptide that IS self returns no wild type",
          idx.one_mismatch("GILGFVFTL") == [])
    check("contains_many agrees with contains",
          list(idx.contains_many(["GILGFVFTL", "WWWWWWWWW"])) == [True, False])


def test_features():
    print("features")
    check("presentation is monotone decreasing in %rank",
          F.f_presentation(0.05) > F.f_presentation(0.5) > F.f_presentation(2.0))
    check("agretopicity rewards mutation-created binding",
          F.f_agretopicity(0.1, 10.0) > F.f_agretopicity(1.0, 1.0) > F.f_agretopicity(10.0, 0.1))
    check("missing wild-type gives a neutral 0.5", F.f_agretopicity(0.1, None) == 0.5)
    check("expression saturates and is monotone",
          0 == F.f_expression(0) < F.f_expression(5) < F.f_expression(50) <= 1.0)
    check("radical substitution scores more dissimilar than conservative",
          F.f_dissimilarity("AAAWAAAAA", "AAAGAAAAA") > F.f_dissimilarity("AAAIAAAAA", "AAAVAAAAA"))
    check("identical peptides are maximally self-like",
          F.f_dissimilarity("AAAAAAAAA", "AAAAAAAAA") == 0.0)
    r = F.tcr_prior_scores(["GILGFVFTL", "WWWWWWWWW"], ["GILGFVFTL"])
    check("TCR prior saturates on an exact reference match", r[0] > 0.99)
    check("TCR prior is near zero for an unrelated peptide", r[1] < 0.01)


def test_score_and_select():
    print("score and selection")
    d = pd.DataFrame({
        "var_id": [f"G{i}:X{i}Y" for i in range(6)],
        "gene": ["G0", "G0", "G0", "G1", "G2", "G3"],
        "allele": ["HLA-A*02:01"] * 4 + ["HLA-B*07:02"] * 2,
        "mut_peptide": [f"AAAAAAAA{c}" for c in "ABCDEF"],
        "mut_rank": [0.01, 0.02, 0.03, 0.04, 5.0, 0.06],
        "tpm": [100] * 6, "ccf": [1.0] * 6, "is_novel_vs_self": [True] * 6,
        "feat_presentation": [0.9, 0.85, 0.8, 0.75, 0.1, 0.7],
        "feat_agretopicity": [0.5] * 6, "feat_expression": [0.8] * 6,
        "feat_clonality": [1.0] * 6, "feat_dissimilarity": [0.5] * 6,
        "feat_tcr_prior": [0.0] * 6, "feat_hydrophobicity": [0.5] * 6,
        "feat_mhc2_support": [0.0] * 6, "binder": ["strong"] * 6,
    })
    gated = S.apply_peptide_gates(d, Gates())
    check("the 5%-rank candidate is gated out", int(gated["passes"].sum()) == 5)
    w = {"presentation": 0.5, "agretopicity": 0.1, "expression": 0.15,
         "clonality": 0.1, "dissimilarity": 0.05, "tcr_prior": 0.05,
         "hydrophobicity": 0.03, "mhc2_support": 0.02}
    sc = S.composite_score(gated[gated["passes"]], w)
    check("score ordering follows presentation here",
          list(sc.sort_values("neo_score", ascending=False)["gene"])[0] == "G0")

    rules = SelectionRules(max_neoantigens=4, max_per_gene=2, min_alleles_covered=2)
    sel = SEL.select_neoantigens(sc, rules, ["HLA-A*02:01", "HLA-B*07:02"])
    check("the per-gene cap is honoured", int((sel["gene"] == "G0").sum()) <= 2)
    check("the payload size cap is honoured", len(sel) <= 4)
    check("both alleles get covered", sel["allele"].nunique() == 2)
    check("every slot records why it was taken", sel["why_selected"].notna().all())

    forced = SelectionRules(max_neoantigens=2, max_per_gene=2,
                            force_include_genes=("G3",))
    sel2 = SEL.select_neoantigens(sc, forced, ["HLA-A*02:01", "HLA-B*07:02"])
    check("a forced gene makes it into the payload", "G3" in set(sel2["gene"]))


def test_construct():
    print("construct")
    a, b = "A" * 25, "C" * 25
    j = C.junction_peptides(a, b, "", 9)
    check("a direct fusion creates 8 novel 9-mers", len(j) == 8, f"got {len(j)}")
    check("no junction peptide exists in either minigene",
          all(p not in a and p not in b for p in j))
    check("a linker creates more junction peptides than a direct fusion",
          len(C.junction_peptides(a, b, "GPGPG", 9)) > len(j))

    prot = "M" * 50 + "K" + "A" * 50
    mg, off = C.minigene(prot, 50, 25)
    check("minigene has the requested length", len(mg) == 25)
    check("mutation is centered", off == 12 and mg[off] == "K")
    check("minigene is a real substring of the mutant protein", mg in prot)
    short, off2 = C.minigene("MKTAYIAKQRQ", 1, 25)
    check("a short protein yields a shorter minigene, not padding", len(short) == 11)

    protein = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK"
    opt = C.codon_optimize(protein, seed=3, avoid=("EcoRI", "BsaI", "BsmBI"))
    check("codon optimization preserves the protein exactly",
          C.translate(opt["dna"]) == protein)
    check("CDS length is 3x the protein", len(opt["dna"]) == 3 * len(protein))
    check("no avoided restriction site survives",
          all(C.RESTRICTION[s] not in opt["dna"] for s in ("EcoRI", "BsaI", "BsmBI")))

    cost = np.array([[0, 9, 1], [9, 0, 9], [1, 9, 0]], float)
    order = C.order_minimizing_junctions(cost, restarts=3)
    tot = sum(cost[order[i], order[i + 1]] for i in range(len(order) - 1))
    check("ordering finds the cheap path", tot <= 10, f"cost {tot}")


def test_auc():
    print("benchmark metrics")
    check("perfect separation gives AUC 1", B.auc([1, 1, 0, 0], [4, 3, 2, 1]) == 1.0)
    check("inverted separation gives AUC 0", B.auc([1, 1, 0, 0], [1, 2, 3, 4]) == 0.0)
    check("all ties give AUC 0.5", B.auc([1, 0, 1, 0], [1, 1, 1, 1]) == 0.5)
    d = pd.DataFrame({"label": [1] * 5 + [0] * 95, "s": list(range(100))[::-1]})
    e = B.topk_enrichment(d, "s", k=5)
    check("top-k enrichment computes correctly",
          e["precision_at_k"] == 1.0 and abs(e["fold_enrichment"] - 20.0) < 1e-9)


def main() -> int:
    for fn in (test_peptides, test_selfindex, test_features, test_score_and_select,
               test_construct, test_auc):
        fn()
    print()
    if FAILS:
        print(f"{len(FAILS)} FAILED: {', '.join(FAILS)}")
        return 1
    print("all self-tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
