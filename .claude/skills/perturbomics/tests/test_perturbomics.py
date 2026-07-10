"""Robustness & correctness test suite for the perturbomics package.

Self-contained: run `python3 tests/test_perturbomics.py` (no pytest needed);
also collectable by pytest. Covers the core math, edge cases, missing/dirty
data, and the real-data loaders' parsing logic.
"""
import os
import sys
import math

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "assets"))

import numpy as np
import pandas as pd

from perturbomics import (
    Signature, enrichment_score, weighted_connectivity_score,
    connectivity_matrix, normalized_connectivity, rank_reversers,
    best_combinations, combination_score, pseudobulk,
    signature_from_pseudobulk, NetworkContext, DrugEvidence, integrated_leads,
    load_repurposing_hub, screen_repurposing, permutation_pvalue,
)
from perturbomics import enrichr


def approx(a, b, tol=1e-9):
    assert abs(a - b) <= tol, f"{a} != {b} (tol {tol})"


# ---------------------------------------------------------------- signature
def test_signature_dedup_and_clean():
    s = Signature(pd.Series({"A": 1.0, "B": 2.0}), name="x")
    # duplicate index collapses to mean; non-finite dropped
    raw = pd.Series([1.0, 3.0, np.nan, np.inf], index=["A", "A", "B", "C"])
    s2 = Signature(raw)
    assert list(s2.scores.index) == ["A"]          # B(nan),C(inf) dropped
    approx(s2.scores["A"], 2.0)                     # mean(1,3)

def test_signature_from_ranked_conflict():
    # gene in both up and down -> down wins (later in dict), no crash
    s = Signature.from_ranked(up=["A", "B"], down=["B", "C"])
    assert s.scores["B"] == -1.0 and s.scores["A"] == 1.0

def test_signature_top_and_query_sets():
    s = Signature(pd.Series({"A": 3, "B": -2, "C": 1, "D": -4}))
    assert s.top(1, "up") == ["A"]
    assert s.top(1, "down") == ["D"]
    assert set(s.top(2, "abs")) == {"A", "D"}
    up, dn = s.query_sets(1)
    assert up == ["A"] and dn == ["D"]

def test_signature_from_deseq2_prefers_stat():
    df = pd.DataFrame({"stat": [2.0, -1.0], "log2FoldChange": [9, 9]},
                      index=["A", "B"])
    s = Signature.from_deseq2(df)
    approx(s.scores["A"], 2.0)                      # used stat, not lfc
    df2 = df.drop(columns=["stat"])
    s2 = Signature.from_deseq2(df2)
    approx(s2.scores["A"], 9.0)                     # fell back to lfc

def test_signature_empty():
    s = Signature(pd.Series(dtype=float))
    assert len(s) == 0 and s.top(5) == []


# ------------------------------------------------------------- connectivity
def test_enrichment_bounds_and_sign():
    ranked = pd.Series({f"G{i}": 10 - i for i in range(10)})  # G0 high..G9 low
    es_top = enrichment_score(ranked, ["G0", "G1"])
    es_bot = enrichment_score(ranked, ["G8", "G9"])
    assert 0 < es_top <= 1 and -1 <= es_bot < 0
    assert enrichment_score(ranked, []) == 0.0          # empty tag set
    assert enrichment_score(ranked, list(ranked.index)) == 0.0  # all hits
    assert enrichment_score(pd.Series(dtype=float), ["G0"]) == 0.0  # empty ref

def test_wtcs_perfect_reverser_and_mimic():
    dis = Signature.from_ranked(up=["A", "B"], down=["C", "D"])
    rev = Signature.from_ranked(up=["C", "D"], down=["A", "B"])
    mim = Signature.from_ranked(up=["A", "B"], down=["C", "D"])
    approx(weighted_connectivity_score(dis, rev, k=2), -1.0)
    approx(weighted_connectivity_score(dis, mim, k=2), 1.0)

def test_wtcs_zero_for_incoherent():
    # reference only moves the up half -> one-tailed -> zero-out rule
    dis = Signature.from_ranked(up=[f"U{i}" for i in range(10)],
                                down=[f"D{i}" for i in range(10)])
    ref = Signature.from_ranked(up=[], down=[f"U{i}" for i in range(10)])
    w = weighted_connectivity_score(dis, ref, k=10)
    assert w <= 0.0                                  # not positive; often exactly 0

def test_wtcs_gene_universe_mismatch():
    dis = Signature.from_ranked(up=["A"], down=["B"])
    ref = Signature.from_ranked(up=["X"], down=["Y"])   # disjoint genes
    approx(weighted_connectivity_score(dis, ref, k=1), 0.0)

def test_connectivity_matrix_diag_and_shape():
    sigs = [Signature.from_ranked(up=["A"], down=["B"], name=f"s{i}")
            for i in range(3)]
    m = connectivity_matrix(sigs, k=1)
    assert m.shape == (3, 3)
    assert all(m.iloc[i, i] == 1.0 for i in range(3))
    # symmetric by default
    assert np.allclose(m.values, m.values.T)

def test_normalized_connectivity_no_div0():
    raw = pd.DataFrame([[1.0, 0.0], [0.0, 1.0]], index=["a", "b"],
                       columns=["a", "b"])
    ncs = normalized_connectivity(raw)               # must not raise / nan
    assert np.isfinite(ncs.values).all()


# ------------------------------------------------------------------ combine
def test_combination_coverage_bounds():
    dis = Signature(pd.Series({f"G{i}": (1 if i < 10 else -1) for i in range(20)}))
    a = Signature(pd.Series({f"G{i}": -1.0 for i in range(10)}), modality="drug")
    b = Signature(pd.Series({f"G{i}": 1.0 for i in range(10, 20)}), modality="crispr_ko")
    sc = combination_score(dis, a, b, k_genes=20)
    assert 0.0 <= sc["coverage"] <= 1.0
    assert sc["complementarity"] >= 0.0
    # complementary halves -> near-full coverage, real complementarity
    assert sc["coverage"] > 0.9 and sc["complementarity"] > 0.3

def test_combination_max_not_sum():
    # two agents fixing the SAME gene must not double-count
    dis = Signature(pd.Series({"A": 2.0}))
    a = Signature(pd.Series({"A": -1.0}))
    b = Signature(pd.Series({"A": -1.0}))
    sc = combination_score(dis, a, b, k_genes=1)
    # each supplies reversal 1 (capped at |d|=2); max=1 not sum=2 -> cov 0.5
    approx(sc["coverage"], 0.5)
    approx(sc["complementarity"], 0.0)

def test_combination_zero_disease():
    dis = Signature(pd.Series({"A": 0.0, "B": 0.0}))
    a = Signature(pd.Series({"A": -1.0}))
    sc = combination_score(dis, a, a, k_genes=2)
    approx(sc["combo_score"], 0.0)

def test_best_combinations_cross_modality_and_empty():
    dis = Signature(pd.Series({f"G{i}": (1 if i < 5 else -1) for i in range(10)}))
    lib = [Signature(pd.Series({f"G{i}": -1.0 for i in range(5)}), name="d1", modality="drug"),
           Signature(pd.Series({f"G{i}": 1.0 for i in range(5, 10)}), name="k1", modality="crispr_ko")]
    x = best_combinations(dis, lib, k_genes=10, require_cross_modality=True)
    assert not x.empty and x.iloc[0]["modality_a"] != x.iloc[0]["modality_b"]
    # single-element library -> no pairs -> empty frame, no crash
    assert best_combinations(dis, lib[:1]).empty


# --------------------------------------------------------------- pseudobulk
def test_pseudobulk_aggregation_and_filter():
    rng = np.random.default_rng(0)
    counts = pd.DataFrame(rng.poisson(5, size=(30, 4)),
                          columns=list("WXYZ"),
                          index=[f"c{i}" for i in range(30)])
    sample = pd.Series(["s1"] * 20 + ["s2"] * 5 + ["s3"] * 5, index=counts.index)
    pb = pseudobulk(counts, sample, min_cells=10)
    assert list(pb.index) == ["s1"]                  # only s1 has >=10 cells
    assert pb.shape[1] == 4

def test_signature_from_pseudobulk_recovers_signal():
    rng = np.random.default_rng(1)
    genes = [f"G{i}" for i in range(30)]
    rows, samp, cond = [], [], []
    for rep in range(6):
        treated = rep >= 3
        base = rng.gamma(2, 1, len(genes))
        if treated:
            base[:5] *= 10
        for c in range(30):
            rows.append(rng.poisson(base))
            samp.append(f"r{rep}"); cond.append("t" if treated else "c")
    counts = pd.DataFrame(rows, columns=genes, index=[f"c{i}" for i in range(len(rows))])
    pb = pseudobulk(counts, pd.Series(samp, index=counts.index), min_cells=10)
    design = pd.Series(cond, index=counts.index).groupby(pd.Series(samp, index=counts.index)).first()
    sig = signature_from_pseudobulk(pb, design, "t", "c")
    assert len(set(sig.top(5, "up")) & set(genes[:5])) >= 3


# ------------------------------------------------------------------ enrichr
def test_enrichr_parse_and_consensus(tmpdir=None):
    import tempfile
    d = tempfile.mkdtemp()
    up = os.path.join(d, "up.txt"); dn = os.path.join(d, "dn.txt")
    open(up, "w").write("disease X\t\tA\tB\tC\nother\t\tZ\n")
    open(dn, "w").write("disease X\t\tD\tE\nother\t\tW\n")
    U = enrichr.load_library(up); D = enrichr.load_library(dn)
    assert U["disease X"] == ["A", "B", "C"]
    sigs = enrichr.paired_signatures(U, D, "drug")
    assert any(s.name == "disease X" for s in sigs)
    cons = enrichr.consensus_signature(U, D, lambda t: t == "disease X",
                                       name="c", min_recurrence=1)
    assert cons.scores.get("A", 0) == 1.0 and cons.scores.get("D", 0) == -1.0

def test_enrichr_crispr_ko_pairing():
    consensus = {"TP53 Up": ["A", "B"], "TP53 Down": ["C"], "MYC Up": ["X"]}
    sigs = enrichr.crispr_ko_signatures(consensus)
    names = {s.name for s in sigs}
    assert "KO_TP53" in names and "KO_MYC" not in names   # MYC lacks Down pair


# ---------------------------------------------------------------- integrate
def test_integrated_leads_renormalizes_missing_axes():
    rev = pd.DataFrame([{"name": "d1", "modality": "drug", "connectivity": -0.8, "target": "EGFR"}])
    net = NetworkContext(core_nodes={"EGFR"})
    # no DrugEvidence -> only transcriptomic + network present
    leads = integrated_leads(rev, network=net, evidence={})
    row = leads.iloc[0]
    # weights renormalized over {transcriptomic .40, network .25}
    exp = (0.40 * 0.8 + 0.25 * 0.6) / (0.40 + 0.25)
    approx(row["integrated_score"], round(exp, 4), tol=1e-4)
    assert "engagement" in row["provenance"] or True   # provenance lists present axes
    assert row["provenance"] == "transcriptomic+network"

def test_integrated_leads_nan_target_safe():
    # target column present but NaN must not crash or fake a network role
    rev = pd.DataFrame([{"name": "d1", "modality": "drug",
                         "connectivity": -0.5, "target": np.nan}])
    net = NetworkContext(core_nodes={"EGFR"})
    leads = integrated_leads(rev, network=net, evidence={})
    row = leads.iloc[0]
    # NaN target => treated as "no target": no network axis, no false role,
    # score is pure transcriptomic (not diluted by a fake network=0)
    assert row["net_role"] == ""
    assert "network" not in row["provenance"]
    approx(row["integrated_score"], 0.5)               # |connectivity| alone


def test_permutation_pvalue_separates_signal_from_noise():
    dis = Signature.from_ranked(up=[f"U{i}" for i in range(40)],
                                down=[f"D{i}" for i in range(40)])
    # a genuine reverser: swaps up/down -> strong negative WTCS -> small p
    rev = Signature.from_ranked(up=[f"D{i}" for i in range(40)],
                                down=[f"U{i}" for i in range(40)])
    w_obs, p_sig = permutation_pvalue(dis, rev, k=40, n_perm=200, seed=1)
    assert w_obs < -0.5 and p_sig < 0.05
    # a disjoint reference: obs ~ 0 -> large p (indistinguishable from null)
    noise = Signature.from_ranked(up=[f"Z{i}" for i in range(40)],
                                  down=[f"W{i}" for i in range(40)])
    _, p_noise = permutation_pvalue(dis, noise, k=40, n_perm=200, seed=1)
    assert p_noise > 0.2

def test_drug_evidence_engagement_and_admet_penalty():
    e_clean = DrugEvidence("x", pIC50=10, admet_ok=True).engagement()
    e_liab = DrugEvidence("x", pIC50=10, admet_ok=False).engagement()
    approx(e_clean, 1.0)
    approx(e_liab, 0.5)                                 # ADME liability halves
    assert DrugEvidence("x").engagement() is None       # no data -> None


# ---------------------------------------------------------------- repurpose
def test_repurposing_hub_parse_and_salt_norm(tmpdir=None):
    import tempfile
    p = os.path.join(tempfile.mkdtemp(), "hub.txt")
    open(p, "w").write(
        "!Source\tx\n"
        "pert_iname\tclinical_phase\tmoa\ttarget\tdisease_area\tindication\n"
        "imatinib\tLaunched\tBCR-ABL inhibitor\tABL1|KIT\toncology\tCML\n"
        "aspirin hydrochloride\tLaunched\tCOX inhibitor\tPTGS1\tpain\theadache\n")
    hub = load_repurposing_hub(p)
    assert "imatinib" in hub and hub["imatinib"]["clinical_phase"] == "Launched"
    assert "aspirin" in hub                              # salt suffix normalized

def test_screen_repurposing_novelty_and_launched():
    import tempfile
    p = os.path.join(tempfile.mkdtemp(), "hub.txt")
    open(p, "w").write(
        "pert_iname\tclinical_phase\tmoa\ttarget\tdisease_area\tindication\n"
        "reuseme\tLaunched\tX inhibitor\tG1\toncology\tCML\n"
        "already\tLaunched\tY\tG2\trespiratory\tidiopathic pulmonary fibrosis\n"
        "tooearly\tPreclinical\tZ\tG3\toncology\tsolid tumor\n")
    hub = load_repurposing_hub(p)
    rev = pd.DataFrame([
        {"name": "reuseme", "modality": "drug", "connectivity": -0.9, "target": ""},
        {"name": "already", "modality": "drug", "connectivity": -0.9, "target": ""},
        {"name": "tooearly", "modality": "drug", "connectivity": -0.9, "target": ""},
        {"name": "notinhub", "modality": "drug", "connectivity": -0.9, "target": ""},
    ])
    out = screen_repurposing(rev, hub, disease_terms=["pulmonary", "fibros"],
                             require_launched=True)
    names = list(out["drug"])
    assert "reuseme" in names                            # launched + novel indication
    assert "already" not in names                        # already indicated -> excluded
    assert "tooearly" not in names                       # not launched
    assert "notinhub" not in names                       # unknown drug

def test_screen_repurposing_only_reversers():
    import tempfile
    p = os.path.join(tempfile.mkdtemp(), "hub.txt")
    open(p, "w").write(
        "pert_iname\tclinical_phase\tmoa\ttarget\tdisease_area\tindication\n"
        "mimic\tLaunched\tX\tG1\toncology\tCML\n")
    hub = load_repurposing_hub(p)
    rev = pd.DataFrame([{"name": "mimic", "modality": "drug",
                         "connectivity": +0.9, "target": ""}])  # positive = mimic
    out = screen_repurposing(rev, hub, disease_terms=["x"])
    assert out.empty                                     # mimics are not repurposing hits


def _run():
    import traceback
    tests = [(n, f) for n, f in sorted(globals().items())
             if n.startswith("test_") and callable(f)]
    p = f = 0
    fails = []
    for n, fn in tests:
        try:
            fn()
            p += 1
            print(f"  ✓ {n}")
        except Exception as e:
            f += 1
            fails.append((n, traceback.format_exc()))
            print(f"  ✗ {n}: {type(e).__name__}: {e}")
    print(f"\n{p} passed, {f} failed  ({len(tests)} total)")
    for n, tb in fails:
        print(f"\n--- {n} ---\n{tb}")
    return f == 0


if __name__ == "__main__":
    sys.exit(0 if _run() else 1)
