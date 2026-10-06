"""Tests for the Boltz measured-negative pilot.

The load-bearing ones:

* `test_variant_sequences_are_distinct_and_mutations_verified` -- every variant
  is built by checking the stated parent residue, so a transcription slip in the
  mutation list cannot silently produce a wrong protein.
* `test_exact_p_matches_the_closed_form` -- the panel's power claim. Hanley-
  McNeil is degenerate at AUC 1.0 with n this small, so the pilot's headline
  ("even a perfect result gives p = 0.1") rests on this enumeration.
* `test_pocket_constraint_is_zero_indexed` -- an off-by-one here would
  constrain the ligand against the wrong residues and quietly invalidate the run.
"""

from __future__ import annotations

import math
import pathlib
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "pylrs"))

from trnaplat import boltz_pocket as bp  # noqa: E402


def test_mc_sequence_is_the_verified_319mer():
    sequence = bp.mc_sequence()
    assert len(sequence) == 319
    assert sequence[32] == "Y", "Y33 in 1-based numbering"
    assert sequence[111] == "Y", "Y112"
    assert sequence[161] == "D", "D162"
    assert sequence[165] == "L", "L166"


def test_apply_mutations_refuses_a_wrong_parent_residue():
    sequence = bp.mc_sequence()
    with pytest.raises(ValueError, match="McTyrRS has"):
        bp.apply_mutations(sequence, "A33G")      # position 33 is Y, not A


def test_variant_sequences_are_distinct_and_mutations_verified():
    rows = bp.variants()
    assert len(rows) == 6, "wild type plus the five measured variants"
    sequences = {v["name"]: v["sequence"] for v in rows}
    assert len(set(sequences.values())) == 6, "two variants share a sequence"
    assert all(len(s) == 319 for s in sequences.values())
    # Mut6+RM is Mut6 plus Y112F, so they differ at exactly one position
    diff = [i for i, (a, b) in enumerate(zip(sequences["Mc Mut6"],
                                             sequences["Mc Mut6+RM"])) if a != b]
    assert diff == [111]
    assert sequences["Mc Mut6"][111] == "Y"
    assert sequences["Mc Mut6+RM"][111] == "F"


def test_ground_truth_covers_exactly_the_measured_variants():
    labels = [v["label"] for v in bp.variants() if v["label"] is not None]
    assert labels.count(1) == 2, "Mut6 and Mut6+RM are the selected pair"
    assert labels.count(0) == 3, "Mut4, Mut5 and Mut7 are the measured rejections"


def test_f_ratios_match_the_published_counts():
    by_name = {v["name"]: v for v in bp.variants()}
    assert by_name["Mc Mut6"]["f_ratio"] == pytest.approx(41 / 37, abs=1e-9)
    assert by_name["Mc Mut6+RM"]["f_ratio"] == pytest.approx(177 / 48, abs=1e-9)
    assert by_name["Mc Mut6+RM"]["f_ratio"] > by_name["Mc Mut6"]["f_ratio"]


def test_pocket_constraint_is_zero_indexed():
    variant = bp.variants()[0]
    spec = bp.job_spec(variant, bp.PAZF_SMILES, "pAzF", bp.MC_POCKET_POSITIONS)
    indices = spec["constraints"][0]["contact_residues"]["A"]
    assert indices == [p - 1 for p in bp.MC_POCKET_POSITIONS]
    assert 32 in indices, "Y33 must appear as index 32"
    # and those indices must really be the pocket residues in the sequence
    sequence = variant["sequence"]
    assert sequence[32] == "Y" and sequence[161] == "D"


def test_job_spec_declares_the_ligand_as_the_binder():
    variant = bp.variants()[0]
    spec = bp.job_spec(variant, bp.PAZF_SMILES, "pAzF", None)
    assert spec["binding"] == {"type": "ligand_protein_binding",
                               "binder_chain_id": "B"}
    assert "constraints" not in spec
    kinds = [e["type"] for e in spec["entities"]]
    assert kinds == ["protein", "ligand_smiles"]


def test_both_ligands_are_scored_for_every_variant():
    """Selection measured discrimination, so Tyr is not optional."""
    assert bp.PAZF_SMILES != bp.TYR_SMILES
    # same backbone, differing only in the para substituent
    assert bp.PAZF_SMILES.startswith("N[C@@H](Cc1ccc(")
    assert bp.TYR_SMILES.startswith("N[C@@H](Cc1ccc(")


def test_auc_counts_ties_as_half():
    assert bp.auc([3.0, 2.0, 1.0], [1, 0, 0]) == 1.0
    assert bp.auc([1.0, 3.0, 2.0], [1, 0, 0]) == 0.0
    assert bp.auc([1.0, 1.0], [1, 0]) == 0.5


def test_auc_is_nan_without_both_classes():
    assert math.isnan(bp.auc([1.0, 2.0], [1, 1]))


def test_exact_p_matches_the_closed_form():
    # A perfect separation is 1 of C(5,2) = 10 arrangements.
    assert bp.exact_p(1.0, 2, 3) == pytest.approx(0.1)
    # And it cannot reach significance, which is the pilot's headline.
    assert bp.exact_p(1.0, 2, 3) > 0.05
    # A bigger panel can: 7 positives against 20 negatives.
    assert bp.exact_p(1.0, 7, 20) < 1e-5


def test_exact_p_is_monotone_in_the_observed_auc():
    values = [bp.exact_p(a, 2, 3) for a in (1.0, 0.833, 0.667, 0.5)]
    assert values == sorted(values), "a worse AUC cannot have a smaller p"


def test_power_note_states_the_limit_numerically():
    note = bp.power_note()
    assert "p = 0.100" in note
    assert "cannot reach p < 0.05" in note


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))


# ------------------------------------------------- the measured run's outcome

RESULTS = ROOT / "trnaplat/data/boltz_results.csv"


def _wide(frame):
    """One row per variant, one column per ligand."""
    return frame.pivot_table(index="name", columns="ligand", values="score")


def test_measured_run_puts_a_rejected_variant_on_top():
    """❌ Pins the negative result, so a later change cannot quietly soften it."""
    if not RESULTS.exists():
        pytest.skip("no Boltz results recorded yet")
    import pandas as pd

    frame = pd.read_csv(RESULTS)
    truth = {v["name"]: v["label"] for v in bp.variants()}
    pazf = frame[(frame["ligand"] == "pAzF")
                 & frame["name"].map(lambda n: truth.get(n) is not None)]
    pazf = pazf.sort_values("score", ascending=False)
    top = pazf.iloc[0]
    assert truth[top["name"]] == 0, (
        "the headline finding is that the top-scoring variant is a measured "
        "rejection; if that changed, rewrite the module docstring"
    )
    assert top["name"] == "Mc Mut7"


def test_measured_run_auc_and_p_are_what_the_docstring_claims():
    if not RESULTS.exists():
        pytest.skip("no Boltz results recorded yet")
    import pandas as pd

    frame = pd.read_csv(RESULTS)
    truth = {v["name"]: v["label"] for v in bp.variants()}
    rows = [(r.score, truth[r.name]) for r in frame.itertuples()
            if r.ligand == "pAzF" and truth.get(r.name) is not None]
    value = bp.auc([s for s, _ in rows], [l for _, l in rows])
    assert value == pytest.approx(2 / 3, abs=1e-9)
    assert bp.exact_p(value, 2, 3) == pytest.approx(0.4, abs=1e-9)


def test_both_boltz_metrics_give_the_same_ranking():
    """So the negative result is not an artifact of which score was picked."""
    if not RESULTS.exists():
        pytest.skip("no Boltz results recorded yet")
    import pandas as pd

    frame = pd.read_csv(RESULTS)
    pazf = frame[frame["ligand"] == "pAzF"]
    by_binding = list(pazf.sort_values("score", ascending=False)["name"])
    by_opt = list(pazf.sort_values("optimization_score", ascending=False)["name"])
    # Mut4 and Mut5 swap at the bottom; the top three are identical, which is
    # what the claim rests on.
    assert by_binding[:3] == by_opt[:3]
    assert by_binding[0] == "Mc Mut7"


def test_structure_confidence_is_uniformly_high_and_uninformative():
    if not RESULTS.exists():
        pytest.skip("no Boltz results recorded yet")
    import pandas as pd

    frame = pd.read_csv(RESULTS)
    conf = frame["structure_confidence"]
    assert conf.min() > 0.93, "every complex folds confidently"
    assert conf.max() - conf.min() < 0.05, (
        "the spread is too small to rank variants, which is why it was not used"
    )


# --------------------------------------------- the pAzF-minus-Tyr comparator


def test_the_tyr_arm_is_complete():
    """The delta needs both ligands for all six proteins, or it is not a delta."""
    if not RESULTS.exists():
        pytest.skip("no Boltz results recorded yet")
    import pandas as pd

    wide = _wide(pd.read_csv(RESULTS))
    assert set(wide.columns) == {"pAzF", "Tyr"}
    assert len(wide) == 6, "five measured variants plus the wild-type control"
    assert not wide.isna().any().any(), "a missing cell silently drops a variant"


def test_wild_type_control_prefers_tyrosine():
    """✅ The one free sanity check: a tyrosyl-tRNA synthetase must prefer Tyr.

    If this flipped, the delta would not be measuring discrimination at all and
    nothing downstream of it could be read.
    """
    if not RESULTS.exists():
        pytest.skip("no Boltz results recorded yet")
    import pandas as pd

    wide = _wide(pd.read_csv(RESULTS))
    row = wide.loc["Mc wild type"]
    assert row["pAzF"] - row["Tyr"] < 0
    # and every pocket mutant must have gone the other way
    mutants = wide.drop(index="Mc wild type")
    assert ((mutants["pAzF"] - mutants["Tyr"]) > 0).all()


def test_the_margin_is_chance_and_loses_to_the_absolute_score():
    """❌ The headline of the Tyr arm.

    The ratio-matched margin was the physically correct correction and the
    obvious rescue for the absolute score's 0.667. It comes out at exactly
    chance. Pinned so a later change cannot quietly turn this into a success.
    """
    if not RESULTS.exists():
        pytest.skip("no Boltz results recorded yet")
    import pandas as pd

    wide = _wide(pd.read_csv(RESULTS))
    truth = {v["name"]: v["label"] for v in bp.variants()}
    rows = [(float(r["pAzF"] - r["Tyr"]), truth[n])
            for n, r in wide.iterrows() if truth.get(n) is not None]
    value = bp.auc([s for s, _ in rows], [l for _, l in rows])
    assert value == pytest.approx(0.5, abs=1e-9)
    assert value < 2 / 3, "the margin must not be reported as beating pAzF"
    assert bp.exact_p(value, 2, 3) == pytest.approx(0.6, abs=1e-9)
    # and it tops the list with a DIFFERENT rejected variant than pAzF did
    top = max(wide.index, key=lambda n: wide.loc[n, "pAzF"] - wide.loc[n, "Tyr"])
    assert top == "Mc Mut5" and truth[top] == 0


def test_tyr_alone_scores_as_well_as_pazf():
    """❌ Why the 0.667 cannot be read as pAzF specificity.

    Ranking the panel by affinity for the substrate it was selected to REJECT
    is exactly as predictive as ranking it by the one it was selected to
    accept. That equality is the strongest single argument in the module
    docstring, so it is pinned numerically.
    """
    if not RESULTS.exists():
        pytest.skip("no Boltz results recorded yet")
    import pandas as pd

    wide = _wide(pd.read_csv(RESULTS))
    truth = {v["name"]: v["label"] for v in bp.variants()}
    named = [n for n in wide.index if truth.get(n) is not None]
    labels = [truth[n] for n in named]
    pazf = bp.auc([float(wide.loc[n, "pAzF"]) for n in named], labels)
    tyr = bp.auc([float(wide.loc[n, "Tyr"]) for n in named], labels)
    assert pazf == pytest.approx(2 / 3, abs=1e-9)
    assert tyr == pytest.approx(pazf, abs=1e-9)


def test_the_resolved_gap_is_an_order_above_the_selected_one():
    """✅ What the model DOES resolve, quantified -- and why it is the wrong axis.

    Wild-type-vs-mutant is a 0.48 effect; the hit-vs-rejection differences the
    selection actually turned on span 0.03. The docstring's "~17x" rests here.
    """
    if not RESULTS.exists():
        pytest.skip("no Boltz results recorded yet")
    import pandas as pd

    wide = _wide(pd.read_csv(RESULTS))
    delta = wide["pAzF"] - wide["Tyr"]
    mutants = delta.drop(index="Mc wild type")
    boundary = mutants.min() - delta["Mc wild type"]
    spread = mutants.max() - mutants.min()
    assert boundary > 0.45 and spread < 0.04
    assert boundary / spread > 15
