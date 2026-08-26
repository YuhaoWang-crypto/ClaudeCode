"""Unit tests. Run with ``python tests/test_vetimmuno.py`` or pytest.

These are fast, offline-once-cached checks on the logic. The scientific
known-answer tests live in ``vetimmuno.validate`` and run as part of every
workflow execution -- they belong in the report, not only in CI.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np

from vetimmuno import epitope, groove, predict, validate
from vetimmuno.insulin import apply_edits, diff, natural_insulin


def test_mature_chain_lengths():
    for species in ("human", "dog", "cat", "bovine", "pig"):
        ins = natural_insulin(species)
        assert len(ins.A) == 21, (species, len(ins.A))
        assert len(ins.B) == 30, (species, len(ins.B))


def test_edit_grammar():
    human = natural_insulin("human")
    glargine = apply_edits(human, ["A21N>G", "B31*>RR"], "glargine")
    assert glargine.A.endswith("G") and len(glargine.A) == 21
    assert glargine.B.endswith("TRR") and len(glargine.B) == 32

    des_b30 = apply_edits(human, ["B30T>-"], "des-B30")
    assert len(des_b30.B) == 29

    for bad in ["A21Q>G", "A99N>G", "nonsense"]:
        try:
            apply_edits(human, [bad], "bad")
        except ValueError:
            continue
        raise AssertionError(f"{bad} should have been rejected")


def test_diff_is_symmetric_in_count():
    a, b = natural_insulin("cat"), natural_insulin("human")
    assert len(diff(a, b)) == len(diff(b, a)) == 4


def test_neo_cores_are_absent_from_self():
    own = natural_insulin("dog")
    drug = natural_insulin("bovine")
    tolerated = epitope.self_core_set(own)
    for nc in epitope.neo_cores(drug, own):
        assert nc.core.sequence not in tolerated


def test_neo_cores_empty_for_identical_donor():
    assert epitope.neo_cores(natural_insulin("pig"), natural_insulin("dog")) == []


def test_tiling_covers_every_residue():
    ins = natural_insulin("cat")
    for chain in ("A", "B"):
        seq = ins.chain(chain)
        covered = set()
        for pep in epitope.tile(chain, seq):
            covered.update(range(pep.start, pep.start + len(pep.sequence)))
        assert covered == set(range(1, len(seq) + 1))


def test_pseudosequence_length_matches_contact_set():
    seq = groove.reference_domain("DRB")
    mol = groove.build_molecule("ref", "human", "DRB", seq)
    assert len(mol.pseudoseq) == len(groove.BETA_CONTACT)
    assert mol.coverage == 1.0
    assert mol.identity_to_ref == 1.0


def test_pseudo_identity_bounds():
    a = "ABCDEFGHIJ"
    assert groove.pseudo_identity(a, a) == 1.0
    assert groove.pseudo_identity(a, "-" * len(a)) == 0.0


def test_background_rank_is_uniform():
    backend = predict.IllustrativeScorer()
    seq = groove.reference_domain("DRB")
    mol = groove.build_molecule("ref", "human", "DRB", seq)
    ranker = predict.BackgroundRank(backend, n_background=8000, seed=1)
    check = validate.check_rank_calibration(ranker, mol, tolerance=1.0)
    assert check.status == "PASS", check.observed


def test_surrogate_prefers_hydrophobic_p1_on_dr1():
    backend = predict.IllustrativeScorer()
    mol = groove.build_molecule("HLA-DRB1*01:01", "human", "DRB",
                                groove.reference_domain("DRB"))
    assert validate.check_surrogate_direction(backend, mol).status == "PASS"


def test_scoring_is_deterministic():
    backend = predict.IllustrativeScorer()
    mol = groove.build_molecule("ref", "human", "DRB", groove.reference_domain("DRB"))
    cores = [c.sequence for c in epitope.all_cores(natural_insulin("human"))]
    assert np.array_equal(backend.score(cores, mol), backend.score(cores, mol))


def test_spearman_matches_known_values():
    assert abs(validate.spearman([1, 2, 3, 4], [1, 2, 3, 4]) - 1.0) < 1e-9
    assert abs(validate.spearman([1, 2, 3, 4], [4, 3, 2, 1]) + 1.0) < 1e-9


def test_netmhciipan_parser():
    text = """
# NetMHCIIpan version 4.3
 Pos MHC Peptide Core Score_EL %Rank_EL
   1 CUSTOM AAAAAAAAA AAAAAAAAA 0.1234 12.00
   2 CUSTOM CCCCCCCCC CCCCCCCCC 0.5678 3.00
"""
    scores = predict.parse_netmhciipan(text, expected=2)
    assert list(scores) == [0.1234, 0.5678]


def main() -> int:
    failures = 0
    for name, fn in sorted(globals().items()):
        if not name.startswith("test_") or not callable(fn):
            continue
        try:
            fn()
            print(f"PASS {name}")
        except Exception as exc:  # noqa: BLE001 - test runner
            failures += 1
            print(f"FAIL {name}: {exc}")
    print(f"\n{failures} failure(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
