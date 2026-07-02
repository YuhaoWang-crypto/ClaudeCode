import pytest

from alphagenome_evo2.types import Sequence, validate_dna


def test_validate_dna_uppercases_and_strips():
    assert validate_dna("  acgt ") == "ACGT"


def test_validate_dna_rejects_non_acgt():
    with pytest.raises(ValueError):
        validate_dna("ACGTX")


def test_sequence_normalises_and_measures_length():
    s = Sequence("s1", "acgtACGT")
    assert s.seq == "ACGTACGT"
    assert len(s) == 8


def test_sequence_with_metadata_is_immutable_copy():
    s = Sequence("s1", "ACGT", metadata={"a": 1})
    s2 = s.with_metadata(b=2)
    assert s.metadata == {"a": 1}
    assert s2.metadata == {"a": 1, "b": 2}
