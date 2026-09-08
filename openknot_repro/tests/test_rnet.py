"""Tests for the RibonanzaNet wrapper and the checkpoint conversion.

The tests that need the checkpoints skip when they are absent; run
`scripts/setup_rnet.sh` to fetch them.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

torch = pytest.importorskip("torch")

from openknot.rnet import PAD_TOKEN, WEIGHTS, encode, mask_diagonal  # noqa: E402

WEIGHTS_PRESENT = (WEIGHTS / "RibonanzaNet.pt").exists() and (
    WEIGHTS / "RibonanzaNet-SS.pt"
).exists()
needs_weights = pytest.mark.skipif(
    not WEIGHTS_PRESENT, reason="RibonanzaNet checkpoints not downloaded"
)


def test_encode_maps_bases_and_pads():
    tokens = encode(["ACGU", "AC"])
    assert tokens.shape == (2, 4)
    assert tokens[0].tolist() == [0, 1, 2, 3]
    assert tokens[1].tolist()[2:] == [PAD_TOKEN, PAD_TOKEN]


def test_encode_accepts_dna_spelling():
    assert encode(["ACGT"])[0].tolist() == encode(["ACGU"])[0].tolist()


def test_mask_diagonal_clears_the_band_only():
    matrix = np.ones((5, 5))
    masked = mask_diagonal(matrix, width=2)
    assert masked[0, 0] == 0 and masked[0, 1] == 0
    assert masked[0, 2] == 1
    assert matrix[0, 0] == 1, "input must not be modified in place"


def test_rename_rule_covers_every_parameter_shape():
    from convert_rnet_weights import original_to_multimolecule  # noqa: PLC0415

    cases = {
        "transformer_encoder.3.self_attn.w_qs.weight":
            "model.encoder.layer.3.attention.self.query.weight",
        "transformer_encoder.0.norm3.bias": "model.encoder.layer.0.conv_norm.bias",
        "transformer_encoder.8.triangle_update_in.to_out.weight":
            "model.encoder.layer.8.pairwise.triangle_mixer_in.out_proj.weight",
        "transformer_encoder.1.pair_transition.3.bias":
            "model.encoder.layer.1.pairwise.output.dense.bias",
        "pos_encoder.linear.weight":
            "model.encoder.pairwise_embeddings.position_embeddings.weight",
        "outer_product_mean.proj_down2.bias":
            "model.encoder.pairwise_embeddings.triangle_proj.out_proj.bias",
    }
    for original, expected in cases.items():
        assert original_to_multimolecule(original) == expected


def test_rename_rule_defers_the_reshaped_parameters():
    from convert_rnet_weights import original_to_multimolecule  # noqa: PLC0415

    for key in ("encoder.weight", "decoder.weight", "decoder.bias"):
        assert original_to_multimolecule(key) is None


def test_rename_rule_rejects_an_unknown_name():
    from convert_rnet_weights import original_to_multimolecule  # noqa: PLC0415

    with pytest.raises(KeyError):
        original_to_multimolecule("transformer_encoder.0.not_a_layer.weight")


@needs_weights
def test_reactivity_has_two_channels_and_matching_length():
    from openknot.rnet import RNet  # noqa: PLC0415

    rnet = RNet(load_ss=False)
    sequences = ["GGGGAAAACCCC", "GGGGAAAACCCCAA"]
    profiles = rnet.reactivity(sequences)
    assert [p.shape for p in profiles] == [(12, 2), (14, 2)]
    assert np.isfinite(profiles[0]).all()


@needs_weights
def test_predicted_structure_pairs_a_hairpin():
    from openknot.rnet import RNet  # noqa: PLC0415

    rnet = RNet(load_reactivity=False)
    structure = rnet.structures(["GGGGGCGAAAGCGCCCCC"])[0]
    assert len(structure) == 18
    assert structure.count("(") == structure.count(")") >= 4


@needs_weights
def test_batching_equal_lengths_is_exact():
    """Sequences of the same length batch together without changing their profiles."""
    from openknot.rnet import RNet  # noqa: PLC0415

    rnet = RNet(load_ss=False)
    alone = rnet.reactivity(["GGGGAAAACCCC"])[0]
    batched = rnet.reactivity(["GGGGAAAACCCC", "GGGGUUUUCCCC"])[0]
    assert np.allclose(alone, batched, atol=1e-6)


@needs_weights
def test_batching_very_different_lengths_shifts_the_profile():
    """Padding is not fully masked out, which is why callers must group by length.

    The model convolves across the padded tail, so a short sequence batched with
    a much longer one gets a different profile than it does alone. Batching
    similar lengths (what scripts/rnet_predict.py does, after sorting) leaves
    predictions unchanged; this pins the behaviour that makes that necessary.
    """
    from openknot.rnet import RNet  # noqa: PLC0415

    rnet = RNet(load_ss=False)
    alone = rnet.reactivity(["GGGGAAAACCCC"])[0]
    with_long = rnet.reactivity(["GGGGAAAACCCC", "GGGGAAAACCCC" * 2])[0]
    assert not np.allclose(alone, with_long, atol=1e-4)
