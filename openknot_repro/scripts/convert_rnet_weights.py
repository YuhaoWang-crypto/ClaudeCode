"""Rebuild the official RibonanzaNet checkpoints from the HuggingFace mirrors.

The official weights (`RibonanzaNet.pt`, `RibonanzaNet-SS.pt`) are distributed
through Kaggle, which needs an account. Two mirrors on HuggingFace serve the
same tensors without one:

  roos23/RibonanzaNet             RibonanzaNet.pt, verbatim
  multimolecule/ribonanzanet      the base model, re-hosted under the
                                  MultiMolecule naming scheme
  multimolecule/ribonanzanet-ss   the secondary-structure fine-tune, same scheme

  chaitjo/gRNAde                  the official RibonanzaNet.pt and
                                  RibonanzaNet-SS.pt, re-hosted with gRNAde

This script converts the MultiMolecule secondary-structure checkpoint into the
layout the original `finetuned_RibonanzaNet` expects, and checks the result.

The renaming rule is not guessed, and the conversion is not taken on trust.
The base model exists in both layouts, so the rule is first required to map
`roos23/RibonanzaNet` onto MultiMolecule's tensors bit-identically, all 576 of
them. The same rule is then applied to the SS checkpoint, and the output is
compared against the official `RibonanzaNet-SS.pt` when that file is present:
all 578 tensors match, which is what licenses using either file.

Usage:  python scripts/convert_rnet_weights.py
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import torch

ROOT = Path(__file__).resolve().parent.parent
WEIGHTS = ROOT / "weights"

# MultiMolecule's vocabulary row for each original token id (A, C, G, U, pad).
# Verified by exact tensor comparison of the embedding rows.
VOCAB_ROWS = [6, 7, 8, 9, 0]

_LAYER_RENAMES = [
    (r"^self_attn\.w_qs\.", "attention.self.query."),
    (r"^self_attn\.w_ks\.", "attention.self.key."),
    (r"^self_attn\.w_vs\.", "attention.self.value."),
    (r"^self_attn\.fc\.", "attention.output.dense."),
    (r"^self_attn\.layer_norm\.", "attention.output.layer_norm."),
    (r"^norm1\.", "attention.output.layer_norm2."),
    (r"^norm2\.", "output.layer_norm."),
    (r"^norm3\.", "conv_norm."),
    (r"^linear1\.", "intermediate.dense."),
    (r"^linear2\.", "output.dense."),
    (r"^pairwise2heads\.", "pairwise_bias_proj."),
    (r"^pairwise_norm\.", "pairwise_bias_norm."),
    (r"^conv\.", "conv."),
    (r"^outer_product_mean\.proj_down1\.", "pairwise.triangle_proj.in_proj."),
    (r"^outer_product_mean\.proj_down2\.", "pairwise.triangle_proj.out_proj."),
    (r"^pair_transition\.0\.", "pairwise.intermediate.layer_norm."),
    (r"^pair_transition\.1\.", "pairwise.intermediate.dense."),
    (r"^pair_transition\.3\.", "pairwise.output.dense."),
]
_TRIANGLE_RENAMES = [
    (r"\.norm\.", ".in_norm."),
    (r"\.to_out_norm\.", ".out_norm."),
    (r"\.to_out\.", ".out_proj."),
]


def original_to_multimolecule(key: str) -> str | None:
    """Map an original RibonanzaNet parameter name to its MultiMolecule name.

    Returns None for the three parameters that need reshaping rather than
    renaming (the token embedding and the two reactivity outputs).
    """
    if key in ("encoder.weight", "decoder.weight", "decoder.bias"):
        return None
    if key.startswith("outer_product_mean."):
        tail = key[len("outer_product_mean."):]
        tail = tail.replace("proj_down1.", "in_proj.").replace("proj_down2.", "out_proj.")
        return f"model.encoder.pairwise_embeddings.triangle_proj.{tail}"
    if key.startswith("pos_encoder.linear."):
        tail = key[len("pos_encoder.linear."):]
        return f"model.encoder.pairwise_embeddings.position_embeddings.{tail}"

    match = re.match(r"^transformer_encoder\.(\d+)\.(.*)$", key)
    if not match:
        raise KeyError(f"unrecognised parameter name: {key}")
    index, tail = match.group(1), match.group(2)

    for pattern, replacement in _LAYER_RENAMES:
        if re.match(pattern, tail):
            tail = re.sub(pattern, replacement, tail)
            return f"model.encoder.layer.{index}.{tail}"

    if tail.startswith("triangle_update_"):
        direction = "out" if tail.startswith("triangle_update_out") else "in"
        rest = tail[len(f"triangle_update_{direction}"):]
        for pattern, replacement in _TRIANGLE_RENAMES:
            rest = re.sub(pattern, replacement, rest)
        return f"model.encoder.layer.{index}.pairwise.triangle_mixer_{direction}{rest}"

    raise KeyError(f"unrecognised parameter name: {key}")


def convert(source: dict, reference_keys) -> dict:
    """Build an original-layout state dict from a MultiMolecule one."""
    out = {}
    for key in reference_keys:
        if key == "encoder.weight":
            out[key] = source["model.embeddings.word_embeddings.weight"][VOCAB_ROWS].clone()
        elif key == "decoder.weight":
            out[key] = torch.cat(
                [source["a3c_head.decoder.weight"], source["dms_head.decoder.weight"]], dim=0
            ).clone()
        elif key == "decoder.bias":
            out[key] = torch.cat(
                [source["a3c_head.decoder.bias"], source["dms_head.decoder.bias"]], dim=0
            ).clone()
        else:
            out[key] = source[original_to_multimolecule(key)].clone()
    return out


def verify_rule(original: dict, multimolecule_base: dict) -> int:
    """Require the rename rule to reproduce the official base checkpoint exactly."""
    rebuilt = convert(multimolecule_base, original.keys())
    if set(rebuilt) != set(original):
        raise AssertionError("converted key set differs from the official checkpoint")
    for key, tensor in original.items():
        if not torch.equal(rebuilt[key], tensor):
            delta = (rebuilt[key].float() - tensor.float()).abs().max().item()
            raise AssertionError(f"{key}: converted tensor differs (max |delta| = {delta})")
    return len(original)


def main() -> int:
    base_official = WEIGHTS / "roos23_RibonanzaNet.pt"
    if not base_official.exists():
        base_official = WEIGHTS / "RibonanzaNet.pt"
    base_mm = WEIGHTS / "multimolecule_ribonanzanet.bin"
    ss_mm = WEIGHTS / "multimolecule_ribonanzanet_ss.bin"
    for path in (base_official, base_mm, ss_mm):
        if not path.exists():
            print(f"missing {path}; run scripts/setup_rnet.sh first", file=sys.stderr)
            return 1

    original = torch.load(base_official, map_location="cpu", weights_only=True)
    mm_base = torch.load(base_mm, map_location="cpu", weights_only=True)
    n = verify_rule(original, mm_base)
    print(f"rename rule verified: all {n} base tensors reproduce bit-exactly")

    mm_ss = torch.load(ss_mm, map_location="cpu", weights_only=True)
    ss_keys = list(original.keys()) + ["ct_predictor.weight", "ct_predictor.bias"]
    converted = convert(mm_ss, [k for k in ss_keys if not k.startswith("ct_predictor")])
    converted["ct_predictor.weight"] = mm_ss["ss_head.decoder.weight"].clone()
    converted["ct_predictor.bias"] = mm_ss["ss_head.decoder.bias"].clone()

    # The fine-tune must have moved the backbone; if it were identical to the
    # base model the conversion has silently picked up the wrong file.
    shared = [k for k in original if not k.startswith("decoder")]
    identical = sum(torch.equal(converted[k], original[k]) for k in shared)
    if identical == len(shared):
        raise AssertionError("SS checkpoint backbone identical to the base model")
    print(
        f"secondary-structure checkpoint: {len(converted)} tensors, "
        f"{len(shared) - identical}/{len(shared)} backbone tensors differ from the base model"
    )

    official_ss = WEIGHTS / "RibonanzaNet-SS.pt"
    if official_ss.exists():
        reference = torch.load(official_ss, map_location="cpu", weights_only=True)
        if set(reference) != set(converted):
            raise AssertionError("converted SS key set differs from the official checkpoint")
        mismatched = [k for k in reference if not torch.equal(reference[k], converted[k])]
        if mismatched:
            raise AssertionError(
                f"converted SS checkpoint differs from the official one at {len(mismatched)} "
                f"tensors, e.g. {mismatched[:3]}"
            )
        print(
            f"cross-check passed: conversion reproduces the official {official_ss.name} "
            f"bit-exactly ({len(reference)} tensors)"
        )
        return 0

    torch.save(converted, official_ss)
    print(f"wrote {official_ss}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
