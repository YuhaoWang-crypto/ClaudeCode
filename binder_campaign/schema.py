"""Frozen campaign vocabularies and the design-sheet schema.

Two files the orchestrator writes ONCE, at kickoff / roster validation, and
freezes:

* ``/state/method_vocab.json`` — the ``structure_method`` and ``seq_method``
  enums (they are **not** free text), the tool-to-family map, and the
  ``alias_map`` block used to canonicalise provenance keys.  A row whose
  canonical token is absent from the enum is rejected.
* ``sheet_schema.json`` — exact column names, dtypes, enums, units and the
  ``mandatory_nonnull`` set.  Per-target sheet writers validate every row
  against it, and the campaign merge script refuses sheets with renamed,
  retyped, missing or duplicated columns.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Any, Mapping

from .scoring import DEFAULT_ARMS

__all__ = [
    "STRUCTURE_METHODS",
    "ASTERISKED_METHODS",
    "SEQ_METHODS",
    "MethodVocab",
    "SheetSchema",
    "default_method_vocab",
    "default_sheet_schema",
]

#: Structure-design / co-design methods named in the prompt roster.  The
#: asterisked ones each owe every target at least 50 backbones.
ASTERISKED_METHODS: tuple[str, ...] = (
    "rfdiffusion",
    "rfdiffusion3",
    "freebindcraft",
    "boltzgen",
    "pxdesign",
    "proteina_complexa",
    "genie3",
)

STRUCTURE_METHODS: tuple[str, ...] = ASTERISKED_METHODS + (
    "mosaic",
    "foldcraft",
    "boltzdesign1",
    "halludesign",
    "protein_hunter",
    "partial_diffusion",  # in silico optimization arm, keeps its parent's root id
)

SEQ_METHODS: tuple[str, ...] = (
    "solublempnn",
    "solublecaliby",
    "solublecaliby_ensemble",  # 32-structure Protpardelle-1c ensemble mode
    "proteinmpnn",             # base variants: backbone search only
    "caliby",
    "codesign_native",         # a co-design model's own sequence
    "mutagenesis",
)

_TOOL_FAMILY = {
    "rfdiffusion": "rfdiffusion",
    "rfdiffusion3": "rfdiffusion",
    "freebindcraft": "bindcraft",
    "boltzgen": "boltz",
    "boltzdesign1": "boltz",
    "pxdesign": "protenix",
    "protein_hunter": "hallucination",
    "halludesign": "hallucination",
    "proteina_complexa": "proteina",
    "genie3": "genie",
    "mosaic": "mosaic",
    "foldcraft": "foldcraft",
    "partial_diffusion": "rfdiffusion",
}

_ALIASES = {
    "RFdiffusion": "rfdiffusion",
    "rf_diffusion": "rfdiffusion",
    "RFdiffusion3": "rfdiffusion3",
    "rfd3": "rfdiffusion3",
    "FreeBindCraft": "freebindcraft",
    "bindcraft": "freebindcraft",
    "BoltzGen": "boltzgen",
    "PXDesign": "pxdesign",
    "Proteina-Complexa": "proteina_complexa",
    "Genie3": "genie3",
    "FoldCraft": "foldcraft",
    "BoltzDesign1": "boltzdesign1",
    "HalluDesign": "halludesign",
    "ProteinHunter": "protein_hunter",
    "SolubleMPNN": "solublempnn",
    "SolubleCaliby": "solublecaliby",
    "ProteinMPNN": "proteinmpnn",
    "Caliby": "caliby",
}


@dataclass(frozen=True)
class MethodVocab:
    """``/state/method_vocab.json``, frozen once at roster validation."""

    structure_methods: tuple[str, ...] = STRUCTURE_METHODS
    asterisked: tuple[str, ...] = ASTERISKED_METHODS
    seq_methods: tuple[str, ...] = SEQ_METHODS
    family: Mapping[str, str] = field(default_factory=lambda: dict(_TOOL_FAMILY))
    alias_map: Mapping[str, str] = field(default_factory=lambda: dict(_ALIASES))

    def canonical_structure_method(self, token: str) -> str:
        """Canonicalise via the alias map; raise if absent from the enum."""
        t = self.alias_map.get(token, token)
        t = self.alias_map.get(t.strip(), t.strip())
        if t not in self.structure_methods:
            raise ValueError(
                f"structure_method {token!r} -> {t!r} is not in the frozen enum"
            )
        return t

    def canonical_seq_method(self, token: str) -> str:
        t = self.alias_map.get(token, token)
        t = self.alias_map.get(t.strip(), t.strip())
        if t not in self.seq_methods:
            raise ValueError(
                f"seq_method {token!r} -> {t!r} is not in the frozen enum"
            )
        return t

    def to_json(self) -> dict:
        return {
            "structure_methods": list(self.structure_methods),
            "asterisked": list(self.asterisked),
            "seq_methods": list(self.seq_methods),
            "family": dict(self.family),
            "alias_map": dict(self.alias_map),
        }


def default_method_vocab() -> MethodVocab:
    return MethodVocab()


# --------------------------------------------------------------------------- #
# sheet schema
# --------------------------------------------------------------------------- #

#: Columns required "at minimum" by the prompt's Design sheet schema item.
_BASE_COLUMNS: dict[str, str] = {
    "design_id": "str",
    "target": "str",
    "sequence": "str",
    "binder_len": "int",
    "rank": "int",
    "rank_zscore": "float",
    "final_score": "float",
    "score_instrument": "str",
    "pose_PASS": "bool",
    "pose_dockq": "float",
    "structure_method": "enum:structure_methods",
    "seq_method": "enum:seq_methods",
    "opt_round": "int",
    "root_backbone_id": "str",
    "parent_design_id": "str",
    "n_seeds": "int",
    "novelty_verdict_path": "str",
    "tm90_cluster_id": "str",
    "fold_class": "enum:fold_class",
    "designed_structure_path": "str",
    # shadow / feasibility metrics
    "monomer_plddt": "float",
    "lcp_score": "float",
    "esmc_ll": "float",
    # disclosure
    "relaxation_step": "str",
    "construct_status": "str",
}

#: Non-null on every ranked row.
_MANDATORY_NONNULL: tuple[str, ...] = (
    "design_id",
    "target",
    "sequence",
    "structure_method",
    "root_backbone_id",
    "rank",
    "rank_zscore",
    "final_score",
    "score_instrument",
    "pose_PASS",
    "pose_dockq",
    "designed_structure_path",
)

#: Required non-null for a row to be *rank-eligible* under the diversity gates.
_DIVERSITY_REQUIRED: tuple[str, ...] = (
    "root_backbone_id",
    "structure_method",
    "tm90_cluster_id",
    "fold_class",
)


@dataclass(frozen=True)
class SheetSchema:
    """``sheet_schema.json``, frozen at campaign kickoff."""

    columns: dict[str, str]
    mandatory_nonnull: tuple[str, ...]
    diversity_required: tuple[str, ...]
    enums: dict[str, tuple[str, ...]]
    arms: tuple[str, ...]
    counter_screened_targets: tuple[str, ...] = ()

    def validate_row(self, row: Mapping[str, Any], ranked: bool = True) -> None:
        """Raise on a renamed, retyped, missing or duplicated column.

        A null in ``mandatory_nonnull`` is fail-loud.  A null in a
        required-but-non-load-bearing column (``seq_method``, a
        ``predicted_structure_path_{arm}``) is *not* fatal here: close-out
        downgrades that row to disclosed-with-deviation instead of dropping it.
        """
        unknown = set(row) - set(self.columns)
        if unknown:
            raise ValueError(f"unknown sheet columns: {sorted(unknown)}")
        missing = set(self.columns) - set(row)
        if missing:
            raise ValueError(f"missing sheet columns: {sorted(missing)}")

        if ranked:
            for col in self.mandatory_nonnull:
                v = row.get(col)
                if v is None or (isinstance(v, float) and v != v) or v == "PENDING":
                    raise ValueError(
                        f"mandatory_nonnull column {col!r} is null on ranked row "
                        f"{row.get('design_id')!r}"
                    )
        for col, dtype in self.columns.items():
            if dtype.startswith("enum:"):
                allowed = self.enums.get(dtype.split(":", 1)[1])
                v = row.get(col)
                if allowed and v not in (None, "") and v not in allowed:
                    raise ValueError(
                        f"column {col!r} value {v!r} not in frozen enum"
                    )

    def to_json(self) -> dict:
        return {
            "columns": dict(self.columns),
            "mandatory_nonnull": list(self.mandatory_nonnull),
            "diversity_required": list(self.diversity_required),
            "enums": {k: list(v) for k, v in self.enums.items()},
            "arms": list(self.arms),
            "counter_screened_targets": list(self.counter_screened_targets),
        }

    def dumps(self) -> str:
        return json.dumps(self.to_json(), indent=2, sort_keys=True)


def default_sheet_schema(
    arms: tuple[str, ...] = DEFAULT_ARMS,
    vocab: MethodVocab | None = None,
    counter_screened_targets: tuple[str, ...] = ("GDF-8",),
) -> SheetSchema:
    """The schema the orchestrator freezes at kickoff.

    One ``predicted_structure_path_{arm}`` column per ranking arm, and one
    column per scoring arm and shadow metric actually computed in the campaign.
    """
    vocab = vocab or default_method_vocab()
    columns = dict(_BASE_COLUMNS)
    for arm in arms:
        columns[f"ipsae_{arm}"] = "float"
        columns[f"sc_DockQ_{arm}"] = "float"
        columns[f"n_seeds_{arm}"] = "int"
        columns[f"predicted_structure_path_{arm}"] = "str"
        columns[f"selectivity_delta_{arm}"] = "float"
        columns[f"ipsae_offtarget_{arm}"] = "float"
        # oligomeric targets: full-occupancy companion columns
        columns[f"ipsae_NN_{arm}"] = "float"
        columns[f"sc_DockQ_NN_{arm}"] = "float"
    columns["binder_binder_clashes_NN"] = "int"

    return SheetSchema(
        columns=columns,
        mandatory_nonnull=_MANDATORY_NONNULL,
        diversity_required=_DIVERSITY_REQUIRED,
        enums={
            "structure_methods": vocab.structure_methods,
            "seq_methods": vocab.seq_methods,
            "fold_class": ("all_alpha", "not_all_alpha"),
        },
        arms=arms,
        counter_screened_targets=counter_screened_targets,
    )
