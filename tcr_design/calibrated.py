"""The calibrated ddG model: load it, and use it to predict.

The model shipped in ``data/calibrated_ddg_model.json`` is a ridge regression
over measured structural descriptors, fitted to pooled ATLAS and SKEMPI 2.0
TCR-pMHC mutation data. Regenerate it with::

    python -m tcr_design benchmark --save-model

Read the performance block in that file before using a prediction for
anything. The headline numbers, validated leave-one-structure-out so the model
never saw the receptor it was scored on, are modest: it explains part of the
variance in ddG and ranks stabilising mutations better than chance, and it is
nowhere near a replacement for flex ddG, MM/PBSA or a binding measurement.

What it is good for: ordering a few hundred candidate mutations so that the
expensive methods are spent on the right dozen.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

DEFAULT_MODEL_PATH = Path(__file__).parent / "data" / "calibrated_ddg_model.json"


@dataclass
class DdgModel:
    """A fitted linear model over the extended structural features."""

    coefficients: dict[str, float]
    intercept: float
    feature_set: str
    performance: dict
    provenance: dict

    @property
    def held_out_spearman(self) -> float | None:
        return self.performance.get("held_out", {}).get("spearman")

    @property
    def held_out_rmse(self) -> float | None:
        return self.performance.get("held_out", {}).get("rmse")

    def summary(self) -> str:
        held = self.performance.get("held_out", {})
        spearman = held.get("spearman")
        rmse = held.get("rmse")
        auc = held.get("auc_find_stabilising")
        return (
            f"calibrated ddG model: held-out Spearman "
            f"{'n/a' if spearman is None else f'{spearman:+.3f}'}, "
            f"RMSE {'n/a' if rmse is None else f'{rmse:.2f}'} kcal/mol, "
            f"AUC for finding stabilising mutations "
            f"{'n/a' if auc is None else f'{auc:.3f}'} "
            f"(n={held.get('n', '?')} mutations, "
            f"{held.get('structures_held_out', '?')} structures)"
        )

    def predict(self, residue, mutation) -> float:
        """Predicted ddG in kcal/mol; positive means weaker binding."""
        from .benchmark import BenchmarkEntry, extended_feature_vector

        entry = BenchmarkEntry(
            pdb_id="", chain=residue.chain, number=residue.number,
            insertion=residue.insertion, wild_type=residue.aa,
            mutant=mutation.mutant, ddg=0.0, side="tcr", location="",
            source="runtime",
        )
        entry.residue = residue
        entry.scored = mutation
        vector = extended_feature_vector(entry)
        return self.intercept + sum(
            coefficient * value
            for coefficient, value in zip(self.coefficients.values(), vector)
        )


def load_model(path: Path | None = None) -> DdgModel | None:
    """Load the shipped model, or return None when it has not been generated."""
    path = path or DEFAULT_MODEL_PATH
    if not path.exists():
        return None
    payload = json.loads(path.read_text())
    return DdgModel(
        coefficients=payload["coefficients"],
        intercept=payload["intercept"],
        feature_set=payload.get("feature_set", "extended"),
        performance=payload.get("performance", {}),
        provenance=payload.get("provenance", {}),
    )


def save_model(model_payload: dict, path: Path | None = None) -> Path:
    path = path or DEFAULT_MODEL_PATH
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(model_payload, indent=2, default=str) + "\n")
    return path
