"""Integration bridge: transcriptomic hits x network dynamics x drug reality.

perturbomics answers WHAT to perturb (which drugs / genes reverse a disease
signature). On its own a top connectivity hit is one line of evidence. This
module fuses three MORE, orthogonal, lines — turning a ranked signature into a
prioritised, multi-scale lead table:

  1. TRANSCRIPTOMIC   |WTCS| reversal of the disease signature      (perturbomics)
  2. NETWORK CONTROL  does the target sit at an irreducible-core /   (network-biomarker
                      bistable-switch node of the disease pathway?    M1/M11 core,
                                                                       M2/M19 switch)
  3. ENGAGEABILITY    can a molecule actually hit it? ChEMBL potency  (ChEMBL + Boltz
                      + Boltz binding/ADME -> engagement E            MCP servers;
                                                                       cf. m6_integrate)
  4. CLINICAL         is it already a plausible/repurposable drug?    (ClinicalTrials
                      trial phase / status                            MCP; cf. m8_clinical)

Each axis is a DIFFERENT experiment; a candidate that passes all four is far
more credible than a connectivity outlier. The axes are combined with explicit,
inspectable weights and full missing-data accounting.

Why these compose at all: every layer speaks the same two words — a TARGET
(gene) and a PERTURBATION MAGNITUDE. perturbomics ranks targets; network-
biomarker maps a target to a node and a perturbation mu and predicts the
dynamical/biomarker consequence; the drug servers ground mu in real binding/PK.

RIGOR
  ✅ rigorous : the |WTCS|, the core/switch membership GIVEN a network, the
                ChEMBL/Boltz numbers, the trial phase, and the weighted sum.
  ⚠️ hypothesis: "sits at a control node => better drug target", the E->effect
                 map (illustrative, per m6_integrate), and any therapeutic call.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping, Sequence

import numpy as np
import pandas as pd


# --------------------------------------------------------------------------
# Axis 2 — network control role (consumes network-biomarker outputs)
# --------------------------------------------------------------------------

@dataclass
class NetworkContext:
    """The disease pathway as summarised by the network-biomarker skill.

    Populate from that skill's ``report()`` dicts:
      * ``core_nodes``   = nodes that survive quotienting (M1 symmetry / M11
        fibration) — the irreducible core.
      * ``switch_nodes`` = nodes on a deficiency>0 / bistable motif (M2 CRNT,
        M19 switch library) — able to flip states.
      * ``biomarker``    = the early-warning readout (M4 DNB / leading-eigenvalue
        gene[s]) to MONITOR a responder, not a target to hit.
    Everything is a plain gene-id set, so perturbomics needs no dependency on
    the grn_pipeline package.
    """
    core_nodes: set[str] = field(default_factory=set)
    switch_nodes: set[str] = field(default_factory=set)
    biomarker: tuple[str, ...] = ()

    def role(self, gene: str) -> str:
        g = str(gene).upper()
        core = g in {x.upper() for x in self.core_nodes}
        sw = g in {x.upper() for x in self.switch_nodes}
        if core and sw:
            return "core+switch"
        if core:
            return "core"
        if sw:
            return "switch"
        return "peripheral"

    def bonus(self, gene: str) -> float:
        """Network-control score in [0,1]. ⚠️ the weighting is a hypothesis."""
        return {"core+switch": 1.0, "switch": 0.7, "core": 0.6,
                "peripheral": 0.0}[self.role(gene)]


# --------------------------------------------------------------------------
# Axes 3 & 4 — molecular engageability + clinical status (MCP drug servers)
# --------------------------------------------------------------------------

@dataclass
class DrugEvidence:
    """Per-compound evidence fetched from the drug MCP servers.

    All fields optional — supply what you have; missing axes are excluded from
    that candidate's score and flagged, never silently zeroed.
      pIC50        : from ChEMBL get_bioactivity (functional potency)   ✅ real
      boltz_bind   : Boltz-2.1 binding confidence in [0,1]              ✅ real
      admet_ok     : bool/None gate from Boltz/Inductive Bio ADME        ✅ real
      clinical_phase: 0..4 from ClinicalTrials / Repurposing Hub         ✅ real
      target       : the gene the compound engages (for network mapping)
    """
    name: str
    target: str | None = None
    pIC50: float | None = None
    boltz_bind: float | None = None
    admet_ok: bool | None = None
    clinical_phase: float | None = None

    def engagement(self) -> float | None:
        """Combine potency + structure into E in [0,1]. ✅ inputs real;
        ⚠️ the squashing map is illustrative (mirror of m6_integrate)."""
        parts = []
        if self.pIC50 is not None:
            parts.append(float(np.clip((self.pIC50 - 4.0) / 6.0, 0, 1)))  # pIC50 4..10
        if self.boltz_bind is not None:
            parts.append(float(np.clip(self.boltz_bind, 0, 1)))
        if not parts:
            return None
        e = float(np.mean(parts))
        if self.admet_ok is False:
            e *= 0.5   # ADME liability discounts, does not kill
        return e

    def clinical_score(self) -> float | None:
        if self.clinical_phase is None:
            return None
        return float(np.clip(self.clinical_phase / 4.0, 0, 1))


# --------------------------------------------------------------------------
# The fusion
# --------------------------------------------------------------------------

DEFAULT_WEIGHTS = {"transcriptomic": 0.40, "network": 0.25,
                   "engagement": 0.20, "clinical": 0.15}


def integrated_leads(
    reversers: pd.DataFrame,
    network: NetworkContext | None = None,
    evidence: Mapping[str, DrugEvidence] | None = None,
    weights: Mapping[str, float] = DEFAULT_WEIGHTS,
    target_col: str = "target",
) -> pd.DataFrame:
    """Fuse perturbomics reversers with network + drug evidence into leads.

    Parameters
    ----------
    reversers : DataFrame from ``perturbomics.rank_reversers`` (needs
        ``name``, ``modality``, ``connectivity`` and ideally a target column).
    network   : NetworkContext from the network-biomarker skill (optional).
    evidence  : name -> DrugEvidence from the MCP drug servers (optional).

    Scoring (each axis normalised to [0,1], only present axes counted, weights
    renormalised over the axes that exist for that row):
      transcriptomic = |WTCS|
      network        = NetworkContext.bonus(target)
      engagement     = DrugEvidence.engagement()
      clinical       = DrugEvidence.clinical_score()

    Returns a DataFrame ranked by ``integrated_score`` with every per-axis value
    and a ``provenance`` string listing which axes actually contributed.
    ✅ the arithmetic; ⚠️ the interpretation (see module docstring).
    """
    ev = dict(evidence or {})
    rows = []
    for _, r in reversers.iterrows():
        name = r["name"]
        tgt = (r.get(target_col) or "") if hasattr(r, "get") else ""
        de = ev.get(name)
        if not tgt and de is not None and de.target:
            tgt = de.target

        axes: dict[str, float] = {}
        axes["transcriptomic"] = float(min(abs(r["connectivity"]), 1.0))
        if network is not None and tgt:
            axes["network"] = network.bonus(tgt)
        if de is not None:
            e = de.engagement()
            if e is not None:
                axes["engagement"] = e
            c = de.clinical_score()
            if c is not None:
                axes["clinical"] = c

        # renormalise weights over the axes present for THIS row
        wsum = sum(weights[a] for a in axes) or 1.0
        score = sum(weights[a] * v for a, v in axes.items()) / wsum

        rows.append({
            "name": name,
            "modality": r.get("modality", "") if hasattr(r, "get") else "",
            "target": tgt,
            "net_role": network.role(tgt) if (network and tgt) else "",
            "integrated_score": round(score, 4),
            "transcriptomic": round(axes.get("transcriptomic", np.nan), 3),
            "network": round(axes["network"], 3) if "network" in axes else np.nan,
            "engagement": round(axes["engagement"], 3) if "engagement" in axes else np.nan,
            "clinical": round(axes["clinical"], 3) if "clinical" in axes else np.nan,
            "provenance": "+".join(axes),
        })
    out = pd.DataFrame(rows).sort_values(
        "integrated_score", ascending=False).reset_index(drop=True)
    return out


def monitoring_biomarker(network: NetworkContext) -> tuple[str, ...]:
    """The early-warning readout to TRACK a responder (from M4 DNB / leading
    eigenvalue). A gene to measure over time, not a target to hit — closes the
    loop from 'what to give' to 'how to tell it's working'. ✅ if the network's
    biomarker was computed rigorously; ⚠️ the toxicity/response mapping is not."""
    return network.biomarker
