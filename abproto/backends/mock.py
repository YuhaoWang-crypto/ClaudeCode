"""Deterministic mock backend so the whole pipeline runs with no GPU/weights.

It is NOT a model. It fabricates plausible-shaped outputs with a built-in
"truth" signal so that Beam Search, self-consistency and clustering have
something real to optimise and we can unit-test the orchestration.

How the fake signal works
-------------------------
Each design gets a hidden ``quality`` in [0, 1] drawn from its own seed. A
higher-quality design is given (a) higher plDDT/iPTM and lower PAE/PDE, and
(b) a refold that lands closer to its design backbone (lower recovery RMSD).
``refine`` nudges children toward higher quality than their parent, so Beam
Search measurably improves mean quality over steps — exactly the behaviour the
real system shows, just synthetic.
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from ..types import ANTIGEN, CDRS, HEAVY, LIGHT, Confidence, Design, DesignSpec, Pose

_AA = "ACDEFGHIKLMNPQRSTVWY"


def _rng(*parts) -> np.random.Generator:
    """A reproducible RNG keyed by arbitrary hashable parts."""
    h = abs(hash(parts)) % (2**32)
    return np.random.default_rng(h)


def _ca_trace(n: int, rng: np.random.Generator, origin=(0.0, 0.0, 0.0)) -> np.ndarray:
    """A self-avoiding-ish random walk at ~3.8 A C-alpha spacing."""
    step = 3.8
    coords = np.zeros((n, 3))
    coords[0] = origin
    direction = rng.normal(size=3)
    direction /= np.linalg.norm(direction) + 1e-9
    for i in range(1, n):
        direction = direction + 0.5 * rng.normal(size=3)
        direction /= np.linalg.norm(direction) + 1e-9
        coords[i] = coords[i - 1] + step * direction
    return coords


def _random_cdr(length: int, rng: np.random.Generator) -> str:
    return "".join(rng.choice(list(_AA), size=length))


class MockGenerator:
    """Fake Generator: grafts random CDRs onto the germline framework and
    attaches a random C-alpha backbone + quality-correlated confidence."""

    def __init__(self, base_quality: float = 0.55):
        self.base_quality = base_quality

    # -- helpers ---------------------------------------------------------
    def _build(self, spec: DesignSpec, idx: int, quality: float,
               parent_id: Optional[str], provenance: str) -> Design:
        rng = _rng(spec.target_name, idx, provenance, round(quality, 4))

        # Assemble sequences: framework with CDR inserts of the sampled lengths.
        heavy, hspans = self._graft(spec.heavy_framework or "EVQLVESGGGLVQPGG", spec,
                                    ("H1", "H2", "H3"), rng)
        light, lspans = self._graft(spec.light_framework or "DIQMTQSPSSLSASVG", spec,
                                    ("L1", "L2", "L3"), rng)
        spans = {**hspans, **lspans}

        backbone = {
            HEAVY: _ca_trace(len(heavy), rng, origin=(0, 0, 0)),
            LIGHT: _ca_trace(len(light), rng, origin=(20, 0, 0)),
        }
        if spec.antigen_coords is not None:
            backbone[ANTIGEN] = spec.antigen_coords
        else:
            backbone[ANTIGEN] = _ca_trace(len(spec.antigen_seq), rng, origin=(-20, 5, 0))

        conf = self._confidence(quality, rng)
        return Design(
            design_id=f"{spec.target_name}_{provenance}_{idx:05d}",
            heavy_seq=heavy, light_seq=light, backbone=backbone,
            cdr_spans=spans, confidence=conf, parent_id=parent_id, provenance=provenance,
        )

    def _graft(self, framework: str, spec: DesignSpec, cdrs, rng):
        """Insert CDRs of the sampled lengths at three points in the framework."""
        seq = framework
        spans: dict[str, tuple[int, int]] = {}
        # Three insertion anchors spread across the framework.
        anchors = [len(seq) // 4, len(seq) // 2, (3 * len(seq)) // 4]
        offset = 0
        for cdr, anchor in zip(cdrs, anchors):
            length = spec.cdr_lengths.get(cdr, {"H1": 5, "H2": 6, "H3": 11,
                                                "L1": 6, "L2": 3, "L3": 9}[cdr])
            loop = _random_cdr(length, rng)
            pos = anchor + offset
            seq = seq[:pos] + loop + seq[pos:]
            spans[cdr] = (pos, pos + length)
            offset += length
        return seq, spans

    def _confidence(self, quality: float, rng: np.random.Generator) -> Confidence:
        q = float(np.clip(quality, 0.0, 1.0))
        plddt = 100.0 * (0.55 + 0.4 * q) + rng.normal(scale=2.0)
        iptm = np.clip(0.3 + 0.6 * q + rng.normal(scale=0.03), 0, 1)
        pae = np.clip(12.0 * (1 - q) + 2.0 + rng.normal(scale=0.5), 0.5, None)
        pde = np.clip(8.0 * (1 - q) + 1.5 + rng.normal(scale=0.4), 0.3, None)
        # Interface heads: more quality-sensitive than the global ones (the
        # interface is where binding is decided), so scoring on them separates
        # good binders better.
        i_pae = np.clip(15.0 * (1 - q) + 1.5 + rng.normal(scale=0.6), 0.3, None)
        i_plddt = np.clip(100.0 * (0.45 + 0.5 * q) + rng.normal(scale=2.5), 0, 100)
        contact = np.clip(0.2 + 0.75 * q + rng.normal(scale=0.03), 0, 1)
        return Confidence(plddt=float(plddt), pae=float(pae), pde=float(pde), iptm=float(iptm),
                          interface_pae=float(i_pae), interface_plddt=float(i_plddt),
                          contact_prob=float(contact))

    # -- Generator protocol ---------------------------------------------
    def generate(self, spec: DesignSpec, n: int, seed: Optional[int] = None) -> list[Design]:
        rng = _rng(spec.target_name, "generate", seed)
        out = []
        for i in range(n):
            quality = float(np.clip(rng.beta(2, 2) * 0.6 + self.base_quality * 0.4, 0, 1))
            out.append(self._build(spec, i if seed is None else seed * 1000 + i,
                                   quality, None, "generator"))
        return out

    def refine(self, parents: list[Design], spec: DesignSpec, n_per_parent: int,
               seed: Optional[int] = None) -> list[Design]:
        out = []
        for p in parents:
            # Parent quality recovered from its confidence (iptm is monotone in q).
            pq = (p.confidence.iptm - 0.3) / 0.6 if p.confidence else 0.5
            rng = _rng(p.design_id, "refine", seed)
            for j in range(n_per_parent):
                # Children drift upward in quality (local search improves).
                child_q = float(np.clip(pq + rng.normal(loc=0.05, scale=0.08), 0, 1))
                child = self._build(spec, abs(hash((p.design_id, j))) % 10**6,
                                    child_q, p.design_id, "generator")
                out.append(child)
        return out


class MockFolder:
    """Fake Folder: returns a refolded backbone whose deviation from the design
    shrinks with design quality, so self-consistency is a meaningful filter."""

    def fold(self, design: Design, spec: DesignSpec, seed: Optional[int] = None) -> Pose:
        q = (design.confidence.iptm - 0.3) / 0.6 if design.confidence else 0.5
        q = float(np.clip(q, 0, 1))
        rng = _rng(design.design_id, "fold", seed)

        # Perturb the design backbone; amplitude ~ (1 - quality).
        noise_amp = 0.4 + 4.5 * (1 - q)
        backbone = {}
        for chain, coords in design.backbone.items():
            if chain == ANTIGEN:
                backbone[chain] = coords  # antigen held fixed
            else:
                backbone[chain] = coords + rng.normal(scale=noise_amp, size=coords.shape)

        conf = MockGenerator()._confidence(q, rng)
        return Pose(backbone=backbone, confidence=conf)
