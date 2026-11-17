r"""Turn parsed structures into IPA training examples.

One example = a CDR design block (CDR residues + a small framework flank, with a
generate_mask marking the CDR positions to be noised/predicted) and an antigen
context block (the K antigen residues nearest the CDRs). Coordinates are centred
on the antigen centroid. Collation pads a minibatch to its max lengths and emits
validity masks so padded positions are excluded from attention and loss.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..data.structures import ComplexStructure
from ..model.ctmc import AA_INDEX

PAD_TOKEN = 0            # padded positions carry a dummy token (masked out anyway)


def _tok(aa: str) -> int:
    return AA_INDEX.get(aa, 0)


def example_from_structure(cx: ComplexStructure, flank: int = 3, ka: int = 32,
                           min_antigen: int = 3, coord_scale: float = 10.0) -> dict | None:
    """Build one numpy example, or None if the complex lacks usable CDRs/antigen.

    Coordinates are centred on the antigen and divided by ``coord_scale`` (~nm)
    so the flow-matching target x1-x0 is O(1), which conditions training well.
    """
    design_ca, design_tok, gen = [], [], []
    for chain in (cx.heavy, cx.light):
        if chain is None or not chain.cdr_spans:
            continue
        cdr_idx = set()
        for lo, hi in chain.cdr_spans.values():
            cdr_idx.update(range(lo, hi))
        if not cdr_idx:
            continue
        crop = sorted({j for i in cdr_idx for j in range(i - flank, i + flank + 1)
                       if 0 <= j < len(chain.seq)})
        for j in crop:
            design_ca.append(chain.ca[j])
            design_tok.append(_tok(chain.seq[j]))
            gen.append(j in cdr_idx)
    if sum(gen) < 3:
        return None
    design_ca = np.asarray(design_ca, float)
    gen = np.asarray(gen, bool)

    # Antigen context: residues nearest to the CDR CAs.
    ag_ca = np.concatenate([c.ca for c in cx.antigen], axis=0)
    ag_seq = "".join(c.seq for c in cx.antigen)
    if len(ag_ca) < min_antigen:
        return None
    cdr_ca = design_ca[gen]
    dmin = np.linalg.norm(ag_ca[:, None, :] - cdr_ca[None, :, :], axis=-1).min(axis=1)
    sel = np.argsort(dmin)[:ka]
    sel = np.sort(sel)                        # keep sequence order
    antigen_ca = ag_ca[sel]
    antigen_tok = np.array([_tok(ag_seq[i]) for i in sel], dtype=np.int64)

    center = antigen_ca.mean(axis=0)
    return {
        "design_coords": ((design_ca - center) / coord_scale).astype("float32"),
        "design_tokens": np.asarray(design_tok, dtype=np.int64),
        "generate_mask": gen,
        "antigen_coords": ((antigen_ca - center) / coord_scale).astype("float32"),
        "antigen_tokens": antigen_tok,
        "pdb_id": cx.pdb_id,
    }


@dataclass
class StructureDataset:
    examples: list[dict]

    @classmethod
    def from_structures(cls, structures: list[ComplexStructure], **kw) -> "StructureDataset":
        ex = [e for e in (example_from_structure(c, **kw) for c in structures) if e]
        return cls(ex)

    def __len__(self):
        return len(self.examples)

    def __getitem__(self, i):
        return self.examples[i]


def collate(batch: list[dict]):
    """Pad a list of examples to the batch max lengths -> torch tensors + masks."""
    import torch

    B = len(batch)
    Ld = max(3, max(len(e["design_tokens"]) for e in batch))
    Ka = max(3, max(len(e["antigen_tokens"]) for e in batch))

    dc = torch.zeros(B, Ld, 3)
    dt = torch.zeros(B, Ld, dtype=torch.long)
    gm = torch.zeros(B, Ld, dtype=torch.bool)
    dp = torch.zeros(B, Ld, dtype=torch.bool)
    ac = torch.zeros(B, Ka, 3)
    at = torch.zeros(B, Ka, dtype=torch.long)
    ap = torch.zeros(B, Ka, dtype=torch.bool)

    for b, e in enumerate(batch):
        ld, ka = len(e["design_tokens"]), len(e["antigen_tokens"])
        dc[b, :ld] = torch.from_numpy(e["design_coords"])
        dt[b, :ld] = torch.from_numpy(e["design_tokens"])
        gm[b, :ld] = torch.from_numpy(e["generate_mask"])
        dp[b, :ld] = True
        ac[b, :ka] = torch.from_numpy(e["antigen_coords"])
        at[b, :ka] = torch.from_numpy(e["antigen_tokens"])
        ap[b, :ka] = True
        # Spread padded CA apart so per-block frames stay finite (not all-coincident).
        if ld < Ld:
            dc[b, ld:] = torch.arange(1, Ld - ld + 1).float()[:, None] * 50.0
        if ka < Ka:
            ac[b, ka:] = torch.arange(1, Ka - ka + 1).float()[:, None] * 50.0

    return {
        "design_coords": dc, "design_tokens": dt, "generate_mask": gm, "design_pad": dp,
        "antigen_coords": ac, "antigen_tokens": at, "antigen_pad": ap,
    }
