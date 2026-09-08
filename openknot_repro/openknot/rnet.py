"""RibonanzaNet (RNet) inference: SHAPE reactivity and secondary structure.

RNet is the model the paper leans on: it replaces 3D structure prediction as
the signal that guides and filters designs. Two checkpoints are used:

  RibonanzaNet.pt      chemical-mapping model, outputs (2A3, DMS) reactivity
  RibonanzaNet-SS.pt   secondary-structure fine-tune, outputs pair logits

The model classes come from the authors' code (Shujun-He/Struct2SeQ,
`Network_test10.py`, MIT), cloned to `third_party/` at a pinned commit by
`scripts/setup_rnet.sh` rather than vendored here. Base-pair decoding uses
arnie's Hungarian matcher, the same one the authors call.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parent.parent
THIRD_PARTY = ROOT / "third_party"
WEIGHTS = ROOT / "weights"

TOKENS = {nt: i for i, nt in enumerate("ACGU")}
PAD_TOKEN = 4


def _ensure_arnie_config() -> None:
    """arnie refuses to import without a config file naming at least one package.

    Only its pure-python Hungarian decoder is used here, so a stub is enough;
    this is what the authors' `make_dummy_arnie.py` writes.
    """
    if os.environ.get("ARNIEFILE") and Path(os.environ["ARNIEFILE"]).exists():
        return
    path = THIRD_PARTY / "arnie_file.txt"
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        path.write_text("linearpartition: .\nTMP: /tmp\n")
    os.environ["ARNIEFILE"] = str(path)


_ensure_arnie_config()

# Config values from the authors' configs/pairwise.yaml. They must match the
# checkpoint: k=5 is the convolution kernel the released weights were trained
# with (Struct2SeQ's own test10_configs/pairwise.yaml says k=9, which does not
# fit these tensors).
MODEL_CONFIG = dict(
    ninp=256,
    nlayers=9,
    nclass=2,
    ntoken=5,
    nhead=8,
    k=5,
    dropout=0.05,
    use_triangular_attention=False,
    pairwise_dimension=64,
    use_grad_checkpoint=False,
)


class _Config:
    def __init__(self, **entries):
        self.__dict__.update(entries)


def _import_network_module():
    """Import the authors' model definitions from the pinned clone."""
    path = THIRD_PARTY / "Struct2SeQ"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found. Run scripts/setup_rnet.sh to fetch the model code."
        )
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))
    import Network_test10  # noqa: PLC0415

    return Network_test10


def encode(sequences: list[str], device: str = "cpu") -> torch.Tensor:
    """Tokenise RNA sequences into a padded (B, L) tensor of token ids."""
    width = max(len(s) for s in sequences)
    out = torch.full((len(sequences), width), PAD_TOKEN, dtype=torch.long)
    for i, sequence in enumerate(sequences):
        for j, base in enumerate(sequence.upper().replace("T", "U")):
            out[i, j] = TOKENS.get(base, PAD_TOKEN)
    return out.to(device)


def mask_diagonal(matrix: np.ndarray, width: int = 3, value: float = 0.0) -> np.ndarray:
    """Zero the band around the diagonal, as the authors do before decoding."""
    out = matrix.copy()
    n = out.shape[0]
    rows, cols = np.indices((n, n))
    out[np.abs(rows - cols) < width] = value
    return out


class RNet:
    """Both RibonanzaNet checkpoints, loaded once and reused.

    >>> rnet = RNet()                             # doctest: +SKIP
    >>> rnet.reactivity(["GGGAAACCC"])[0][:, 0]   # doctest: +SKIP
    """

    def __init__(
        self,
        reactivity_weights: Path | None = None,
        ss_weights: Path | None = None,
        device: str = "cpu",
        load_reactivity: bool = True,
        load_ss: bool = True,
    ):
        module = _import_network_module()
        self.device = device
        self.reactivity_model = None
        self.ss_model = None

        if load_reactivity:
            path = reactivity_weights or WEIGHTS / "RibonanzaNet.pt"
            model = module.RibonanzaNet(_Config(**MODEL_CONFIG))
            model.load_state_dict(torch.load(path, map_location="cpu", weights_only=True))
            self.reactivity_model = model.to(device).eval()

        if load_ss:
            path = ss_weights or WEIGHTS / "RibonanzaNet-SS.pt"
            model = module.finetuned_RibonanzaNet(_Config(**MODEL_CONFIG))
            model.load_state_dict(torch.load(path, map_location="cpu", weights_only=True))
            self.ss_model = model.to(device).eval()

    @torch.no_grad()
    def reactivity(self, sequences: list[str]) -> list[np.ndarray]:
        """Predicted (2A3, DMS) reactivity per sequence, shape (L, 2).

        Batch sequences of similar length. The convolutions run across the
        padded tail, so a short sequence batched with a much longer one comes
        back with a different profile than it does alone; equal lengths are
        exact (see tests/test_rnet.py).
        """
        if self.reactivity_model is None:
            raise RuntimeError("reactivity model not loaded")
        tokens = encode(sequences, self.device)
        mask = (tokens != PAD_TOKEN).long()
        output = self.reactivity_model(tokens, src_mask=mask).cpu().numpy()
        return [output[i, : len(s)] for i, s in enumerate(sequences)]

    @torch.no_grad()
    def pair_probabilities(self, sequences: list[str]) -> list[np.ndarray]:
        """Predicted base-pair probability matrix per sequence, shape (L, L)."""
        if self.ss_model is None:
            raise RuntimeError("secondary-structure model not loaded")
        tokens = encode(sequences, self.device)
        logits = self.ss_model(tokens).sigmoid().cpu().numpy()
        return [logits[i, : len(s), : len(s)] for i, s in enumerate(sequences)]

    def structures(
        self, sequences: list[str], theta: float = 0.5, min_len_helix: int = 1
    ) -> list[str]:
        """Predicted secondary structures in dot-bracket notation.

        Decoding matches the authors' `DQN_env.get_structure`: Hungarian
        matching on the pair probabilities with the diagonal band masked out.
        """
        from arnie.pk_predictors import _hungarian  # noqa: PLC0415

        out = []
        for matrix in self.pair_probabilities(sequences):
            structure, _ = _hungarian(
                mask_diagonal(matrix), theta=theta, min_len_helix=min_len_helix
            )
            out.append(structure)
        return out
