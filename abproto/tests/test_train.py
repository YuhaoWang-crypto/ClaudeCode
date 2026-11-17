"""Offline tests for the structure dataset + IPA training on real-shaped data.

No network: builds synthetic ComplexStructures. Torch parts skip if absent.
"""

from __future__ import annotations

import numpy as np
import pytest

from abproto.data.structures import Chain, ComplexStructure
from abproto.train.dataset import StructureDataset, collate, example_from_structure


def _synth_complex(pdb="TEST", seed=0):
    rng = np.random.default_rng(seed)
    h = Chain("H", "E" * 110, rng.normal(size=(110, 3)) * 10, "heavy",
              {"1": (25, 34), "2": (50, 58), "3": (97, 105)})
    l = Chain("L", "D" * 105, rng.normal(size=(105, 3)) * 10 + 20, "light",
              {"1": (24, 34), "2": (50, 53), "3": (89, 98)})
    ag = Chain("A", "A" * 60, rng.normal(size=(60, 3)) * 10, "antigen")
    return ComplexStructure(pdb, heavy=h, light=l, antigen=[ag])


def test_example_shapes_and_masks():
    e = example_from_structure(_synth_complex(), flank=3, ka=32)
    assert e is not None
    Ld = len(e["design_tokens"])
    assert e["design_coords"].shape == (Ld, 3)
    assert e["generate_mask"].shape == (Ld,)
    assert e["generate_mask"].sum() >= 3           # CDR positions present
    assert not e["generate_mask"].all()            # framework flank present too
    assert e["antigen_coords"].shape[0] <= 32
    # centred on antigen centroid (before scaling it was zero-mean)
    assert np.abs(e["antigen_coords"].mean(axis=0)).max() < 1e-5


def test_collate_pads_and_masks():
    torch = pytest.importorskip("torch")
    e1 = example_from_structure(_synth_complex("A", 1))
    e2 = example_from_structure(_synth_complex("B", 2))
    batch = collate([e1, e2])
    B, Ld = batch["design_tokens"].shape
    assert B == 2
    assert batch["design_coords"].shape == (2, Ld, 3)
    assert batch["design_pad"].dtype == torch.bool
    # the shorter example has padding; generate positions are a subset of valid
    assert (batch["generate_mask"] & ~batch["design_pad"]).sum() == 0


def test_ipa_trains_on_structure_batches():
    torch = pytest.importorskip("torch")
    from abproto.model.net import CoDesignModel
    torch.manual_seed(0)
    ds = StructureDataset([example_from_structure(_synth_complex(f"C{i}", i)) for i in range(6)])
    model = CoDesignModel(d_model=48, n_layers=2, lr=2e-3, backbone="ipa")
    losses = []
    rng = np.random.default_rng(0)
    for _ in range(40):
        idx = rng.integers(0, len(ds), size=3)
        losses.append(model.training_step(collate([ds[i] for i in idx]))["loss"])
    assert np.mean(losses[-10:]) < np.mean(losses[:10])


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(f"PASS {name}")
    print("all train tests passed")
