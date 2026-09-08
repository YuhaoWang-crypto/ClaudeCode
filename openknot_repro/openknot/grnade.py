"""gRNAde inference: sample designs for a target structure.

gRNAde is one of the three deep-learning frameworks the paper found competitive
with expert human designers. Its code (chaitjo/geometric-rna-design, MIT) is
cloned to `third_party/` at a pinned commit by `scripts/setup_grnade.sh`, which
also fetches the checkpoint the paper used from HuggingFace.

Two design modes, both from the authors' `design.py`:

  "2d"  conditioned on the target secondary structure only. This is the
        gRNAde-no3d variant of the paper -- no 3D backbone is needed.
  "3d"  conditioned on a 3D backbone as well, read from a PDB file.

The authors' code imports two compiled PyTorch Geometric extensions,
torch_scatter and torch_cluster, for one function each. Building them from
source takes far longer than this pipeline needs, so when they are absent a
pure-torch stand-in for those two functions is registered instead; the real
packages take precedence whenever they are importable.
"""

from __future__ import annotations

import importlib.machinery
import sys
import types
from pathlib import Path

import torch

ROOT = Path(__file__).resolve().parent.parent
GRNADE_PATH = ROOT / "third_party" / "geometric-rna-design"
CHECKPOINTS = ROOT / "weights" / "grnade"


def _scatter_add(src, index, dim: int = -1, out=None, dim_size: int | None = None):
    """torch_scatter.scatter_add, in plain torch."""
    if dim < 0:
        dim = src.dim() + dim
    if out is None:
        size = list(src.shape)
        size[dim] = dim_size if dim_size is not None else int(index.max()) + 1
        out = torch.zeros(size, dtype=src.dtype, device=src.device)
    shape = [1] * src.dim()
    shape[dim] = -1
    expanded = index.view(shape).expand_as(src)
    return out.scatter_add_(dim, expanded, src)


def _knn_graph(x, k: int, batch=None, loop: bool = False, flow: str = "source_to_target"):
    """torch_cluster.knn_graph, in plain torch.

    Returns edge_index [2, N*k]; with the default flow the first row holds the
    neighbour and the second the query node, matching torch_cluster.
    """
    distance = torch.cdist(x, x)
    if batch is not None:
        different = batch.view(-1, 1) != batch.view(1, -1)
        distance = distance.masked_fill(different, float("inf"))
    if not loop:
        distance.fill_diagonal_(float("inf"))
    k = min(k, x.size(0) - (0 if loop else 1))
    neighbours = distance.topk(k, dim=1, largest=False).indices
    targets = torch.arange(x.size(0), device=x.device).view(-1, 1).expand(-1, k).reshape(-1)
    sources = neighbours.reshape(-1)
    if flow == "source_to_target":
        return torch.stack([sources, targets])
    return torch.stack([targets, sources])


def _register_shim(name: str, **functions) -> None:
    module = types.ModuleType(name)
    # importlib.util.find_spec() raises on a module without a spec, and
    # torch_geometric probes for both of these at import time.
    module.__spec__ = importlib.machinery.ModuleSpec(name, None)
    for attribute, value in functions.items():
        setattr(module, attribute, value)
    sys.modules[name] = module


def _install_geometric_shims() -> list[str]:
    """Provide torch_scatter/torch_cluster stand-ins only if the real ones are missing."""
    installed = []
    try:
        import torch_scatter  # noqa: F401, PLC0415
    except ImportError:
        _register_shim("torch_scatter", scatter_add=_scatter_add, scatter=_scatter_add)
        installed.append("torch_scatter")
    try:
        import torch_cluster  # noqa: F401, PLC0415
    except ImportError:
        _register_shim("torch_cluster", knn_graph=_knn_graph)
        installed.append("torch_cluster")
    return installed


def _set_grnade_env() -> None:
    """gRNAde reads its paths from the environment via a .env file.

    Only PROJECT_PATH is load-bearing here. X3DNA and EternaFold are optional
    external tools used to derive secondary structures from coordinates, which
    this pipeline does not do -- but `src.constants` joins their paths at import
    time, so they need to be set to something.
    """
    import os  # noqa: PLC0415

    os.environ.setdefault("PROJECT_PATH", str(GRNADE_PATH))
    os.environ.setdefault("DATA_PATH", str(GRNADE_PATH / "data"))
    os.environ.setdefault("X3DNA", str(GRNADE_PATH / "tools" / "x3dna-v2.4"))
    os.environ.setdefault("ETERNAFOLD", str(GRNADE_PATH / "tools" / "EternaFold"))
    os.environ.setdefault("WANDB_MODE", "disabled")


def _prepare_imports():
    if not GRNADE_PATH.exists():
        raise FileNotFoundError(
            f"{GRNADE_PATH} not found. Run scripts/setup_grnade.sh to fetch gRNAde."
        )
    _set_grnade_env()
    shims = _install_geometric_shims()
    if str(GRNADE_PATH) not in sys.path:
        sys.path.insert(0, str(GRNADE_PATH))
    return shims


# Model shape from the authors' configs/design.yaml; it must match the checkpoint.
MODEL_CONFIG = dict(
    node_in_dim=(15, 4),
    node_h_dim=(128, 16),
    edge_in_dim=(132, 3),
    edge_h_dim=(64, 4),
    num_layers=4,
    drop_rate=0.5,
    out_dim=4,
)
FEATURIZER_CONFIG = dict(
    radius=0.0,
    top_k=32,
    num_rbf=32,
    num_posenc=32,
    max_num_conformers=1,
    noise_scale=0.1,
    drop_prob_3d=0.75,
)


class GRNAde:
    """The gRNAde model plus its featurizer, loaded once and reused."""

    def __init__(self, checkpoint: Path | None = None, mode: str = "2d", device: str = "cpu"):
        self.shims = _prepare_imports()
        from src.data.featurizer import RNAGraphFeaturizer  # noqa: PLC0415
        from src.models import gRNAde  # noqa: PLC0415

        self.mode = mode
        self.device = device
        self.featurizer = RNAGraphFeaturizer(
            split="test" if mode == "3d" else "test_2d", **FEATURIZER_CONFIG
        )
        path = checkpoint or CHECKPOINTS / "gRNAde_drop3d@0.75_maxlen@500.h5"
        model = gRNAde(**MODEL_CONFIG)
        model.load_state_dict(torch.load(path, map_location="cpu"))
        self.model = model.to(device).eval()

    def featurize_2d(self, sequence: str, sec_struct: str):
        """Featurise a target from its secondary structure alone (3D masked out)."""
        from src.constants import FILL_VALUE  # noqa: PLC0415

        raw = {
            "sequence": sequence,
            "coords_list": [torch.ones(len(sec_struct), 3, 3) * FILL_VALUE],
            "sec_struct_list": [sec_struct],
        }
        return self.featurizer(raw).to(self.device)

    def featurize_pdb(self, pdb_path: Path, sec_struct: str | None = None):
        """Featurise a target from a PDB backbone."""
        _, raw = self.featurizer.featurize_from_pdb_file(str(pdb_path))
        if sec_struct is not None:
            raw["sec_struct_list"] = [sec_struct]
        return self.featurizer(raw).to(self.device), raw

    @torch.no_grad()
    def sample(self, data, n_samples: int = 8, temperature: float = 0.5) -> list[str]:
        """Sample sequences for a featurised target."""
        from src.constants import NUM_TO_LETTER  # noqa: PLC0415

        samples, _ = self.model.sample(data, n_samples, temperature, return_logits=True)
        mask = data.mask_seq.cpu().numpy()
        array = samples.cpu().numpy()
        out = []
        for row in array:
            out.append("".join(NUM_TO_LETTER[int(i)] for i in row[mask]))
        return out


def load_targets() -> "list[dict]":
    """The 40 Round 3 and Round 4 targets, from gRNAde's own benchmark metadata."""
    import pandas as pd  # noqa: PLC0415

    base = GRNADE_PATH / "projects" / "openknot_benchmark"
    frames = []
    for name, round_number in (("metadata_7a.csv", 3), ("metadata_7b.csv", 4)):
        frame = pd.read_csv(base / name)
        frame["round"] = round_number
        frames.append(frame)
    targets = pd.concat(frames, ignore_index=True)
    targets["structure_dir"] = targets["puzzleID"].map(lambda p: base / "structures" / str(p))
    return targets.to_dict("records")
