"""Download and load the trained AF2BIND logistic-regression heads.

The released checkpoints are pickled jax arrays, but the payload underneath is
plain numpy. `load_head` therefore uses a restricted unpickler that swaps
`jax._src.array._reconstruct_array` for a numpy-only equivalent, so the head can
be loaded (and evaluated) without jax installed.
"""

from __future__ import annotations

import io
import os
import pickle
import urllib.request
import zipfile
from pathlib import Path

import numpy as np

WEIGHTS_URL = (
    "https://raw.githubusercontent.com/sokrypton/af2bind/main/"
    "attempt_7_2k_lam0-03.zip"
)
WEIGHTS_ZIP_NAME = "attempt_7_2k_lam0-03.zip"
WEIGHTS_ROOT = "attempt_7_2k_lam0-03"

#: The two AF2-pair-only heads from the paper. The `_nosc` variant is trained on
#: features from a target whose side chains were masked, and *must* be paired
#: with ``mask_sidechains=True`` at feature-extraction time.
MODEL_TYPES = {
    True: "split_nosc_pair_A_split_nosc_pair_B",   # mask_sidechains=True
    False: "split_pair_A_split_pair_B",            # mask_sidechains=False
}

#: Checkpoints ship as 10 independently-trained folds, seed 0..9.
N_SEEDS = 10


def default_cache_dir() -> Path:
    return Path(os.environ.get("AF2BIND_CACHE", Path.home() / ".cache" / "af2bind"))


def _reconstruct_array(fun, args, arr_state, aval_state):
    """numpy-only stand-in for ``jax._src.array._reconstruct_array``."""
    value = fun(*args)
    value.__setstate__(arr_state)
    return value


class _NumpyUnpickler(pickle.Unpickler):
    """Unpickler that refuses anything but numpy reconstruction + the jax shim."""

    _ALLOWED = {
        ("numpy", "ndarray"),
        ("numpy", "dtype"),
        ("numpy.core.multiarray", "_reconstruct"),
        ("numpy._core.multiarray", "_reconstruct"),
    }

    def find_class(self, module, name):
        if module == "jax._src.array" and name == "_reconstruct_array":
            return _reconstruct_array
        if (module, name) in self._ALLOWED:
            return super().find_class(module, name)
        raise pickle.UnpicklingError(f"blocked global {module}.{name}")


def ensure_weights(cache_dir: Path | None = None) -> Path:
    """Download the weights zip if absent. Returns the path to the zip."""
    cache_dir = Path(cache_dir or default_cache_dir())
    cache_dir.mkdir(parents=True, exist_ok=True)
    zip_path = cache_dir / WEIGHTS_ZIP_NAME
    if zip_path.exists() and zip_path.stat().st_size > 1_000_000:
        return zip_path
    tmp = zip_path.with_suffix(".zip.part")
    with urllib.request.urlopen(WEIGHTS_URL, timeout=300) as r, open(tmp, "wb") as f:
        f.write(r.read())
    if not zipfile.is_zipfile(tmp):
        tmp.unlink(missing_ok=True)
        raise RuntimeError(
            f"downloaded {WEIGHTS_URL} but it is not a zip file "
            "(a proxy or auth wall probably intercepted the request)"
        )
    tmp.replace(zip_path)
    return zip_path


def load_head(
    mask_sidechains: bool = True,
    seed: int = 0,
    cache_dir: Path | None = None,
) -> dict[str, np.ndarray]:
    """Load one trained head.

    Returns a dict with ``mean``/``std`` (5120,) standardisation stats and
    ``w`` (5120, 1) / ``b`` (1,) logistic-regression coefficients.
    """
    if not 0 <= seed < N_SEEDS:
        raise ValueError(f"seed must be in 0..{N_SEEDS - 1}, got {seed}")
    zip_path = ensure_weights(cache_dir)
    model_type = MODEL_TYPES[bool(mask_sidechains)]
    member = f"{WEIGHTS_ROOT}/{model_type}_{seed}.pickle"
    with zipfile.ZipFile(zip_path) as z:
        raw = z.read(member)
    obj = _NumpyUnpickler(io.BytesIO(raw)).load()
    flat = {**obj["~"], **obj["linear"]}
    return {k: np.asarray(v, dtype=np.float32) for k, v in flat.items()}


def load_heads(
    mask_sidechains: bool = True,
    seeds=(0,),
    cache_dir: Path | None = None,
) -> list[dict[str, np.ndarray]]:
    return [load_head(mask_sidechains, s, cache_dir) for s in seeds]
