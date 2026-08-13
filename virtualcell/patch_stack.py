"""Patch a crash in ``arc-stack`` 0.1.3's iterative generation.

``stack-generation`` fails partway through with::

    ValueError: Quantiles must be in the range [0, 1]
    stack/models/core/inference.py:700 in get_incontext_prediction

reproducible on the project's own tutorial (``notebooks/tutorial-predict.ipynb``,
OpenProblems donor 1, ``--num-steps 5``).

**Cause.** Generation unmasks cells over ``num_steps`` rounds on a schedule
``mr_list = 1 - t``.  Each round computes how far to unmask as a quantile of the
predicted logits::

    f = is_masked.sum() / is_masked.shape[0]          # fraction still masked
    unmask_rate = (f - mask_rate) * is_masked.shape[0] / is_masked.sum()
                = 1 - mask_rate / f

That is only in ``[0, 1]`` while ``f >= mask_rate``.  The schedule would hold it
there, except the very next lines re-mask on the model's own output::

    is_masked[~cell_indices_to_keep] = False
    is_masked[new_test_logit > 0] = True

so ``f`` is data-dependent, not schedule-determined.  Whenever the model unmasks
faster than the schedule anticipates, ``f`` drops below ``mask_rate``,
``unmask_rate`` goes negative, and ``np.quantile`` raises.  ``f == 0`` divides by
zero as well.

**Fix.** Clamp to ``[0, 1]`` and guard the empty case.  ``unmask_rate = 0`` takes
the quantile at the minimum, so nothing further is unmasked this round — the
correct reading of "already at or past the scheduled mask rate".

This monkey-patches at import time rather than editing the installed package, so
the change is visible in version control and travels with the analysis.  Import
it before calling any Stack generation entry point::

    python -m virtualcell.patch_stack -- stack-generation --checkpoint ...
"""

from __future__ import annotations

import numpy as np


def apply() -> bool:
    """Wrap ``np.quantile`` so out-of-range rates from Stack are clamped.

    Patching the numpy call rather than rewriting Stack's method keeps the fix
    to one line of behaviour and survives point releases that move code around.
    Only values that would already have raised are altered.
    """
    if getattr(np.quantile, "_stack_patched", False):
        return False

    original = np.quantile

    def quantile(a, q, *args, **kwargs):
        arr = np.asarray(a)
        if np.isscalar(q) and (not np.isfinite(q) or q < 0.0 or q > 1.0):
            if arr.size == 0:
                return np.float64("nan")
            q = float(np.clip(np.nan_to_num(q, nan=0.0), 0.0, 1.0))
        return original(arr, q, *args, **kwargs)

    quantile._stack_patched = True          # type: ignore[attr-defined]
    quantile._stack_original = original     # type: ignore[attr-defined]
    np.quantile = quantile                  # type: ignore[assignment]
    return True


def main() -> None:
    """Apply the patch, then hand off to a Stack console entry point."""
    import sys

    apply()
    argv = sys.argv[1:]
    if argv and argv[0] == "--":
        argv = argv[1:]
    if not argv:
        print("patched numpy.quantile; nothing to run")
        return

    entry = argv[0]
    sys.argv = argv
    mapping = {
        "stack-generation": ("stack.cli.generation", "main"),
        "stack-embedding": ("stack.cli.embedding", "main"),
        "stack-train": ("stack.cli.launch_training", "main"),
        "stack-finetune": ("stack.cli.launch_finetuning", "main"),
    }
    if entry not in mapping:
        raise SystemExit(f"unknown entry point {entry!r}; expected one of "
                         f"{', '.join(sorted(mapping))}")
    module, func = mapping[entry]
    import importlib

    getattr(importlib.import_module(module), func)()


if __name__ == "__main__":
    main()
