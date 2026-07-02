"""Oracle backends (the AlphaGenome-shaped judge)."""

from .base import Oracle
from .mock import MockOracle

__all__ = ["Oracle", "MockOracle", "AlphaGenomeOracle"]


def __getattr__(name: str):
    # Lazily expose the real adapter without importing the heavy dependency.
    if name == "AlphaGenomeOracle":
        from .alphagenome import AlphaGenomeOracle

        return AlphaGenomeOracle
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
