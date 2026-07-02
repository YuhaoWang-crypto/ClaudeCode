"""Generator backends (the Evo 2-shaped designer)."""

from .base import Generator
from .mock import MockGenerator

__all__ = ["Generator", "MockGenerator", "Evo2Generator"]


def __getattr__(name: str):
    # Lazily expose the real adapter without importing the heavy dependency.
    if name == "Evo2Generator":
        from .evo2 import Evo2Generator

        return Evo2Generator
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
