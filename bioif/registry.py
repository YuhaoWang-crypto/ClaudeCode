"""
bioif.registry -- routing.

The payoff of typing the seams is that composition becomes a search problem
rather than an integration project. You register an adapter once; every chain
that needs that hop finds it. Swapping a splice model for a better one, or
adding a model that reports Ki instead of IC50, is a registration, not a
rewrite of the pipeline.

Routing prefers short paths, and among equal-length paths prefers the one
with fewer assumption-bearing bridges -- i.e. it will take arithmetic over a
cross-level assumption whenever both reach the goal.
"""
from __future__ import annotations

from collections import deque

from .adapter import Adapter, BRIDGE
from .core import Quantity
from .quantities import LOSSLESS


class Registry:
    def __init__(self, include_lossless: bool = True):
        self._by_input: dict[str, list[Adapter]] = {}
        self.adapters: list[Adapter] = []
        if include_lossless:
            for c in LOSSLESS:
                self.register(c)

    def register(self, adapter: Adapter) -> Adapter:
        self._by_input.setdefault(adapter.consumes.key, []).append(adapter)
        self.adapters.append(adapter)
        return adapter

    def route(self, src: Quantity, dst: Quantity,
              max_len: int = 8) -> list[Adapter] | None:
        """Breadth-first shortest path, tie-broken toward fewer bridges."""
        if src.key == dst.key:
            return []
        best: list[Adapter] | None = None
        q: deque[tuple[str, list[Adapter]]] = deque([(src.key, [])])
        seen_cost: dict[str, int] = {src.key: 0}
        while q:
            node, path = q.popleft()
            if len(path) >= max_len:
                continue
            for a in self._by_input.get(node, []):
                if any(a is p for p in path):
                    continue
                npath = path + [a]
                nk = a.produces.key
                if nk == dst.key:
                    if best is None or _cost(npath) < _cost(best):
                        best = npath
                    continue
                cost = len(npath)
                if seen_cost.get(nk, 10 ** 9) < cost:
                    continue
                seen_cost[nk] = cost
                q.append((nk, npath))
        return best

    def describe(self) -> str:
        lines = ["registered edges:"]
        for a in self.adapters:
            lines.append(f"  [{a.kind:9s}] {a.name:28s} "
                         f"{a.consumes.name} -> {a.produces.name}")
        return "\n".join(lines)


def _cost(path: list[Adapter]) -> tuple[int, int]:
    return (len(path), sum(1 for a in path if a.kind == BRIDGE))
