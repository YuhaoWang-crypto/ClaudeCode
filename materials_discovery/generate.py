"""
generate.py -- GNoME's two candidate pipelines.

  structural()   -- substitute ions into a known prototype (the dominant GNoME
                    pipeline). Enumerate role->element assignments, keep only
                    charge-neutral (and, for perovskite, tolerance-factor-formable)
                    compositions. This is the cheap physical filter that turns a
                    combinatorial explosion into a screenable candidate set.

  compositional() -- enumerate stoichiometries over a fixed element palette within
                    small integer bounds, keep charge-neutral ones. The randomized/
                    template-free arm of GNoME.

Both are deterministic (sorted enumeration) so runs are reproducible.
"""
from __future__ import annotations
from itertools import product
from dataclasses import dataclass
from .composition import Composition
from .prototypes import (PROTOTYPES, OXIDATION_STATES, SHANNON_RADII,
                         tolerance_factor)


@dataclass(frozen=True)
class Candidate:
    composition: Composition
    prototype: str
    assignment: tuple        # ((role, element, oxidation), ...)
    source: str              # "structural" | "compositional"


def _charge(assignment_counts) -> float:
    """assignment_counts: list of (element, oxidation, count) -> net charge."""
    return sum(ox * n for _, ox, n in assignment_counts)


def structural(prototype_name: str, role_candidates: dict,
               tol_range=(0.8, 1.02)) -> list:
    """Enumerate charge-neutral substitutions of `role_candidates` into a prototype.

    role_candidates: {role: [element, ...]} allowed elements per site role.
    """
    proto = PROTOTYPES[prototype_name]
    out = []
    roles = proto.roles
    choices = [role_candidates[r] for r in roles]
    for combo in product(*choices):
        # assign an oxidation state per role element consistent with role sign
        ox_options = []
        ok = True
        for r, el in zip(roles, combo):
            sign = proto.role_sign[r]
            oxs = [o for o in OXIDATION_STATES.get(el, []) if (o > 0) == (sign > 0)]
            if not oxs:
                ok = False
                break
            ox_options.append(oxs)
        if not ok:
            continue
        for oxset in product(*ox_options):
            counts = []
            comp = {}
            for r, el, ox in zip(roles, combo, oxset):
                n = proto.counts[r]
                counts.append((el, ox, n))
                comp[el] = comp.get(el, 0.0) + n
            if abs(_charge(counts)) > 1e-9:
                continue
            # no element playing both cation and anion role
            if len({el for el, _, _ in counts}) < len(counts):
                continue
            # perovskite: Goldschmidt tolerance-factor filter
            if prototype_name == "perovskite":
                try:
                    rA = SHANNON_RADII[(combo[0], oxset[0])]
                    rB = SHANNON_RADII[(combo[1], oxset[1])]
                    rX = SHANNON_RADII[(combo[2], oxset[2])]
                except KeyError:
                    continue
                t = tolerance_factor(rA, rB, rX)
                if not (tol_range[0] <= t <= tol_range[1]):
                    continue
            assignment = tuple((r, el, ox) for r, el, ox in zip(roles, combo, oxset))
            out.append(Candidate(Composition.from_dict(comp), prototype_name,
                                 assignment, "structural"))
    # dedupe by reduced formula, keep first (deterministic)
    seen = {}
    for c in out:
        k = c.composition.reduced_formula()
        if k not in seen:
            seen[k] = c
    return sorted(seen.values(), key=lambda c: c.composition.reduced_formula())


def compositional(elements: list, max_count: int = 4) -> list:
    """Enumerate charge-neutral integer stoichiometries over `elements`."""
    ox_lists = []
    for el in elements:
        oxs = OXIDATION_STATES.get(el, [])
        if not oxs:
            return []
        ox_lists.append(oxs)
    out = {}
    ranges = [range(0, max_count + 1) for _ in elements]
    for counts in product(*ranges):
        if sum(counts) < 2 or sum(1 for c in counts if c > 0) < 2:
            continue
        # must contain at least one cation and one anion element present
        for oxset in product(*ox_lists):
            charge = sum(o * c for o, c in zip(oxset, counts))
            if abs(charge) > 1e-9:
                continue
            has_cat = any(o > 0 and c > 0 for o, c in zip(oxset, counts))
            has_an = any(o < 0 and c > 0 for o, c in zip(oxset, counts))
            if not (has_cat and has_an):
                continue
            comp = {el: c for el, c in zip(elements, counts) if c > 0}
            cand = Composition.from_dict(comp)
            k = cand.reduced_formula()
            if k not in out:
                out[k] = Candidate(cand, "", (), "compositional")
            break
    return sorted(out.values(), key=lambda c: c.composition.reduced_formula())
