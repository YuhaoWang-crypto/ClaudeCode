#!/usr/bin/env python3
"""
T cell exhaustion <-> effector: a minimal Boolean-network demo.

What this shows, end to end:
  1. A ~14-node core regulatory circuit of CD8 T cell exhaustion vs effector,
     written once as a PyBoolNet `bnet` model (single source of truth).
  2. Its stable cell states (attractors / steady states) -- the "homeostatic
     fixed points" we discussed: one effector/memory-like, one exhausted.
  3. A single- and double-node perturbation scan that forces nodes ON/OFF
     (gene overexpression / knockdown) and asks: which minimal combination
     flips the cell OUT of the exhausted attractor INTO the effector attractor?

Biology encoded (deliberately simplified, but literature-motivated):
  - Chronic antigen -> sustained NFAT WITHOUT AP-1 partner drives the
    exhaustion program (TOX, NR4A, PDCD1).            (Martinez 2015; Chen 2019)
  - TOX is the master exhaustion TF, self-sustaining once on.   (Khan/Scott 2019)
  - TCF7 (TCF1) maintains stem/memory identity and antagonizes TOX. (Siddiqui 2019)
  - AP-1 / c-Jun restores effector function (TBX21, IFNG).        (Lynn 2019)

Run:  python3 tcell_exhaustion_boolean_demo.py
Deps: pyboolnet  (pip install pyboolnet)  -- used for parsing + steady states.
      The perturbation engine below is plain Python and needs nothing else.
"""

from itertools import combinations, product
from pyboolnet.file_exchange import bnet2primes
from pyboolnet.trap_spaces import compute_steady_states

# --------------------------------------------------------------------------
# 1. The model. PyBoolNet bnet syntax:  node, boolean_expression
#    !  = NOT     &  = AND     |  = OR
#    NFAT is the environmental input (chronic stimulation) -> held ON.
#    AP1 is the therapeutically perturbable effector drive (c-Jun).
# --------------------------------------------------------------------------
BNET = """
TCF7,   !TOX & !PRDM1
TOX,    !TCF7 & (NFAT & !AP1 | TOX)
AP1,    AP1
NFAT,   NFAT
TBX21,  AP1 & !TOX & !NR4A
PRDM1,  TBX21 & !TCF7
NR4A,   TOX | (NFAT & !AP1)
BATF,   TOX
IRF4,   BATF & NFAT
PDCD1,  TOX | NR4A
HAVCR2, TOX & NR4A
LAG3,   TOX
EOMES,  TBX21 | TOX
IFNG,   TBX21 & AP1 & !TOX
"""

# Markers used to *classify* an attractor (not forced during the scan).
EFFECTOR_SIG  = {"TCF7": 1, "TOX": 0, "TBX21": 1, "NR4A": 0, "PDCD1": 0, "IFNG": 1}
EXHAUSTED_SIG = {"TCF7": 0, "TOX": 1, "TBX21": 0, "NR4A": 1, "PDCD1": 1, "IFNG": 0}

# Environmental context we analyse: chronic antigen present.
CONTEXT = {"NFAT": 1}

# Nodes we allow a "drug" to perturb. NFAT is environment (not a drug);
# pure readout markers are excluded so a hit means we moved the *program*,
# not just painted the readout.
PERTURBABLE = ["TCF7", "TOX", "AP1", "TBX21", "PRDM1", "NR4A", "BATF", "IRF4"]


# --------------------------------------------------------------------------
# 2. A tiny synchronous Boolean engine driven directly by PyBoolNet primes,
#    so the dynamics can never drift from the published model above.
# --------------------------------------------------------------------------
def _node_next(primes, node, state):
    """Next value of `node`: 1 iff state satisfies any of its '1' prime implicants."""
    for implicant in primes[node][1]:            # primes[node][1] = implicants -> 1
        if all(state[k] == v for k, v in implicant.items()):
            return 1
    return 0


def step(primes, state, overrides):
    """One synchronous update; overridden nodes are pinned (a drug held constant)."""
    nxt = {}
    for node in primes:
        nxt[node] = overrides[node] if node in overrides else _node_next(primes, node, state)
    return nxt


def run_to_attractor(primes, start, overrides, max_iter=200):
    """Iterate from `start` under `overrides` to a fixed point or limit cycle.
    Returns (kind, states) where kind in {'point','cycle'}."""
    seen = {}
    state = dict(start)
    state.update(overrides)
    trajectory = []
    for i in range(max_iter):
        key = tuple(sorted(state.items()))
        if key in seen:
            cycle = trajectory[seen[key]:]
            return ("point" if len(cycle) == 1 else "cycle", cycle)
        seen[key] = i
        trajectory.append(state)
        state = step(primes, state, overrides)
    return ("cycle", [state])  # did not settle; treat tail as attractor


def classify(state, sig):
    return all(state[k] == v for k, v in sig.items())


def matches_effector(attractor_states, free_only):
    """An attractor is 'effector' if some state in it matches the effector
    signature, ignoring nodes that were forced (so we credit the program, not
    the override itself)."""
    sig = {k: v for k, v in EFFECTOR_SIG.items() if k in free_only}
    return any(all(s[k] == v for k, v in sig.items()) for s in attractor_states)


# --------------------------------------------------------------------------
# 3. Find the resting attractors of the unperturbed network (context fixed).
# --------------------------------------------------------------------------
def all_attractors(primes, context):
    """Exhaustive synchronous attractor search over the free-node state space."""
    free = [n for n in primes if n not in context]
    found = {}      # frozen attractor -> representative states
    for bits in product((0, 1), repeat=len(free)):
        start = dict(context)
        start.update(dict(zip(free, bits)))
        kind, states = run_to_attractor(primes, start, context)
        key = frozenset(tuple(sorted(s.items())) for s in states)
        found[key] = states
    return list(found.values())


# --------------------------------------------------------------------------
# 4. Perturbation scan: which forced combos move the exhausted state -> effector?
# --------------------------------------------------------------------------
def scan(primes, exhausted_state, context, max_combo=2):
    candidates = [(n, v) for n in PERTURBABLE for v in (0, 1)]
    hits = []  # (size, dict_of_overrides)
    for size in range(1, max_combo + 1):
        for combo in combinations(candidates, size):
            nodes = [n for n, _ in combo]
            if len(set(nodes)) != size:          # no node forced two ways at once
                continue
            overrides = dict(context)
            overrides.update(dict(combo))
            free_only = [n for n in primes if n not in overrides]
            kind, states = run_to_attractor(primes, exhausted_state, overrides)
            if matches_effector(states, free_only):
                hits.append((size, dict(combo)))
    return hits


def fmt(overrides):
    arrow = {1: "↑ON ", 0: "↓OFF"}
    return " + ".join(f"{n} {arrow[v]}" for n, v in sorted(overrides.items()))


# --------------------------------------------------------------------------
# 5. Drive it all and print a readable report.
# --------------------------------------------------------------------------
def main():
    primes = bnet2primes(BNET)
    order = list(primes.keys())

    print("=" * 70)
    print("T CELL EXHAUSTION <-> EFFECTOR  |  Boolean network demo")
    print("=" * 70)
    print(f"Nodes ({len(order)}): {', '.join(order)}")
    print(f"Context (environmental input held fixed): NFAT = 1  (chronic antigen)\n")

    # --- Resting attractors (steady states), cross-checked two ways ---
    print("-" * 70)
    print("STEP 1 | Stable cell states (attractors = homeostatic fixed points)")
    print("-" * 70)

    pbn_steady = [s for s in compute_steady_states(primes) if s["NFAT"] == 1]
    attractors = all_attractors(primes, CONTEXT)
    print(f"PyBoolNet steady states (NFAT=1): {len(pbn_steady)}")
    print(f"Synchronous attractors found (NFAT=1): {len(attractors)}  "
          f"(engine cross-check)\n")

    exhausted_states = []
    effector_states = []
    for att in attractors:
        s = att[0]
        tag = ("EXHAUSTED" if classify(s, EXHAUSTED_SIG)
               else "EFFECTOR" if classify(s, EFFECTOR_SIG)
               else "other")
        on = [n for n in order if s[n] == 1]
        kind = "fixed point" if len(att) == 1 else f"cycle(len={len(att)})"
        print(f"  [{tag:9}] {kind:12} ON: {', '.join(on) if on else '(none)'}")
        if tag == "EXHAUSTED":
            exhausted_states.append(s)
        elif tag == "EFFECTOR":
            effector_states.append(s)

    # Canonical chronic exhaustion = sustained NFAT WITHOUT an AP-1 partner.
    # Start the scan from the AP1-OFF exhausted state (the one therapy must escape).
    exhausted_state = next((s for s in exhausted_states if s["AP1"] == 0),
                           exhausted_states[0] if exhausted_states else None)
    effector_state = next((s for s in effector_states if s["AP1"] == 1),
                          effector_states[0] if effector_states else None)

    if not (exhausted_state and effector_state):
        print("\n!! Expected both attractors; check the model. Aborting scan.")
        return
    print(f"\nScan start  = canonical EXHAUSTED state (AP1 OFF, chronic NFAT).")
    print(f"Scan target = EFFECTOR program (AP1 ON).")

    # --- Perturbation scan ---
    print()
    print("-" * 70)
    print("STEP 2 | Perturbation scan: minimal combo to flip EXHAUSTED -> EFFECTOR")
    print("-" * 70)
    print("Forcing nodes ON (overexpression) / OFF (knockdown), starting from")
    print("the exhausted attractor, then letting the network re-settle.\n")

    hits = scan(primes, exhausted_state, CONTEXT, max_combo=2)

    singles = sorted({frozenset(h[1].items()) for h in hits if h[0] == 1})
    doubles = sorted({frozenset(h[1].items()) for h in hits if h[0] == 2})

    print(f"Single-node interventions that flip the state: {len(singles)}")
    for h in singles:
        print(f"    - {fmt(dict(h))}")
    if not singles:
        print("    (none -- no monotherapy escapes the exhausted attractor)")

    # Report only doubles that are not supersets of a working single.
    minimal_doubles = [d for d in doubles
                       if not any(s.issubset(d) for s in singles)]
    print(f"\nMinimal two-node combinations that flip the state: {len(minimal_doubles)}")
    for d in minimal_doubles:
        print(f"    - {fmt(dict(d))}")

    # --- Interpretation ---
    print()
    print("-" * 70)
    print("STEP 3 | Read-out")
    print("-" * 70)
    if not singles and minimal_doubles:
        print("No single node flips the switch: the TOX self-sustaining loop and the")
        print("NFAT-without-AP1 drive must BOTH be overcome. Every working pair couples")
        print("an effector drive (AP1/c-Jun UP) with an exhaustion brake (TOX DOWN or")
        print("TCF7 UP) -- mirroring c-Jun-overexpression + TOX-deletion rescue in vivo.")
    print("\nThis is the concrete 'complex / combination needed to lever a signal':")
    print("the minimal forced-node sets above are exactly that, read straight off")
    print("the attractor structure of the regulatory network.")


if __name__ == "__main__":
    main()
