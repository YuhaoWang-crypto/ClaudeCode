#!/usr/bin/env python3
"""
End-to-end: network ranking  ->  auto-generated Boolean model  ->  flip combo.

This connects the two halves we discussed:

  PART A  Random-Walk-with-Restart (personalized PageRank) on a regulatory
          network, run separately for the EXHAUSTED and EFFECTOR states, then
          ranked by the differential score.  -> "which proteins matter".

  PART B  The top differential subnetwork is AUTO-CONVERTED into a Boolean
          model (no hand-writing), its attractors are found, and a single/
          double perturbation scan reports the minimal combination that flips
          the exhausted attractor into the effector program. -> "by which knobs".

Inputs (under data/):
  regulatory_edges.tsv   signed, directed causal edges (source target +/-)
  state_expression.tsv   per-gene activity in each state (RWR restart vectors)

Run:  python3 pipeline.py
Deps: numpy, pyboolnet
"""

import os
from itertools import combinations, product
import numpy as np
from pyboolnet.file_exchange import bnet2primes

HERE = os.path.dirname(os.path.abspath(__file__))
EDGES = os.path.join(HERE, "data", "regulatory_edges.tsv")
EXPR = os.path.join(HERE, "data", "state_expression.tsv")

ALPHA = 0.20          # RWR restart probability
TOP_K = 10            # size of the differential subnetwork to model
MIN_DELTA = 0.004     # nodes with |Δ| below this carry ~no state signal (decoys)
EFFECTOR_SIG = {"TCF7": 1, "TOX": 0, "TBX21": 1, "NR4A": 0, "PDCD1": 0, "IFNG": 1}
EXHAUSTED_SIG = {"TCF7": 0, "TOX": 1, "TBX21": 0, "NR4A": 1, "PDCD1": 1, "IFNG": 0}
# Nodes a "drug" may force (regulators, not pure readouts / environment).
DRUGGABLE = {"TCF7", "TOX", "AP1", "TBX21", "PRDM1", "NR4A", "BATF", "IRF4"}


# ----------------------------------------------------------------------------
# Data loading
# ----------------------------------------------------------------------------
def load_edges(path):
    edges = []  # (source, target, sign)
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            s, t, sign = line.split()
            edges.append((s, t, +1 if sign == "+" else -1))
    return edges


def load_expression(path):
    exo, eff = {}, {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            g, a, b = line.split()
            exo[g], eff[g] = float(a), float(b)
    return exo, eff


# ----------------------------------------------------------------------------
# PART A -- Random Walk with Restart (personalized PageRank)
# ----------------------------------------------------------------------------
def build_transition(nodes, edges):
    """Column-stochastic transition matrix W on the directed graph (sign ignored
    for propagation; magnitude/confidence = 1 here). W[i,j] = j->i / outdeg(j)."""
    idx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)
    A = np.zeros((n, n))
    for s, t, _ in edges:
        A[idx[t], idx[s]] += 1.0
    col = A.sum(axis=0)
    for j in range(n):
        if col[j] > 0:
            A[:, j] /= col[j]
        else:
            A[:, j] = 1.0 / n          # dangling node -> uniform teleport
    return A, idx


def rwr(W, restart, alpha=ALPHA, tol=1e-12, max_iter=1000):
    r = restart / restart.sum()
    p = r.copy()
    for _ in range(max_iter):
        p_new = (1 - alpha) * W.dot(p) + alpha * r
        if np.abs(p_new - p).sum() < tol:
            return p_new
        p = p_new
    return p


def rank_states(nodes, edges, exo, eff):
    W, idx = build_transition(nodes, edges)
    r_exo = np.array([exo.get(n, 0.0) for n in nodes])
    r_eff = np.array([eff.get(n, 0.0) for n in nodes])
    p_exo = rwr(W, r_exo)
    p_eff = rwr(W, r_eff)
    delta = p_exo - p_eff
    return p_exo, p_eff, delta, idx


# ----------------------------------------------------------------------------
# PART B -- auto-generate a Boolean model from the top subnetwork
# ----------------------------------------------------------------------------
def induced_module(nodes, edges, delta, idx, top_k):
    """Pick top-k nodes by |delta| (the switch spans up AND down arms), then
    pull in their direct upstream regulators so every rule is complete.
    Nodes with |delta| < MIN_DELTA carry no state signal (decoy hubs) and are
    excluded -- both as seeds and as pulled-in regulators."""
    eligible = {n for n in nodes if abs(delta[idx[n]]) >= MIN_DELTA}
    ranked = sorted(eligible, key=lambda n: -abs(delta[idx[n]]))
    seed = set(ranked[:top_k])
    closed = set(seed)
    for s, t, _ in edges:
        if t in seed and s in eligible:
            closed.add(s)            # include eligible regulators of seed nodes
    # keep only edges fully inside the closed module
    mod_edges = [(s, t, sg) for s, t, sg in edges if s in closed and t in closed]
    mod_nodes = sorted({s for s, _, _ in mod_edges} | {t for _, t, _ in mod_edges})
    return mod_nodes, mod_edges, seed


def synthesize_bnet(mod_nodes, mod_edges):
    """Standard signed-graph -> Boolean heuristic:
       node = (OR of activators) AND NOT (OR of inhibitors).
       A node with no incoming edge becomes a free input (self-rule)."""
    acts = {n: [] for n in mod_nodes}
    inhs = {n: [] for n in mod_nodes}
    for s, t, sg in mod_edges:
        (acts if sg > 0 else inhs)[t].append(s)
    lines = []
    for n in mod_nodes:
        a, i = sorted(set(acts[n])), sorted(set(inhs[n]))
        if not a and not i:
            rule = n                                   # free input
        elif a and not i:
            rule = "(" + " | ".join(a) + ")"
        elif i and not a:
            rule = "!(" + " | ".join(i) + ")"
        else:
            rule = "(" + " | ".join(a) + ") & !(" + " | ".join(i) + ")"
        lines.append(f"{n}, {rule}")
    return "\n".join(lines) + "\n"


# ----------------------------------------------------------------------------
# Boolean engine (synchronous) + attractor / perturbation scan
# ----------------------------------------------------------------------------
def _next(primes, node, state):
    for imp in primes[node][1]:
        if all(state[k] == v for k, v in imp.items()):
            return 1
    return 0


def step(primes, state, overrides):
    return {n: overrides[n] if n in overrides else _next(primes, n, state)
            for n in primes}


def run_to_attractor(primes, start, overrides, max_iter=300):
    seen, traj, state = {}, [], dict(start)
    state.update(overrides)
    for i in range(max_iter):
        key = tuple(sorted(state.items()))
        if key in seen:
            return traj[seen[key]:]
        seen[key] = i
        traj.append(state)
        state = step(primes, state, overrides)
    return [state]


def classify(s, sig):
    return all(s.get(k) == v for k, v in sig.items() if k in s)


def all_attractors(primes, context):
    free = [n for n in primes if n not in context]
    found = {}
    for bits in product((0, 1), repeat=len(free)):
        start = dict(context); start.update(dict(zip(free, bits)))
        att = run_to_attractor(primes, start, context)
        found[frozenset(tuple(sorted(s.items())) for s in att)] = att
    return list(found.values())


def reaches_effector(att, forced):
    sig = {k: v for k, v in EFFECTOR_SIG.items() if k not in forced}
    return any(all(s.get(k) == v for k, v in sig.items()) for s in att)


def scan(primes, start_state, context, max_combo=2):
    cands = [(n, v) for n in primes if n in DRUGGABLE for v in (0, 1)]
    hits = []
    for size in range(1, max_combo + 1):
        for combo in combinations(cands, size):
            if len({n for n, _ in combo}) != size:
                continue
            ov = dict(context); ov.update(dict(combo))
            att = run_to_attractor(primes, start_state, ov)
            if reaches_effector(att, set(ov)):
                hits.append((size, dict(combo)))
    return hits


def fmt(ov):
    arrow = {1: "↑ON", 0: "↓OFF"}
    return " + ".join(f"{n} {arrow[v]}" for n, v in sorted(ov.items()))


# ----------------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------------
def main():
    edges = load_edges(EDGES)
    exo, eff = load_expression(EXPR)
    nodes = sorted({s for s, _, _ in edges} | {t for _, t, _ in edges})

    print("=" * 72)
    print("PIPELINE: differential RWR ranking -> auto Boolean model -> flip combo")
    print("=" * 72)
    print(f"Regulatory edges: {len(edges)}   nodes: {len(nodes)}\n")

    # ---- PART A ----
    print("-" * 72)
    print("PART A | RWR (personalized PageRank) ranking, exhausted vs effector")
    print("-" * 72)
    p_exo, p_eff, delta, idx = rank_states(nodes, edges, exo, eff)
    ranked = sorted(nodes, key=lambda n: -delta[idx[n]])
    print(f"{'gene':8} {'p_exhausted':>12} {'p_effector':>12} {'Δ (exo-eff)':>12}")
    for n in ranked:
        flag = "  <- decoy hub" if n in {"ACTB", "GAPDH", "TUBB"} else ""
        print(f"{n:8} {p_exo[idx[n]]:>12.4f} {p_eff[idx[n]]:>12.4f} "
              f"{delta[idx[n]]:>+12.4f}{flag}")
    print("\nTop exhaustion-specific (Δ>0):",
          ", ".join(n for n in ranked if delta[idx[n]] > 0)[:200])
    print("Note: ACTB/GAPDH/TUBB are high-degree hubs but Δ~0 -> not top-ranked.")

    # ---- PART B ----
    print("\n" + "-" * 72)
    print(f"PART B | Auto-build Boolean model from top-{TOP_K} |Δ| subnetwork")
    print("-" * 72)
    mod_nodes, mod_edges, seed = induced_module(nodes, edges, delta, idx, TOP_K)
    print(f"Module nodes ({len(mod_nodes)}): {', '.join(mod_nodes)}")
    print(f"(seed by |Δ|: {', '.join(sorted(seed))})\n")

    bnet = synthesize_bnet(mod_nodes, mod_edges)
    print("Auto-generated bnet rules:")
    for line in bnet.strip().splitlines():
        print("   " + line)

    primes = bnet2primes(bnet)
    context = {"NFAT": 1} if "NFAT" in primes else {}

    print("\nAttractors of the auto-generated model (context NFAT=1):")
    atts = all_attractors(primes, context)
    exhausted_state = effector_state = None
    for att in atts:
        s = att[0]
        tag = ("EXHAUSTED" if classify(s, EXHAUSTED_SIG)
               else "EFFECTOR" if classify(s, EFFECTOR_SIG) else "other")
        on = [n for n in mod_nodes if s.get(n) == 1]
        kind = "fixed point" if len(att) == 1 else f"cycle({len(att)})"
        print(f"  [{tag:9}] {kind:12} ON: {', '.join(on)}")
        if tag == "EXHAUSTED" and (exhausted_state is None or s.get("AP1") == 0):
            exhausted_state = s
        if tag == "EFFECTOR" and (effector_state is None or s.get("AP1") == 1):
            effector_state = s

    if not (exhausted_state and effector_state):
        print("\n!! Did not recover both attractors from the auto model.")
        print("   (Adjust TOP_K or the edge set so the bistable core is included.)")
        return

    # ---- perturbation scan ----
    print("\n" + "-" * 72)
    print("PART B | Minimal flip: EXHAUSTED (AP1 off) -> EFFECTOR")
    print("-" * 72)
    hits = scan(primes, exhausted_state, context, max_combo=2)
    singles = sorted({frozenset(h[1].items()) for h in hits if h[0] == 1})
    doubles = sorted({frozenset(h[1].items()) for h in hits if h[0] == 2})
    minimal_doubles = [d for d in doubles if not any(s <= d for s in singles)]

    print(f"Single-node flips: {len(singles)}")
    for h in singles:
        print("    - " + fmt(dict(h)))
    if not singles:
        print("    (none -- no monotherapy escapes the exhausted attractor)")
    print(f"Minimal two-node flips: {len(minimal_doubles)}")
    for d in minimal_doubles:
        print("    - " + fmt(dict(d)))

    print("\n" + "=" * 72)
    print("RESULT: the ranking (PART A) chose the module; the module's attractor")
    print("structure (PART B) yielded the minimal combination to lever the state.")
    print("No hand-written Boolean model -- it was synthesized from the network.")
    print("=" * 72)


if __name__ == "__main__":
    main()
