#!/usr/bin/env python3
"""
End-to-end: network ranking -> auto-generated Boolean model -> flip combo.

PART A  Random-Walk-with-Restart (personalized PageRank) on a signed regulatory
        network, run for the EXHAUSTED and EFFECTOR states, ranked by the
        differential score. -> "which proteins matter".
PART B  The top differential subnetwork is auto-converted to a Boolean model
        (no hand-writing), its attractors are found, and a perturbation scan
        reports the minimal combination that flips exhausted -> effector.

This version is configurable, so you can watch the flip-set change with the
modelling assumptions:
    --edges / --expr   input files (default: toy; use *_real.tsv for real data)
    --update           sync | async        (synchronous vs asynchronous dynamics)
    --synthesis        veto | majority      (hard AND-NOT vs soft threshold logic)
    --max-combo        N                    (search up to N-node perturbations)

Examples:
    python3 pipeline.py
    python3 pipeline.py --edges data/signed_edges_real.tsv --expr data/state_expression_real.tsv
    python3 pipeline.py --synthesis majority --max-combo 3
    python3 pipeline.py --update async

Deps: numpy, networkx  (pyboolnet not required by this script).
"""

import argparse
import os
import numpy as np

import bnlib

HERE = os.path.dirname(os.path.abspath(__file__))
ALPHA = 0.20
TOP_K = 10
MIN_DELTA = 0.004
EFFECTOR_SIG = {"TCF7": 1, "TOX": 0, "TBX21": 1, "NR4A": 0, "PDCD1": 0, "IFNG": 1}
EXHAUSTED_SIG = {"TCF7": 0, "TOX": 1, "TBX21": 0, "NR4A": 1, "PDCD1": 1, "IFNG": 0}
DRUGGABLE = {"TCF7", "TOX", "AP1", "TBX21", "PRDM1", "NR4A", "BATF", "IRF4"}


def load_edges(path):
    edges = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            s, t, sign = parts[0], parts[1], parts[2]
            edges.append((s, t, +1 if sign == "+" else -1))
    return edges


def load_expression(path):
    exo, eff = {}, {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            g, a, b = line.split()[:3]
            exo[g], eff[g] = float(a), float(b)
    return exo, eff


# ---------------- PART A: RWR -----------------------------------------------
def rwr_rank(nodes, edges, exo, eff, alpha=ALPHA):
    idx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)
    A = np.zeros((n, n))
    for s, t, _ in edges:
        A[idx[t], idx[s]] += 1.0
    col = A.sum(axis=0)
    for j in range(n):
        A[:, j] = A[:, j] / col[j] if col[j] > 0 else 1.0 / n

    def walk(restart):
        r = restart / restart.sum()
        p = r.copy()
        for _ in range(2000):
            p2 = (1 - alpha) * A.dot(p) + alpha * r
            if np.abs(p2 - p).sum() < 1e-12:
                return p2
            p = p2
        return p

    r_exo = np.array([exo.get(x, 0.0) for x in nodes])
    r_eff = np.array([eff.get(x, 0.0) for x in nodes])
    return walk(r_exo), walk(r_eff), idx


# ---------------- PART B: module -------------------------------------------
def induced_module(nodes, edges, delta, idx, top_k):
    eligible = {x for x in nodes if abs(delta[idx[x]]) >= MIN_DELTA}
    ranked = sorted(eligible, key=lambda x: -abs(delta[idx[x]]))
    seed = set(ranked[:top_k])
    closed = set(seed)
    for s, t, _ in edges:
        if t in seed and s in eligible:
            closed.add(s)
    mod_edges = [(s, t, sg) for s, t, sg in edges if s in closed and t in closed]
    mod_nodes = sorted({s for s, _, _ in mod_edges} | {t for _, t, _ in mod_edges})
    return mod_nodes, mod_edges, seed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--edges", default=os.path.join(HERE, "data", "regulatory_edges.tsv"))
    ap.add_argument("--expr", default=os.path.join(HERE, "data", "state_expression.tsv"))
    ap.add_argument("--update", choices=["sync", "async"], default="sync")
    ap.add_argument("--synthesis", choices=["veto", "majority"], default="veto")
    ap.add_argument("--max-combo", type=int, default=2)
    ap.add_argument("--top-k", type=int, default=TOP_K)
    args = ap.parse_args()

    edges = load_edges(args.edges)
    exo, eff = load_expression(args.expr)
    nodes = sorted({s for s, _, _ in edges} | {t for _, t, _ in edges})

    print("=" * 74)
    print("PIPELINE | differential RWR -> auto Boolean model -> flip combo")
    print(f"edges={os.path.basename(args.edges)}  expr={os.path.basename(args.expr)}  "
          f"update={args.update}  synthesis={args.synthesis}  max_combo={args.max_combo}")
    print("=" * 74)
    print(f"regulatory edges: {len(edges)}   nodes: {len(nodes)}")

    # PART A
    p_exo, p_eff, idx = rwr_rank(nodes, edges, exo, eff)
    delta = p_exo - p_eff
    ranked = sorted(nodes, key=lambda x: -delta[idx[x]])
    print("\n-- PART A | RWR ranking (Δ = p_exhausted - p_effector) --")
    decoys = {"ACTB", "GAPDH", "TUBB"}
    for x in ranked:
        tag = "  <- decoy hub" if x in decoys else ""
        print(f"  {x:8} p_exo={p_exo[idx[x]]:.4f}  p_eff={p_eff[idx[x]]:.4f}  "
              f"Δ={delta[idx[x]]:+.4f}{tag}")

    # PART B
    mod_nodes, mod_edges, seed = induced_module(nodes, edges, delta, idx, args.top_k)
    rules = bnlib.synthesize_rules(mod_nodes, mod_edges, synthesis=args.synthesis)
    print(f"\n-- PART B | auto Boolean model from top-{args.top_k} |Δ| module "
          f"({len(mod_nodes)} nodes) --")
    for n in mod_nodes:
        print("   " + bnlib.rule_text(n, rules[n]))

    context = {"NFAT": 1} if "NFAT" in rules else {}
    atts = bnlib.all_attractors(rules, context, update=args.update)
    counts = {"EXHAUSTED": 0, "EFFECTOR": 0, "other": 0}
    exhausted_state = effector_state = None
    for att in atts:
        s = att[0]
        tag = ("EXHAUSTED" if bnlib.matches(att, EXHAUSTED_SIG)
               else "EFFECTOR" if bnlib.matches(att, EFFECTOR_SIG) else "other")
        counts[tag] += 1
        if tag == "EXHAUSTED" and (exhausted_state is None or s.get("AP1") == 0):
            exhausted_state = s
        if tag == "EFFECTOR" and (effector_state is None or s.get("AP1") == 1):
            effector_state = s
    print(f"\nAttractors ({args.update}, context={context}): {len(atts)} total  "
          f"[EXHAUSTED={counts['EXHAUSTED']}, EFFECTOR={counts['EFFECTOR']}, "
          f"other={counts['other']}]")
    for label, st in [("EXHAUSTED", exhausted_state), ("EFFECTOR", effector_state)]:
        if st:
            on = [n for n in mod_nodes if st.get(n) == 1]
            print(f"  canonical [{label:9}] ON: {', '.join(on)}")

    if not (exhausted_state and effector_state):
        print("\n!! Did not recover both attractors under these settings; no scan.")
        print("   (Real-data note: public regulon DBs are incomplete; try a")
        print("    different --synthesis/--update or enrich the edge set.)")
        return

    # perturbation scan
    print(f"\n-- PART B | minimal flip EXHAUSTED -> EFFECTOR "
          f"({args.update}, {args.synthesis}, ≤{args.max_combo} nodes) --")
    hits = bnlib.scan_flips(rules, exhausted_state, context, DRUGGABLE,
                            EFFECTOR_SIG, update=args.update, max_combo=args.max_combo)
    by_size = {}
    for size, ov in hits:
        by_size.setdefault(size, set()).add(frozenset(ov.items()))
    minimal = []
    for size in sorted(by_size):
        for combo in sorted(by_size[size]):
            if not any(sm <= combo for s2 in range(1, size) for sm in by_size.get(s2, [])):
                minimal.append(combo)
    if not minimal:
        print(f"  none up to {args.max_combo} nodes.")
    else:
        seen_sizes = sorted({len(c) for c in minimal})
        for size in seen_sizes:
            combos = [c for c in minimal if len(c) == size]
            print(f"  {size}-node:")
            for c in combos:
                print("    - " + bnlib.fmt_combo(dict(c)))


if __name__ == "__main__":
    main()
