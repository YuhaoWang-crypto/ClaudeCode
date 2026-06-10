"""
Boolean-network library: rule synthesis (veto vs majority), synchronous and
asynchronous dynamics, attractor finding, and a perturbation scan.

Used by pipeline.py. Kept dependency-light: only networkx (for async SCCs).

State = dict {node: 0/1}. A "rule" for a node is one of:
    ("input", None)          free input -> keeps its value
    ("veto",  (acts, inhs))  next = OR(acts) AND NOT OR(inhs)   (inhibitor = hard veto)
    ("maj",   (acts, inhs))  next = sign(sum(acts) - sum(inhs)) (soft threshold; ties hold)
"""

from itertools import combinations, product
import networkx as nx


# --------------------------------------------------------------------------
# Rule synthesis from a signed edge list
# --------------------------------------------------------------------------
def synthesize_rules(nodes, edges, synthesis="veto"):
    """edges: list of (source, target, sign in {+1,-1}). Returns {node: rule}."""
    acts = {n: [] for n in nodes}
    inhs = {n: [] for n in nodes}
    for s, t, sg in edges:
        if s in nodes and t in nodes:
            (acts if sg > 0 else inhs)[t].append(s)
    rules = {}
    for n in nodes:
        a, i = sorted(set(acts[n])), sorted(set(inhs[n]))
        if not a and not i:
            rules[n] = ("input", None)
        else:
            rules[n] = (("veto" if synthesis == "veto" else "maj"), (a, i))
    return rules


def rule_text(node, rule):
    kind, data = rule
    if kind == "input":
        return f"{node}  :=  input (held)"
    a, i = data
    if kind == "veto":
        if a and i:
            return f"{node}  :=  ({' | '.join(a)}) & !({' | '.join(i)})"
        if a:
            return f"{node}  :=  ({' | '.join(a)})"
        return f"{node}  :=  !({' | '.join(i)})"
    # majority
    act = " + ".join(a) if a else "0"
    inh = " + ".join(i) if i else "0"
    base = " + 1(baseline)" if (i and not a) else ""
    return f"{node}  :=  sign[({act}){base} - ({inh})] > 0   (soft)"


# --------------------------------------------------------------------------
# Update functions
# --------------------------------------------------------------------------
def node_next(rules, node, state, overrides):
    if node in overrides:
        return overrides[node]
    kind, data = rules[node]
    if kind == "input":
        return state[node]
    a, i = data
    if kind == "veto":
        if a and i:
            return int(any(state[x] for x in a) and not any(state[x] for x in i))
        if a:
            return int(any(state[x] for x in a))
        return int(not any(state[x] for x in i))
    # majority / soft threshold; inhibitor-only nodes get a +1 constitutive baseline
    net = sum(state[x] for x in a) - sum(state[x] for x in i)
    if i and not a:
        net += 1
    if net > 0:
        return 1
    if net < 0:
        return 0
    # tie: hysteresis (hold) only for genuinely bistable mixed regulation;
    # pure-activator / inhibitor-only nodes resolve deterministically to OFF,
    # otherwise readout nodes latch and the landscape fragments artefactually.
    return state[node] if (a and i) else 0


def sync_successor(rules, state, overrides):
    return {n: node_next(rules, n, state, overrides) for n in rules}


def async_successors(rules, state, overrides):
    """All states reachable by updating exactly one node that wants to change."""
    out = []
    for n in rules:
        v = node_next(rules, n, state, overrides)
        if v != state[n]:
            nxt = dict(state)
            nxt[n] = v
            out.append(nxt)
    return out


# --------------------------------------------------------------------------
# Attractor finding
# --------------------------------------------------------------------------
def _key(state, order):
    return tuple(state[n] for n in order)


def sync_attractor_from(rules, start, overrides, max_iter=500):
    """Deterministic: follow synchronous trajectory to its cycle/fixed point."""
    order = list(rules)
    seen, traj, state = {}, [], dict(start)
    for n, v in overrides.items():
        state[n] = v
    for i in range(max_iter):
        k = _key(state, order)
        if k in seen:
            return traj[seen[k]:]
        seen[k] = i
        traj.append(state)
        state = sync_successor(rules, state, overrides)
    return [state]


def build_async_stg(rules, context, overrides):
    """Async state-transition graph over the free-node space (context+overrides pinned)."""
    order = list(rules)
    fixed = dict(context)
    fixed.update(overrides)
    free = [n for n in order if n not in fixed]
    g = nx.DiGraph()
    for bits in product((0, 1), repeat=len(free)):
        state = dict(fixed)
        state.update(dict(zip(free, bits)))
        k = _key(state, order)
        g.add_node(k)
        succ = async_successors(rules, state, overrides)
        if not succ:
            g.add_edge(k, k)                       # fixed point
        for s in succ:
            g.add_edge(k, _key(s, order))
    return g, order


def async_attractors(rules, context, overrides):
    """Attractors = terminal SCCs of the async STG."""
    g, order = build_async_stg(rules, context, overrides)
    cond = nx.condensation(g)
    atts = []
    for scc_id in cond.nodes:
        if cond.out_degree(scc_id) == 0:
            states = [dict(zip(order, k)) for k in cond.nodes[scc_id]["members"]]
            atts.append(states)
    return atts


def all_attractors(rules, context, update="sync"):
    if update == "async":
        return async_attractors(rules, context, {})
    order = list(rules)
    free = [n for n in order if n not in context]
    found = {}
    for bits in product((0, 1), repeat=len(free)):
        start = dict(context)
        start.update(dict(zip(free, bits)))
        att = sync_attractor_from(rules, start, context)
        found[frozenset(_key(s, order) for s in att)] = att
    return list(found.values())


# --------------------------------------------------------------------------
# Reachability + perturbation scan
# --------------------------------------------------------------------------
def reachable_attractors(rules, start, overrides, update="sync"):
    """Attractors reachable from `start` under `overrides`."""
    if update == "sync":
        return [sync_attractor_from(rules, start, overrides)]
    # async: BFS the STG from start, return terminal SCCs inside the reachable set
    order = list(rules)
    g, _ = build_async_stg(rules, {}, overrides)
    s0 = _key({**start, **overrides}, order)
    if s0 not in g:
        return []
    reach = nx.descendants(g, s0) | {s0}
    sub = g.subgraph(reach)
    cond = nx.condensation(sub)
    atts = []
    for scc_id in cond.nodes:
        if cond.out_degree(scc_id) == 0:
            atts.append([dict(zip(order, k)) for k in cond.nodes[scc_id]["members"]])
    return atts


def matches(attractor, signature):
    """True if some state in the attractor matches all signature entries present."""
    return any(all(s.get(k) == v for k, v in signature.items() if k in s)
               for s in attractor)


def all_reach_effector(attractors, effector_sig, forced):
    """All reachable attractors commit to the effector program (ignoring forced nodes)."""
    sig = {k: v for k, v in effector_sig.items() if k not in forced}
    return bool(attractors) and all(matches(a, sig) for a in attractors)


def scan_flips(rules, start_state, context, druggable, effector_sig,
               update="sync", max_combo=2):
    """Smallest forced-node sets that drive start_state to commit to effector."""
    cands = [(n, v) for n in rules if n in druggable for v in (0, 1)]
    hits = []
    for size in range(1, max_combo + 1):
        for combo in combinations(cands, size):
            if len({n for n, _ in combo}) != size:
                continue
            ov = dict(context)
            ov.update(dict(combo))
            atts = reachable_attractors(rules, start_state, ov, update=update)
            if all_reach_effector(atts, effector_sig, set(ov)):
                hits.append((size, dict(combo)))
    return hits


def fmt_combo(ov):
    arrow = {1: "↑ON", 0: "↓OFF"}
    return " + ".join(f"{n} {arrow[v]}" for n, v in sorted(ov.items()))
