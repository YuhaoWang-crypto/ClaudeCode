"""
Module 2 — Chemical Reaction Network Theory (CRNT) deficiency and switching.

Idea (rigorous):
  For a mass-action reaction network the DEFICIENCY is a non-negative
  integer topological invariant
        delta = n - l - s
  n = number of complexes, l = number of linkage classes,
  s = rank of the stoichiometric matrix.
  Feinberg's Deficiency-Zero Theorem: if delta = 0 and the network is
  weakly reversible, then for EVERY choice of positive rate constants the
  system has exactly one positive steady state and it is stable -> the
  network CANNOT be a switch (no bistability), no matter the kinetics.
  Bistability (a genuine on/off switch) therefore requires delta >= 1.

We compute delta mechanically for two canonical mass-action networks and
then CONFIRM the prediction by numerically integrating the mass-action ODEs
from many initial conditions and counting stable attractors.
"""
import numpy as np
from scipy.integrate import solve_ivp


# ---------- deficiency from a complex-based specification ----------------
# A complex is a frozenset-free dict {species: stoichiometric coefficient}.
# A reaction is (reactant_complex, product_complex).

def _complex_vec(cx, species):
    return np.array([cx.get(s, 0) for s in species], dtype=float)


def deficiency(reactions, species, chemostat=()):
    """Deficiency of the EFFECTIVE network. Chemostatted (reservoir) species
    are held constant, so they fold into the rate constants and must be
    stripped from every complex before counting complexes/linkage classes."""
    chemostat = set(chemostat)

    def strip(cx):
        return {s: c for s, c in cx.items() if s not in chemostat}

    reactions = [(strip(r), strip(p)) for r, p in reactions]
    species = [s for s in species if s not in chemostat]

    complexes = []
    def cx_key(cx):
        return tuple(sorted(cx.items()))
    index = {}
    for r, p in reactions:
        for cx in (r, p):
            k = cx_key(cx)
            if k not in index:
                index[k] = len(complexes)
                complexes.append(cx)
    n = len(complexes)

    # linkage classes: connected components of the complex graph
    import networkx as nx
    Cg = nx.Graph()
    Cg.add_nodes_from(range(n))
    for r, p in reactions:
        Cg.add_edge(index[cx_key(r)], index[cx_key(p)])
    l = nx.number_connected_components(Cg)

    # stoichiometric matrix rank
    cols = []
    for r, p in reactions:
        cols.append(_complex_vec(p, species) - _complex_vec(r, species))
    S = np.array(cols).T if cols else np.zeros((len(species), 0))
    s = int(np.linalg.matrix_rank(S)) if S.size else 0

    return {"n": n, "l": l, "s": s, "delta": n - l - s,
            "complexes": complexes, "S": S}


# ---------- mass-action ODE simulation to count attractors ---------------

def massaction_rhs(reactions, species, rates, chemostat):
    """Build dc/dt for dynamic (non-chemostatted) species."""
    dyn = [sp for sp in species if sp not in chemostat]
    dyn_idx = {sp: i for i, sp in enumerate(dyn)}

    def rhs(t, y):
        conc = dict(chemostat)
        for sp, i in dyn_idx.items():
            conc[sp] = max(y[i], 0.0)
        dydt = np.zeros(len(dyn))
        for (r, p), k in zip(reactions, rates):
            rate = k
            for sp, coef in r.items():
                rate *= conc.get(sp, 0.0) ** coef
            for sp, coef in {**r, **{}}.items():
                pass
            net = {}
            for sp, coef in p.items():
                net[sp] = net.get(sp, 0) + coef
            for sp, coef in r.items():
                net[sp] = net.get(sp, 0) - coef
            for sp, coef in net.items():
                if sp in dyn_idx:
                    dydt[dyn_idx[sp]] += coef * rate
        return dydt

    return rhs, dyn


def count_attractors(reactions, species, rates, chemostat,
                     n_starts=60, cmax=6.0, seed=0):
    rhs, dyn = massaction_rhs(reactions, species, rates, chemostat)
    rng = np.random.default_rng(seed)
    fixed = []
    for _ in range(n_starts):
        y0 = rng.uniform(0, cmax, size=len(dyn))
        sol = solve_ivp(rhs, (0, 4000), y0, method="LSODA",
                        rtol=1e-8, atol=1e-10)
        yf = np.clip(sol.y[:, -1], 0, None)
        if not any(np.allclose(yf, f, atol=1e-3) for f in fixed):
            fixed.append(yf)
    # merge near-duplicates
    merged = []
    for f in fixed:
        if not any(np.linalg.norm(f - m) < 1e-2 for m in merged):
            merged.append(f)
    return dyn, merged


# ------------------------------ networks ---------------------------------

def net_deficiency_zero():
    """Weakly reversible, delta = 0: A<->B<->C<->A. Must be monostable."""
    A, B, C = ({"A": 1}, {"B": 1}, {"C": 1})
    reactions = [(A, B), (B, A), (B, C), (C, B), (C, A), (A, C)]
    species = ["A", "B", "C"]
    rates = [1.0, 1.0, 1.2, 0.8, 0.9, 1.1]
    chemostat = {}
    # conservation: total A+B+C fixed; count attractors on the simplex
    return reactions, species, rates, chemostat


def net_schlogl():
    """Schlogl model — canonical 1-species bistable mass-action network.
        A + 2X -> 3X ,  3X -> A + 2X ,  X -> B ,  B -> X
    A and B are chemostatted (reservoirs)."""
    reactions = [
        ({"A": 1, "X": 2}, {"X": 3}),   # A+2X -> 3X
        ({"X": 3}, {"A": 1, "X": 2}),   # 3X -> A+2X
        ({"X": 1}, {"B": 1}),           # X -> B
        ({"B": 1}, {"X": 1}),           # B -> X
    ]
    species = ["A", "B", "X"]
    # rate constants tuned into the classic bistable regime
    rates = [3.0, 0.6, 3.5, 1.0]
    chemostat = {"A": 1.0, "B": 1.0}
    return reactions, species, rates, chemostat


def report():
    print("=" * 68)
    print("MODULE 2 — CRNT deficiency  ->  switching capacity (bistability)")
    print("=" * 68)

    for name, factory, conserve in [
        ("Deficiency-zero  A<->B<->C<->A", net_deficiency_zero, True),
        ("Schlogl bistable network       ", net_schlogl, False),
    ]:
        reactions, species, rates, chemostat = factory()
        d = deficiency(reactions, species, chemostat=chemostat)
        dyn = species  # default
        if name.startswith("Deficiency-zero"):
            # simulate on a fixed conservation class (total = 3)
            dyn, attr = _count_conserved(reactions, species, rates, total=3.0)
        else:
            dyn, attr = count_attractors(reactions, species, rates, chemostat)
        print(f"\n{name}")
        print(f"  complexes n={d['n']}  linkage classes l={d['l']}  "
              f"rank s={d['s']}  =>  DEFICIENCY delta = {d['delta']}")
        print(f"  stable attractors found by simulation: {len(attr)} "
              f"({'MONOSTABLE — no switch' if len(attr)==1 else 'BISTABLE — genuine switch'})")
        for a in attr:
            print("      steady state  " +
                  ", ".join(f"{s}={v:.3f}" for s, v in zip(dyn, a)))
    print("\nInterpretation: delta=0 (weakly reversible) forbids bistability")
    print("for ALL rate constants (Feinberg). A switch REQUIRES delta>=1,")
    print("which the Schlogl network has — topology gates switching capacity.")


def _count_conserved(reactions, species, rates, total):
    """Simulate the closed A/B/C network on the simplex A+B+C=total."""
    rhs, dyn = massaction_rhs(reactions, species, rates, chemostat={})
    rng = np.random.default_rng(1)
    attr = []
    for _ in range(40):
        w = rng.uniform(0, 1, size=3)
        y0 = total * w / w.sum()
        sol = solve_ivp(rhs, (0, 3000), y0, method="LSODA",
                        rtol=1e-9, atol=1e-11)
        yf = np.clip(sol.y[:, -1], 0, None)
        if not any(np.linalg.norm(yf - a) < 1e-3 for a in attr):
            attr.append(yf)
    return dyn, attr


if __name__ == "__main__":
    report()
