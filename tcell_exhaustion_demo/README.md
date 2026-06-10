# T cell exhaustion ↔ effector: minimal Boolean-network demo

A small, runnable illustration of the idea we discussed: a cell state is a
**stable attractor (steady state) of a regulatory network**, and the
"complex / combination needed to lever a signal" is the **minimal set of
forced nodes that moves the system from one attractor to another**.

## What it does

1. Encodes a ~14-node core CD8 T-cell circuit (TCF7, TOX, PDCD1, NR4A, TBX21,
   AP1/c-Jun, NFAT, BATF/IRF4, HAVCR2, LAG3, EOMES, IFNG, PRDM1) once as a
   PyBoolNet `bnet` model.
2. Computes the network's **attractors** — the homeostatic fixed points. Two are
   biologically labelled: an **effector/memory** state and an **exhausted** state.
3. Runs a **single- and double-node perturbation scan** (force genes ON =
   overexpression, OFF = knockdown) and reports the **minimal combination** that
   flips the canonical exhausted state into the effector program.

## Run

```bash
pip install pyboolnet          # parsing + steady-state cross-check
python3 tcell_exhaustion_boolean_demo.py
```

No NuSMV/graphviz needed — steady states use the bundled ASP solver, and the
perturbation engine is plain Python.

## Result (current model)

```
Single-node interventions that flip the state: 0   (no monotherapy works)
Minimal two-node combinations that flip the state:
    - AP1 ↑ON + TOX ↓OFF
    - AP1 ↑ON + TCF7 ↑ON
```

**Why two nodes, not one.** The exhausted attractor is held by two things at
once: the *self-sustaining TOX loop* and the *NFAT-without-AP-1* drive on NR4A.
A flip must (a) supply an effector drive — AP1/c-Jun up — **and** (b) release an
exhaustion brake — TOX down or TCF7 up. This mirrors published in-vivo rescue:
c-Jun overexpression (Lynn 2019) and TOX deletion (Khan/Scott 2019).

## How this maps to the earlier discussion

- **Steady-state matrix / eigenvalues** → here the attractors *are* the steady
  states; in a Boolean model stability is read off the state-transition
  structure rather than a Jacobian's eigenvalues, but the concept is the same:
  a basin you get pulled back into.
- **"What complex combination perturbs the signal"** → exactly the minimal
  forced-node sets printed by the scan. This is the discrete analogue of
  controllability / driver-node analysis.

## Full pipeline: ranking → auto-generated model (`pipeline.py`)

`tcell_exhaustion_boolean_demo.py` hand-writes the Boolean model. `pipeline.py`
instead **derives it from network ranking**, joining the two halves we discussed:

- **Part A — differential RWR.** Random-Walk-with-Restart (personalized
  PageRank) is run on the signed regulatory graph
  (`data/regulatory_edges.tsv`) with restart vectors from per-state activity
  (`data/state_expression.tsv`), once for the exhausted and once for the
  effector state. Proteins are ranked by Δ = p(exhausted) − p(effector).
  Decoy housekeeping hubs (ACTB/GAPDH/TUBB) are high-degree but Δ≈0, so they
  rank low — the demo's illustration of why differential ranking beats raw
  degree/PageRank.
- **Part B — auto-Boolean.** The top-|Δ| subnetwork (plus its direct upstream
  regulators) is converted to a Boolean model with no hand-writing
  (`node = OR(activators) AND NOT OR(inhibitors)`), its attractors are found,
  and the same flip-scan runs.

```bash
pip install numpy pyboolnet
python3 pipeline.py
```

Result: the model recovers both attractors, and the scan reports **AP1 ↑ON**
as a single-node flip.

### Why the auto-model says one node but the hand-model said two

This difference is the point, not a bug. The auto-synthesis turns *every*
inhibitory edge into a hard veto (`AND NOT`), so forcing AP1 on simultaneously
(a) drives TBX21/IFNG, (b) vetoes TOX, and (c) vetoes NR4A — one node releases
all three brakes. The hand model encoded a softer, self-sustaining TOX loop, so
it needed a second hit. **Take-away:** the ranking robustly picks the right
*module*, but the minimal flip-set depends on the rule-synthesis assumptions,
which must be validated against biology (dose, kinetics, partial knockdown)
before trusting them — this is why you don't fully automate the logic layer.
Both answers, incidentally, are literature-consistent: c-Jun (AP-1)
overexpression alone substantially rescues exhaustion (Lynn 2019).

## Extending it

- Swap in your own state's circuit (tumor EMT: ZEB1/SNAI1/miR-200/CDH1;
  inflammatory macrophage: NFKB/IRF5/STAT1 vs STAT6/KLF4). Only `BNET`,
  the marker signatures, and `PERTURBABLE` change.
- Raise `max_combo` to 3 to find triple combinations when no pair works.
- Replace synchronous update with asynchronous (PyBoolNet supports it) for a
  more realistic — but slower — attractor analysis.
- Feed the resting node values from scRNA-seq pseudobulk to *seed* the start
  state instead of hand-picking it.
