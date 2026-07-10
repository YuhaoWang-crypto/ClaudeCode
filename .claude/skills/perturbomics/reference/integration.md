# Integration — perturbomics × network-biomarker × drug-MCP

How this skill composes with the **pathway/dynamics analysis** and **drug-
discovery** capabilities into one multi-scale funnel. The three already share a
vocabulary — a **target (gene)** and a **perturbation magnitude** — so they plug
together with a thin bridge (`perturbomics/integrate.py`), no rewrites.

## The funnel

```
                perturbomics                    network-biomarker              drug-MCP servers
                (transcriptome)                 (pathway dynamics)             (structure / clinic)
 disease  ─►  rank_reversers / best_combinations
 signature     │  ranked drugs + gene targets
               │  + drug+CRISPR combinations
               ▼
          ┌───────────────────────────────────────────────────────────────────────────┐
          │  integrate.integrated_leads(reversers, network, evidence)                   │
          │                                                                             │
          │  axis 1 TRANSCRIPTOMIC  |WTCS| reversal                 ← perturbomics       │
          │  axis 2 NETWORK CONTROL core / switch node?             ← M1/M11, M2/M19     │
          │  axis 3 ENGAGEABILITY   ChEMBL pIC50 + Boltz bind/ADME  ← ChEMBL/Boltz  (m6) │
          │  axis 4 CLINICAL        trial phase / repurposing       ← ClinicalTrials (m8)│
          └───────────────────────────────────────────────────────────────────────────┘
               │  integrated_score (weighted, per-row renormalised)
               ▼
          prioritised leads  +  monitoring biomarker (M4 DNB)  ─►  what to give AND how to tell it works
```

Each axis is an **independent experiment**. A connectivity outlier that is also
a network control node, druggable, and clinically plausible is far more credible
than a top-|WTCS| hit alone — and, as the demo shows, the fusion routinely
*demotes* a strong-but-promiscuous transcriptomic hit and *promotes* a moderate
one with a good target.

## Wiring each axis to the real skills

### Axis 2 — network control (network-biomarker)
That skill's modules emit exactly what `NetworkContext` needs:
- **`m1_symmetry` / `m11_fibration`** → `quotient`/orbits ⇒ the **irreducible
  core** nodes → `NetworkContext.core_nodes`.
- **`m2_crnt` (deficiency>0) / `m19_switch_library`** → **bistable / switch**
  nodes → `NetworkContext.switch_nodes`.
- **`m4_dnb_lyapunov`** → the **early-warning biomarker** (leading-eigenvalue /
  DNB readout) → `NetworkContext.biomarker` — a gene to *monitor a responder*,
  not a target to hit. Retrieve with `monitoring_biomarker(net)`.

```python
from grn_pipeline import m1_symmetry, m19_switch_library, m4_dnb_lyapunov
core  = set(m1_symmetry.report()["core_nodes"])
switch= set(m19_switch_library.report()["switch_nodes"])
bmk   = tuple(m4_dnb_lyapunov.report()["biomarker_genes"])
net = NetworkContext(core_nodes=core, switch_nodes=switch, biomarker=bmk)
```

(The exact result keys vary by module; read each module's `report()` — the point
is that core/switch/biomarker sets are what cross the boundary, as plain gene
ids, so perturbomics never imports grn_pipeline.)

### Axes 3 & 4 — engageability + clinical (the drug-discovery pipeline, MCP)
Populate one `DrugEvidence` per nominated compound from the same MCP servers the
network-biomarker drug chain (m6–m9) already uses:
- **ChEMBL** `get_bioactivity` → `pIC50`; `get_mechanism`/`target_search` →
  `target`; `get_admet` → `admet_ok`.
- **Boltz** `start_structure_and_binding` → `boltz_bind` (binding confidence)
  and its ADME estimates. (This is m6/m7's exact input.)
- **ClinicalTrials** `search_trials` → `clinical_phase`; or the Drug Repurposing
  Hub `clinical_phase` column for repurposing candidates.

```python
ev = {}
for name, target in nominated:               # top perturbomics reversers
    act  = chembl.get_bioactivity(target=target)   # -> pIC50
    bind = boltz.get_structure_and_binding(...)     # -> boltz_bind, admet
    tri  = clinicaltrials.search_trials(intervention=name)  # -> phase
    ev[name] = DrugEvidence(name, target=target, pIC50=..., boltz_bind=...,
                            admet_ok=..., clinical_phase=...)
leads = integrated_leads(reversers, network=net, evidence=ev)
```

Missing axes are **excluded and flagged** (`provenance` column), never zeroed —
a genetic CRISPR hit with no compound simply scores on the axes it has.

## The reverse direction (this makes the whole loop close)
network-biomarker's **m6_integrate** takes a drug's engagement E, maps it to a
network perturbation μ, and reads the **M4 stability biomarker**. So a compound
nominated by perturbomics can be pushed *back through* the dynamics to predict
**how far it moves the disease network toward its tipping point** — and which
early-warning signal should change if it works. perturbomics says *what*,
network-biomarker says *where/how much + how to measure*, the drug servers say
*whether the molecule is real*.

## Run it

```bash
python3 -m perturbomics.demo_integrate      # offline 4-axis funnel, deterministic
python3 examples/real_leads_ipf.py          # LIVE ChEMBL+ClinicalTrials worked run
```

Then swap the synthetic `NetworkContext` for real network-biomarker output and
the illustrative `DrugEvidence` for live ChEMBL/Boltz/ClinicalTrials calls.

## A real worked run (IPF, live ChEMBL + ClinicalTrials)

`examples/real_leads_ipf.py` takes the real IPF drug reversers, fetches their
targets/potency/phase live, and fuses. The fetched evidence (this session):

| drug | ChEMBL target | best IC50 → pIC50 | Ro5 viol | phase | CT.gov |
|---|---|---|---|---|---|
| canertinib | **EGFR** (+ERBB2/4) | 1.5 nM → 8.82 | 0 | 3 | 3 |
| pevonedistat (MLN4924) | NAE / NEDD8 | 4.7 nM → 8.33 | 0 | 3 | 41 |
| trichostatin A | HDAC | 1.4 nM → 8.84 | 0 | 1 | — |
| procaterol | ADRB2 (β2) | agonist | 0 | 3 | — |

Fused result — the funnel's whole point, on real data:

```
name             target  net_role   transcriptomic  network  engagement  clinical  INTEGRATED
canertinib       EGFR    core            0.950        0.60      0.803      0.75      0.803  ◄ #1
pevonedistat     NAE1    peripheral      0.991        0.00      0.722      0.75      0.653
procaterol       ADRB2   peripheral      1.000        0.00        —        0.75      0.641
trichostatin A   HDAC1   peripheral      0.996        0.00      0.807      0.25      0.597
```

`procaterol` **wins the transcriptomic axis alone** (−1.000) but is an off-
pathway β2-agonist → 3rd; `trichostatin`/`pevonedistat` are the most potent but
their targets aren't network-control nodes → mid-pack; **`canertinib` rises to
#1 despite weaker connectivity** because EGFR is on-pathway (EGFR/ErbB signalling
is genuinely implicated in IPF) AND it is potent + phase-3. No single axis ranks
it first.

Caveats kept explicit in the script: the network axis is an **⚠️ illustrative**
IPF core/switch set (not a computed M1 quotient — swap in `m1_symmetry.report()`
on an IPF GRN to make it rigorous), and the **Boltz** binding job (async, minutes)
was not run this session, so engagement uses ChEMBL potency alone.

## Rigor discipline (unchanged, carried across the boundary)
- ✅ **rigorous**: every |WTCS|, core/switch membership *given a network*, the
  ChEMBL/Boltz numbers, the trial phase, and the weighted sum.
- ⚠️ **hypothesis**: "sits at a control node ⇒ better target", the E→μ→effect
  map (explicitly illustrative in m6), and any therapeutic/synergy claim.
Keep the two labels separate in any integrated report, exactly as each
contributing skill already does.
