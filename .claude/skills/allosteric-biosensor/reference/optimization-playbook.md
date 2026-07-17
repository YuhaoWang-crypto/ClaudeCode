# Optimization playbook — how to maximize biosensor design success

The *way of thinking* distilled from the biosensor-engineering literature
(Guo/Baker Nat.Biotechnol. 2026; Fan 2008 & Binkowski 2011 cpFluc; Dixon 2016
NanoBiT; the chimera-design skill), turned into an explicit decision workflow
(`campaign.py`). The goal: given a **new analyte**, choose the design that has
the **highest probability of switching**, and triage in silico before the bench.

## 9 principles the papers actually optimize on

1. **Match the linear range to the analyte, not to the tightest binder**
   (Binkowski 2011). A sensor whose Kd sits at the *low end* of the analyte's
   physiological span saturates; one at the *high end* is insensitive. Target
   **Kd ≈ the operating concentration** so the linear range brackets it. The
   "better" receptor is the affinity-matched one — a lesson that flips the naive
   "maximize affinity" instinct.
2. **Small, rigid receptors win** (Guo/Baker). <150 aa ML/de-novo binders with
   no global conformational change → a **<10-variant** circular-permutation
   library suffices, and the single-component entropic switch works.
3. **Never permute the pocket** (chimera-design). Circular-permutation sites must
   be loops *away* from the ligand-contact residues, or binding dies.
4. **The insertion-site × linker grid IS the dynamic-range knob**
   (Fan/Binkowski/Guo). Longer Gly linker → higher basal activity & kcat but
   **lower** DR; Ser rigidification reverts. Don't guess — **scan** a small grid.
5. **Pick the architecture from analyte type + readout**, not habit:
   - small molecule, colorimetric → **domain insertion / TEM-1**
   - electrochemical (wearable, point-of-care) → **PQQ-GDH insertion**
   - luminescent / in-cell / low-abundance → **NanoLuc CP-fusion** or cpFluc
   - protein / two-epitope / PPI → **split complementation (NanoBiT)**
   - multiple inputs / conditional → **logic gate (ligand-gated modules)**
6. **Engineer fragment affinity weak** for split designs (Dixon). LgBiT/SmBiT
   ~190 µM so the *target* interaction, not the tags, drives signal; optimize
   stability on the big fragment, affinity on the small one — separately.
7. **Offset the circular-permutation affinity penalty** (Guo). CP costs ~4× Kd;
   recover it with an **auxiliary binding domain** (raises local ligand conc.)
   or avidity — especially when the receptor is already borderline-weak.
8. **Triage in silico before the bench.** Boltz holo/apo/control → fold +
   **binding retention** (chimera vs native ligand_iptm) + **active-site
   integrity** → rank; `coupling.py`/`md_entropy.py` for the entropic mechanism.
   **DR itself is never predicted** — it gates which few go to the bench.
9. **Focused, not exhaustive.** The small-binder payoff is a tiny library; spend
   the budget on the site×linker grid and on 2–3 replicas, not on hundreds of
   constructs.

## The workflow (`campaign.py`)

```
analyte spec ─▶ (1) triage: class, physiological range → target-Kd window,
                              candidate architectures
             ─▶ (2) receptor sourcing: mine PDB (discover.py) / de-novo / given;
                              check affinity match → auxiliary-domain flag
             ─▶ (3) architecture + reporter选择 (readout × analyte type)
             ─▶ (4) library: pocket-safe CP sites × linker grid  (<10 core)
             ─▶ (5) in-silico plan: Boltz holo/apo/control → retention/active-site
             ─▶ (6) optimize: site×linker grid on marginal variants
             ─▶ (7) SUCCESS SCORECARD: weighted readiness + specific fixes
```

## The success scorecard (⚠️ heuristic, transparent)

Each axis scored 0–1, weighted, → a **readiness** score with *specific*
recommendations (not a black box). It encodes the principles above:

| axis | good (→1) | fires recommendation |
|---|---|---|
| receptor size/rigidity | small, single-domain | "receptor large — expect harder switch" |
| **affinity match** | Kd within/near physiological range | "too tight → saturation, use weaker variant"; "too weak → add auxiliary domain" |
| pocket-safe sites | ≥3 loop sites away from pocket | "few safe sites — widen loop search / de-novo redesign" |
| binding retention (Boltz) | ≥0.85 vs native | "retention low — scan sites/linker or graft elsewhere" |
| active-site integrity | intact constellation | "active site perturbed — try another insertion loop" |
| architecture fit | matches readout+class | "readout mismatch — switch reporter/topology" |
| DR-tunability | site×linker grid available | (always available here) |

Readiness is **not** a probability of success and **not** a dynamic range — it
is a transparent prioritization that tells you *what to fix first* and *whether
to spend bench effort*. Ground truth remains the wet-lab kobs/luminescence
titration.
