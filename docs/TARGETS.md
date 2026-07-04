# Target systems — rationale, structures, wet-lab notes

All PDB files are in `data/pdb/`; canonical sequences in `data/fasta/all_targets.fasta`.
UniProt IDs and references are in `data/targets.json`.

---

## 1. GCK — human glucokinase (P35557) · monomer · **positive control #1**
- **Pair:** `1V4T` (super-open / inactive, apo) → `1V4S` (closed / active, glucose + synthetic allosteric activator).
- **Why:** single-domain-motion positive control with clear open/closed transition and a
  documented small-molecule **allosteric activator site**. Ground truth is explicit (Kamata et al. 2004).
- **Demo:** align inactive vs active, compute domain motion / Cα displacement / pocket opening,
  check whether the workflow ranks the activator pocket in the top regions. **It does** (regions
  65–77 and 446–461 bracket the activator site; small-domain 85–147 is the dominant motion).
- **Wet-lab:** human GCK expression/stability + coupled kinase (G6PDH) assay need optimization —
  harder than *E. coli* enzymes. Best used as the computational flagship.

## 2. PTP1B — protein tyrosine phosphatase 1B (P18031) · monomer · **positive control #2**
- **Pair:** `1SUG` (apo, WPD open) → `1T49` (allosteric-inhibitor bound); `1PTY` (orthosteric / WPD closed) as a third state.
- **Why:** known **distal allosteric inhibitor site ~20 Å** from catalytic Cys215; inhibitor restricts
  WPD/catalytic-loop mobility (Wiesmann et al. 2004). Drug-discovery relevant.
- **Demo:** recover the allosteric site from active-site state + WPD-loop motion + distal pocket coupling.
  Displacements are subtle (overall RMSD < 1 Å) → weight **loop dynamics / contact-network changes**, not raw max displacement.
- **Wet-lab:** monomer, simple **pNPP** or fluorogenic (DiFMUP) assay. Excellent allosteric-pocket drug-discovery demo.

## 3. AdK — *E. coli* adenylate kinase (P69441) · monomer · **engineering target #1 (Zn-switch)**
- **Pair:** `4AKE` (apo / open) → `1AKE` (Ap5A closed / transition-state mimic).
- **Why:** large open/closed motion, clear hinge, small protein, easy expression. Ap5A mimics ATP+AMP
  (transition-state model). The paper describes a distal region whose mobility increases on substrate
  binding — a "counterweight" balancing binding energy: exactly a de-novo Zn-switch target.
- **Demo:** apo vs TS-mimic alignment → distal LID (117–163) / NMP (31–70) / hinge hotspots → design
  His-pair / His-triplet there → test Zn²⁺/Ni²⁺/Co²⁺ on kcat / Km.
- **Wet-lab:** small, easy to express, big conformational change → best de-novo metal-switch platform.
  This is an **engineering ("R&D") demo**, not a known-site recovery control.

## 4. PFK — *E. coli* phosphofructokinase (P0A796) · homotetramer · oligomeric benchmark
- **Pair:** `2PFK` (unliganded) → `1PFK` (products + allosteric activator **ADP/Mg²⁺**).
- **Why:** classic natural allosteric enzyme; conformational change around the effector site.
  Tests whether the workflow handles **oligomeric** allostery.
- **Demo:** go beyond residue displacement — compare **tetramer interface / effector site / active site
  contact maps**. Tetramer motion complicates single-chain alignment interpretation → not the first monomeric demo.
- **Wet-lab:** mature assay (coupled to aldolase/GDH), strong natural allostery.

## 5. TEM-1 β-lactamase (P62593) · monomer · **cpTEM-1 replication scaffold**
- **Structure:** `1BTL` (native TEM-1). Circular-permutant constructs (cpTEM-1-His3, cpTEM-1-3M-His2) have
  no direct PDB → use the **sequence** to design/rebuild constructs.
- **Why:** most literature-faithful replication target. Rationally designed cpTEM-1-His3 shows weak Zn²⁺
  activation; evolved cpTEM-1-3M-His2 is **down-regulated** by Zn²⁺/Ni²⁺/Co²⁺ (Ni²⁺/Co²⁺ mainly hit kcat;
  Zn²⁺ ~80% loss in kcat/Km). Same metal binding, different site geometry → **opposite regulatory direction**.
- **Wet-lab:** keep to purified enzyme + safe chromogenic substrate (nitrocefin/CENTA); do **not** center
  the project on antibiotic-resistance screening.

## 6. ATCase — *E. coli* aspartate transcarbamoylase (P0A786) · A6B6 12-mer · strong-allostery benchmark
- **Pair:** `6AT1` (T-state) → `1D09` (R-state, PALA bisubstrate analogue, tetrahedral-intermediate-like).
- **Why:** canonical strong natural allostery; tests recovery of **regulatory chains / subunit interfaces / T↔R transition**.
- **Caveat:** very large, complex oligomer → computational benchmark only, **not** a first Zn-engineering target.

## 7. GP — human liver glycogen phosphorylase (P06737) · homodimer · drug-site benchmark
- **Structure:** `3CEH` (tense/inactive state + allosteric inhibitor AVE5688 at the AMP site).
- **Why:** drug-type allosteric-site benchmark — can the workflow surface the **AMP / allosteric inhibitor site**?
- **Caveat:** large, dimeric, human → expression + assay costlier than bacterial enzymes.

---

### Notes on active/inactive pairing
- Displacement is computed after a **robust core superposition** (iterative outlier rejection) so the
  reference frame sits on the rigid domain — essential for hinge-type motions (GCK, AdK).
- For subtle-motion systems (PTP1B), raw Cα displacement under-reports coupling; complement with
  loop RMSF, contact-map differences, or an elastic-network / normal-mode analysis.
- Oligomeric systems (PFK, ATCase, GP) need interface-aware comparison across chains, not just chain A.
