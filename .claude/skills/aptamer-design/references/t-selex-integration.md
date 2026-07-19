# Upgraded workflow: T-SELEX × Boltz consensus in-silico SELEX

Integrates the **T-SELEX** toolchain (CMCDD/T_SELEX; iScience 2025) with this skill's
Boltz-2.1 co-folding + LM + counter-screen pipeline into a single **iterative in-silico
SELEX** with **two orthogonal 3D validators**.

## Why combine them (the core idea)
Two independent ways to judge an aptamer–target 3D interaction:

| | T-SELEX path | Boltz path |
|---|---|---|
| method | classic biophysics: 2D→3D→rigid dock | deep-learning co-folding |
| tools | ViennaRNA → RNAComposer → HDOCK | Boltz-2.1 |
| score | HDOCK score (−, lower better) + confidence % | ipTM / pLDDT |
| strengths | interpretable, template/physics-based, fast | joint fold, induced fit, **handles DNA** |
| weaknesses | rigid-body, **RNA-only (RNAComposer)** | "optimistic" for flexible ssNA |

They fail differently → a candidate that **both** rank highly is much more trustworthy.
**Consensus (rank-agreement) of HDOCK + Boltz is the central upgrade.** Neither score is
an affinity; consensus is a *relative confidence* signal.

## T-SELEX API (verified from the repo)
```python
from T_SELEX_program import (gen_aptamers, fold_and_composition,
                             tertiary_structure, intarna, Mol_docking_calc,
                             DMBA, BMA, PDCA)
RNA        = gen_aptamers(seed=1, length="randomize", aptamers_num=100)   # pool
df         = fold_and_composition(RNA)                                    # ViennaRNA 2D + MFE
structures = tertiary_structure(df["Aptamer"], df["MFE structure"])       # RNAComposer 3D (Selenium)
inter      = intarna('aptamers.csv','Aptamer', target_seq,'name', True)   # IntaRNA RNA-RNA
dock       = Mol_docking_calc(data_frame=df, MFE_column='Minimum free Energy',
                receptor_name="6zsl", receptor='6zsl.pdb',
                ligands_directory='ligands', directory_path='software/HDOCKlite-v1.1',
                Ap_folded=True)                                           # HDOCK dock+score
DMBA(...); BMA(...); PDCA(...)                                            # post-dock analytics
```

## The optimized pipeline (8 stages)

```
 round r ──────────────────────────────────────────────────────────────────┐
 0. POOL      gen_aptamers()  ∪  LM generate_lm.py (motif-biased)  ∪ AptamerBase
 1. 2D FILTER fold_and_composition()  → keep low-MFE, well-structured (cheap first cut)
 2. LIABILITY intarna()  → drop self-dimerizing / pool-cross-hybridizing / off-target RNA
 3. 3D MODEL  tertiary_structure()  (RNAComposer; RNA only — DNA skips to Boltz)
 4. DUAL 3D SCORING  (+ scramble decoy in both):
        A) Mol_docking_calc()  → HDOCK score + confidence     ┐
        B) Boltz-2.1 co-fold   → ipTM / pLDDT                 �it� consensus_rank.py
                                → rank-consensus (Borda)       ┘
 5. SPECIFICITY  repeat stage 4 vs paralogs → require on-target ≫ off-target in BOTH
 6. SELECT+EVOLVE  keep top-k by consensus → LM mutates/recombines winners (doped pool)
 └────────────  if not converged, feed new pool back to stage 1 (next SELEX round) ──┘
 7. REPORT   DMBA/BMA/PDCA analytics + rank_candidates.py (use-case weighted) + report_template
```

**What changed vs the old linear workflow**
- Added a real **iterative SELEX loop** (evolve winners with the LM instead of one-shot design).
- Added **IntaRNA liability screening** (self-dimerization / pool competition) — a failure mode Boltz alone misses.
- Added a **second, orthogonal 3D scorer (HDOCK)** and require **consensus** with Boltz.
- Kept the **scramble decoy** calibration and **paralog counter-screen** (decisive for diagnostics).

## Division of labour (who scores what)
- **RNA aptamers** → full T-SELEX 3D (RNAComposer→HDOCK) **and** Boltz → consensus.
- **DNA aptamers** → RNAComposer can't model them; use Boltz (native DNA) as primary, HDOCK
  only if a DNA 3D model is supplied by another tool. Report which scorers were available.

## Local install notes (verified on a bare Linux x86_64 sandbox)
- **HDOCKlite runs locally** once its one missing dep is satisfied. The `hdock`/`createpl`
  ELF binaries need `libfftw3.so.3`; if apt can't provide it, extract it from the **pyfftw
  wheel**: `pip download pyfftw --no-deps`, unzip, copy
  `pyfftw.libs/libfftw3-*.so.3.*` → `libs/libfftw3.so.3`, and run with
  `LD_LIBRARY_PATH=libs`. Verified working (`ldd hdock` resolves fftw).
- Usage (README order = **receptor first**): `hdock receptor.pdb ligand.pdb -out H.out`
  then `createpl H.out top10.pdb -nmax 10 -complex -models`. The docking score is in each
  model header as `REMARK Score: -NNN.NN` (more negative = better). A ~120-aa receptor ×
  33-nt RNA takes ~400 s single-thread.
- **HDOCK confidence** = `1/(1+exp(0.02*(score+150)))` (HDOCK server formula).
- `T_SELEX/Docking.py` is hardwired to the author's paths (`/home/s1800206/...`) with the
  docking loop commented out — **drive `hdock`/`createpl` directly**, don't rely on the wrapper.
- **ViennaRNA** installs cleanly via `pip install ViennaRNA` (the 2D stage runs natively;
  you don't need the T_SELEX wrapper for it).
- **RNAComposer / IntaRNA** remain the web/bioconda pieces to provision separately; if you
  already have an aptamer 3D PDB (e.g. from RNAComposer), HDOCK docks it directly.
- **Generating a decoy/new 3D without RNAComposer:** fold the sequence alone with Boltz,
  extract the RNA chain as PDB, and feed it to HDOCK — lets you calibrate the HDOCK side too.

## Honest limitations
- **RNAComposer & HDOCK are web servers**; T-SELEX drives RNAComposer with **Selenium**
  (brittle, rate-limited) and needs Linux + a **manual HDOCK download**. Not runnable in a
  bare sandbox — provision a Linux env per the repo before using stages 2–5's T-SELEX half.
- **HDOCK score and Boltz ipTM are different scales** → combine by **rank** (Borda), never raw sum.
- HDOCK is **rigid-body** (limited induced fit); Boltz is **learned** (optimistic). Consensus
  mitigates each but neither is a Kd.
- RNAComposer is **RNA-only**; DNA needs Boltz or a DNA-specific 3D method.
- Everything here is still a **prioritization for wet-lab SELEX**, not a validated binder.
