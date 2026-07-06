# Hard-won lessons (read before running — each cost real debugging time)

These are the failure modes actually hit while building this pipeline. Most are environment/setup
traps that produce *plausible-looking wrong numbers* or hard crashes.

## 1. Modal + sandbox proxy
- The agent sandbox routes outbound through a SOCKS proxy. `modal run` fails with
  `'python-socks' not installed` until you `pip install 'modal[api-proxy-support]' python-socks`.
- Module-scope file reads execute **inside the Modal container too** (the container imports your
  script). `open("data/x")` at module top → `FileNotFoundError` remotely. **Read inputs only inside
  `local_entrypoint()`** and pass contents as function args.
- Long jobs: `modal run --detach` + write results to a `modal.Volume` so a dropped local client
  doesn't lose them. Poll the volume with a background bash job (harness re-invokes you on exit);
  do NOT foreground-`sleep`.
- `tied_featurize` (ProteinMPNN) return order matters: `residue_idx` is index **12**, not 11
  (11 is `omit_AA_mask`). Wrong index → `gather_edges` dimension error.

## 2. MM-GBSA (OpenMM + AmberTools)
- **tleap renames terminal nucleotides** (DA→DA5/DA3, DC→DC5/DC3, …). A `{DA,DT,DG,DC}` name
  filter then miscounts the protein/DNA boundary → receptor and ligand get mismatched coordinates
  → energies of ~−10⁷ kcal/mol. **Split complex into receptor/ligand by ATOM INDEX** (protein is
  the first N atoms), never by residue name. Sanity: component energies should be ~−10⁴ kcal/mol.
- Let **tleap add hydrogens**, not pdbfixer — pdbfixer's HIS H-naming (HD1 on HIE) clashes with
  tleap templates. Strip H before tleap.
- Strip **5′-terminal DNA phosphates** (P/OP1/OP2 on the first residue of each DNA chain) — the
  tleap DX5 template has none, else "atom does not have a type".
- Set `PBRadii mbondi2` in tleap for OBC2 GB. Native OBC2 works via
  `AmberParm.createSystem(implicitSolvent=OBC2)` — the prmtop route, not the XML route.
- Single-trajectory MM-GBSA **overestimates charge-adding mutations** (bare Coulomb to phosphates
  not offset by approximate GB desolvation). Report **direction only**; use FEP for magnitude.

## 3. QM cluster (xtb + Psi4)
- Automated `obabel -p 7` protonation of a metal cluster often yields an **odd electron count** at
  the naïve charge → xtb runs an unphysical open-shell doublet / Psi4 fails chg/mult. Compute
  electron parity and choose an even-electron charge in the **anionic** direction (0/−1/−2 for an
  acidic active site), not +1. Verify closed-shell singlet.
- Psi4 `b3lyp-d3bj` needs the external `s-dftd3`/`dftd3` binary (else `ResourceError`). Use plain
  `b3lyp` or conda-install `dftd3`.
- Set charge/mult on the molecule object (`mol.set_molecular_charge`; `mol.update_geometry()`), not
  only in the geometry string, to avoid qcelemental chg/mult reconciliation errors.

## 4. QM/MM reaction barrier
- **Distance audit first**: a barrier is only informative if the mutation is near the *scissile
  atom*, not just the catalytic residues. (Here V257 was 17 Å from the DEDD tetrad but ~20 Å from
  the scissile phosphate → uninformative; A13 was 6.7 Å → the right target.)
- A **product-state** cryo-EM structure is not a reactant. Build a pre-cleavage Michaelis complex:
  re-ligate the scissile bond, place the 2nd catalytic metal at its canonical site (near a free
  non-bridging phosphate O + a second carboxylate, ~3.8 Å from metal 1), add an in-line
  metal-activated nucleophile (~180° to the leaving group). See `build_reactant.py`.
- **A bare gas-phase cluster cannot hold a two-metal geometry** — without the protein scaffold's
  electrostatics the two Mg²⁺ collapse together (saw Mg–Mg 2.9 Å vs correct ~3.8). This forces the
  conclusion: a real barrier needs **QM/MM with the protein environment**, not a cluster.
- xtb SCF for metal clusters: raise `--etemp` (400–1000) and `--iterations`; if a hand-built
  geometry has clashes, **GFN-FF pre-optimize first** (`--gfnff`) — but note it can distort metal
  coordination, so re-check geometry after.
- Enforce a **quality gate**: reject the profile unless it is smooth (no >~60 kcal/mol single
  steps), the barrier is physical (~3–80 kcal/mol), and the TS is interior. DFT single-points on
  non-DFT-optimized scan geometries are NOT a valid barrier (they gave an impossible negative here).

## 5. FEP/TI (perses)
- The FEP *engine* (OpenMM + openmmtools + pymbar on CUDA) installs and runs fine.
- Adding `pdbfixer`/`mdtraj` to the perses conda solve pulls **incompatible openff** versions
  (`ImportError: cannot import name 'Unit' from 'openff.units'`). Use the minimal working solve and
  add pdbfixer via `pip_install(..., extra_options="--no-deps")`.
- perses `PointMutationExecutor.__init__` **hard-imports `openeye.oechem`** (proprietary). The
  conda package installs without a license, but **runtime segfaults** without `OE_LICENSE` (Modal
  aborts after 8 crash-retries; a native segfault is NOT catchable by try/except).
- Therefore, in a bare sandbox: rigorous FEP needs a **licensed OpenEye** env, or the OpenEye-free
  **GROMACS + pmx** pipeline, or a custom single-topology TI in openmmtools. Validate with a
  **5-iteration smoke run** before committing GPU-hours — it catches all of the above cheaply.
- perses charge correction for charge-changing mutations: `transform_waters_into_ions_for_charge_changes=True`.

## 6. Artifacts / figures
- Artifact CSP blocks external scripts/fonts/CDNs. **Inline** `3Dmol-min.js` and embed the PDB as
  text; use system-font stacks (no webfont URLs). Validate the HTML headless with Playwright
  (`/opt/pw-browsers/chromium`) and capture console errors before publishing.

## 7. General honesty rules that paid off
- Validate cheap before expensive (5-iter FEP smoke; xtb before DFT; a distance audit before a barrier).
- Gate every physics number; report negatives as negatives; discard artifacts explicitly.
- Confirm the input sequence/numbering with the user when it came from an unreadable upload.
