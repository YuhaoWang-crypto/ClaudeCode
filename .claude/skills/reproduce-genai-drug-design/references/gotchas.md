# Gotchas — hard-won fixes (each cost a debugging round)

## PDF / environment
- **`pdftoppm not installed`** → can't render PDF pages. Extract text with pypdf:
  `PdfReader(path).pages[i].extract_text()`. If it dies with
  `ModuleNotFoundError: _cffi_backend` (broken `cryptography`), run
  `pip install --force-reinstall cffi`.
- **apt is often unavailable** in these containers; rely on `pip`/conda only.

## Structure resolution (Part 1)
- **OpenBabel's CDX reader can hard-crash (SIGABRT, exit 134)** on some files and
  kill the whole process. Convert **one file per `subprocess.run(["obabel", f,
  "-osmi"])`** with a timeout so a crash only loses that file. Do NOT loop with
  the Python `pybel` API — one bad file aborts everything.
- **`.cdx` filenames are inconsistent** (`TCA1.cdx` vs `ardp.202200168_3a.cdx`).
  Match to activity rows by `(DOI-folder = DOI without "/", trailing label)`.
- **OPSIN fails on semi-systematic names** (~30%). Use it only as a fallback to
  CDX. `py2opsin` needs Java on PATH.

## Docking (Part 2)
- **meeko `mk_prepare_receptor.py` fatals on altlocs / cofactor bonds.** Strip to
  chain A, keep only altloc A and blank the altloc column, separate the cofactor,
  and run with `--allow_bad_res --default_altloc A`. Append the cofactor PDBQT
  (obabel `-xr`) to the protein PDBQT to keep it rigid in the receptor.
- **Vina pose → ligand for MD:** `obabel pose.pdbqt -osdf` **drops nonpolar H**
  (you get ~3 H). antechamber then fails. Reconstruct the full ligand (all H +
  bond orders) with `meeko.PDBQTMolecule.from_file` + `RDKitMolCreate.from_pdbqt_mol`
  — the pdbqt carries the SMILES for exactly this.

## MD parametrisation / tleap (Part 3) — the big one
Sequence of failures hit while getting `modal_md_app.py` to run, in order:
1. **`openff Molecule.to_file(mol2)` fails** — no MOL2 writer without OpenEye in
   a conda image. **Fix:** antechamber reads SDF natively — pass `-fi sdf` and
   take the net charge from RDKit `GetFormalCharge` (with a `sanitize=False`
   fallback for FAD, whose bond perception can fail).
2. **tleap `FATAL: Atom HIE.HD1 does not have a type`** — PDBFixer/OpenMM hydrogen
   names clash with ff14SB templates. **Fix:** don't add H in PDBFixer; add only
   missing heavy atoms, then `pdb4amber -y` (strips H, adds TER at chain breaks,
   HIS→HIE) and let **tleap** protonate with Amber naming.
3. **20–26 Å bonds between consecutive residues** — chain breaks / built loops.
   **Fix:** `fixer.missingResidues = {}` (don't build long loops) and rely on
   pdb4amber TER cards so tleap doesn't bond across gaps.
4. **`addIonsRand complex Na+ 0 Cl- 0` → FATAL "'0' not allowed for second ion"**
   — and `addIons` rejects it too. **Fix:** two single-ion calls
   `addIons complex Na+ 0` then `addIons complex Cl- 0`; tleap adds only the
   counterion matching the net-charge sign and no-ops the other (sign-agnostic).
5. **MMPBSA `Unrecognized namelist`** — `&general ... /` on one line doesn't
   parse. **Fix:** put the namelist name and the closing `/` on their **own
   lines**.
6. **MMPBSA `Could not open com.prmtop`** — `ante-MMPBSA.py` **skips writing the
   complex topology when there's no solvent to strip** (the dry complex == input).
   **Fix:** use `complex_dry.prmtop` directly as `-cp`; take rec/lig from
   ante-MMPBSA; pass `-sp SYS.prmtop` so MMPBSA strips water from the solvated
   trajectory to match. Set `set default PBRadii mbondi2` in tleap so complex
   radii match the igb=5 GB model.

## Cofactors
- **FAD** (no metal) parametrises fine with GAFF2 — but AM1-BCC on 88 atoms via
  `sqm` takes ~15 min. Use Gasteiger (`-c gas`) for fast test iterations (its
  internal charges cancel in the ligand MM-GBSA ΔG); AM1-BCC for production.
- **HEM/Fe is NOT supported** by GAFF2 (needs bonded-metal parameters). Keep heme
  only in the rigid docking receptor, not in the MD.

## Modal
- Set `MODAL_TOKEN_ID`/`MODAL_TOKEN_SECRET`; `pip install
  'modal[api-proxy-support]'` (proxy needs `python-socks`) or auth fails with a
  proxy error. First image build downloads ~1 GB (ambertools+openmm) then caches.
- Stream logs to a file and wait on `pgrep -f "modal run"` — the CLI's own exit
  can be masked by wrappers. Fetch per-run outputs with `modal volume get`.
