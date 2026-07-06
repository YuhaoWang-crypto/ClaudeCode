---
name: md-screen
description: Run making-it-rain-style molecular dynamics (MD) on Modal serverless GPUs to screen small molecules for activity, binding affinity (MM/GBSA), and allosteric effects. Use when the user wants to run protein MD, small-molecule MD, protein-ligand MD, estimate binding free energy, rank ligands, judge whether a compound is an active binder / enzyme inhibitor / activator, or detect allosteric modulation — using their Modal GPU account. Covers the full loop: prepare (PDBFixer + OpenFF + GAFF2), simulate (OpenMM), analyze (RMSD/RMSF/Rg/contacts), score (MM/GBSA), and allostery (DCCM + apo/holo).
---

# md-screen: MD-based small-molecule screening on Modal GPUs

This skill runs the **mdscreen** system (in this repo: `mdscreen/` + `modal_app.py`)
which reimplements the [making-it-rain](https://github.com/pablo-arantes/making-it-rain)
MD notebooks (OpenMM + GAFF2 + ff14SB + TIP3P) on **serverless Modal GPUs**, to
judge whether small molecules are active binders, enzyme inhibitors/activators,
or allosteric modulators.

## When to use
- "Run MD for this protein / ligand / complex on my Modal GPU"
- "Estimate the binding free energy / rank these compounds"
- "Is this molecule likely active / an inhibitor / allosteric?"
- "Screen this ligand series against this target"

## Quick reference — the four workflows

```bash
# protein-only MD (apo dynamics, stability)
modal run modal_app.py::protein --pdb 1ubq --ns 5

# small-molecule-only MD (conformational stability)
modal run modal_app.py::ligand --ligand "CCO" --ns 5          # SMILES
modal run modal_app.py::ligand --ligand ligand.sdf --ns 5     # file

# protein-ligand complex MD + MM/GBSA binding free energy
modal run modal_app.py::complex --pdb 3poz.pdb --ligand ligand.sdf --ns 2

# full activity + allostery screen (parallel across GPUs)
modal run modal_app.py::screen --config examples/screen_config.yaml

# pull results from the persistent volume
modal volume get mdscreen-outputs <run_id> ./results
```

Pick GPU via env var (default `A10G`): `MDSCREEN_GPU=A100 modal run ...`

---

## Step-by-step procedure (follow this order)

### 1. One-time environment setup
The Modal client needs proxy support in this sandboxed environment, or auth
hangs with `ImportError: python-socks not installed`.

```bash
pip install --quiet modal pyyaml 'python-socks' 'modal[api-proxy-support]'
```

Auth uses the `MODAL_TOKEN_ID` / `MODAL_TOKEN_SECRET` env vars already set in the
session. Verify with a trivial remote call before doing real work:

```bash
python3 - <<'PY'
import modal
app = modal.App('auth-check')
@app.function()
def ping(): return 'pong'
with modal.enable_output():
    with app.run():
        print('AUTH OK:', ping.remote())
PY
```

If `modal_app.py` imports cleanly (`python3 -c "import modal_app"`) the app is
wired up.

### 2. Prepare inputs
- **Receptor**: a local `.pdb` path OR a 4-letter RCSB id (auto-downloaded &
  cleaned by PDBFixer).
- **Ligand**: a **SMILES** string OR an `.sdf`/`.mol2` file. AM1-BCC charges and
  a 3D conformer are generated automatically.
- ⚠️ **Ligand pose**: `::complex` and `::screen` place the ligand at its own
  coordinates and do NOT dock. For a meaningful complex the ligand must already
  be a sensible bound pose (co-crystal / docked SDF in the pocket). If the user
  only has SMILES, add a docking step first (Vina/gnina) or confirm they have
  poses. Say this explicitly rather than producing a meaningless floating-ligand
  run.

### 3. Launch as a background job + monitor (runs are minutes-to-hours)
`modal run` streams live; run it in the background and watch the log so you
catch prep/force-field errors early (they surface right after the image build,
before the long MD).

```bash
MDSCREEN_GPU=A10G modal run modal_app.py::protein --pdb 1ubq --ns 1 \
  > /tmp/run.log 2>&1 &     # use the Bash run_in_background flag, not a raw &
```

Then attach a Monitor that emits progress AND every failure signature:
```
tail -f /tmp/run.log | grep -E --line-buffered \
  "minimisation|equilibration|production|Analysing|Analysis summary|Pull results|Error|Traceback|Exception|UnicodeEncode|Failed|App completed"
```
See `scripts/run_and_watch.sh` for a ready-made launcher.

### 4. Pull and interpret results
Directory `modal volume get` can error `Is a directory`; fetch files
individually into an existing dir:
```bash
modal volume ls mdscreen-outputs <run_id>/<name>
mkdir -p out
for f in apo_rmsd_backbone.png apo_rmsf.png apo_rgyr.png apo_analysis.json; do
  modal volume get mdscreen-outputs "<run_id>/<name>/$f" "out/$f"
done
```
Read the PNGs and the `*_analysis.json` summary; surface the plots to the user
with SendUserFile.

---

## Interpreting the output

| Metric | Meaning | Healthy / signal |
|---|---|---|
| `rmsd_backbone_*_nm` | protein stability vs. start | plateaus at ~0.1–0.3 nm |
| `rmsf_*_ang` | per-residue flexibility | rigid core ~0.5 Å, flexible loops/termini spike |
| `rgyr_mean_ang` | compactness | stable ⇒ folded |
| `ligand_rmsd_*_nm` | pose stability | ≤ 0.25 nm ⇒ stable pose |
| `mmgbsa.dg_bind_kcal` | binding score (no entropy) | more negative ⇒ stronger; ranking only |
| `allosteric_verdict` | apo→holo distal dynamics | distal RMSF/DCCM change ⇒ allosteric |

**Activity call**: strong ΔG + stable pose ⇒ *likely active* (orthosteric ⇒
competitive inhibitor candidate). Distal allosteric signal ⇒ *allosteric
modulator* (activator or inhibitor — confirm mechanism by assay). MM/GBSA is a
*ranking* score for a congeneric series, not an absolute affinity.

---

## Verified working environment (as of last run)
Image built from conda-forge, all force-field strings resolved:
- openmm 8.1.2, cudatoolkit 11.8, ambertools 24.8, openff-toolkit 0.18,
  openmmforcefields 0.15.1, mdanalysis 2.10, rdkit 2025.03, pdbfixer 1.9
- Force fields that resolve: `amber/ff14SB.xml`, `amber/tip3p_standard.xml`,
  `gaff-2.11`, `implicit/gbn2.xml` (MM/GBSA)
- Timing (A10G): ubiquitin ~17k atoms, 1 ns production ≈ 100 s; +~30 s min/equil;
  first run adds ~150 s image build (cached after).

## Gotchas learned (don't rediscover these)
1. **Proxy**: install `python-socks` + `modal[api-proxy-support]` or Modal auth
   errors in the sandbox.
2. **UTF-8**: the container defaults to ASCII; writing CSV headers containing
   "Å" crashes with `UnicodeEncodeError`. Fixed in code via `encoding="utf-8"`
   on every CSV open + `PYTHONUTF8=1`/`PYTHONIOENCODING=utf-8` in the image env.
   If you add new file writes, keep them UTF-8.
3. **`.env()` on the image** creates a new cheap layer — it does NOT re-run the
   150 s conda install (that layer is content-cached).
4. **`modal volume get <dir>`** may fail on a directory target; get files
   individually.
5. **Atom ordering** for MM/GBSA assumes complex = protein-atoms-then-ligand
   (as built by `prepare.build_complex`). Don't reorder.
6. Non-standard residues, metals, cofactors, covalent/charged-metal ligands, and
   membranes are not auto-handled — flag to the user.

## Scaling a screen
`::screen` runs the apo reference once, then fans every ligand to its own GPU
container via `complex_md.starmap`, then runs allostery comparisons. N ligands ≈
wall-clock of 1. Raise `production_ns` (and add replicas) for production rigor;
2 ns is a screening default.

## Files
- `modal_app.py` — Modal image, GPU functions, `local_entrypoint`s
- `mdscreen/` — engine (config/prepare/simulate/analyze/binding/allostery/pipeline/cli)
- `examples/screen_config.yaml` — worked screen config
- `scripts/setup_modal.sh`, `scripts/run_and_watch.sh` — helpers in this skill
