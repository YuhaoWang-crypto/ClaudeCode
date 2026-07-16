---
name: materials-sim
description: >-
  Computational materials design & validation on cloud GPU (Modal): build
  crystal structures with pymatgen, relax + phonon + bulk-modulus with FairChem
  UMA, DFT DOS / band-filling with GPAW, electron-phonon coupling (lambda,
  Allen-Dynes Tc) with Quantum ESPRESSO DFPT, and generative inverse design
  with MatterGen. Use for: designing a material for a target property (quantum
  dot emission, upconversion host, superconductor); element substitution
  (isovalent Ca->Sr->Ba, Pb->Sn) or doping (Nb->Ti, Er->Ca) property studies;
  high-throughput descriptor screening to push a property (higher quantum
  efficiency, higher Tc); generating candidates by a target property and
  verifying them; or overlaying simulation vs experimental data. Enforces
  ✅-verified / ⚠️-limitation labeling and knows which descriptor predicts which
  property (and when a single descriptor fails, e.g. phonons for Tc).
---

# Computational materials design & validation (UMA / GPAW / QE / MatterGen on Modal)

A reusable pipeline that turns a materials-design question into **computed
numbers on concrete crystals**, validated against experiment, always labeling
what is rigorous vs. limited. All heavy compute runs on **Modal GPU/CPU**; the
local side just builds structures and plots.

## What each engine computes (and cannot)

| Engine | Gives | Does NOT give | Script |
|---|---|---|---|
| **FairChem UMA** (`uma-s-1p1`) | relaxed structure, energy, forces, Γ phonons, bulk modulus (EOS) | band gap, Tc, carriers | `fairchem_relax_modal.py`, `subst_dope_modal.py` |
| **GPAW/PBE** | DOS, N(E_F), band gap, projected DOS, n/p-type | Tc; PBE underestimates gap | `gpaw_dos_modal.py`, `gpaw_diboride_modal.py` |
| **Quantum ESPRESSO DFPT** | phonon spectrum, electron-phonon λ, Allen-Dynes Tc | — (λ needs dense k/q grids to converge) | `qe_epc_modal.py` |
| **MatterGen** | generate crystals conditioned on a target property | property values (verify separately with UMA/DFT) | `mattergen_modal.py` |

**Property → right tool** (the key judgment):
- Quantum-dot emission / band gap → Brus equation (analytic) + GPAW gap; UMA only confirms structure.
- Upconversion efficiency → **max phonon energy** (UMA Γ phonon); lower = less multiphonon quenching. Fluorides < oxides.
- Superconductor Tc → **NOT phonon frequency alone** (fails: MgB₂ ranks low). Needs electronic descriptor (N(E_F), σ-band character via GPAW) AND ultimately electron-phonon λ (QE DFPT + Allen-Dynes).
- Hardness / incompressibility → bulk modulus (UMA EOS).

## Modal setup (do this first every session)

```bash
# tokens are usually in env: MODAL_TOKEN_ID / MODAL_TOKEN_SECRET / HF_TOKEN
uv venv /root/.modal-venv --python 3.11 && source /root/.modal-venv/bin/activate
uv pip install modal numpy "python-socks[asyncio]"     # python-socks: gRPC over the CONNECT proxy
export SSL_CERT_FILE=/root/.ccr/ca-bundle.crt          # trust the agent-proxy CA
export MODAL_TOKEN_ID MODAL_TOKEN_SECRET HF_TOKEN
modal app list                                         # verify auth
```
Behind the agent proxy, Modal's gRPC fails unless `python-socks` is installed and
`SSL_CERT_FILE` points at the CA bundle. UMA weights are gated → needs `HF_TOKEN`.
A second venv (`/root/.egnome-venv`, python 3.10, `pymatgen ase matplotlib`) builds
structures and plots locally. Modal function results contain numpy → the local
venv running `modal run` must have `numpy` (else DeserializationError, silent stale output).

## Workflow

1. **Build structures** with pymatgen (`Structure(Lattice..., species, coords)`), write CIFs to a JSON `{name: cif}`.
2. **Run the engine** on Modal: `modal run scripts/<engine>.py --candidates cands.json --out results.json`. Each run is background + monitored; the first run builds the image (cached after).
3. **Analyze + plot**: read results JSON, overlay simulation vs experimental values, save PNG. Use CJK font `WenQuanYi Zen Hei`, Okabe–Ito colorblind-safe palette, simulation=line/bar, experiment=points.
4. **Verify generative output**: MatterGen generate-by-target → feed the generated CIFs back through UMA to confirm the property (closed loop).
5. **Report**: build a self-contained HTML (base64-embed PNGs), `chromium --headless --print-to-pdf` for PDF.

## Patterns that worked (and pitfalls)

- **Isovalent substitution**: same-group swap (Ca→Sr→Ba, Pb→Sn, Nb₃Sn→V₃Si) → lattice/phonon/bulk-modulus trend, validate vs experiment. Bigger/heavier → lattice↑, phonon↓, modulus↓.
- **Doping**: build a 2×2×2 supercell, substitute 1 site (+ interstitial for charge compensation, e.g. Er³⁺→Ca²⁺ + F⁻). Look at Vegard lattice expansion and local bonds; use GPAW for the n/p-type DOS shift.
- **Screening**: enumerate → descriptor → filter stability → rank. Be honest when a single descriptor fails; a monotonic property↔structure link (upconversion↔phonon) screens well, an electronic-structure property (Tc) does not.
- **Pitfalls**: GPAW single-point `raw_dos` at E_F is unreliable for metals → compute N(E_F) by manual eigenvalue Gaussian histogram; give metals dense k + `nbands="200%"`; QE `electron_phonon='interpolated'` needs a dense nscf grid to converge λ (coarse = under-converged); Brus overestimates gap below ~3 nm; container restarts kill local `modal run` — results only persist when the local entrypoint's json.dump completes (or use a Modal Volume).

## Honesty labeling (required on every claim)
- ✅ **verified**: structure is a stable UMA minimum; lattice ≈ experiment; a descriptor trend matches experiment; a generated structure's property confirmed by an independent engine.
- ⚠️ **limitation**: property not directly computed (proxy/literature used); method known to be off (PBE gap, harmonic vs anharmonic phonon, under-converged λ); idealized/metastable structure.

Scripts are in `scripts/`. See the companion `materials_design/` project for a full worked example (quantum dots, upconversion, superconductors; substitution & doping; screening; inverse design; DFPT) with figures and a report.
