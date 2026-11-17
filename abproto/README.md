# abproto — Pythia-1-style antibody design pipeline (minimal runnable prototype)

A dependency-light scaffold of the design pipeline from the technical plan. The
**orchestration is real and runs on CPU with just numpy**; the two heavy learned
models sit behind swappable backends, with a deterministic mock so the whole
thing runs today with no GPU, no weights, no network.

```
input spec
  └─(1) all-atom joint generation ───────┐  Generator backend (Multiflow / La-Proteina)
     (2) confidence-guided Beam Search ◄──┘  test-time compute scaling
  └─(3) CTMC-style sequence redesign        liability-constrained variants
  └─(4) self-consistency refolding          pose recovery, Folder backend (Boltz-2 / AF3)
  └─(5) physics/sequence filter + 3Å cluster
  └─ ≤48 wet-lab candidates
```

## Run the dry run (no GPU)

```bash
pip install numpy pyyaml
python -m abproto.run --spec abproto/examples/il11_spec.json --backend mock --out results.json
```

You'll see the Beam Search mean confidence climb each step (test-time scaling),
then `generated → passed_filters → clusters → selected`.

## Tests

```bash
python -m pytest abproto/tests -q      # or: python abproto/tests/test_pipeline.py
```

Covers Kabsch RMSD (identity & rotation invariance), Beam Search monotonic
improvement, CTMC liability removal, clustering, and the end-to-end pipeline.

## What is real vs stubbed

| Component | Status | File |
| --- | --- | --- |
| Beam Search (test-time scaling) | **real** | `beam_search.py` |
| Kabsch superposition + RMSD | **real** | `geometry.py` |
| Self-consistency / pose recovery | **real** | `selfconsistency.py` |
| Sequence-liability filters (N-glyc, deamidation, isomerization, Cys) | **real** | `filters.py` |
| CTMC-style liability-constrained redesign | **real logic**, learned model stubbed | `ctmc.py` |
| 3Å clustering (avg-linkage via scipy, else single-linkage) | **real** | `cluster.py` |
| Generator (joint seq+structure) | **adapter + mock** | `backends/multiflow.py`, `backends/mock.py` |
| Folder (refold + pLDDT/PAE/PDE) | **adapter + mock** | `backends/boltz.py`, `backends/mock.py` |

## Wiring in the real backends

- **Folder → Boltz-2** (MIT): `pip install boltz -U`, then use `--backend boltz-multiflow`.
  `backends/boltz.py` writes the Boltz YAML, runs `boltz predict`, and parses the
  confidence JSON + mmCIF C-alpha coords. Works out of the box once Boltz is installed.
- **Generator → Multiflow** (MIT): needs a checkpoint **fine-tuned on SAbDab** for
  antibody CDR co-design conditioned on the antigen + germline framework (route B
  in the plan). `backends/multiflow.py` marks the exact integration points; until
  that checkpoint exists it raises a clear `NotImplementedError`.

Both backends satisfy the `Generator` / `Folder` protocols in `backends/base.py`,
so anything implementing those (La-Proteina, AF3, …) drops straight in.

## Layout

```
abproto/
  types.py          DesignSpec, Design, Pose, Confidence, ScoredDesign
  geometry.py       Kabsch RMSD, clash proxy
  confidence.py     scalar objective from pLDDT/PAE/PDE/iPTM
  beam_search.py    (2) test-time compute scaling
  ctmc.py           (3) liability-constrained sequence redesign
  selfconsistency.py(4) refold + pose recovery
  filters.py        (5) sequence liabilities + physics gate
  cluster.py        (5) 3Å agglomerative clustering + representatives
  pipeline.py       end-to-end orchestration
  run.py            CLI
  backends/         base protocols + mock + boltz + multiflow adapters
  examples/         il11_spec.json (toy)
  tests/            unit + smoke tests
```

The mock is **not a model** — it fabricates plausibly-shaped outputs with a hidden
quality signal so Beam Search, self-consistency and clustering have something real
to optimise. Replace the backends and the same orchestration drives real designs.
