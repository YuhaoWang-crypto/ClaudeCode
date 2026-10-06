# AlphaProtein Novo: audit of the Biomni run package, and a runnable replacement

```bash
# audit any AP Novo package (no GPU, no jax)
python3 -m trnaplat.apnovo.audit_apnovo_package <package_dir>
# ... and cross-check against the real upstream validator and sampler
python3 -m trnaplat.apnovo.audit_apnovo_package <package_dir> --repo <alphaprotein-novo>
# rebuild a runnable package from a broken one
python3 -m trnaplat.apnovo.build_fixed_package <package_dir> --out fixed/
```

## Verdict

| | |
|---|---|
| ✅ | `google-deepmind/alphaprotein-novo` **exists** — Apache 2.0, `run_pipeline.py` / `run_generator.py` / `run_ligandmpnn.py` / `run_alphafold.py` / `evaluate_design.py`, 5 evaluation suites exactly as the package claims |
| ✅ | Both weight URLs are live: `generator.bin.zst` **547,140,096 B**, `af3_leaving_atom.bin.zst` **1,020,518,661 B**, HTTP 200 |
| ✅ | The install commands match the repo's own README verbatim |
| ✅ | The motif CIF is well-formed, and **19/19** residues carry wild-type MmPylRS's own residue at their stated position (UniProt Q8PWY1) |
| ✅ | CCD `YLY` matches the ligand's atom set exactly (C22 H35 N8 O9 P = 22 C, 8 N, 9 O, 1 P) |
| ❌ | **Both design jobs fail before the first diffusion step.** Demonstrated with AP Novo's own `sample_designable_lengths` |
| ❌ | `run_commands.sh` step 3 raises `app.UsageError` immediately, and never installs LigandMPNN |
| ⚠️ | The ligand is the **adenylate intermediate**, not the free ncAA: 10 of the 19 "ncAA pocket" residues are nearer the nucleotide |

**So: not as shipped. After three fixes, yes** — `fixed/` passes the upstream
`Manifest.from_file` and `sample_designable_lengths` for all four jobs.

---

## Defect 1 — 18 of 19 pocket residues are silently discarded (fatal)

The rule, from `alphaprotein_novo/data/residue_mapping.py`:

```python
def _is_segment_designable(segment_str: str) -> bool:
  return segment_str[0].isnumeric()
```

A comma-separated segment starting with a **letter** is a fixed motif residue;
one starting with a **digit** is a designable length. The package writes

```
A300,302,305,306,330,332,338,339,342,344,346,348,384,396,399,401,417,426,437,60-120,40-80,40-80,40-80/B1
```

so AP Novo fixes **`A300` only**, and reads the other eighteen as designable
segments of exactly 302, 305, 306, … 437 residues — **6,492 residues** of
requested design. With the four real ranges the designable total is pinned to
[6672, 6852] against a declared `seq_length` of 260–360.

✅ Run against the upstream sampler:

```
ValueError: Cannot sample values in range [260, 360]:
            the range of possible sample values is [6672, 6853]
```

The official grammar interleaves them, every fixed residue carrying its chain:

```
"20-100,A56,20-70,A100,20-70,A125,20-100/A2001/A3001"
```

or lifts the motif out with `|` and `{}` placeholders, or — the mode documented
for exactly this problem — puts the residues in `unindexed_motif_residues` and
leaves `motif_str` a plain length range.

## Defect 2 — the two length fields contradict each other (fatal)

The second job's motif syntax is correct (residues in
`unindexed_motif_residues`), but `motif_str: "60-250"` caps the designable chain
at 250 while `seq_length: "260-360"` demands at least 260:

```
ValueError: Cannot sample values in range [260, 360]:
            the range of possible sample values is [60, 250]
```

The non-obvious part, and the likely cause: **`seq_length` bounds the sum of the
*designable segments*, not the finished protein.**
`sample_designable_lengths` passes only the designable ranges into
`sample_values_with_sum_in_range` and the fixed motif residues never enter that
sum.

## Defect 3 — the launch command cannot start (fatal)

The manifest sets `resequence.enabled: true`, and `run_pipeline.py` enforces:

> `The manifest sets 'resequence.enabled', but --ligandmpnn_dir …` → `app.UsageError`

Step 3 of `run_commands.sh` passes neither `--ligandmpnn_dir` nor
`--ligandmpnn_python`, and the script never clones LigandMPNN or fetches its
weights — the repo's README covers that in a separate section. Step 3 is also
labelled a "quick validation" while running the unmodified manifest at 1000
sampling steps × 50 designs × 2 jobs.

## ⚠️ The schema passes anyway

```
Manifest.from_file : PASSED
❌ 2/2 job(s) cannot run as written.
```

Pydantic validates types; the motif grammar is only checked when the sampler
runs — on a GPU node, after the weights have loaded. That gap is the reason
`audit_apnovo_package.py` exists.

---

## The scientific finding: this is not the ncAA pocket

PDB **2Q7H** is titled *"Pyrrolysyl-tRNA synthetase bound to **adenylated
pyrrolysine** and pyrophosphate"*, and CCD **`YLY`** is
`C22 H35 N8 O9 P` — pyrrolysyl-**adenylate**. The CIF's chain B carries
`O5' C5' C4' O4' C3' O3' C2' O2' C1' N9 N7 C2 N6` and a phosphate: ribose and
adenine. Free pyrrolysine has 18 heavy atoms; this ligand has 40.

✅ Splitting the ligand at the ester oxygen and taking each motif residue's
minimum heavy-atom distance to each half:

| nearer the amino-acid half (9) | nearer the nucleotide half (10) |
|---|---|
| M300, A302, L305, Y306, **N346**, **C348**, Y384, V401, W417 | R330, E332, H338, L339, F342, M344, E396, S399, R426, I437 |

**10 of 19 residues in a motif described as "the ncAA binding pocket" sit closer
to AMP than to the amino acid.** Conditioning on it asks the model to rebuild the
amino-acid pocket *and* the ATP site at once. For a synthetase that is arguably
correct — it must bind both — but it is not what the package says it is doing,
and it roughly doubles the problem's difficulty.

The 9 amino-acid-side residues include exactly **N346 and C348**, the two
positions that distinguish the FPFORCOM parent IFRS from wild-type MmPylRS
(`trnaplat/demo1_multisite.py`). The pocket this motif really constrains and the
pocket the directed-evolution literature mutates are the same pocket — an
independent consistency check across two unrelated sources.

⚠️ The ligand is **not** trimmed in `fixed/`. Removing the AMP half would leave
an atom set that no longer matches CCD `YLY`, and `folding.states` builds its
ligand from that code — so a trimmed ligand needs a custom CCD entry or a
supplied component CIF. Until then, even the 9-residue job designs around the
adenylate.

---

## Two corrections to the HTML report

**1. IFRS is `N346S/C348I`, not `N346I/C348S`.** The report's Demo 2 table
(§2.2) states `IFRS (N346I/C348S)`. Three independent sources say otherwise, and
one of them is the report's own motif file:

| source | residue 346 | residue 348 |
|---|---|---|
| UniProt Q8PWY1 (wild-type MmPylRS) | **N** | **C** |
| `pylrs_pyl_motif.cif` in this very package | **ASN**346 | **CYS**348 |
| the released IFRS sequence vs that wild type | N→**S** | C→**I** |

So the wild type carries N346 and C348, and IFRS mutates them to S and I. The
report has the two substituted residues swapped.

**2. `D76N` makes Com2 worse, so a top prediction containing it is a miss.** The
report (§1.3) states the top prediction `N7Y/H63L/K67N/V74W/D76N` is *"完全一致"*
with the published Com2-IFRS. In the paper's own measured panel:

| variant | measured fitness |
|---|---|
| `N7Y/H63L/K67N/V74W` | **2.752** ← the maximum, 4 mutations |
| `N7Y/H63L/K67N/V74W/D76H` | 2.597 |
| `N7Y/H63L/K67N/V74W/D76N` | **2.519** |

Adding `D76N` costs 8.5% of activity. The published Com2-IFRS is the
**quadruple**; a quintuple top-1 is a measurable miss, not agreement — and the
report labels that row `literature-validated`.

---

## Capability boundary

The package's own §4 is right, and the audit sharpens it:

- ✅ **No RNA conditioning.** The generator is motif + ligand conditioned and
  folds with AF3; nothing in the pipeline conditions on an RNA chain. The
  TBD–tRNA^Pyl interface — where 6 of the 7 COM1 activity mutations sit, per
  `trnaplat/demo1_multisite.py` — is **out of scope**, and that is the half of
  aaRS engineering with the published activity gains.
- ✅ **No matching evaluation suite.** The five built-in suites
  (`kemp_eliminase`, `serine_esterase`, `dehp_esterase`, `carbene_transfer`,
  `nitrene_transfer`) have no aaRS analogue, so scoring falls back to AF3
  confidence plus a custom pocket-RMSD script.
- ✅ **Partial diffusion is the better formulation here.** `examples/` ships
  `dehp_esterase_partial_diffusion`, and `partial_diffusion_input_file` /
  `partial_diffusion_num_steps` exist. Starting from the real 2Q7H fold and
  perturbing it preserves the measured pocket geometry, whereas de novo
  scaffolding asks a new fold to rediscover it. `manifest_partial_diffusion.json`
  sets this up; point it at a full 2Q7H chain A, not at the motif fragment.
- ⚠️ **Where it does fit the platform:** generating alternative scaffolds that
  hold a given ncAA pocket geometry — a *hypothesis generator* upstream of the
  Demo 2 scaffold comparison, not a replacement for it. Any output still needs
  an orthogonality argument built from scratch, which is the cost a new scaffold
  carries and an existing cross-species scaffold does not.
- ⚠️ **The GPU estimates are unverified.** "~minutes/design on an A100,
  1–2 days for 200 designs/job" is the package's extrapolation from the repo's
  kemp-eliminase benchmark, not a measurement. There is no GPU in this
  environment either, so neither of us has timed it.

---

## What `fixed/` contains

| file | |
|---|---|
| `motif_active_site_full.cif` | all 19 residues + YLY — the whole active site |
| `motif_ncaa_pocket.cif` | the 9 amino-acid-side residues + YLY, 120 atoms |
| `manifest_fixed.json` | 3 jobs: ncAA pocket unindexed, full site unindexed, ncAA pocket indexed (correctly interleaved) |
| `manifest_partial_diffusion.json` | the formulation that keeps the real geometry |
| `run_commands_fixed.sh` | installs LigandMPNN, passes its two flags, and smoke-tests with `--only_stage=generation` at 50 steps × 1 design |

✅ All four jobs pass upstream `Manifest.from_file` and
`sample_designable_lengths`. Sampled example from the indexed job:

```
39,A300,33,A302,42,A305,28,A306,31,A346,44,A348,36,A384,32,A401,29,A417,25/B1
```

⚠️ "Passes validation and sampling" is not "produces a good design". Nothing
here has been near a GPU, and the evaluation gap above is unfixed: there is no
aaRS suite to score the output against.
