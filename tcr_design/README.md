# `tcr_design` — computational TCR engineering: the landscape, and a runnable demo

Two things in here:

1. **[A survey of what actually works](#part-1--the-landscape)** in computational TCR
   design, affinity optimisation and mutation analysis, as of late 2026, with
   the honest scorecard of what is mature and what is near-random.
2. **[A working demo](#part-2--the-demo)** — a CPU-only pipeline that measures a
   real TCR–pMHC interface, triages single mutants, gates them on specificity
   and developability risk, screens the epitope against the human proteome, and
   writes the input files for the tools that can actually score the survivors.

The demo includes a retrospective benchmark **that its own scoring function
fails**. That result is reported prominently rather than buried, because it is
the most useful thing the demo has to say.

---

## Part 1 — the landscape

### 1.1 The short version

| Layer | Status | Use |
|---|---|---|
| TCR-alone structure | ✅ mature | TCRBuilder2+ (seconds, CPU) |
| TCR–pMHC complex structure, class I | ✅ good enough to design on | AlphaFold3 best; TCRdock / TCRmodel2 open alternatives |
| Repertoire distance / clustering | ✅ mature | `tcrdist3`, SCEPTR embeddings |
| Databases | ✅ mature | TCR3d 2.0, STCRDab (both weekly), ATLAS for mutant affinities |
| Construct assembly | ✅ mature | Stitchr |
| **De novo pMHC binder / TCR-mimic design** | ✅ **working in the wet lab** | RFdiffusion + ProteinMPNN + AF2 |
| Affinity maturation from a crystal structure | ⚠️ works, but slow physics + wet-lab iteration | Rosetta flex ddG, MM/PBSA |
| Epitope specificity from sequence, **unseen epitopes** | ❌ near random | — |
| Effect of a point mutation in a peptide | ❌ **chance-level** | — |
| Affinity from co-folding confidence | ❌ not affinity | — |
| A TCR-specific ΔΔG predictor | ❌ **does not exist** (this package fits a weak one, §2) | — |
| An "RF-TCR" diffusion model | ❌ **does not exist** | — |

### 1.2 Structure prediction — the mature layer

Two independent 2025 benchmarks, both including structures deposited *after*
the models' training cutoffs, agree that **AlphaFold3 leads** on TCR–pMHC, with
roughly 90% of models at acceptable quality or better and ~20% high-quality.
TCRmodel2 and Chai-1 are close behind.

The finding that matters most for engineering: **Boltz-2 and tFold-TCR produced
zero medium- or high-quality models across a 20-complex benchmark.** They
recover the fold (TCR interface RMSD < 2 Å) but place the receptor in the wrong
docking geometry — which is exactly the thing interface design depends on. Open
weights do not rescue a wrong docking mode.

| Tool | Repo | License | Notes |
|---|---|---|---|
| **TCRdock** | [phbradley/TCRdock](https://github.com/phbradley/TCRdock) | MIT + Apache-2.0 | AF2 with native TCR docking-geometry templates; also TCRdist + docking-RMSD utilities |
| **TCRmodel2** | [piercelab/tcrmodel2](https://github.com/piercelab/tcrmodel2) · [server](https://tcrmodel.ibbr.umd.edu) | check repo | End-to-end TCR / TCR–pMHC, >10× faster than default AF |
| **ImmuneBuilder / TCRBuilder2+** | [oxpig/ImmuneBuilder](https://github.com/oxpig/ImmuneBuilder) | see repo | TCR **alone**, no pMHC. Seconds, CPU-feasible, ships weights |
| **AlphaFold3** | DeepMind | code open, **weights gated by application**, non-commercial | Best for this problem |
| **Boltz-1 / Boltz-2** | [jwohlwend/boltz](https://github.com/jwohlwend/boltz) | **MIT** (code + weights) | Open, but see the docking-geometry caveat above |
| **Chai-1** | [chaidiscovery/chai-lab](https://github.com/chaidiscovery/chai-lab) | Apache-2.0 (repo) | Roughly AF3-comparable in one benchmark |
| **NetTCRfold-1.0** | [mnielLab/NetTCRfold-1.0](https://github.com/mnielLab/NetTCRfold-1.0) | not stated | TCR–pMHC-specialised AF3 pipeline (Aug 2026). New; requires SLURM |
| ESMFold | Meta ESM | MIT | Single-chain monomer predictor — **not** a complex predictor. Don't use it for interfaces |

**Second caveat: structure confidence is not affinity.** This is worth stating
with numbers, because using ΔipTM between a wild type and its point mutants as
an affinity signal is a common instinct and it has been directly tested:

- **AlphaFold3 ipTM delta vs measured ΔΔG on SKEMPI v2** (7,085 mutations, 348
  complexes): **Pearson R = 0.102.** The stated reason is mechanical — ipTM
  operates on a coarse, near-binary scale and cannot express a continuous ΔΔG
  ([ICLR 2025 GEM workshop](https://openreview.net/pdf/bce908983c69f44ef567607cf4a37711cdaa9a2e.pdf)).
- **Boltz-2 fine-tuned specifically for protein–protein affinity**, evaluated on
  TCR3d 2.0 (~251 TCR complexes): **r = 0.153**, and feeding it true
  experimental structures instead of predicted ones barely moves it (r = 0.159),
  so structure quality is not the bottleneck. Plain **ESM2-650M sequence
  embeddings reach r = 0.239** on the same split — a sequence-only model beats
  the co-folding model ([arXiv:2512.06592](https://arxiv.org/abs/2512.06592)).
- Separately, AF3 *structures* are excellent ΔΔG substrate when fed to a
  downstream model (r = 0.86 vs 0.88 with crystal structures), but **ipTM failed
  to flag the structures with large errors**
  ([JCIM 2024](https://pubs.acs.org/doi/10.1021/acs.jcim.4c00976)).

The pattern is consistent: the **structures** are useful, the **confidence
scores** are not the signal. Use co-folding to reject mutants that break the
interface, not to rank the ones that don't. If you want to use a co-folding
model quantitatively, take deltas on its *internal embeddings* rather than its
confidence scores — that is what
[PreFold-dG](https://academic.oup.com/bioinformatics/article/42/Supplement_2/btag489/8767279)
does.

### 1.3 Sequence-based specificity prediction — the unreliable layer

This is where the gap between published AUCs and usable performance is widest.

- **[Nature Methods 2026;23:248–259](https://doi.org/10.1038/s41592-025-02910-0)**
  — 50 models, 21 datasets, 762 epitopes. Conclusion in the authors' own words:
  "all struggle to generalize to unseen epitopes." Also: the *source of negative
  examples* substantially changes apparent accuracy.
- **[ePytope-TCR, Cell Genomics 2025;5:100946](https://doi.org/10.1016/j.xgen.2025.100946)**
  — 21 predictors on a controlled benchmark. **Only 6 averaged AUC > 0.6**; best
  was MixTCRpred at **0.63 ± 0.17**. Performance by how well-represented the
  epitope is in public databases: **AUC 0.85 (high) / 0.57 (medium) / 0.59
  (low)** — that gradient *is* the memorisation signature.
- **[arXiv:2606.04994](https://arxiv.org/abs/2606.04994)** (2026) — the sharpest
  result. Deduplication allowing ≤3 CDR3β substitutions removed **40–70% of test
  sequences**. Key diagnosis: *performance degradation is driven by **epitope**
  novelty, not TCR novelty.* When only the TCR was novel, performance held. The
  models memorised epitopes.

**The experiment that matters most for engineering:** 21 TCRs × 172
single-amino-acid mutants of one epitope. **All eight tested models were
indistinguishable from random**, with no significant correlation for any model.
Scoring the effect of a point mutation in a peptide *is* the off-target
screening task — and it is the task these models fail hardest.

Still worth using, for what they are good at:

| Tool | Repo | Good for |
|---|---|---|
| **tcrdist3** | [kmayerb/tcrdist3](https://github.com/kmayerb/tcrdist3) | Distance/clustering/metaclonotypes. Mature, actively maintained, makes no generalisation claim that can fail |
| **SCEPTR** | [yutanagano/sceptr](https://github.com/yutanagano/sceptr) | **MIT**, ~153k params, **runs on CPU**. Embeddings, not epitope lookup |
| **MixTCRpred** | [GfellerLab/MixTCRpred](https://github.com/GfellerLab/MixTCRpred) | Best general method in ePytope-TCR. Architecturally **per-epitope** — cannot address unseen epitopes at all |
| **ePytope-TCR** | [SchubertLab/epytope](https://github.com/SchubertLab/epytope) | **Use this instead of installing 21 models.** Peer-reviewed, one interface over all of them |

### 1.4 Affinity optimisation — what has actually worked

**The reference success is still from 2014.**
[Pierce et al., PLoS Comput Biol 10:e1003478](https://doi.org/10.1371/journal.pcbi.1003478)
redesigned the DMF5 TCR against MART-1/HLA-A2: 12 designed mutations tested,
**best variant ~400-fold tighter, with no detectable binding to non-cognate
ligands**, and a crystal structure of the design ([PDB 4L3E](https://www.rcsb.org/structure/4L3E)).
Twelve years on, nothing ML-based has beaten that *with measured affinities*.
It was Rosetta-style physics on a crystal structure plus wet-lab iteration, and
it treated **specificity as an explicit design constraint**.

The ML picture:

- **[PLoS Comput Biol 2024;20:e1012489](https://doi.org/10.1371/journal.pcbi.1012489)**
  benchmarked ProteinMPNN, ESM-IF1 and Rosetta Design on TCR interfaces.
  Native CDR3-interface recovery: **ESM-IF1 50.1%, ProteinMPNN 43.9%, Rosetta
  31.7%**. Removing the pMHC from the input cost ~20% — *always keep the pMHC as
  context.* **No wet-lab validation**, and the paper found that randomly chosen
  dissimilar sequences sometimes still scored high model confidence, so
  confidence alone is not a filter. Its MM/PBSA step is the most credible
  TCR-specific energetics number available: **r = 0.9 across 7 mutants with
  solved structures, r = 0.74 across 35 modelled mutants** (vs ATLAS).
- **TCRAD** (preprint, Jan 2026) generated de novo CDR3β specificity: 29 designs
  grafted into the 1G4 framework, **8/29 bound tetramer, 5/29 (17%) gave
  antigen-induced activation**. No KD values.
- **Orthogonal TCR α/β interface design** (preprint, Jan 2026) — Rosetta
  multistate design, **~95% correct α/β pairing** validated by MS and
  crystallography. A genuine win, but on *pairing fidelity*, not affinity.

**ΔΔG tooling, ranked by how much TCR-specific validation exists:**

| Method | TCR-specific validation | Notes |
|---|---|---|
| MM/PBSA vs ATLAS | **r = 0.74–0.9**, but n = 7 and n = 35 | Best available. Omitted explicit waters and entropy |
| Rosetta flex ddG | **none found** | Validated on 1,240 general PPI mutations; largest gains on small→large and antibody–antigen. TCR interfaces are flat and low-affinity — the regime where ΔΔG methods do worst |
| FoldX | none found | Fast, free academic licence, closed source. Weakest on charged and buried-polar substitutions |
| SSIPe | not TCR-specific | Pearson 0.61 / RMSE 1.93 kcal/mol on 2,204 SKEMPI mutations |
| ThermoMPNN | ⚠️ **stability, not binding** | Do not use for interface affinity |
| **A TCR-specific ΔΔG ML predictor** | — | **Does not exist.** Closest is a 2021 random forest on ATLAS classifying above/below an affinity threshold at ~75% accuracy. This is an open gap |

### 1.5 The 2025 strategic shift: don't affinity-mature a TCR, design a binder

Three papers in the **same 24 July 2025 issue of *Science*** are the real state
of the art in "bind a pMHC with designed specificity":

- **[Liu et al., Science 389:386](https://doi.org/10.1126/science.adv0185)** —
  RFdiffusion-designed pMHC-I binders; specific binders for **11 target pMHCs**;
  built into CARs, designs for 8 targets conferred peptide-specific T cell
  activation.
- **[Johansen et al., Science 389:380](https://doi.org/10.1126/science.adv0422)** —
  RFdiffusion + ProteinMPNN + AF2 minibinders against NY-ESO-1/HLA-A*02.
  Notable for **"in silico cross-panning"** as a computational specificity
  prescreen, and for MD improving in-vitro predictability. Cryo-EM confirmed;
  killed NY-ESO-1⁺ melanoma as a CAR.
- **[Householder et al., Science 389:375](https://doi.org/10.1126/science.adv3813)** —
  de novo α-helical **TCR mimic**, **KD = 9.5 nM**, 2.05 Å structure, rigid
  TCR-like docking focused on up-facing peptide side chains. And: a
  structure-informed in-silico screen of **14,363 HLA-A*02 peptides correctly
  predicted two off-target peptides.**

**There is no released TCR-specific diffusion model.** RFantibody
([RosettaCommons/RFantibody](https://github.com/RosettaCommons/RFantibody),
*Nature* 2025) is antibodies and nanobodies only — no TCR mode. If your goal is
potent, specific pMHC recognition rather than a natural αβ TCR per se, the
designed-binder route is where the tooling is mature and the wet-lab evidence
is strongest.

### 1.6 Off-target screening — the step that is not optional

The affinity-enhanced MAGE-A3 "A3A" TCR caused **fatal cardiac toxicity in two
patients** via cross-reactivity with a titin peptide
([Linette et al., Blood 2013](https://pubmed.ncbi.nlm.nih.gov/23926201/)). The
details that should shape any screening workflow:

- The titin peptide was found **only on retrospective review**.
- Titin was **undetectable in cardiac-derived primary cell lines** and showed up
  only in a **beating** myocyte culture.
- **Mouse models lack the equivalent reactive peptide** — animal models cannot
  catch this class of problem.
- What did find it: a **Gly/Ala scan** of the receptor → a **degenerate
  recognition motif** → a **directed in-silico proteome search** → functional
  confirmation against beating iPSC-derived cardiomyocytes.

That three-step pattern remains best practice. **No formal regulatory standard
or consensus guideline for in-silico TCR off-target screening appears to exist.**

| Tool | What it does | Reported performance |
|---|---|---|
| **CrossDome** | [AntunesLab/crossdome](https://github.com/AntunesLab/crossdome) — R package; **TCR-centred mode weights positions by contact hotspots** | Moved the titin peptide from **rank 27 → rank 6 of 36,000**. Across 16 known cases, top-50 enrichment **63% → 82%** with the TCR-centred mode |
| **Yeast-display pHLA libraries** | ~10⁸ 9-mers | **Identified titin reactivity for A3A** — would have caught the fatality prospectively |
| **Structure-informed in-silico screen** | Screen the designed binder against the HLA-matched peptide universe | 14,363 peptides → **correctly predicted both real off-targets** (Householder 2025) |
| pMHC-pseudotyped virus + yeast display | Scales to many TCRs | 47M mutant TCRs × 92 pMHCs (preprint) |

**Do not use a sequence-based TCR–epitope predictor for off-target screening.**
See §1.3: on peptide point mutations these models are at chance.

### 1.7 End-to-end frameworks

Nothing stitches the whole stack together. The partial ones worth knowing:
**ePytope-TCR** (specificity layer, peer-reviewed, best-engineered thing in this
survey), **NetTCRfold-1.0** (structure layer at scale, new/immature),
**TCRAD** (design layer, preprint), **[Aseamann/TRain](https://github.com/Aseamann/TRain)**
(sequence → predicted complexes, *BMC Bioinformatics* 2025), and **Stitchr**
(V/J/CDR3 → spliced cDNA; mature, narrow).

### 1.8 Data resources

| Resource | Content | Note |
|---|---|---|
| **ATLAS** ([atlas.ibbr.umd.edu](https://atlas.ibbr.umd.edu)) | **WT + mutant TCR–pMHC affinities tied to structures** | **The key resource for ΔΔG work.** Last updated 2023-08-15; two of its own pages report different entry counts (988 vs 1,174) — verify live |
| **TCR3d 2.0** ([tcr3d.ibbr.umd.edu](https://tcr3d.ibbr.umd.edu)) | >1,330 pMHC structures (1,098 class I + 235 class II) + affinities + docking tools | Weekly updates. [NAR 2025](https://doi.org/10.1093/nar/gkae840) |
| **STCRDab** ([OPIG](https://opig.stats.ox.ac.uk/webapps/stcrdab-stcrpred/About)) | Auto-curated TCR structures, IMGT-numbered | Weekly |
| **VDJdb** | TCR→antigen specificity | 61,049 records in the Jul-2019 release; **current size unverified**. CC BY-NC 4.0 per the 2019 paper |
| **IEDB**, **McPAS-TCR** | Epitopes; TCR↔pathology | McPAS size/licence not verified |
| **10x dextramer dataset** | Dextramer-labelled CD8⁺ cells | ⚠️ **Both ICON and ITRAP flag ~75% of the raw data as noise**; only 6 pMHCs have ≥500 receptors. Treat as hazardous |

### 1.9 Can you computationally design a higher-affinity TCR and have it work?

**Yes — via structure-based physics with experimental iteration, and
specificity does not come for free.** Pierce 2014's ~400-fold DMF5 improvement
is the standing proof. ML inverse folding reaches 44–50% native recovery on TCR
interfaces with no wet-lab check; the one ML design paper with wet-lab data
reports a 17% functional hit rate and no KDs. And the field's own 2025 answer
was to sidestep native TCR maturation entirely in favour of de novo binders.

Carry forward Pierce 2014's warning: **very high affinity can reduce T cell
activity**, and affinity maturation is exactly what produced the titin fatality.
Design the §1.6 screening cascade in from the start.

---

## Part 2 — the demo

A CPU-only pipeline over a deposited TCR–pMHC structure. No GPU, no model
weights, ~14 seconds end to end.

```bash
pip install biopython matplotlib numpy pandas
python -m tcr_design run 1AO7 --outdir results/1ao7
```

### What it does

| Stage | What it computes | Honesty label |
|---|---|---|
| 1. Interface | Chain typing, CDR loops, contact maps, Shrake–Rupley burial, side-chain void probing | **measured** — deterministic geometry of the input structure |
| 2. Scan | Saturation mutagenesis over contacting CDR positions, six-term component score | **heuristic** — uncalibrated, arbitrary units |
| 2b. Calibrated ΔΔG | Ridge model over measured structural descriptors, fitted to ATLAS + SKEMPI | **calibrated but weak** — held-out Spearman +0.43, RMSE 0.88 kcal/mol |
| 3. Gates | Peptide-vs-MHC contact partitioning; developability liability motifs | **heuristic** — flags risk, does not measure it |
| 4. Off-target | Contact-weighted epitope motif → HLA anchor filter → reviewed human proteome | **heuristic** — a similarity shortlist, not a prediction |
| 5. Handoff | AlphaFold3 / Boltz / Rosetta resfiles / FoldX / ProteinMPNN inputs + exact commands | **needs external tool** |

Chain typing uses constant-domain sequence fingerprints, falling back to the
conserved diagonal docking polarity (TCRα over the MHC α2 helix, β over α1) when
a construct has no constant domain. CDR3 uses the IMGT junction definition
(conserved Cys → J-region FGXG); CDR1 is anchored on Cys23/Trp41. On 1AO7 this
recovers the published A6 loops exactly: CDR3α `CAVTTDSWGKLQF`, CDR3β
`CASRPGLAGGRPEQYF`, CDR1α `DRGSQS`, CDR1β `MNHEY`.

### The two design decisions worth explaining

**1. The specificity gate is kept out of the score.** A TCR buries a lot of MHC
surface, so most of the "best" single mutations buy affinity at MHC-contacting
positions — which is peptide-independent recognition, the documented route to
off-target toxicity. On 1AO7, **195 of 418 scanned mutants are rejected on this
basis**, and that is the main result of the scan. Blending specificity into a
single score would let a large predicted affinity gain average it away.

**2. Peptide-reading positions are ranked against each other.** Because MHC
contacts dominate a global ranking, mutations at positions that touch the
peptide get a separate z-score computed within that subset.

### The off-target screen

Rather than comparing whole peptides, positions are weighted by **measured TCR
contact**, so the screen asks "what else in the proteome looks like this *at the
positions this receptor reads*" — the same principle as CrossDome's TCR-centred
mode. On 1AO7 the contact profile correctly identifies **Tyr5 of the Tax
peptide (92 contacts)** as the dominant read-out position and the C-terminal
Val9 (0 contacts) as an MHC anchor the TCR never touches.

Requiring identity at the hottest position and HLA-A*02:01 anchors at P2/P9
reduces 11.25M proteome windows to a ranked shortlist; the top human mimic of
`LLFGYPVYV` is `LLNNYPVYI` (BEND4), matching 6 of 8 weighted positions including
both tyrosines.

### Calibration against measured ΔΔG

```bash
python -m tcr_design benchmark --save-model
```

This pulls both public sources of measured TCR–pMHC mutation affinities,
scores every mutation structurally, and fits a model — validated
**leave-one-structure-out**, so it never sees the receptor it is scored on.

| source | entries |
|---|---|
| [ATLAS](https://atlas.ibbr.umd.edu) — TCR-side single mutants with ΔΔG | 273 |
| [SKEMPI 2.0](https://life.bsc.es/pid/skempi2/) — `TCR/pMHC` subset, single mutants | 189 |
| overlapping between the two | 94 |
| **pooled unique** | **368** |
| mapped onto structures and scored | 261 (25 structures) |

ΔΔG = RT·ln(Kd_mut/Kd_wt), positive meaning weaker binding.

**Result 1 — the hand-weighted score is worthless.**

| | Spearman vs ΔΔG | AUC, finding stabilising |
|---|---|---|
| shipped heuristic | −0.082 | 0.542 |

(A working score gives a *negative* Spearman; AUC above 0.5 is correct.) Per
structure it ranges from −0.78 to +0.75 — i.e. noise.

**Result 2 — refitting the six weights does not help.** Held-out Spearman
+0.125, RMSE 0.980 kcal/mol against 0.979 for predicting the training mean. The
six hand-crafted terms carry essentially no transferable signal, and no choice
of weights rescues them.

**Result 3 — the raw structural measurements do carry signal.** Fitting ridge
regression on 20 measured descriptors (burial, ΔSASA, contact counts, void,
partner charge/polarity, Δvolume, Δhydropathy, and so on) instead of the six
derived terms:

| | held out | in sample |
|---|---|---|
| Spearman | **+0.429** | +0.550 |
| Pearson | **+0.438** | +0.587 |
| RMSE (kcal/mol) | **0.884** | 0.795 |
| AUC, finding stabilising | **0.718** | — |

versus 0.979 RMSE for predicting the mean. The fitted coefficients are
physically sensible — larger side chain → tighter; mutating a buried,
peptide-contacting position → weaker; mutating away from glycine → weaker —
which is a reasonable sign it learned energetics rather than noise.

For context: Rosetta flex ddG reaches R ≈ 0.6–0.7 on *general* protein
interfaces (no TCR-specific benchmark exists), and MM/PBSA reached r = 0.74–0.9
against ATLAS on n = 7 and n = 35. This model is weaker than both — and it runs
in CPU-seconds, on 25 structures, with honest leave-one-structure-out numbers.

The model ships at `tcr_design/data/calibrated_ddg_model.json` and the pipeline
uses it automatically, reporting predicted ΔΔG in kcal/mol. `--no-calibration`
falls back to the raw heuristic.

### The benchmark both models still fail

```bash
python -m tcr_design validate 3QDG 4L3E --calibrated
```

`3QDG` is wild-type DMF5; `4L3E` is the Pierce 2014 designed variant at
~400-fold improved affinity. The pipeline diffs the structures, recovers the
two substitutions (**αD26Y** in CDR1α, **βL98W** in CDR3β), scans the wild type
blind, and looks up where they ranked. With `--calibrated` the model is
**refitted without 3QDG and 4L3E** (243 mutations from 23 other structures), so
there is no leakage.

| mutation | uncalibrated rank | calibrated rank | rank at its own position | predicted ΔΔG | measured ΔΔG |
|---|---|---|---|---|---|
| `D:D26Y` | 237/456 (52%) | **112/456 (25%)** | 15/19 → **5/19** | +0.14 | **−1.79** |
| `E:L98W` | 174/456 (38%) | **117/456 (26%)** | 5/19 → **1/19** | +0.16 | — |

(The DMF5 double mutant αD26Y + βL98W is recorded in ATLAS at **−3.28 kcal/mol**;
SKEMPI records the αD26Y single at −1.79.)

Calibration roughly **halves** the percentile, and βL98W becomes the
**top-ranked substitution of 19** at its position. That is real improvement in
the question a designer actually faces: *given this position, which residue?*

But both are still predicted to **weaken** binding (+0.14, +0.16 kcal/mol) when
αD26Y in fact tightens it by −1.79. Neither clears the threshold, so neither
would be shortlisted. **The verdict is still `FAILED`.**

The reason is visible in the middle panel of `results/calibration/calibration.png`:
predictions are compressed into roughly 0 to +2 kcal/mol while the measurements
span −2 to +8. The benchmark is dominated by alanine scanning, so the model
mostly learned *deleting a buried side chain weakens binding*. Adding aromatic
volume to gain affinity is the rarer direction, and it is exactly the direction
design cares about.

The compression is easy to quantify on the A6 CDR1α position-26 series, where
ATLAS has four measurements. **1AO7 is in the training set, so this is in-sample
and not evidence of generalisation** — but it shows the failure mode cleanly:

| mutation | predicted ΔΔG | measured ΔΔG |
|---|---|---|
| `D26W` | −0.14 | **−1.08** |
| `D26M` | +0.11 | +0.19 |
| `D26V` | +0.23 | +1.45 |
| `D26A` | +0.49 | +1.31 |

The *ordering* is essentially right. The *range* is about four-fold too narrow,
which is precisely why nothing ever clears a threshold set in real kcal/mol.

One encouraging sign, on an unseen receptor: the calibrated model's best picks
for the A6 TCR are aromatic insertions at peptide-contacting CDR positions
(`T93W`, `G28W`, `G101W`, `T93F/Y`, `R95W`) — the same chemistry as both
validated designs in this field. It is pointing at the right kind of change
while being unable to say how big the change is.

So, honestly:

- **Calibration bought real ranking skill within a position, and no ability to
  recognise an affinity-improving mutation as such.**
- **The specificity gate's contact-partitioning proxy is crude.** αD26Y contacts
  only MHC in the wild-type structure and is rejected on that basis, yet the
  published variant retained its specificity.

This is why the package remains a funnel into AlphaFold3 / Rosetta / FoldX /
MM-PBSA rather than a predictor. To extend the benchmark with your own
measurements, `validate.benchmark_from_csv` takes a CSV of chain/position/mutant/ΔΔG.

### Does stacking mutations add up?

```bash
python -m tcr_design epistasis
```

Greedy combination — take the best singles, stack them, assume the free
energies add — is the standard move in computational affinity maturation, and
it is an assumption that ATLAS can test. Wherever a measured multi-mutant
decomposes into disjoint subsets that were *also* measured, the sum of the
parts can be compared with the measured whole. 51 such decompositions exist
across six TCR systems; 13 have every component measured in a single study,
which controls for the fact that wild-type affinities differ between papers
(DMF5 against ELAGIGILTV is reported at 5.6, 9.5, 11 and 36.4 µM — about
1 kcal/mol of baseline noise).

| mutations stacked | n (within one study) | mean error | within 0.5 kcal/mol |
|---|---|---|---|
| **2** | 6 | **+0.00** | **6/6** |
| 4 | 2 | −0.73 | 1/2 |
| 6 | 1 | −0.32 | 1/1 |
| 8 | 3 | **+1.60** | 1/3 |
| 9 | 1 | **+3.38** | 0/1 |

Error is measured minus predicted, so positive means the stack underdelivered.

**Pairs add up; deep stacks do not.** Every two-mutation combination measured
within one study lands inside 0.5 kcal/mol of the sum of its parts. By eight or
nine simultaneous mutations the measured gain falls short by 1.6 to 3.4
kcal/mol — the diminishing return of an interface that is already optimised.
The size-dependence survives the within-study restriction, so it is not an
artefact of baseline noise.

This is why the published deep variants — 1G4 reaches 14 simultaneous
mutations and 26 pM — were found by directed evolution rather than by
computational accumulation. Practically: predicting a double by summing two
singles is sound, triples are borderline and should be checked, and beyond
about four mutations summation is worthless.

#### A note on reaching ATLAS

ATLAS serves a valid certificate but omits its InCommon intermediate, so a
default trust store cannot build the chain and the download fails. The package
repairs this properly — it fetches the intermediate from the address named in
the certificate's own Authority Information Access extension and appends it to
the bundle (`benchmark.build_atlas_ca_bundle`). TLS verification stays on; no
`-k`, no `verify=False`.

### Commands

```bash
python -m tcr_design run 1AO7 --outdir results/1ao7          # full pipeline
python -m tcr_design interface 1AO7                          # geometry only
python -m tcr_design scan 1AO7 --loops A-CDR3,B-CDR3         # restrict to CDR3s
python -m tcr_design run 1AO7 --skip-offtarget               # skip the proteome download
python -m tcr_design run 1AO7 --no-calibration               # raw heuristic only
python -m tcr_design benchmark --save-model                  # calibrate against ATLAS + SKEMPI
python -m tcr_design epistasis                               # does stacking mutations add up?
python -m tcr_design validate 3QDG 4L3E --calibrated         # leakage-free retrospective test
python -m tcr_design run 5XOV --chains mhc=A,peptide=C,tcr_alpha=D,tcr_beta=E
```

Verified to run on 1AO7 (A6/Tax/HLA-A2), 3QDG and 4L3E (DMF5/MART-1).

### Outputs

```
results/calibration/
  CALIBRATION.md             # the full calibration write-up
  calibration_data.csv       # every benchmark mutation, features, held-out prediction
  calibration.png            # unfitted vs fitted vs per-structure
results/1ao7/
  REPORT.md                  # labelled markdown report
  summary.json
  interface_residues.csv     # per-residue measurements
  mutation_scan.csv          # all mutants, all score components, verdicts
  offtarget_candidates.csv
  interface_profile.png      # peptide read-out + paratope contact split
  mutation_heatmap.png       # position × amino acid
  gate_outcomes.png
  tier2/
    RUN_TIER2.md             # exact commands + caveats for each tool
    alphafold3/*.json  boltz/*.yaml  rosetta/*.resfile
    foldx/individual_list.txt  proteinmpnn_*.jsonl
```

### Limits

- Single point mutations only; real maturation explores combinations.
- One static structure — no side-chain or loop rearrangement on binding, and
  CDR3 loops are flexible.
- No expression, pairing, stability or immunogenicity assessment.
- The proteome screen covers only the canonical reviewed proteome: no isoforms,
  modifications, non-canonical products or tissue-restricted variants. The HLA
  anchor filter is a coarse motif, not NetMHCpan.
- Higher affinity is not automatically better.

### Tests

```bash
python -m pytest tests/test_tcr_design.py -q      # 36 tests
```

Unit tests run offline. Integration tests skip cleanly without cached
structures, and assert against published ground truth: the A6 CDR loops, the
TCR–pMHC buried-surface range, Tyr5 as the dominant Tax contact, and recovery
of the two DMF5 design mutations.
