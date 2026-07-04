# Library Construction

Turn the selected spacer list into a cloned, QC'd, pooled lentiviral plasmid
library.

## 1. Vector choice

| Strategy | Vector(s) | When |
|----------|-----------|------|
| **Two-vector (recommended)** | Cell line stably expressing **Cas9** (e.g. lentiCas9-Blast) + **lentiGuide-Puro** (sgRNA only) | Higher titer, higher/typically more uniform editing, best for screens |
| One-vector | **lentiCRISPRv2** (Cas9 + sgRNA + Puro) | When you can't pre-make a Cas9 line; lower titer |

Recommended: build a **Cas9-expressing derivative of your screen cell line**
first, validate Cas9 activity (e.g., a reporter or a control guide against a
known essential), then clone the kinome guides into **lentiGuide-Puro**.

## 2. Oligo pool design

Each library member is a synthesized oligo:

```
[5' PCR/amplification adapter] - [BsmBI site] - [G] - [20-nt spacer] - [BsmBI site] - [3' adapter]
```

- BsmBI (Esp3I) Golden-Gate sites are oriented to drop the 20-mer into the
  vector's sgRNA scaffold and are removed during cloning.
- Add a leading `G` if not present (U6 transcription).
- Keep all oligos the same length; the variable region is only the 20-mer.
- Order as an **array-synthesized oligo pool** (Twist, Agilent/SurePrint,
  GenScript). For ~5,000 members this is a single small pool.

## 3. Cloning workflow

1. **Amplify** the pool by limited-cycle PCR (minimize cycles to avoid skew and
   recombination; use a high-fidelity polymerase).
2. **Golden-Gate assembly** (BsmBI/Esp3I) into the digested vector — one-pot
   digest+ligate. (Gibson is an alternative if using a linearized backbone.)
3. **Electroporate** into high-efficiency recombination-deficient *E. coli*
   (Endura / Lucigen). Plate dilutions to **count colonies**.
4. **Coverage during cloning:** collect **≥ 50–100× colonies per guide**
   (e.g., ≥ 250k–500k colonies for a 5k library) to preserve representation.
   Scrape all colonies together; **do not pick individually**.
5. **Maxiprep** the pooled plasmid.

## 4. Plasmid library QC (before making virus)

- **NGS the plasmid pool** (same amplicon readout as the screen). Compute:
  - **Representation:** ≥ 95% of designed guides detected.
  - **Skew:** top-10% : bottom-10% read ratio **< 10×** (ideally < 6×).
  - **Gini index** of guide counts **< 0.10** (tight distribution).
  - No large dropout of any gene's full guide set.
- Sanger-spot-check a few colonies for correct scaffold junctions.

## 5. Lentivirus production and titering

1. Co-transfect the plasmid library with packaging plasmids (psPAX2 + pMD2.G)
   into HEK293T; harvest and concentrate virus.
2. **Titer on the actual screen cell line** by transducing a dilution series
   and measuring puromycin-resistant fraction → compute functional MOI.
3. Screen at **MOI ≈ 0.3** (so ~1 integrant/cell by Poisson; ~26% of cells
   infected, and among infected cells the vast majority carry a single guide).

## Deliverables from this stage

- Cloned, sequence-verified plasmid library (maxiprep).
- Plasmid-QC NGS run (**FASTQ referenced by name**, e.g.
  `plasmid_qc_R1.fastq.gz` — large, not linked in summary) + a small guide-count
  table and QC-metrics report you can link.
- Titered, aliquoted lentiviral stock with measured functional titer.
