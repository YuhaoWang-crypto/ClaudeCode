# Screen Scale, Coverage & Sequencing Depth

The single most important experimental parameter is **representation
(coverage)** — cells per guide — maintained at every step. Under-sampling at
any point (transduction, selection, passaging, gDNA, or reads) permanently
loses guides to bottleneck drift and creates false hits. Use
`scripts/coverage_calculator.py` to get exact numbers for your library.

## 1. Library size (worked default)

```
518 protein kinases × 6 guides           = 3,108
+ ~20 lipid kinases × 6                   =   120
+ 250 non-targeting controls             =   250
+ 50 intergenic cutting controls         =    50
+ ~100 essential/non-essential benchmark =   100
--------------------------------------------------
≈ 3,600–3,700 guides  (round planning to ~4,000; ~5,000 if 8 guides/gene)
```

## 2. Coverage (cells per guide)

- **Minimum 500×**, **1000× recommended** for dropout/essentiality screens
  (kinome screens are usually dropout-oriented, which are more bottleneck-
  sensitive than strong positive selections).
- Coverage must be held at **every** stage: number of cells transduced (as
  integrants), cells carried through selection, cells seeded at each passage,
  cells harvested per timepoint, genome copies in the gDNA PCR, and reads.

For a 4,000-guide library at 1000×:
- Integrants needed = 4,000 × 1000 = **4.0 × 10⁶ cells** carrying a guide.
- At MOI 0.3 (~26% infected), cells to plate = 4.0e6 / 0.26 ≈ **1.5 × 10⁷**.
- Maintain **≥ 4.0 × 10⁶ cells at every passage and per timepoint**.

## 3. MOI and infection

- **MOI ≈ 0.3**, low enough that most infected cells carry a single guide
  (avoids multi-guide confounding).
- Puromycin-select 2–3 days post-transduction; the surviving fraction should
  match the predicted MOI (a titer sanity check).

## 4. Timepoints and replicates

| Sample | Timing | Purpose |
|--------|--------|---------|
| **Plasmid** | pre-virus | Library QC / input reference |
| **T0** | just after selection (day ~3–7) | Reference for computing depletion/enrichment |
| **T_end** | ~day 14–21 (≈ 8–12 population doublings) | Dropout of essential kinases |
| **Condition arms** | drug vs. DMSO, etc. | Chemogenomic / synthetic-lethal readout |

- **≥ 3 biological replicates** (independent infections) per arm. Replicate
  correlation is a primary QC.
- Let the screen run **enough doublings** (≥ ~8–10) for essential-gene dropout
  to separate from noise, without so many that the library over-drifts.

## 5. gDNA readout and how much you need

You must PCR from **enough genome copies to preserve coverage** — this is a
common failure point. Human diploid gDNA ≈ **6.6 pg/cell**.

- Genome copies needed = library_size × coverage
  (= 4,000 × 1000 = 4.0 × 10⁶ genomes at 1000×).
- gDNA mass = genomes × 6.6 pg = 4.0e6 × 6.6 pg ≈ **26 µg** per sample.
- Use **all** of that gDNA as template, split across **multiple parallel PCR
  reactions** (~10 µg per 100 µL reaction is a typical ceiling), then pool.
- Use a two-step PCR: PCR1 amplifies the integrated sgRNA cassette; PCR2 adds
  Illumina adapters + sample **indices** (with staggered/variable-length primers
  to add nucleotide diversity for the sequencer).

## 6. Sequencing depth (reads)

- Target **~500–1000 reads per guide** (1000× recommended).
  - 4,000 guides × 1000 reads = **4 × 10⁶ reads/sample** → round up to
    ~5M reads/sample to absorb non-uniformity.
  - 5,000-guide library → ~5–6M reads/sample.
- **Read length:** single-end **1×75 bp** is sufficient to read the 20-nt
  spacer; add index reads for demultiplexing.
- **Instrument:** a single NextSeq/NovaSeq lane yields hundreds of millions of
  reads → you can multiplex **many** samples. Example: 3 arms × 3 reps × 2
  timepoints + plasmid + T0 ≈ 20 samples × 5M ≈ 100M reads — a small fraction
  of one flow cell.
- **Add a PhiX spike-in (~10–20%)** to offset the low base diversity of
  amplicon libraries.

## 7. Quick reference table (1000× coverage)

| Library size | Integrants (1000×) | Cells to plate (MOI 0.3) | gDNA/sample | Reads/sample (1000/guide) |
|--------------|--------------------|--------------------------|-------------|---------------------------|
| 4,000 | 4.0 × 10⁶ | ~1.5 × 10⁷ | ~26 µg | ~4–5 M |
| 5,000 | 5.0 × 10⁶ | ~1.9 × 10⁷ | ~33 µg | ~5–6 M |

Run `scripts/coverage_calculator.py` to regenerate these for your exact
parameters (coverage, MOI, ploidy, cell mass, replicates, timepoints).
