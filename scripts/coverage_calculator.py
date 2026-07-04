#!/usr/bin/env python3
"""Scale & sequencing-depth calculator for a pooled CRISPR knockout screen.

Computes, from a handful of biological parameters, everything you need to size
a kinome (or any) pooled screen so that library *representation* (cells/guide)
is preserved at every bottleneck: transduction, selection, passaging, the gDNA
PCR template, and sequencing reads.

Example:
    python3 scripts/coverage_calculator.py --genes 518 --guides-per-gene 6 \
        --controls 400 --coverage 1000 --moi 0.3 --replicates 3 --arms 2 \
        --timepoints 2

All model assumptions are printed so the numbers are auditable.
"""
import argparse
import math


def poisson_infected_fraction(moi: float) -> float:
    """Fraction of cells receiving >=1 virion under Poisson statistics."""
    return 1.0 - math.exp(-moi)


def human_number(x: float) -> str:
    """Format a large number compactly (e.g. 1.5e7 -> '1.50e+07  (15.0 M)')."""
    if x >= 1e6:
        return f"{x:.2e}  ({x/1e6:.1f} M)"
    if x >= 1e3:
        return f"{x:.2e}  ({x/1e3:.1f} k)"
    return f"{x:.0f}"


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--genes", type=int, default=518,
                   help="number of target genes (default 518 protein kinases)")
    p.add_argument("--guides-per-gene", type=int, default=6)
    p.add_argument("--controls", type=int, default=400,
                   help="total control guides (NTC + intergenic + benchmarks)")
    p.add_argument("--coverage", type=int, default=1000,
                   help="cells (and reads) per guide to maintain (default 1000x)")
    p.add_argument("--moi", type=float, default=0.3,
                   help="functional MOI used for transduction (default 0.3)")
    p.add_argument("--reads-per-guide", type=int, default=1000,
                   help="target sequencing reads per guide (default 1000)")
    p.add_argument("--read-overshoot", type=float, default=1.25,
                   help="multiplier on reads to absorb non-uniformity (default 1.25)")
    p.add_argument("--ploidy-pg", type=float, default=6.6,
                   help="gDNA mass per cell in pg (human diploid ~6.6)")
    p.add_argument("--max-gdna-per-pcr-ug", type=float, default=10.0,
                   help="max gDNA (ug) per PCR reaction (default 10)")
    p.add_argument("--replicates", type=int, default=3)
    p.add_argument("--arms", type=int, default=1,
                   help="experimental arms (e.g. drug vs DMSO = 2)")
    p.add_argument("--timepoints", type=int, default=2,
                   help="sequenced timepoints per arm (e.g. T0 + Tend = 2)")
    p.add_argument("--include-plasmid-t0", action="store_true", default=True,
                   help="add a plasmid QC sample to the sample count")
    args = p.parse_args()

    lib = args.genes * args.guides_per_gene + args.controls
    integrants = lib * args.coverage                       # cells carrying a guide
    infected_frac = poisson_infected_fraction(args.moi)
    cells_to_plate = integrants / infected_frac            # before selection
    maintain = integrants                                  # per passage / timepoint

    # gDNA per sequenced sample: enough genome copies to preserve coverage.
    genomes_needed = lib * args.coverage
    gdna_ug = genomes_needed * args.ploidy_pg / 1e6        # pg -> ug
    n_pcr = math.ceil(gdna_ug / args.max_gdna_per_pcr_ug)

    reads_per_sample = lib * args.reads_per_guide * args.read_overshoot

    # Sample accounting.
    screen_samples = args.arms * args.timepoints * args.replicates
    total_samples = screen_samples + (1 if args.include_plasmid_t0 else 0)
    total_reads = reads_per_sample * total_samples

    print("=" * 68)
    print(" POOLED CRISPR SCREEN — SCALE & DEPTH")
    print("=" * 68)
    print("\n-- Library --")
    print(f"  genes ................... {args.genes}")
    print(f"  guides/gene ............. {args.guides_per_gene}")
    print(f"  control guides .......... {args.controls}")
    print(f"  LIBRARY SIZE ............ {lib:,} guides")

    print("\n-- Per-sample cell requirements (coverage = "
          f"{args.coverage}x) --")
    print(f"  integrants needed ....... {human_number(integrants)}")
    print(f"  infected fraction @ MOI {args.moi} = {infected_frac:.2%}")
    print(f"  cells to transduce ...... {human_number(cells_to_plate)}")
    print(f"  MAINTAIN >= ............. {human_number(maintain)} cells "
          "at every passage & timepoint")

    print("\n-- gDNA readout (per sequenced sample) --")
    print(f"  genome copies to sample . {human_number(genomes_needed)}")
    print(f"  gDNA required ........... {gdna_ug:.1f} ug "
          f"(@ {args.ploidy_pg} pg/cell)")
    print(f"  PCR1 reactions .......... {n_pcr} "
          f"(<= {args.max_gdna_per_pcr_ug:.0f} ug each), pooled")

    print("\n-- Sequencing --")
    print(f"  reads/sample ............ {human_number(reads_per_sample)} "
          f"({args.reads_per_guide}/guide x {args.read_overshoot} overshoot)")
    print(f"  samples ................. {total_samples} "
          f"({args.arms} arms x {args.timepoints} timepoints x "
          f"{args.replicates} reps"
          f"{' + plasmid' if args.include_plasmid_t0 else ''})")
    print(f"  TOTAL READS ............. {human_number(total_reads)}")
    print(f"  read length ............. 1x75 bp SE + index; +10-20% PhiX spike")

    print("\n-- Assumptions --")
    print("  * Poisson single-integrant model; MOI ~0.3 keeps most infected")
    print("    cells at one guide.")
    print("  * Human diploid gDNA ~6.6 pg/cell; use ALL gDNA as PCR template.")
    print("  * Coverage must hold at EVERY bottleneck, not just the average.")
    print("=" * 68)


if __name__ == "__main__":
    main()
