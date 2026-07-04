"""Latch workflow: kinome guide off-target + on-target scoring (FlashFry).

Wraps scripts/flashfry_pipeline.py so the two pending guide-QC fields
(on-target + CFD off-target specificity) can be filled on Latch. Inputs are the
30-mer context table produced by the design step and a genome FASTA; output is a
scored, filtered guide table.

Register:  latch register latch_workflow --remote
Then it appears in workspace 42942 and can be launched by name/id.
"""
from pathlib import Path
import subprocess

from latch import large_task, workflow
from latch.types import (
    LatchFile,
    LatchDir,
    LatchOutputDir,
    LatchMetadata,
    LatchParameter,
    LatchAuthor,
)

metadata = LatchMetadata(
    display_name="Kinome Guide Off-Target Scoring (FlashFry)",
    documentation="https://github.com/mckenna-lab/FlashFry",
    author=LatchAuthor(name="Kinome CRISPR screen design"),
    parameters={
        "context_tsv": LatchParameter(
            display_name="Guide context table (rs3_context.tsv)",
            description="TSV with columns id, gene, spacer, pam, context_30mer, source. "
            "Produced by the design step (data/rs3_context.tsv).",
        ),
        "genome_fasta": LatchParameter(
            display_name="Genome FASTA (hg38)",
            description="Reference FASTA used to build the FlashFry off-target database.",
        ),
        "cfd_min": LatchParameter(
            display_name="Min CFD specificity to PASS",
            description="Guides with CFD specificity below this are flagged FAIL (default 0.2).",
        ),
        "output_dir": LatchParameter(
            display_name="Output directory",
            description="Latch Data folder for guides_scored.tsv.",
        ),
    },
)


@large_task
def score_task(
    context_tsv: LatchFile,
    genome_fasta: LatchFile,
    output_dir: LatchOutputDir,
    cfd_min: float = 0.2,
) -> LatchFile:
    local_out = Path("/root/guides_scored.tsv")
    subprocess.run(
        [
            "python3", "/root/scripts/flashfry_pipeline.py", "run",
            "--genome", genome_fasta.local_path,
            "--context", context_tsv.local_path,
            "--out", str(local_out),
            "--ff", "/opt/FlashFry.jar",
            "--cfd-min", str(cfd_min),
            "--workdir", "/root/ff_work",
        ],
        check=True,
    )
    return LatchFile(str(local_out), f"{output_dir.remote_path}/guides_scored.tsv")


@workflow(metadata)
def kinome_offtarget_wf(
    context_tsv: LatchFile,
    genome_fasta: LatchFile,
    output_dir: LatchOutputDir,
    cfd_min: float = 0.2,
) -> LatchFile:
    """Score a fixed kinome guide list for off-target (CFD) and on-target (Doench 2016).

    For each guide in `context_tsv`, FlashFry builds an index from `genome_fasta`,
    computes on-target and CFD off-target specificity, and the workflow writes a
    merged, filtered table (PASS if CFD specificity >= `cfd_min`).
    """
    return score_task(
        context_tsv=context_tsv,
        genome_fasta=genome_fasta,
        output_dir=output_dir,
        cfd_min=cfd_min,
    )
