"""Species MHC-II panel construction.

Dog: IPD-MHC ships a curated DLA class II set (DRB1, DQA1, DQB1), maintained by
the ISAG nomenclature committee, so the canine panel is an actual allele panel.

Cat: IPD-MHC has *no* feline entries. The FLA panel has to be assembled from
GenBank/RefSeq records, which means no official allele names, no curation, and
no way to know what fraction of feline diversity is covered. That asymmetry is
a result of this module, not a footnote -- it is reported in the panel summary.

Neither species has an allele *frequency* database, so a panel here is a
diversity sample, never a "% of the population covered" claim.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Sequence

from . import data, groove
from .groove import Molecule

# Headers that must never enter a class II beta panel.
_EXCLUDE = re.compile(
    r"class I antigen|class I histocompat|transactivator|invariant chain|"
    r"alpha chain|DRA|DQA|DMA|DOA|CD74|beta-2|B2M",
    re.I,
)


@dataclass
class LocusPanel:
    locus: str
    species: str
    molecules: List[Molecule] = field(default_factory=list)
    n_records: int = 0
    n_rejected: int = 0
    curated: bool = True
    source: str = ""

    @property
    def n_unique_grooves(self) -> int:
        return len({m.pseudoseq for m in self.molecules})

    def summary(self) -> Dict[str, object]:
        return {
            "locus": self.locus,
            "species": self.species,
            "records_considered": self.n_records,
            "records_rejected": self.n_rejected,
            "molecules": len(self.molecules),
            "unique_grooves": self.n_unique_grooves,
            "curated_nomenclature": self.curated,
            "source": self.source,
        }


def _accept(mol: Molecule, min_identity: float, min_coverage: float) -> bool:
    return mol.coverage >= min_coverage and mol.identity_to_ref >= min_identity


def from_ipd(species: str, locus: str, prefixes: Sequence[str],
             min_identity: float = 0.40, min_coverage: float = 0.90) -> LocusPanel:
    alleles = data.ipd_alleles(prefixes)
    panel = LocusPanel(locus=locus, species=species, curated=True,
                       source=f"IPD-MHC {', '.join(prefixes)}")
    for name, seq in sorted(alleles.items()):
        panel.n_records += 1
        mol = groove.build_molecule(name, species, locus, seq, source="IPD-MHC")
        if _accept(mol, min_identity, min_coverage):
            panel.molecules.append(mol)
        else:
            panel.n_rejected += 1
    return panel


def from_ncbi_feline(locus: str, include: Sequence[str],
                     min_identity: float = 0.40, min_coverage: float = 0.90) -> LocusPanel:
    """Assemble a feline panel from GenBank/RefSeq protein records."""
    panel = LocusPanel(locus=locus, species="cat", curated=False,
                       source="NCBI Protein (Felis catus); IPD-MHC has no feline entries")
    include_re = re.compile("|".join(re.escape(x) for x in include), re.I)
    for header, seq in data.ncbi_feline_class2():
        if not include_re.search(header) or _EXCLUDE.search(header):
            continue
        panel.n_records += 1
        acc = header.split()[0]
        mol = groove.build_molecule(f"FLA-{locus}:{acc}", "cat", locus, seq,
                                    source=header)
        if _accept(mol, min_identity, min_coverage):
            panel.molecules.append(mol)
        else:
            panel.n_rejected += 1
    return panel


def deduplicate(panel: LocusPanel) -> LocusPanel:
    """Collapse molecules with an identical groove -- one point in model space."""
    seen: Dict[str, Molecule] = {}
    aliases: Dict[str, List[str]] = {}
    for mol in panel.molecules:
        if mol.pseudoseq in seen:
            aliases.setdefault(mol.pseudoseq, []).append(mol.name)
        else:
            seen[mol.pseudoseq] = mol
    for ps, names in aliases.items():
        seen[ps].name = f"{seen[ps].name} (+{len(names)} identical-groove)"
    out = LocusPanel(panel.locus, panel.species, list(seen.values()),
                     panel.n_records, panel.n_rejected, panel.curated, panel.source)
    return out


def build(species: str, spec: Sequence[dict]) -> Dict[str, LocusPanel]:
    panels: Dict[str, LocusPanel] = {}
    for item in spec:
        locus = item["locus"]
        if item["source"] == "ipd":
            panel = from_ipd(species, locus, item["prefixes"])
        elif item["source"] == "ncbi_feline":
            panel = from_ncbi_feline(locus, item["include"])
        else:
            raise ValueError(f"unknown panel source {item['source']!r}")
        if item.get("deduplicate", True):
            panel = deduplicate(panel)
        limit = item.get("max_molecules")
        if limit and len(panel.molecules) > limit:
            panel.molecules = _diverse_subset(panel.molecules, limit)
        panels[locus] = panel
    return panels


def _diverse_subset(molecules: Sequence[Molecule], k: int) -> List[Molecule]:
    """Max-min diversity pick over groove pseudosequences (greedy farthest-first).

    When a panel has to be capped, keeping the *most different* grooves samples
    the model's input space better than keeping the first k alphabetically.
    """
    chosen = [molecules[0]]
    remaining = list(molecules[1:])
    while len(chosen) < k and remaining:
        best, best_d = None, -1.0
        for mol in remaining:
            d = min(1.0 - groove.pseudo_identity(mol.pseudoseq, c.pseudoseq) for c in chosen)
            if d > best_d:
                best, best_d = mol, d
        chosen.append(best)
        remaining.remove(best)
    return chosen


# --------------------------------------------------------------------------
# NetMHCIIpan custom-molecule inputs
# --------------------------------------------------------------------------

NETMHCIIPAN_NOTE = """\
# NetMHCIIpan-4.3 custom-molecule inputs generated by vetimmuno.
#
# DR-type molecules need the beta chain only; DQ/DP-type molecules need both
# the alpha and the beta chain. Three caveats change how you read the output:
#
#   1. In custom-molecule (-mhcfsa) mode NetMHCIIpan does not emit %Rank -- the
#      rank reference distributions only exist for its built-in alleles. Use the
#      raw EL/BA score and convert it with vetimmuno's background-rank module
#      (`vetimmuno.predict.BackgroundRank`), which builds a per-molecule
#      reference distribution from random peptides.
#   2. DLA/FLA molecules are outside the model's training species. Check
#      applicability_domain.csv before reading any score as a risk estimate.
#   3. DQ alpha/beta pairs below are COMBINATORIAL. Real DLA class II
#      haplotypes are linked, so most of these pairs do not exist in any dog.
#      Supply real haplotypes in the species config to replace this.
#
# Verify the exact -mhcfsa chain-pairing convention against your local
# NetMHCIIpan-4.3 release notes before running the DQ commands.
"""


def write_netmhciipan_inputs(panels: Dict[str, LocusPanel], outdir: Path,
                             peptide_fasta: Path) -> Path:
    """Write per-molecule FASTA files plus a ready-to-run driver script."""
    outdir.mkdir(parents=True, exist_ok=True)
    mol_dir = outdir / "molecules"
    mol_dir.mkdir(exist_ok=True)

    lines = [NETMHCIIPAN_NOTE, "set -euo pipefail", 'NETMHCIIPAN="${NETMHCIIPAN:-netMHCIIpan}"',
             f'PEPTIDES="{peptide_fasta}"', 'OUT="$(dirname "$0")/raw"', 'mkdir -p "$OUT"', ""]

    dq_alpha = panels.get("DQA")
    for locus, panel in panels.items():
        if locus == "DQA":
            continue
        for mol in panel.molecules:
            safe = re.sub(r"[^A-Za-z0-9._-]", "_", mol.name.split(" ")[0])
            if locus == "DRB":
                path = mol_dir / f"{safe}.fasta"
                data.write_fasta(path, [(safe, mol.sequence)])
                lines.append(
                    f'"$NETMHCIIPAN" -f "$PEPTIDES" -inptype 0 -length 15 '
                    f'-mhcfsa "{path}" -xls -xlsfile "$OUT/{safe}.xls" > "$OUT/{safe}.txt"'
                )
            elif locus == "DQB" and dq_alpha is not None:
                for alpha in dq_alpha.molecules:
                    a_safe = re.sub(r"[^A-Za-z0-9._-]", "_", alpha.name.split(" ")[0])
                    pair = f"{a_safe}__{safe}"
                    path = mol_dir / f"{pair}.fasta"
                    data.write_fasta(path, [(f"{pair}_alpha", alpha.sequence),
                                            (f"{pair}_beta", mol.sequence)])
                    lines.append(
                        f'"$NETMHCIIPAN" -f "$PEPTIDES" -inptype 0 -length 15 '
                        f'-mhcfsa "{path}" -xls -xlsfile "$OUT/{pair}.xls" > "$OUT/{pair}.txt"'
                    )
    script = outdir / "run_netmhciipan.sh"
    script.write_text("\n".join(lines) + "\n")
    script.chmod(0o755)
    return script
