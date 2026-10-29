"""Gene lookup + rule-based interpretation (paper Table 1 and Discussion).

Works offline from the committed score table. If gnomAD v2.1.1 has been
fetched (python -m gismo fetch), LOEUF is joined in and the cross-metric
rules (recessive-like, X-linked LOEUF discordance) are applied too.
"""
import numpy as np
import pandas as pd

from gismo import paths
from gismo.core import load_scores

LOEUF_INTOLERANT = 0.35   # gnomAD v2.1.1 convention used by the paper (~top 15%)
LOEUF_TOLERANT = 0.6      # clearly tolerant of heterozygous LoF
SHORT_CDS = 1000          # bp; LOEUF has too few expected LoF sites below this


def _loeuf():
    p = paths.REF / "gnomad_v2.1.1_lof.txt.bgz"
    if not p.exists():
        return None
    v2 = pd.read_csv(p, sep="\t", compression="gzip",
                     usecols=["gene", "oe_lof_upper", "cds_length"])
    return (v2.dropna(subset=["oe_lof_upper"]).sort_values("oe_lof_upper")
              .drop_duplicates("gene").set_index("gene"))


def classify(row):
    """Return list of interpretation tags (⚠️ heuristic, decile cut-offs from Table 1)."""
    tags = []
    gd, md = row.get("gismo_decile"), row.get("gismo_mis_decile")
    lo = row.get("loeuf")
    if pd.notna(gd) and gd <= 2 or pd.notna(md) and md <= 2:
        tags.append("deep-mammalian constraint (deciles 1-2): Mendelian/essential candidate")
    cds = row.get("cds_length")
    short = pd.notna(cds) and cds < SHORT_CDS
    if short and pd.notna(lo) and lo >= LOEUF_INTOLERANT and pd.notna(gd) and gd <= 2:
        tags.append(f"short CDS ({int(cds)} bp): LOEUF underpowered here; "
                    "GISMO/GISMO-mis are the more reliable constraint readout")
    elif pd.notna(gd) and gd <= 4 and pd.notna(lo) and lo >= LOEUF_TOLERANT:
        tags.append("recessive-like candidate: retained across mammals but tolerates "
                    "heterozygous LoF in humans (homozygous selection)")
    if pd.notna(lo) and lo < LOEUF_INTOLERANT and pd.notna(gd) and gd >= 9:
        if str(row.get("chrom")) in ("X", "chrX"):
            tags.append("LOEUF-constrained but often lost in mammals, X-linked: "
                        "LOEUF may be inflated (hemizygosity / mapping)")
        else:
            tags.append("human-specific constraint or lineage-specific dispensability")
    hqd = row.get("gismo_hq_decile")
    if pd.notna(gd) and pd.notna(hqd) and gd - hqd >= 3:
        tags.append(f"assembly-sensitive: losses sit in fragmented genomes; GISMO decile "
                    f"{int(gd)} -> {int(hqd)} using high-quality assemblies only")
    gcp = row.get("gc_cds_pct")
    gc_risk = pd.notna(gd) and gd >= 8 and pd.notna(gcp) and gcp >= 0.8
    if gc_risk:
        tags.append(f"GC-rich CDS (GC {row.get('gc_cds', float('nan')):.0%}, top 20%): many 'losses' "
                    "may be assembly dropout - such d8-10 genes are mouse-lethal 3x as often "
                    "(16% vs 5.5%); do not read as dispensable without other evidence")
    if pd.notna(gd) and gd >= 9 and (pd.isna(md) or md >= 8) and not gc_risk:
        tags.append("frequently lost / fast-evolving across mammals (olfactory/immune-like)")
    if pd.notna(row.get("n_loss")) and row["n_loss"] == 0:
        tags.append("never lost in any of 462 mammals")
    return tags


def lookup(genes):
    d = load_scores()
    d["symbol_u"] = d.symbol.astype(str).str.upper()
    lo = _loeuf()
    out = []
    for g in genes:
        hit = d[(d.symbol_u == g.upper()) | (d.gene_id == g)]
        if hit.empty:
            out.append({"query": g, "found": False})
            continue
        r = hit.iloc[0].to_dict()
        if lo is not None and r["symbol"] in lo.index:
            r["loeuf"] = float(lo.loc[r["symbol"], "oe_lof_upper"])
            r["cds_length"] = float(lo.loc[r["symbol"], "cds_length"])
        r["query"], r["found"] = g, True
        r["tags"] = classify(r)
        out.append(r)
    return out


def fmt(rows):
    lines = []
    for r in rows:
        if not r["found"]:
            lines.append(f"{r['query']}: not in GISMO table")
            continue
        def f(x, p=3):
            return "NA" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:.{p}g}"
        lines.append(
            f"{r['symbol']} ({r['gene_id']}, chr{r['chrom']}): "
            f"GISMO={f(r['gismo'])} d{f(r['gismo_decile'])} (lost in {f(r['n_loss'])}/462) | "
            f"GISMO-mis={f(r['gismo_mis'])} d{f(r['gismo_mis_decile'])} | "
            f"LOEUF={f(r.get('loeuf'))}")
        for t in r["tags"]:
            lines.append(f"    - {t}")
    return "\n".join(lines)
