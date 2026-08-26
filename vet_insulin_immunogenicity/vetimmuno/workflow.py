"""End-to-end per-species workflow driver."""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import yaml

from . import data, epitope, groove, panel, predict, report, validate
from .groove import DomainCall, Molecule
from .insulin import Insulin, ProductSpec, build_product, natural_insulin

PKG_ROOT = Path(__file__).resolve().parent.parent
CONFIG_DIR = PKG_ROOT / "config"


def load_config(species_or_path: str) -> dict:
    path = Path(species_or_path)
    if not path.exists():
        path = CONFIG_DIR / f"{species_or_path}.yaml"
    return yaml.safe_load(path.read_text())


def _products(cfg: dict) -> List[Insulin]:
    out: List[Insulin] = []
    for item in cfg["products"]:
        out.append(build_product(ProductSpec(
            name=item["name"], parent=item["parent"],
            edits=tuple(item.get("edits", ())), note=item.get("note", ""),
        )))
    return out


def _domain_summary(species: str, per_locus: Dict[str, dict]) -> dict:
    """Per-locus applicability rows plus one overall verdict.

    Kept per locus on purpose: each locus is compared against its own training
    pool, so pooling the numbers across loci would compare distances measured on
    different scales.
    """
    rows = []
    calls: List[DomainCall] = []
    for locus, blob in sorted(per_locus.items()):
        lc = blob["calls"]
        calls.extend(lc)
        ident = [c.identity for c in lc]
        rows.append({
            "locus": locus,
            "training_molecules": len(blob["pool"]),
            "train_median": float(np.median(blob["loo"])),
            "train_p5": float(np.percentile(blob["loo"], 5)),
            "panel_median": float(np.median(ident)) if ident else float("nan"),
            "panel_min": float(min(ident)) if ident else float("nan"),
            "panel_max": float(max(ident)) if ident else float("nan"),
            "n_out": sum(c.verdict.startswith("OUT-OF-DOMAIN") for c in lc),
            "n_marginal": sum(c.verdict.startswith("marginal") for c in lc),
            "n_total": len(lc),
        })
    ident = [c.identity for c in calls]
    n_out = sum(c.verdict.startswith("OUT-OF-DOMAIN") for c in calls)
    n_marg = sum(c.verdict.startswith("marginal") for c in calls)
    n_near = len(calls) - n_out - n_marg
    n_training = sum(r["training_molecules"] for r in rows)
    train_p5 = float(np.mean([r["train_p5"] for r in rows]))
    med = float(np.median(ident))
    if n_out == len(calls):
        verdict = (
            f"**Every {species} molecule in this panel is outside the model's training "
            f"space.** Its nearest training neighbour is further away than the nearest "
            f"neighbour of essentially any molecule the model was actually trained on. "
            f"Binding scores for this species are extrapolation, and there is no public "
            f"{species} benchmark that would tell you how far off they are. Report them "
            f"as research-use hypotheses, never as a risk quantification.")
    elif med < train_p5:
        verdict = (
            f"**The {species} panel sits below the training set's own 5th percentile "
            f"({med:.2f} vs {train_p5:.2f}).** {n_out}/{len(calls)} molecules are frank "
            f"extrapolation and {n_near}/{len(calls)} are close enough to training "
            f"molecules that predictions are defensible as directional. Stratify the "
            f"output by this column rather than treating the panel as uniform.")
    else:
        verdict = (
            f"**The {species} panel largely falls inside the training set's own "
            f"nearest-neighbour range** (median {med:.2f} vs 5th percentile "
            f"{train_p5:.2f}). Predictions are still cross-species extrapolation with no "
            f"species-specific benchmark, but they are not geometric extrapolation.")
    return {
        "rows": rows, "n_training": n_training, "train_p5": train_p5,
        "panel_median": med, "panel_min": float(min(ident)),
        "panel_max": float(max(ident)), "n_out": n_out, "n_marginal": n_marg,
        "n_near": n_near, "n_total": len(calls), "verdict_text": verdict,
    }


def run(species: str, outdir: Optional[Path] = None, backend_name: str = "surrogate",
        n_background: int = 20000, offline: bool = False,
        score_top_n: int = 15) -> dict:
    cfg = load_config(species)
    species = cfg["species"]
    data.set_offline(offline)

    outdir = Path(outdir or PKG_ROOT / "results" / species)
    figdir = outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    # -- 1. sequences -------------------------------------------------------
    self_ins = natural_insulin(cfg["self_insulin"])
    products = _products(cfg)
    burdens = [epitope.foreign_burden(p, self_ins) for p in products]
    report.write_csv(outdir / "foreign_burden.csv", [
        {k: (", ".join(v) if isinstance(v, list) else v) for k, v in b.items()}
        for b in burdens])

    # -- 2. panel -----------------------------------------------------------
    panels = panel.build(species, cfg["panel"])
    all_molecules: List[Molecule] = [m for p in panels.values() for m in p.molecules]
    report.write_csv(outdir / "panel_summary.csv", [p.summary() for p in panels.values()])
    report.write_csv(outdir / "panel_molecules.csv", [
        {"locus": m.locus, "chain": m.chain, "molecule": m.name,
         "groove_pseudosequence": m.pseudoseq, "contact_coverage": round(m.coverage, 3),
         "identity_to_human_reference": round(m.identity_to_ref, 3), "source": m.source}
        for m in all_molecules])

    # -- 3. applicability domain -------------------------------------------
    calls: List[DomainCall] = []
    loo_all: List[float] = []
    per_locus: Dict[str, dict] = {}
    for locus, p in panels.items():
        pool = groove.training_space(locus)
        loo = groove.loo_identities(pool)
        locus_calls = groove.applicability(p.molecules, pool, loo)
        per_locus[locus] = {"pool": pool, "loo": loo, "calls": locus_calls}
        calls.extend(locus_calls)
        if locus == "DRB" or not loo_all:
            loo_all = loo
    report.write_csv(outdir / "applicability_domain.csv", [
        {"molecule": c.molecule, "locus": c.locus, "nearest_training_molecule": c.nearest,
         "groove_identity": round(c.identity, 3), "groove_blosum": round(c.blosum, 3),
         "percentile_vs_training_LOO": round(c.percentile_vs_training, 2),
         "verdict": c.verdict}
        for c in calls])
    domain_summary = _domain_summary(species, per_locus)

    # -- 4. binding predictions --------------------------------------------
    if backend_name == "netmhciipan":
        backend = predict.NetMHCIIpanBackend()
        if not backend.available:
            raise SystemExit(
                "NetMHCIIpan not found. Inputs were still generated under "
                f"{outdir / 'netmhciipan'}; run them on a licensed host."
            )
    else:
        backend = predict.IllustrativeScorer()
    ranker = predict.BackgroundRank(backend, n_background=n_background)

    scoring_molecules = [m for m in panels[cfg.get("scoring_locus", "DRB")].molecules]
    hits: List[predict.Hit] = []
    for product in products:
        peptides = epitope.neo_peptides(product, self_ins)
        product_hits = predict.score_peptides(peptides, scoring_molecules, backend, ranker)
        for h in product_hits:
            hits.append(h)
        report.write_csv(outdir / f"hits_{_slug(product.name)}.csv", [
            {"product": product.name, **{k: v for k, v in asdict(h).items()}}
            for h in sorted(product_hits, key=lambda x: x.rank)[:2000]])

    # Non-self cores only: the rows a reviewer would actually look at.
    by_core: Dict[str, List[predict.Hit]] = {}
    for h in hits:
        by_core.setdefault(h.core, []).append(h)
    neo_rows = []
    for product in products:
        for nc in epitope.neo_cores(product, self_ins):
            best = by_core.get(nc.core.sequence, [])
            best_rank = min((h.rank for h in best), default=float("nan"))
            best_mol = min(best, key=lambda h: h.rank).molecule if best else ""
            neo_rows.append({
                "product": product.name, "chain": nc.core.chain,
                "core": nc.core.sequence, "position": nc.core.label,
                "foreign_residues": ", ".join(nc.foreign_residues) or "—",
                "best_rank_any_molecule": round(best_rank, 3) if best else "",
                "best_molecule": best_mol,
                "class": predict.classify(best_rank) if best else "",
                "score_provenance": backend.name,
            })
    report.write_csv(outdir / "non_self_cores.csv", neo_rows)

    # One row per distinct non-self core, listing which products carry it --
    # the same core recurs across products and repeating it adds no information.
    grouped: Dict[tuple, dict] = {}
    for row in neo_rows:
        if row["best_rank_any_molecule"] == "":
            continue
        key = (row["chain"], row["core"], row["position"])
        entry = grouped.setdefault(key, {**row, "products": []})
        entry["products"].append(row["product"])
    top_hits = []
    for row in sorted(grouped.values(), key=lambda r: r["best_rank_any_molecule"])[:score_top_n]:
        top_hits.append([
            f"{len(row['products'])} product(s): " + "; ".join(row["products"]),
            row["chain"], row["core"], row["position"], row["foreign_residues"],
            row["best_molecule"], row["best_rank_any_molecule"], row["class"]])

    # -- 5. NetMHCIIpan inputs ---------------------------------------------
    pep_fasta = outdir / "netmhciipan" / "peptides.fasta"
    records = []
    for product in products:
        for ch in ("A", "B"):
            records.append((f"{_slug(product.name)}|{ch}", product.chain(ch)))
    data.write_fasta(pep_fasta, records)
    panel.write_netmhciipan_inputs(panels, outdir / "netmhciipan", pep_fasta)

    # -- 6. validation ------------------------------------------------------
    reference_mol = _human_reference_molecule()
    vreport = validate.run_all(
        species=species, products=products, molecules=all_molecules,
        backend=backend, ranker=ranker, domain_calls=calls, loo=loo_all,
        identity_donor=cfg.get("identity_control_donor"),
        reference_molecule=reference_mol,
    )
    report.write_csv(outdir / "validation.csv", vreport.to_rows())

    # -- 7. figures + report -----------------------------------------------
    figs = []
    f1 = report.figure_difference_map(self_ins, products, figdir / "difference_map.png")
    figs.append(("Residues of each product that are foreign to this species", f"figures/{f1.name}"))
    f2 = report.figure_applicability(loo_all, [c for c in calls if c.locus == "DRB"],
                                     species, figdir / "applicability_domain.png")
    figs.append(("Distance from the predictor's training space", f"figures/{f2.name}"))
    worst = max(burdens, key=lambda b: b["n_neo_cores"])
    worst_product = next(p for p in products if p.name == worst["product"])
    cov = epitope.residue_coverage(worst_product, epitope.neo_cores(worst_product, self_ins))
    f3 = report.figure_core_landscape(worst_product, cov, figdir / "core_landscape.png")
    figs.append((f"Non-self core landscape for {worst_product.name}", f"figures/{f3.name}"))

    ctx = {
        "species": species, "common_name": cfg["common_name"], "self_insulin": self_ins,
        "burdens": burdens, "panel_summaries": [p.summary() for p in panels.values()],
        "domain_summary": domain_summary, "backend": backend.name,
        "backend_trained": getattr(backend, "trained", True), "top_hits": top_hits,
        "validation": vreport, "figures": figs,
        "interpretation": cfg["interpretation"].strip(),
    }
    (outdir / "report.md").write_text(report.build_report(ctx))
    (outdir / "summary.json").write_text(json.dumps({
        "species": species, "backend": backend.name,
        "burdens": burdens, "domain_summary": domain_summary,
        "validation": {"passed": vreport.passed, "total": len(vreport.checks),
                       "failed": [c.id for c in vreport.failed]},
    }, indent=2))
    return {"outdir": outdir, "validation": vreport, "domain": domain_summary,
            "burdens": burdens}


def _slug(name: str) -> str:
    return "".join(c if c.isalnum() else "_" for c in name).strip("_").lower()


def _human_reference_molecule() -> Optional[Molecule]:
    """HLA-DRB1*01:01 -- the molecule the surrogate's direction check uses."""
    try:
        table = data.imgt_hla("DRB")
        seq = table.get("DRB1*01:01:01:01")
        if not seq:
            return None
        return groove.build_molecule("HLA-DRB1*01:01", "human", "DRB", seq)
    except Exception:
        return None
