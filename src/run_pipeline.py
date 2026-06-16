"""
run_pipeline.py
---------------
End-to-end driver:

    python src/run_pipeline.py            # all targets in targets.yaml
    python src/run_pipeline.py AdK PTP1B  # selected targets

For each target it resolves the conformer pair (offline ProDy data or RCSB
mirrors), aligns them, ranks distal candidate regions + hinge residues, designs a
His-pair/His-triplet library for the top regions, and writes a Markdown report to
results/<TARGET>.md. Targets whose structures can't be fetched (egress blocked)
are skipped with a clear note so the rest still run.
"""
from __future__ import annotations
import os
import sys
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import fetch_structures as fetch          # noqa: E402
import analyze_conformers as ana          # noqa: E402
import design_his_sites as design         # noqa: E402

RESULTS = os.path.join(ROOT, "results")


def load_targets():
    with open(os.path.join(ROOT, "targets.yaml")) as fh:
        doc = yaml.safe_load(fh)
    return doc["defaults"], doc["targets"]


def run_target(key, spec, defaults):
    cfg = dict(defaults)
    pair = spec.get("pair") or {}
    ref, alt = pair.get("ref"), pair.get("alt")
    if not ref or not alt:
        return {"key": key, "status": "skipped", "reason": "no conformer pair (single structure / sequence-only target)"}

    offline = spec.get("offline_files") or {}
    try:
        ref_path = fetch.resolve_structure(ref["pdb"], offline.get("ref"))
        alt_path = fetch.resolve_structure(alt["pdb"], offline.get("alt"))
    except Exception as e:                 # network blocked etc.
        return {"key": key, "status": "skipped", "reason": str(e)}

    result = ana.analyze(ref_path, ref["chain"], alt_path, alt["chain"],
                         spec.get("active_site", []), cfg)
    regions = result["regions"]
    top_regions = regions[:max(1, min(20, len(regions)))]
    library = design.design_library(result, regions, cfg, focus_regions=top_regions)
    recovery = ana.recovery_stats(result, spec.get("known_allosteric", []))

    write_report(key, spec, ref, alt, result, top_regions, library, recovery)
    return {"key": key, "status": "ok", "n_regions": len(regions),
            "n_pairs": len(library["pairs"]), "n_triplets": len(library["triplets"]),
            "recovery": recovery}


def write_report(key, spec, ref, alt, result, regions, library, recovery):
    os.makedirs(RESULTS, exist_ok=True)
    L = []
    a = L.append
    a(f"# {key} — {spec['name']}\n")
    a(f"**Conformer pair:** `{ref['pdb']}:{ref['chain']}` ({ref['state']})  ⇄  "
      f"`{alt['pdb']}:{alt['chain']}` ({alt['state']})\n")
    a(f"- matched residues: **{result['n_matched']}**")
    a(f"- rigid-core residues: **{int(result['core_mask'].sum())}**  "
      f"(core RMSD {result['rmsd_core']:.2f} Å, all-CA RMSD {result['rmsd_all']:.2f} Å)")
    a(f"- max Cα displacement after core fit: **{result['max_disp']:.1f} Å**")
    a(f"- active-center = centroid of {len(spec.get('active_site', []))} active-site CA; "
      f"distal cutoff {result['distal_cutoff']:.0f} Å; "
      f"responsive ≥ {result['disp_threshold']:.2f} Å (P{int(result['displacement_pctl'])})\n")
    if spec.get("expected_signal"):
        a(f"> _Expected signal:_ {spec['expected_signal']}\n")

    a("## Distal candidate regions (ranked)\n")
    a("| # | residues | len | mean Δ (Å) | max Δ (Å) | min dist→active (Å) |")
    a("|---|----------|-----|------------|-----------|---------------------|")
    for i, r in enumerate(regions, 1):
        a(f"| {i} | {r['start']}–{r['end']} | {r['length']} | {r['mean_disp']:.1f} | "
          f"{r['max_disp']:.1f} | {r['min_dist_active']:.1f} |")
    a("")

    # top hinge residues
    hs = result["hinge_score"]
    order = hs.argsort()[::-1]
    a("## Top hinge residues (clamp points)\n")
    a("Distal residues with the steepest local displacement gradient — where a "
      "moving domain meets the rigid core. Best single-metal clamp points.\n")
    a("| residue | resname | hinge score | Δ (Å) | dist→active (Å) |")
    a("|---------|---------|-------------|-------|-----------------|")
    shown = 0
    for k in order:
        if hs[k] <= 0:
            break
        a(f"| {int(result['resnums'][k])} | {result['resnames'][k]} | {hs[k]:.2f} | "
          f"{result['displacement'][k]:.1f} | {result['r_to_active'][k]:.1f} |")
        shown += 1
        if shown >= 12:
            break
    a("")

    a("## Suggested His-pair / His-triplet library\n")
    a("_Backbone-geometry heuristics — validate with rotamer-level metal modelling "
      "before synthesis. Screen each under Zn²⁺/Ni²⁺/Co²⁺ for kcat, Km, kcat/Km._\n")
    a("**His-pairs (bidentate clamp):**\n")
    a("| mutations | Cα–Cα (Å) | score |")
    a("|-----------|-----------|-------|")
    for p in library["pairs"][:15]:
        muts = " / ".join(f"{x}H" for x in p["positions"])
        a(f"| {muts} | {p['ca_dist']} | {p['score']} |")
    a("")
    if library["triplets"]:
        a("**His-triplets (facial site):**\n")
        a("| mutations | Cα–Cα pairwise (Å) | score |")
        a("|-----------|--------------------|-------|")
        for t in library["triplets"][:10]:
            muts = " / ".join(f"{x}H" for x in t["positions"])
            a(f"| {muts} | {t['ca_dists']} | {t['score']} |")
        a("")

    if recovery is not None:
        a("## Ground-truth recovery (known allosteric site)\n")
        a(f"- known allosteric residues: {recovery['n_known']}")
        a(f"- recovered inside top regions: **{recovery['n_recovered']}** "
          f"({100*recovery['fraction']:.0f}%) → {recovery['recovered']}\n")

    with open(os.path.join(RESULTS, f"{key}.md"), "w") as fh:
        fh.write("\n".join(L))


def main(argv):
    defaults, targets = load_targets()
    keys = argv[1:] if len(argv) > 1 else list(targets)
    summary = []
    for k in keys:
        if k not in targets:
            print(f"!! unknown target {k}")
            continue
        res = run_target(k, targets[k], defaults)
        summary.append(res)
        if res["status"] == "ok":
            rec = res.get("recovery")
            rstr = ""
            if rec:
                rstr = f", recovered {rec['n_recovered']}/{rec['n_known']} known-allosteric"
            print(f"[ok]   {k}: {res['n_regions']} regions, "
                  f"{res['n_pairs']} His-pairs, {res['n_triplets']} His-triplets{rstr}  -> results/{k}.md")
        else:
            print(f"[skip] {k}: {res['reason']}")
    write_summary(summary, targets)
    return summary


def write_summary(summary, targets):
    os.makedirs(RESULTS, exist_ok=True)
    L = ["# Run summary\n",
         "| target | status | distal regions | His-pairs | His-triplets | known-allosteric recovered |",
         "|--------|--------|----------------|-----------|--------------|-----------------------------|"]
    for r in summary:
        k = r["key"]
        if r["status"] == "ok":
            rec = r.get("recovery")
            rstr = f"{rec['n_recovered']}/{rec['n_known']}" if rec else "— (design target)"
            L.append(f"| {k} | ok | {r['n_regions']} | {r['n_pairs']} | {r['n_triplets']} | {rstr} |")
        else:
            reason = "egress blocked" if "egress" in r["reason"] else r["reason"][:40]
            L.append(f"| {k} | skipped | — | — | — | {reason} |")
    L.append("\n_Skipped targets need a structural-data host in the network egress "
             "allowlist (e.g. `files.rcsb.org`), or their PDBs dropped into `data/raw/`._")
    with open(os.path.join(RESULTS, "SUMMARY.md"), "w") as fh:
        fh.write("\n".join(L))


if __name__ == "__main__":
    main(sys.argv)
