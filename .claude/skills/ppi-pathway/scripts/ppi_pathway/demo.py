"""
demo.py — end-to-end smoke/integration test for the ppi-pathway skill.

Exercises every data source the skill touches and prints a readable report with
PASS/FAIL checks, so you can confirm the whole pipeline works after install or
after an environment change.

    python -m ppi_pathway.demo                 # full run (downloads humanPPI ~14MB)
    python -m ppi_pathway.demo --online-only   # skip the humanPPI download
    python -m ppi_pathway.demo --with-string   # also download+use bulk STRING (~83MB)

Sources covered:
  1. STRING REST API  — PPI coherence + pathway enrichment (KEGG/Reactome/WikiP.)
  2. humanPPI (Cong)  — download structural predictions, build network, extract
                        the module around the demo genes, rank hubs
  3. InterPro         — annotate the top hub with domains/families/GO
  4. Reactome         — pathway enrichment (or graceful Cloudflare fallback)
  5. (optional) STRING bulk + BioGRID overlay via network.overlay()

Exit code is non-zero if any required check fails.
"""
from __future__ import annotations

import argparse
import sys

# A coherent DNA-damage / cell-cycle module — should light up cleanly.
DEMO_GENES = ["TP53", "BRCA1", "BRCA2", "EGFR", "MYC", "CDK2", "CDK4", "CDK6",
              "RB1", "MDM2", "ATM", "ATR", "CHEK1", "CHEK2", "CCNE1", "E2F1"]

_checks: list[tuple[str, bool, str]] = []


def _check(name: str, ok: bool, detail: str = ""):
    _checks.append((name, ok, detail))
    mark = "PASS" if ok else "FAIL"
    print(f"  [{mark}] {name}" + (f" — {detail}" if detail else ""))


def step_string_online():
    print("\n(1) STRING REST API — coherence + pathway enrichment")
    from . import string_api
    coh = string_api.ppi_enrichment(DEMO_GENES)
    p = float(coh.get("p_value", 1)) if coh else 1.0
    _check("STRING PPI-coherence p-value < 0.05", p < 0.05,
           f"p={p:.2e}, observed edges={coh.get('number_of_edges')}, "
           f"expected≈{coh.get('expected_number_of_edges')}")
    paths = string_api.top_pathways(DEMO_GENES, n=8)
    _check("STRING returns enriched pathways (KEGG/RCTM/WikiPathways)", len(paths) > 0,
           f"{len(paths)} pathways; top: {paths[0]['description'] if paths else '—'}")
    cats = {p["category"] for p in paths}
    print("      pathway categories seen:", ", ".join(sorted(cats)) or "—")
    for r in paths[:5]:
        print(f"        {r['category']:13s} fdr={r['fdr']:.1e}  {r['description']}")
    return paths


def step_humanppi(online_only: bool):
    print("\n(2) humanPPI (Cong lab, structural predictions) — download + network")
    if online_only:
        print("      skipped (--online-only)")
        return None, None
    from . import humanppi, subnetwork
    humanppi.download(precision=80)
    H = humanppi.load(precision=80)
    _check("humanPPI graph built", H.number_of_edges() > 10000,
           f"{H.number_of_nodes():,} nodes / {H.number_of_edges():,} edges")
    sub = subnetwork.induced(H, DEMO_GENES, expand=1)
    s = subnetwork.summary(sub)
    _check("module extracted around demo genes", s["seeds_found"] >= 5,
           f"{s['seeds_found']} seeds found, module {s['nodes']} nodes / {s['edges']} edges")
    hubs = subnetwork.rank_genes(sub)[:8]
    print("      top hubs (degree):")
    for r in hubs:
        tag = "seed" if r["is_seed"] else "connector"
        print(f"        {r['gene']:10s} deg={r['degree']:3d} btw={r['betweenness']:.3f} ({tag})")
    return H, sub


def step_interpro(H):
    print("\n(3) InterPro — annotate a hub with domains / families / GO")
    from . import interpro
    # find a UniProt accession for TP53 from the humanPPI edge metadata (no extra service)
    acc = None
    if H is not None:
        for u, v, d in H.edges("TP53", data=True):
            if u == "TP53" and d.get("up1"):
                acc = d["up1"]; break
            if v == "TP53" and d.get("up2"):
                acc = d["up2"]; break
    acc = acc or "P04637"  # TP53
    entries = interpro.entries_for_protein(acc)
    _check("InterPro returns entries for TP53", len(entries) > 0,
           f"{len(entries)} entries for {acc}")
    for e in entries[:4]:
        gos = ",".join(g["id"] for g in e["go_terms"][:3])
        print(f"        {e['accession']} {e['type']:22s} {e['name']}  [{gos}]")


def step_reactome():
    print("\n(4) Reactome — pathway enrichment (or graceful fallback)")
    from . import reactome
    try:
        rows = reactome.top_pathways(DEMO_GENES, n=5)
        _check("Reactome AnalysisService reachable", len(rows) > 0,
               f"{len(rows)} pathways; top: {rows[0]['name'] if rows else '—'}")
        for r in rows:
            print(f"        fdr={r['fdr']:.1e}  {r['found']}/{r['total']}  {r['name']}")
    except reactome.ReactomeBlocked as e:
        # Not a failure — expected in Cloudflare-blocked datacenter environments.
        _check("Reactome blocked → detected & fallback available", True,
               "Cloudflare challenge; STRING RCTM enrichment used instead")
        from . import string_api
        for r in string_api.top_pathways(DEMO_GENES, categories=("RCTM",), n=5):
            print(f"        (STRING/RCTM) fdr={r['fdr']:.1e}  {r['description']}")


def step_overlay(H, with_string: bool):
    if not with_string:
        return
    print("\n(5) Multi-source overlay — STRING(bulk) ∪ BioGRID ∪ humanPPI")
    from . import network, download
    download.string_files()
    S = network.load_string(min_score=700)
    try:
        download.biogrid_file()
        B = network.load_biogrid()
    except Exception as e:
        print("      BioGRID skipped:", e); B = None
    ov = network.overlay(string=S, biogrid=B, humanppi=H)
    multi = sum(1 for _, _, d in ov.edges(data=True) if d.get("support", 0) >= 2)
    _check("overlay built with multi-source consensus edges", multi > 0,
           f"{ov.number_of_edges():,} edges, {multi:,} supported by ≥2 sources")


def main(argv=None):
    ap = argparse.ArgumentParser(description="ppi-pathway end-to-end demo test")
    ap.add_argument("--online-only", action="store_true",
                    help="skip the humanPPI download (online checks only)")
    ap.add_argument("--with-string", action="store_true",
                    help="also download bulk STRING + BioGRID and test the overlay")
    a = ap.parse_args(argv)

    print("=" * 68)
    print("ppi-pathway skill — end-to-end demo test")
    print(f"demo gene set ({len(DEMO_GENES)}): {', '.join(DEMO_GENES)}")
    print("=" * 68)

    try:
        step_string_online()
    except Exception as e:
        _check("STRING online step", False, f"{type(e).__name__}: {e}")
    H, _sub = (None, None)
    try:
        H, _sub = step_humanppi(a.online_only)
    except Exception as e:
        _check("humanPPI step", False, f"{type(e).__name__}: {e}")
    try:
        step_interpro(H)
    except Exception as e:
        _check("InterPro step", False, f"{type(e).__name__}: {e}")
    try:
        step_reactome()
    except Exception as e:
        _check("Reactome step", False, f"{type(e).__name__}: {e}")
    try:
        step_overlay(H, a.with_string)
    except Exception as e:
        _check("overlay step", False, f"{type(e).__name__}: {e}")

    print("\n" + "=" * 68)
    passed = sum(1 for _, ok, _ in _checks if ok)
    total = len(_checks)
    print(f"RESULT: {passed}/{total} checks passed")
    print("=" * 68)
    failed = [n for n, ok, _ in _checks if not ok]
    if failed:
        print("Failed:", "; ".join(failed))
        return 1
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except BrokenPipeError:
        sys.exit(0)
