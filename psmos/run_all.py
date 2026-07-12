"""
run_all.py — end-to-end PSMOS pipeline for a pathway.

    python3 -m psmos.run_all Notch            # full pipeline incl. live Evo2
    python3 -m psmos.run_all Notch --no-evo2  # skip the Modal/Evo2 step

Steps: UniProt ortholog search (gate) → Ensembl CDS → Evo2 scoring on Modal →
scoring engine → dashboard. Each step caches under psmos/data/, so re-runs are
cheap and the dashboard rebuilds offline.
"""

import os
import sys

os.environ.setdefault("SSL_CERT_FILE", "/root/.ccr/ca-bundle.crt")


def run(pathway_key="Notch", do_evo2=True):
    from .pathways import get_pathway
    from .orthologs import fetch_pathway_gates
    from .cds import fetch_pathway_cds

    pw = get_pathway(pathway_key)
    print(f"\n=== 1. UniProt ortholog search (gate) — {pw.label} ===")
    orths = fetch_pathway_gates(pw, use_cache=True)

    print(f"\n=== 2. Ensembl CDS ===")
    fetch_pathway_cds(pw, orths, use_cache=True)

    print(f"\n=== 2b. Ensembl Compara paralogue counts (redundancy) ===")
    try:
        from .paralogs import fetch_pathway_paralogs
        fetch_pathway_paralogs(pw, use_cache=True)
    except Exception as e:
        print(f"  [skipped] Compara step: {e}")

    print(f"\n=== 2c. AlphaGenome R-layer (human/mouse) ===")
    try:
        from .alphagenome_r import compute_R
        status, _ = compute_R(pw)
        print(f"  AlphaGenome: {status}")
    except Exception as e:
        print(f"  [skipped] AlphaGenome step: {e}")

    if do_evo2:
        print(f"\n=== 3. Evo2-7B scoring on Modal H100 (live) ===")
        try:
            from .run_evo2_scoring import main as evo2_main
            evo2_main(pathway_key)
        except Exception as e:
            print(f"  [skipped] Evo2 step failed: {e}")
            print("  (dashboard will render curated + gate; rerun the Evo2 step later)")
    else:
        print("\n=== 3. Evo2 scoring — SKIPPED (--no-evo2) ===")

    print(f"\n=== 4-5. Scoring + dashboard ===")
    if pathway_key == "Hippo":
        from .build_hippo_dashboard import build as build_hippo
        build_hippo()
    else:
        from .scoring import compute, evo2_vs_curated_correlation
        _, scores = compute(pathway_key)
        live = [s for s in scores if s.evo2]
        r, n = evo2_vs_curated_correlation(scores)
        print(f"  gate pass: {sum(s.gate_ok for s in scores)}/{len(scores)}  |  Evo2 live: {len(live)}")
        if r is not None:
            print(f"  Evo2 constraint vs curated fidelity: r={r:+.2f} (n={n})")
        from .build_dashboard import build
        build(pathway_key)
    print("\ndone.")


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    pk = args[0] if args else "Notch"
    run(pk, do_evo2="--no-evo2" not in sys.argv)
