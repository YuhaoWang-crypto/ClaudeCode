"""
run_evo2_scoring.py — score pathway ortholog CDS with Evo2 on Modal (LIVE).

Reads the resolved CDS cache (psmos/data/cds/<pathway>_gates_cds.json), sends the
real coding sequences to the Evo2 scoring service on Modal, and writes the
per-(species,family) mean log-likelihoods to psmos/data/evo2/<pathway>_scores.json.

This is the step that turns the "fidelity / sequence-constraint" layer from a
divergence-time proxy into a model-measured quantity.

    pip install modal python-socks
    export SSL_CERT_FILE=/root/.ccr/ca-bundle.crt
    python3 -m psmos.run_evo2_scoring Notch

Cost: one H100 cold-start + a few seconds per sequence. Keep the CDS set small
(the gate families across ~8 species = ~15 sequences).
"""

import json
import os
import sys

import modal

from .pathways import get_pathway

DATA = os.path.join(os.path.dirname(__file__), "data")


def main(pathway_key="Notch", max_nt=None):
    os.environ.setdefault("SSL_CERT_FILE", "/root/.ccr/ca-bundle.crt")
    pw = get_pathway(pathway_key)

    cds_path = os.path.join(DATA, "cds", f"{pathway_key}_gates_cds.json")
    cds = json.load(open(cds_path))
    live = [c for c in cds if c["found"] and c["sequence"]]
    print(f"{len(live)} CDS to score with Evo2 ({pathway_key})")

    seqs, keys = [], []
    for c in live:
        s = c["sequence"]
        if max_nt:
            s = s[:max_nt]
        seqs.append(s)
        keys.append((c["species_key"], c["family"], c["nt_length"]))

    # Import the Modal app lazily so this module also imports without modal.
    from .evo2_modal import app, Evo2Scorer

    results = {}
    with modal.enable_output(), app.run():
        scorer = Evo2Scorer()
        print("model:", scorer.info.remote())
        # score in one batch (small set); Evo2.score_sequences handles the list
        lls = scorer.score.remote(seqs)
        for (sk, fam, ntlen), ll in zip(keys, lls):
            results.setdefault(sk, {})[fam] = round(float(ll), 6)
            print(f"  {sk:12s} {fam:12s} nt={ntlen:6d}  meanLL={ll:.5f}")

    os.makedirs(os.path.join(DATA, "evo2"), exist_ok=True)
    out = os.path.join(DATA, "evo2", f"{pathway_key}_scores.json")
    json.dump(results, open(out, "w"), indent=1)
    print(f"\nwrote {out}  ({sum(len(v) for v in results.values())} scores, "
          f"{len(results)} species)")


if __name__ == "__main__":
    pk = sys.argv[1] if len(sys.argv) > 1 else "Notch"
    main(pk)
