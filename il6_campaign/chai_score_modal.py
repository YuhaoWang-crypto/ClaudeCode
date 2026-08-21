"""Score designed binder-target complexes with Chai-1 on a Modal GPU.

Chai-1 did not participate in the design, so it is an orthogonal judge: it re-folds the
delivered sequence pair from scratch and returns its own PAE, from which the same
ipSAE_min used for Boltz is computed.

    modal run chai_score_modal.py --input shortlist.json --out chai_scores.json

`shortlist.json` is a list of {"design_id", "binder", "target"} records.
Structures and PAE come back inline so the caller can run ipsae.py on them locally.
"""

import json
import os

import modal

CHAI = "chai_lab==0.6.1"

image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("git")
    .pip_install("torch==2.5.1", extra_index_url="https://download.pytorch.org/whl/cu121")
    .pip_install(CHAI, "numpy<2")
    .env({"CHAI_DOWNLOADS_DIR": "/weights"})
)

app = modal.App("chai1-binder-scoring")
weights = modal.Volume.from_name("chai1-weights", create_if_missing=True)


@app.function(image=image, gpu="A10G", volumes={"/weights": weights},
              timeout=60 * 60, max_containers=6)
def score(job: dict) -> dict:
    """Fold one binder-target pair with Chai-1 and return scores + PAE + structure."""
    import tempfile
    from pathlib import Path

    import numpy as np
    import torch
    from chai_lab.chai1 import run_inference

    tmp = Path(tempfile.mkdtemp())
    fasta = tmp / "in.fasta"
    fasta.write_text(
        f">protein|target\n{job['target']}\n>protein|binder\n{job['binder']}\n"
    )
    out = tmp / "out"

    cands = run_inference(
        fasta_file=fasta,
        output_dir=out,
        num_trunk_recycles=3,
        num_diffn_timesteps=200,
        seed=job.get("seed", 0),
        device="cuda:0",
        use_esm_embeddings=True,
    )
    cands = cands.sorted()

    rank = cands.ranking_data[0]
    pae = cands.pae[0].numpy().astype(np.float32)
    cif = Path(cands.cif_paths[0]).read_text()

    def scalar(x):
        return float(x.item()) if hasattr(x, "item") else float(x)

    weights.commit()
    return {
        "design_id": job["design_id"],
        "seed": job.get("seed", 0),
        "aggregate_score": scalar(rank.aggregate_score),
        "ptm": scalar(rank.ptm_scores.complex_ptm),
        "iptm": scalar(rank.ptm_scores.interface_ptm),
        "plddt": scalar(rank.plddt_scores.complex_plddt),
        "has_inter_chain_clashes": bool(
            torch.any(rank.clash_scores.has_inter_chain_clashes).item()),
        "pae": pae.tolist(),
        "cif": cif,
    }


@app.local_entrypoint()
def main(input: str, out: str, seeds: int = 1):
    jobs = json.loads(open(input).read())
    todo = [dict(j, seed=s) for j in jobs for s in range(seeds)]
    print(f"scoring {len(todo)} complexes with Chai-1 on Modal")
    results = []
    for r in score.map(todo, order_outputs=False, return_exceptions=True):
        if isinstance(r, Exception):
            print("  FAILED:", repr(r)[:200])
            continue
        results.append(r)
        print(f"  {r['design_id']} seed{r['seed']}: iptm={r['iptm']:.3f} "
              f"plddt={r['plddt']:.3f} agg={r['aggregate_score']:.3f}")
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    with open(out, "w") as fh:
        json.dump(results, fh)
    print(f"wrote {out} ({len(results)}/{len(todo)} succeeded)")
