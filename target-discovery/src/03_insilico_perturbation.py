"""
Step 3 — In-silico perturbation (the hypothesis engine). GPU required.

For each gene, virtually DELETE it from pathological cells and measure how far
the cell's Geneformer embedding shifts toward the HEALTHY reference state.
Genes whose knockout most strongly moves 'aberrant_basaloid' /
'pathologic_fibroblast' cells toward the healthy state are the candidate targets.

This is counterfactual, not correlational: it asks "what does removing this gene
DO to the cell state", which is what makes it capable of surfacing *novel*
targets that differential-expression alone misses.

Uses geneformer.InSilicoPerturber with a start_state -> goal_state model, the
standard idiom for this analysis. See bionemo-recipes/recipes/geneformer for
obtaining/fine-tuning the checkpoint.
"""
from __future__ import annotations
import argparse
import os

import yaml


def load_config(path: str) -> dict:
    with open(path) as fh:
        return yaml.safe_load(fh)


def run_state(cfg: dict, start_state: str) -> str:
    """Run in-silico deletion for one pathological start state.

    Returns the path to the per-gene shift table for this state.
    """
    from geneformer import InSilicoPerturber, InSilicoPerturberStats

    p = cfg["perturbation"]
    healthy = cfg["healthy_reference_state"]
    out_dir = os.path.join(p["out_dir"], start_state)
    os.makedirs(out_dir, exist_ok=True)

    # Model the transition pathological(start) -> healthy(goal) in embedding space.
    cell_states_to_model = {
        "state_key": "cell_state",
        "start_state": start_state,
        "goal_state": healthy,
        "alt_states": [],
    }

    isp = InSilicoPerturber(
        perturb_type=p["perturb_type"],          # "delete"
        genes_to_perturb="all",                   # scan the whole expressed genome
        model_type=p["model_type"],
        emb_mode=p["emb_mode"],                    # "cell"
        cell_states_to_model=cell_states_to_model,
        max_ncells=p["max_cells_per_state"],
        forward_batch_size=p["batch_size"],
        nproc=p["num_workers"],
    )
    isp.perturb_data(
        model_directory=p["model_dir"],
        input_data_file=cfg["tokenize"]["out_dataset"],
        output_directory=out_dir,
        output_prefix=f"insilico_{start_state}",
    )

    # Aggregate into a per-gene "shift toward goal" statistic.
    stats = InSilicoPerturberStats(
        mode="goal_state_shift",
        cell_states_to_model=cell_states_to_model,
    )
    stats_out = os.path.join(out_dir, f"stats_{start_state}.csv")
    stats.get_stats(
        input_data_directory=out_dir,
        null_dist_data_directory=None,
        output_directory=out_dir,
        output_prefix=f"stats_{start_state}",
    )
    print(f"[{start_state}] shift table -> {stats_out}")
    return stats_out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.yaml")
    args = ap.parse_args()
    cfg = load_config(args.config)

    os.makedirs(cfg["perturbation"]["out_dir"], exist_ok=True)
    produced = []
    for state in cfg["target_cell_states"]:
        produced.append(run_state(cfg, state["name"]))

    manifest = os.path.join(cfg["perturbation"]["out_dir"], "manifest.txt")
    with open(manifest, "w") as fh:
        fh.write("\n".join(produced) + "\n")
    print(f"[done] {len(produced)} state(s) perturbed; manifest -> {manifest}")


if __name__ == "__main__":
    main()
