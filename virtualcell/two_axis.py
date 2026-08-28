"""The two transferability axes, measured apart, with Arc's protein embeddings.

Everything this project has reported so far pools two questions that behave
differently:

* **cross-context** -- the knockdown was measured in a source line, but *this*
  cell line is new.  Every benchmark here so far sits on this axis.
* **unseen target** -- the knockdown was never perturbed in any source line, so
  there is no measured effect to transfer at all.  22 of the official 300 are in
  this position, and so are 10,852 of the 18,533 official genes.

Two independent sources say to separate them.  ``transferability-prior-eval``
measured that they respond oppositely -- representations help the unseen-target
axis, while even an *oracle* context feature scored worse than none on the
cross-context axis -- and CellFluxV2 (bioRxiv 2026.01.19.696785) crosses its
plate and perturbation splits for the same reason.  A single pooled number hides
which axis a change acted on.

The intervention tested here is Arc's **SE-600M protein embeddings**
(``protein_embeddings.pt``, ESM-derived, 5,120-d, 19,790 symbols).  They are a
poor candidate for the cross-context axis -- a protein embedding is identical in
every cell line, so it cannot express a context term, the same structural
objection that rules out a DNA language model.  They are a *good* candidate for
the unseen-target axis, and for one reason specific to this model: the existing
route for an unseen knockdown needs the silenced gene to be **measured**, so that
its response profile across the rest of the panel can identify functional
neighbours.  A protein embedding needs only a sequence.  It therefore reaches
genes the current route cannot reach at all, which is where the headroom is.

    python -m virtualcell.two_axis
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from . import metrics as M
from .benchmark import de_truth, evaluate
from .data import load_all, shared_perturbations
from .model import (ContextTransferModel, GlobalMeanBaseline, Hyper, SourceBank)

EMBED = Path("/home/user/vcc_data/state/protein_embeddings.pt")
# An .npz beside it wins: reading it needs no torch, which is the
# difference between working and not working on a lean remote image.
EMBED_NPZ = EMBED.with_suffix(".npz")
MIXES = (0.0, 0.25, 0.5, 0.75, 1.0)


class BlindToGene(ContextTransferModel):
    """Same model, but the silenced gene's own measured column is hidden.

    This is the third axis, and it is the one that matters for the official
    panel.  Withholding a knockdown from the source still leaves the *gene*
    measured, so the response-profile route can find its functional neighbours
    and a protein embedding has nothing to add -- which is exactly what the
    first run of this benchmark measured, and why that run could not have
    detected the effect it was built to look for.

    26 of the official 300 are neither perturbed nor measured in the source.
    For those the response-profile route has no input at all and the model falls
    back to the generic response.  Forcing ``gene=None`` reproduces that
    situation on data where the truth is known, so the real question becomes
    answerable: **is embedding-routed prediction better than giving up?**
    """

    def _from_neighbours(self, gene: int | None,
                         name: str | None = None) -> np.ndarray:
        return super()._from_neighbours(None, name)


def load_embeddings(symbols: np.ndarray, required: bool = False
                    ) -> np.ndarray | None:
    """Unit-norm protein embeddings on this gene axis; zero rows where absent.

    ``required`` turns a missing file into an error.  Returning None quietly is
    right for a benchmark arm that can be skipped, and wrong for a build that
    asked for embedding routing: on a fresh machine the file is simply absent
    and the submission silently loses the only route it has to the 26 panel
    targets with no measurement at all.
    """
    if EMBED_NPZ.exists():
        z = np.load(EMBED_NPZ, allow_pickle=True)
        blob = dict(zip(z["symbols"].astype(str), z["embeddings"]))
    elif EMBED.exists():
        import torch
        blob = torch.load(EMBED, map_location="cpu", weights_only=False)
    elif required:
        raise SystemExit(
            f"esm_mix > 0 but neither {EMBED_NPZ} nor {EMBED} exists. Fetch "
            f"with:\n  curl -sSL -o {EMBED} https://huggingface.co/arcinstitute/"
            f"SE-600M/resolve/main/protein_embeddings.pt")
    else:
        return None
    dim = len(next(iter(blob.values())))
    out = np.zeros((symbols.size, dim), dtype=np.float32)
    hit = 0
    for i, s in enumerate(symbols):
        v = blob.get(str(s))
        if v is None:
            continue
        e = np.asarray(v, dtype=np.float32)
        n = np.linalg.norm(e)
        if n > 0:
            out[i] = e / n
            hit += 1
    print(f"  protein embeddings: {hit:,}/{symbols.size:,} genes covered "
          f"({hit / symbols.size:.1%}), {dim}-d")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-eval", type=int, default=250)
    ap.add_argument("--n-held-out-perts", type=int, default=250)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, default=Path("results/two_axis.json"))
    args = ap.parse_args()

    lines, _, symbols = load_all()
    common = np.array(shared_perturbations(lines))
    embed = load_embeddings(symbols)
    rng = np.random.default_rng(args.seed)

    # Split the shared panel once, for every fold: knockdowns in ``withheld``
    # are removed from every source bank, so predicting them exercises the
    # unseen-target route while ``kept`` exercises the cross-context one.  Both
    # are scored in the same held-out line, so the only difference between the
    # two axes is whether the source ever saw the perturbation.
    order = rng.permutation(common.size)
    withheld = np.sort(common[order[:args.n_held_out_perts]])
    kept = np.sort(common[order[args.n_held_out_perts:]])
    eval_kept = np.sort(rng.choice(kept, min(args.n_eval, kept.size),
                                   replace=False))
    eval_unseen = np.sort(rng.choice(withheld, min(args.n_eval, withheld.size),
                                     replace=False))
    print(f"{common.size} shared knockdowns -> {kept.size} kept in the source, "
          f"{withheld.size} withheld from it")

    results: dict[str, dict] = {}
    for target in lines:
        sources = [cl for cl in lines if cl.name != target.name]
        # Strip the withheld knockdowns out of every source line.
        from dataclasses import replace as dc_replace
        stripped = []
        for s in sources:
            keep = np.array([n not in set(withheld) for n in s.names])
            stripped.append(dc_replace(s, names=s.names[keep],
                                       pert=s.pert[keep], ncells=s.ncells[keep],
                                       _delta=None))
        bank = SourceBank.build(stripped, np.ones(len(stripped)))

        fold: dict[str, dict] = {}
        for axis, perts in (("cross-context", eval_kept),
                            ("unseen target", eval_unseen),
                            ("unseen target, gene unmeasured", eval_unseen)):
            blind = axis.endswith("unmeasured")
            truth = de_truth(target, perts)
            base = evaluate(GlobalMeanBaseline().fit(stripped, target.mu, symbols)
                            .predict(perts), target, perts, truth, symbols)
            row = {"global mean [baseline]": base}
            for mix in MIXES:
                if mix > 0 and embed is None:
                    continue
                hp = Hyper(temp=0.1, esm_mix=mix)
                cls = BlindToGene if blind else ContextTransferModel
                m = cls(hp).fit(stripped, target.mu, symbols,
                                bank=bank, gene_embed=embed)
                row[f"ContextTransfer esm_mix={mix:g}"] = evaluate(
                    m.predict(perts), target, perts, truth, symbols)
            fold[axis] = row
        results[target.name] = fold

        print(f"\n--- held out {target.name} ---")
        for axis, row in fold.items():
            print(f"  {axis}")
            base = row["global mean [baseline]"]
            for name, m in row.items():
                print(f"    {name:<34}PDS={m['discrimination_score_l1']:.3f}  "
                      f"ovl@100={m['overlap_at_100']:.3f}  "
                      f"MAE={m['mae']:.4f}  "
                      f"score={M.vcc_score(m, base)['avg_score']:+.3f}", flush=True)

    print("\n" + "=" * 74)
    print("mean over folds")
    for axis in ("cross-context", "unseen target",
                 "unseen target, gene unmeasured"):
        print(f"\n  {axis}")
        names = list(results[lines[0].name][axis])
        for name in names:
            ms = [results[t][axis][name] for t in results]
            bases = [results[t][axis]["global mean [baseline]"] for t in results]
            sc = np.mean([M.vcc_score(m, b)["avg_score"]
                          for m, b in zip(ms, bases)])
            print(f"    {name:<34}"
                  f"PDS={np.mean([m['discrimination_score_l1'] for m in ms]):.3f}  "
                  f"ovl@100={np.mean([m['overlap_at_100'] for m in ms]):.3f}  "
                  f"MAE={np.mean([m['mae'] for m in ms]):.4f}  "
                  f"score={sc:+.3f}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(results, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
