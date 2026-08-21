"""Can a language model supply the context-specific part this model cannot?

The measurements say where the headroom is, and it is not everywhere:

* **Effect shape across contexts is nearly noise-limited.**  Cross-line effect
  correlation has a median of 0.11-0.26, and two runs of the *same* K562 screen
  agree at only 0.319 (``gwps --ceiling``).  Nothing -- language model, sequence
  model, or larger transformer -- recovers signal that the assay does not
  reproduce.
* **Effect magnitude is already transferred**, conserved between lines at
  Spearman 0.62-0.69 and read straight off the source measurement.
* **Genes with no measurement anywhere have unlimited headroom**, because the
  current answer for them is "predict no change".  22 of the official 300
  knockdowns are in that position, as are 10,852 of the 18,533 official genes.

That leaves one honest question for a language model, and this module is built
to answer it rather than to assume it:

    **Does the model's prediction change when the cell context changes?**

If two contexts produce the same gene list, the model cannot improve
*transferability* however good its biology is -- it is a context-free prior, and
a context-free prior is exactly what this project already has in the source
measurement.  The test is therefore not "is the language model any good at
biology" but "is its answer a function of the context at all", scored against
measured differential expression in two lines that really do differ.

Deliberately the *easy* version: the cell lines are named, not anonymised as the
challenge anonymises them.  A negative result under the most favourable
conditions is worth more than a negative result under hard ones.

**What it measured.**  The premise turned out to be wrong in an informative way.
The model is *not* context-blind -- its two answers share only a quarter of their
genes, and its stated reasoning is specifically lineage-aware -- and it still
scores **below the context-free baseline in both contexts** (0.105 and 0.089
against 0.123 and 0.138), and below this project's model by a factor of 2.5.
Blending it into that model's ranking at any weight makes it monotonically
worse.  So the failure is not "it ignores the context"; it is that its
context-specific variation is confidently wrong, which is the harder failure to
use.

    python -m virtualcell.llm_prior --prepare        # writes the task file
    python -m virtualcell.llm_prior --score preds/   # scores whatever came back
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from . import metrics as M
from .benchmark import de_truth
from .data import load_all, shared_perturbations
from .model import ContextTransferModel, GlobalMeanBaseline, NaiveTransferBaseline
from .predict import DEFAULT_HYPER

TASK = Path("results/llm_task.json")

# Jurkat and RPE1 are the two atlas lines that map onto official contexts A
# (T lymphoid) and B (RPE-like epithelium), which makes any context effect found
# here directly relevant rather than an analogy.
CONTEXTS = {
    "Jurkat": "Jurkat — human T-lymphoblast leukaemia line, T-lymphoid lineage, "
              "suspension, CD3/TCR-positive, immature T cell",
    "RPE1": "RPE1 (hTERT RPE-1) — human retinal pigment epithelial line, "
            "non-cancerous, near-diploid, adherent epithelium",
}
N_GENES = 50          # per direction, so 100 predicted genes -> overlap@100


def prepare(n_knockdowns: int = 24, seed: int = 0) -> dict:
    """Choose knockdowns and write the prompt payload plus the answer key."""
    lines, _, symbols = load_all()
    by_name = {cl.name: cl for cl in lines}
    common = shared_perturbations(lines)

    # Bias toward knockdowns with a real measured response in both lines: a
    # knockdown that does nothing cannot separate a good prior from a bad one.
    a, b = by_name["Jurkat"], by_name["RPE1"]
    ia = {n: i for i, n in enumerate(a.names)}
    ib = {n: i for i, n in enumerate(b.names)}
    strength = np.array([min(np.linalg.norm(a.delta[ia[p]]),
                             np.linalg.norm(b.delta[ib[p]])) for p in common])
    order = np.argsort(-strength)
    pool = np.array(common)[order[:400]]
    rng = np.random.default_rng(seed)
    perts = np.sort(rng.choice(pool, n_knockdowns, replace=False))

    payload = {
        "contexts": CONTEXTS,
        "knockdowns": perts.tolist(),
        "n_genes_per_direction": N_GENES,
        "gene_universe_size": int(symbols.size),
        "instructions": (
            "For each knockdown in each context, list the gene symbols you "
            "expect to change most. CRISPRi knockdown, steady state, "
            "pseudobulk mRNA. Use HGNC symbols."),
    }
    TASK.parent.mkdir(exist_ok=True)
    TASK.write_text(json.dumps(payload, indent=2))
    print(f"wrote {TASK}: {perts.size} knockdowns x {len(CONTEXTS)} contexts")
    print(f"  {', '.join(perts[:8])}...")
    return payload


def _pseudo_lfc(pred: dict[str, dict[str, list[str]]], perts: np.ndarray,
                symbols: np.ndarray) -> tuple[np.ndarray, float]:
    """Turn ranked up/down symbol lists into a signed effect matrix.

    Rank is preserved as magnitude so that ``overlap_at_k`` sees the model's own
    ordering; genes it did not name score zero.  Also returns the fraction of
    named symbols that exist in the measured gene universe, which bounds what
    the score can be.
    """
    index = {str(s): i for i, s in enumerate(symbols)}
    out = np.zeros((perts.size, symbols.size))
    named = hit = 0
    for i, p in enumerate(perts):
        entry = pred.get(str(p)) or {}
        for key, sign in (("down", -1.0), ("up", 1.0)):
            genes = entry.get(key) or []
            for r, g in enumerate(genes[:N_GENES]):
                named += 1
                j = index.get(str(g).strip().upper())
                if j is None:
                    continue
                hit += 1
                out[i, j] = sign * (N_GENES - r)
    return out, (hit / named if named else 0.0)


def score(pred_dir: Path, task: Path = TASK) -> dict:
    """Score the returned predictions against measured DE, with the baselines."""
    payload = json.loads(task.read_text())
    perts = np.array(payload["knockdowns"])
    lines, _, symbols = load_all()
    by_name = {cl.name: cl for cl in lines}

    preds = {}
    for name in payload["contexts"]:
        f = pred_dir / f"{name}.json"
        if not f.exists():
            raise SystemExit(f"missing {f}")
        preds[name] = json.loads(f.read_text())

    results: dict[str, dict] = {}
    lfcs: dict[str, np.ndarray] = {}
    for name in payload["contexts"]:
        target = by_name[name]
        sources = [cl for cl in lines if cl.name != name]
        truth = de_truth(target, perts)
        idx = {n: i for i, n in enumerate(target.names)}
        true_delta = target.delta[[idx[p] for p in perts]]

        llm, coverage = _pseudo_lfc(preds[name], perts, symbols)
        lfcs[name] = llm

        row = {"symbol coverage": coverage}
        arms = {
            "language model": llm,
            "global mean [challenge baseline]":
                GlobalMeanBaseline().fit(sources, target.mu, symbols).predict(perts),
            "naive transfer":
                NaiveTransferBaseline().fit(sources, target.mu, symbols).predict(perts),
            "ContextTransfer (ours)":
                ContextTransferModel(DEFAULT_HYPER).fit(
                    sources, target.mu, symbols).predict(perts),
        }
        for arm, delta in arms.items():
            row[arm] = float(np.mean(
                M.overlap_at_k(delta, true_delta, truth.sig, k=100)))

        # A weak standalone arm can still carry signal the strong arm lacks, and
        # the DE metrics only read the *ranking*, so blend on rank rather than
        # on value -- the two arms are not on a common scale.
        def rank_signed(a: np.ndarray) -> np.ndarray:
            order = np.argsort(np.argsort(-np.abs(a), axis=1), axis=1)
            return np.sign(a) * (a.shape[1] - order)

        ours = rank_signed(arms["ContextTransfer (ours)"])
        llm_r = rank_signed(llm)
        row["blend"] = {}
        for w in (0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0):
            mix = (1 - w) * ours / max(np.abs(ours).max(), 1) + \
                  w * llm_r / max(np.abs(llm_r).max(), 1)
            row["blend"][f"{w:g}"] = float(np.mean(
                M.overlap_at_k(mix, true_delta, truth.sig, k=100)))
        results[name] = row

    # The decisive number: does the language model's answer depend on context?
    a, b = list(payload["contexts"])
    same = []
    for i in range(perts.size):
        sa = set(np.flatnonzero(lfcs[a][i] != 0))
        sb = set(np.flatnonzero(lfcs[b][i] != 0))
        if sa or sb:
            same.append(len(sa & sb) / max(len(sa | sb), 1))
    ta, tb = by_name[a], by_name[b]
    ida = {n: i for i, n in enumerate(ta.names)}
    idb = {n: i for i, n in enumerate(tb.names)}
    measured = M.pearson_delta(ta.delta[[ida[p] for p in perts]],
                               tb.delta[[idb[p] for p in perts]])
    results["_context"] = {
        "jaccard_between_contexts": float(np.mean(same)),
        "measured_cross_context_r": float(np.median(measured)),
    }
    Path("results/llm_prior.json").write_text(json.dumps(results, indent=2))

    print(f"\nDE overlap@100 against measured differential expression, "
          f"{perts.size} knockdowns\n")
    head = f"{'arm':<36}" + "".join(f"{n:>12}" for n in payload["contexts"])
    print(head); print("-" * len(head))
    for arm in ("global mean [challenge baseline]", "naive transfer",
                "ContextTransfer (ours)", "language model"):
        print(f"{arm:<36}" + "".join(f"{results[n][arm]:>12.4f}"
                                     for n in payload["contexts"]))
    print(f"\n{'symbols found in the measured universe':<36}"
          + "".join(f"{results[n]['symbol coverage']:>12.1%}"
                    for n in payload["contexts"]))

    print("\nBlending the language model into this model's ranking "
          "(weight 0 = ours alone):")
    ws = list(results[list(payload["contexts"])[0]]["blend"])
    print(f"  {'weight':<10}" + "".join(f"{w:>12}" for w in ws))
    for n in payload["contexts"]:
        print(f"  {n:<10}" + "".join(f"{results[n]['blend'][w]:>12.4f}"
                                     for w in ws))
    c = results["_context"]
    print(f"\nOverlap between the two contexts' predicted gene lists: "
          f"{c['jaccard_between_contexts']:.3f} Jaccard")
    print(f"Measured cross-context effect correlation, same knockdowns: "
          f"{c['measured_cross_context_r']:.3f}")
    print("\nRead these two numbers together. The Jaccard is what the language\n"
          "model's context-sensitivity looks like; the overlap@100 above is what\n"
          "that sensitivity was worth. A model that differentiates strongly and\n"
          "still scores below a context-free baseline is not ignoring context --\n"
          "it is confidently varying its answer in the wrong direction.")
    return results


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--prepare", action="store_true")
    ap.add_argument("--score", type=Path)
    ap.add_argument("--n-knockdowns", type=int, default=24)
    args = ap.parse_args()
    if args.prepare:
        prepare(args.n_knockdowns)
    elif args.score:
        score(args.score)
    else:
        ap.error("pass --prepare or --score")


if __name__ == "__main__":
    main()
