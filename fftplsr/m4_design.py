"""M4 -- designing a new enzyme variant: would this pipeline have found Com2-IFRS?

M1-M3 establish what the model can and cannot do. This module is the part you
actually reuse: the `design.design_round` driver, exercised on the one question a
design method has to answer -- *given only what was known before the final round,
does it put the eventual winner on the order list?*

The setup reproduces the paper's position at the start of its last round. Com1-IFRS
is the parent; saturation mutagenesis has returned 27 measured single mutants across
6 positions (7, 63, 67, 68, 74, 76); recombining the best of them spans ~1.5k
variants, of which the paper went on to measure 21. The answer turns out to be
Com2-IFRS = N7Y/H63L/K67N/V74W, the best variant in the whole released panel
(2.75x Com1-IFRS, 30.8x the original IFRS).

Two rounds are run, because the difference between them is the real lesson:

  A. singles only (28 measurements)   -- the cheap case, no epistasis observed yet
  B. singles + 92 doubles (120)       -- what the paper actually trained on

Run: ``python3 -m fftplsr.m4_design``
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import datasets
from .design import design_round, improved_sites
from .variants import parse_variant

COM2 = "N7Y/H63L/K67N/V74W"


def _singles_only(table) -> dict[str, float]:
    return {
        label: float(value)
        for label, value in zip(table["Variants"], table["Fitness"])
        if len(parse_variant(label)) <= 1
    }


def _rank_of(ranking: pd.DataFrame, variant: str):
    hit = ranking.index[ranking["variant"] == variant]
    return int(hit[0]) + 1 if len(hit) else None


def _score_picks(picks: pd.DataFrame, measured: dict[str, float]) -> str:
    known = [(row.variant, measured[row.variant]) for row in picks.itertuples()
             if row.variant in measured]
    if not known:
        return "none of the picks were among the 21 variants the paper went on to measure"
    best = max(known, key=lambda kv: kv[1])
    return (
        f"{len(known)}/{len(picks)} picks were later measured; "
        f"mean {np.mean([v for _, v in known]):.2f}x, best {best[0]} at {best[1]:.2f}x"
    )


def design(measured: dict[str, float], parent: str, label: str, pick: int = 8):
    # No per-position pruning: see `improved_sites`. Three of Com2's four
    # substitutions are in the bottom half of their position's singles, so
    # ranking-based pruning would delete the answer before fitting anything.
    sites = improved_sites(measured, threshold=1.05, per_position=None)
    print(f"\n{'=' * 88}\n{label}\n{'=' * 88}")
    print(f"training measurements: {len(measured)}   beneficial sites kept: {len(sites)}")
    print(f"  {', '.join(sites)}")
    report = design_round(
        parent=parent,
        measured=measured,
        sites=sites,
        n_rounds=1,
        cv=None,
        pick=pick,
        verbose=True,
    )
    return report


def main() -> None:
    pd.set_option("display.width", 200)
    table4, parent = datasets.trainset(4)
    full = dict(zip(table4["Variants"], table4["Fitness"]))
    singles = _singles_only(table4)

    panel, _ = datasets.measured_panel("Com1-IFRS")
    truth = dict(zip(panel["Variants"], panel["Fitness"]))
    later = {k: v for k, v in truth.items() if k not in full}

    print("Target: rediscover Com2-IFRS = " + COM2)
    print(f"  measured at {truth[COM2]:.3f}x Com1-IFRS -- the best of the {len(panel)} released variants")
    print(f"  {len(later)} variants were measured after this training set; the rest of the space is unknown")

    results = {}
    for tag, measured in [
        ("ROUND A -- 27 saturation singles only (no epistasis observed yet)", singles),
        ("ROUND B -- singles + 92 measured doubles (as the paper trained)", full),
    ]:
        report = design(measured, parent, tag)
        rank = _rank_of(report.ranking, COM2)
        results[tag] = (report, rank)

        print(f"\n{report.headline}")
        print("\nBaselines (leave-one-out on the training set):")
        print(report.baselines.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
        print(f"\nCom2 ({COM2}) rank: "
              + (f"{rank} of {report.n_space}" if rank else "not in the enumerated space"))
        print(f"\nPick-list ({len(report.picks)} variants):")
        shown = report.picks.copy()
        shown["later_measured"] = [truth.get(v, np.nan) for v in shown["variant"]]
        print(shown.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
        print(f"  -> {_score_picks(report.picks, later)}")
        for note in report.notes:
            print(f"  note: {note}")

    print(f"\n{'=' * 88}\nSUMMARY\n{'=' * 88}")
    rows = []
    for tag, (report, rank) in results.items():
        picks = report.picks["variant"].tolist()
        known = [later[v] for v in picks if v in later]
        rows.append({
            "round": tag.split(" -- ")[0],
            "n_train": report.n_train,
            "descriptor": report.selection.label,
            "cvR2": round(report.selection.cv_r2, 3),
            "best baseline cvR2": round(report.baselines["cvR2"].max(), 3),
            "Com2 rank": f"{rank}/{report.n_space}" if rank else "absent",
            "picks later measured": len(known),
            "best pick measured": round(max(known), 3) if known else None,
        })
    print(pd.DataFrame(rows).to_string(index=False))
    print(f"\n  for reference, the best variant reachable in this space measured "
          f"{max(later.values()):.3f}x and chance over the later-measured set is "
          f"{np.mean(list(later.values())):.3f}x")


if __name__ == "__main__":
    main()
