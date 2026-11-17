"""Leakage-controlled train/test splits.

Two complementary controls, both used by Pythia-1 / Chai-2:

  - time split  : train on complexes deposited before a cutoff date, test on
                  those after (mimics prospective design on novel structures).
  - homology split: cluster by ANTIGEN sequence and keep whole clusters on one
                  side, so no antigen homolog of a test target is ever in train.

The strict split is the intersection: a test item must be both after the cutoff
AND in an antigen cluster with no training member.
"""

from __future__ import annotations

from dataclasses import dataclass

from .classify import AbAgComplex
from .dedup import cluster_sequences


@dataclass
class Split:
    train: list[AbAgComplex]
    test: list[AbAgComplex]
    stats: dict


def make_splits(complexes: list[AbAgComplex], date_cutoff: str = "2021-09-30",
                antigen_homology: float = 0.4,
                dedup_identity: float = 0.95) -> Split:
    """Return a Split with dedup + time + antigen-homology controls applied."""
    # 1. Dedup near-identical complexes (same Ab against same Ag re-deposits).
    keys = [c.antibody_seq + "|" + c.antigen_seq for c in complexes]
    dd_labels = cluster_sequences(keys, threshold=dedup_identity)
    seen: dict[int, AbAgComplex] = {}
    for c, lab in zip(complexes, dd_labels):
        # keep the highest-resolution (lowest value) representative
        prev = seen.get(int(lab))
        if prev is None or (c.resolution or 1e9) < (prev.resolution or 1e9):
            seen[int(lab)] = c
    uniq = list(seen.values())

    # 2. Antigen-homology clusters.
    ag_labels = cluster_sequences([c.antigen_seq for c in uniq], threshold=antigen_homology)
    cluster_of = {id(c): int(lab) for c, lab in zip(uniq, ag_labels)}

    # 3. Time split, then enforce cluster-disjointness: any antigen cluster that
    #    has a training member is forced entirely into train.
    after = [c for c in uniq if c.deposit_date and c.deposit_date > date_cutoff]
    before = [c for c in uniq if not (c.deposit_date and c.deposit_date > date_cutoff)]
    train_clusters = {cluster_of[id(c)] for c in before}

    train, test = list(before), []
    for c in after:
        if cluster_of[id(c)] in train_clusters:
            train.append(c)          # antigen homolog seen in train -> leak, keep in train
        else:
            test.append(c)

    stats = {
        "input": len(complexes),
        "after_dedup": len(uniq),
        "antigen_clusters": len(set(ag_labels.tolist())),
        "date_cutoff": date_cutoff,
        "train": len(train),
        "test": len(test),
        "test_gpcr": sum(1 for c in test if c.is_gpcr),
        "leak_forced_to_train": len(after) - len(test),
    }
    return Split(train=train, test=test, stats=stats)
