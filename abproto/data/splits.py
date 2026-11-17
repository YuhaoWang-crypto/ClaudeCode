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

from .anarci import concat_cdrs, has_anarci
from .classify import AbAgComplex
from .cluster_seqs import cluster as seq_cluster, has_mmseqs
from .dedup import cluster_sequences


@dataclass
class Split:
    train: list[AbAgComplex]
    test: list[AbAgComplex]
    stats: dict


def make_splits(complexes: list[AbAgComplex], date_cutoff: str = "2021-09-30",
                antigen_homology: float = 0.4,
                dedup_identity: float = 0.95,
                cdr_dedup: bool = True) -> Split:
    """Return a Split with dedup + time + antigen-homology controls applied.

    Uses mmseqs2 for antigen-homology clustering when available (else k-mer
    proxy). With ``cdr_dedup`` and ANARCI present, near-duplicate antibodies are
    detected on their CONCATENATED CDRs (what actually defines the binder)
    instead of the whole chain.
    """
    # 1. Dedup near-identical complexes (same Ab against same Ag re-deposits).
    use_cdr = cdr_dedup and has_anarci()
    keys = []
    for c in complexes:
        hseq = c.heavy[0].sequence if c.heavy else ""
        lseq = c.light[0].sequence if c.light else ""
        ab_key = ""
        if use_cdr:
            ab_key = concat_cdrs(hseq, lseq)
        if not ab_key:                      # ANARCI off or no CDRs found
            ab_key = c.antibody_seq
        keys.append(ab_key + "|" + c.antigen_seq)
    dd_labels = cluster_sequences(keys, threshold=dedup_identity)
    seen: dict[int, AbAgComplex] = {}
    for c, lab in zip(complexes, dd_labels):
        # keep the highest-resolution (lowest value) representative
        prev = seen.get(int(lab))
        if prev is None or (c.resolution or 1e9) < (prev.resolution or 1e9):
            seen[int(lab)] = c
    uniq = list(seen.values())

    # 2. Antigen-homology clusters (mmseqs2 when available).
    ag_labels = seq_cluster([c.antigen_seq for c in uniq], min_seq_id=antigen_homology)
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
        "identity_backend": "mmseqs2" if has_mmseqs() else "kmer-jaccard(fallback)",
        "cdr_annotation": "anarci" if (cdr_dedup and has_anarci()) else "whole-chain(fallback)",
    }
    return Split(train=train, test=test, stats=stats)
