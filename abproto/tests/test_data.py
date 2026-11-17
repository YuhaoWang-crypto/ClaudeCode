"""Offline tests for the Stage 1 data pipeline (no network)."""

from __future__ import annotations

from abproto.data.classify import to_complexes
from abproto.data.dedup import cluster_sequences, jaccard, dedup_representatives
from abproto.data.rcsb import Entity, Entry
from abproto.data.splits import make_splits


def _entry(pdb, date, *ents):
    return Entry(pdb_id=pdb, deposit_date=date, release_date=date, resolution=2.0,
                 entities=[Entity(f"{pdb}_{i}", d, s, "Protein")
                           for i, (d, s) in enumerate(ents, 1)])


LYS = "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"
HER2 = "STQVCTGTDMKLRLPASPETHLDMLRHLYQGCQVVQGNLELTYLPTNASLSFLQDIQEVQGYVLIAHNQVRQVPLQRLRIVRGTQLFEDNYALAVLDNGDPLNNTTPVTGASPG"
HEAVY = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAK"
LIGHT = "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPL"


def test_classify_detects_complex_and_skips_non_complex():
    entries = [
        _entry("AAAA", "2020-01-01",
               ("Herceptin Fab heavy chain", HEAVY),
               ("Herceptin Fab light chain", LIGHT),
               ("Receptor tyrosine kinase erbB-2", HER2)),
        _entry("BBBB", "2020-01-01", ("some enzyme", LYS)),  # no antibody -> skip
    ]
    good, skipped = to_complexes(entries)
    assert len(good) == 1
    assert good[0].pdb_id == "AAAA"
    assert good[0].heavy and good[0].light and good[0].antigen
    assert "BBBB" in skipped


def test_jaccard_and_clustering():
    assert jaccard(set(), {"a"}) == 0.0
    labels = cluster_sequences([LYS, LYS, HER2], threshold=0.9)
    assert labels[0] == labels[1]          # identical cluster
    assert labels[0] != labels[2]          # different antigen separate
    reps = dedup_representatives([LYS, LYS, HER2], threshold=0.9)
    assert len(reps) == 2


def test_time_split_holds_out_post_cutoff_novel_antigen():
    entries = [
        _entry("OLD1", "2018-05-01",
               ("Fab heavy chain", HEAVY), ("Fab light chain", LIGHT),
               ("lysozyme", LYS)),
        _entry("NEW1", "2023-01-01",
               ("Fab heavy chain", HEAVY), ("Fab light chain", LIGHT),
               ("CXC chemokine receptor", HER2)),  # novel antigen, post-cutoff
    ]
    good, _ = to_complexes(entries)
    split = make_splits(good, date_cutoff="2021-09-30", antigen_homology=0.4)
    test_ids = {c.pdb_id for c in split.test}
    train_ids = {c.pdb_id for c in split.train}
    assert "NEW1" in test_ids and "OLD1" in train_ids


def test_homology_leak_forced_to_train():
    # Same antigen (lysozyme) bound by two DIFFERENT antibodies, before AND after
    # the cutoff. The antibodies differ enough to survive dedup, but the shared
    # antigen cluster must pull the post-cutoff one back into train (no leak).
    heavy2 = "QVQLQESGPGLVKPSETLSLTCTVSGGSISSYYWSWIRQPPGKGLEWIGYIYYSGSTNYNPSLKS"
    light2 = "EIVLTQSPGTLSLSPGERATLSCRASQSVSSSYLAWYQQKPGQAPRLLIYGASSRATGIPDRFSG"
    entries = [
        _entry("OLD", "2018-05-01",
               ("Fab heavy chain", HEAVY), ("Fab light chain", LIGHT), ("lysozyme", LYS)),
        _entry("NEW", "2023-01-01",
               ("Fab heavy chain", heavy2), ("Fab light chain", light2), ("lysozyme", LYS)),
    ]
    good, _ = to_complexes(entries)
    split = make_splits(good, date_cutoff="2021-09-30", antigen_homology=0.4)
    assert all(c.pdb_id != "NEW" for c in split.test)   # leaked homolog not in test
    assert split.stats["leak_forced_to_train"] >= 1


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(f"PASS {name}")
    print("all data tests passed")
