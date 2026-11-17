"""Smoke + unit tests. Run: python -m pytest abproto/tests -q  (or see run_tests.py)."""

from __future__ import annotations

import numpy as np

from abproto.backends.mock import MockFolder, MockGenerator
from abproto.beam_search import BeamConfig, beam_search
from abproto.cluster import cluster
from abproto.ctmc import redesign
from abproto.filters import cdr_liabilities
from abproto.geometry import kabsch_rmsd
from abproto.pipeline import DesignPipeline, PipelineConfig
from abproto.types import DesignSpec


def make_spec() -> DesignSpec:
    return DesignSpec(
        target_name="TEST",
        antigen_seq="ACDEFGHIKLMNPQRSTVWY" * 3,
        heavy_framework="EVQLVESGGGLVQPGGSLRLSCAAS",
        light_framework="DIQMTQSPSSLSASVGDRVTITCRAS",
        cdr_lengths={"H1": 5, "H2": 6, "H3": 11, "L1": 6, "L2": 3, "L3": 9},
    )


def test_kabsch_identity_is_zero():
    rng = np.random.default_rng(0)
    P = rng.normal(size=(30, 3))
    assert kabsch_rmsd(P, P.copy()) < 1e-9


def test_kabsch_rotation_invariant():
    rng = np.random.default_rng(1)
    P = rng.normal(size=(40, 3))
    theta = 0.7
    R = np.array([[np.cos(theta), -np.sin(theta), 0],
                  [np.sin(theta), np.cos(theta), 0],
                  [0, 0, 1]])
    Q = P @ R.T + np.array([5.0, -2.0, 1.0])
    assert kabsch_rmsd(P, Q) < 1e-6


def test_beam_search_improves_mean_score():
    spec = make_spec()
    gen = MockGenerator()
    trace: list[dict] = []
    beam_search(spec, gen, BeamConfig(beam_width=8, branching=4, steps=4, init_pool=48, seed=3),
                trace=trace)
    # Mean beam score should be non-decreasing across steps (elitism + refine).
    means = [row["mean_score"] for row in trace]
    assert means[-1] >= means[0]
    assert all(b >= a - 1e-9 for a, b in zip(means, means[1:]))


def test_ctmc_removes_nglyc_and_cys():
    spec = make_spec()
    gen = MockGenerator()
    d = gen.generate(spec, 1, seed=5)[0]
    variants = redesign(d, 5)
    for v in variants:
        liab = cdr_liabilities(v)
        assert liab["nglyc"] == 0
        assert liab["cys"] == 0
        assert v.provenance == "ctmc"


def test_cluster_singletons_and_dups():
    spec = make_spec()
    gen = MockGenerator()
    ds = gen.generate(spec, 4, seed=7)
    # identical backbones should merge; distinct should separate at small cutoff
    ds[1].backbone = {k: v.copy() for k, v in ds[0].backbone.items()}
    ds[1].heavy_seq, ds[1].light_seq = ds[0].heavy_seq, ds[0].light_seq
    ds[1].cdr_spans = dict(ds[0].cdr_spans)
    labels = cluster(ds, cutoff=0.5)
    assert labels[0] == labels[1]


def test_pipeline_end_to_end():
    spec = make_spec()
    cfg = PipelineConfig(
        beam=BeamConfig(beam_width=12, branching=4, steps=3, init_pool=48, seed=0),
        ctmc_variants=3, max_candidates=48, seed=0,
    )
    pipe = DesignPipeline(MockGenerator(), MockFolder(), cfg)
    results = pipe.run(spec)
    assert 0 < len(results) <= 48
    # representatives are flagged and scores are sorted-ish at the top
    assert any(s.is_representative for s in results)
    assert all(s.passed_filters for s in results)
    assert pipe.trace["n_clusters"] >= 1
    # every selected design has zero CDR glycosylation sequons (hard constraint)
    for s in results:
        assert s.liabilities["nglyc"] == 0


if __name__ == "__main__":
    # Allow running without pytest.
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"PASS {name}")
    print("all tests passed")
