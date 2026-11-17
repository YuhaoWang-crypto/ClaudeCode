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
from abproto.selfconsistency import pose_recovery
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


def test_beam_search_best_is_monotone():
    spec = make_spec()
    gen = MockGenerator()
    trace: list[dict] = []
    # diverse=False => pure elitist beam; best AND mean are non-decreasing.
    beam_search(spec, gen, BeamConfig(beam_width=8, branching=4, steps=4, init_pool=48,
                                      seed=3, diverse=False), trace=trace)
    best = [row["best_score"] for row in trace]
    means = [row["mean_score"] for row in trace]
    assert all(b >= a - 1e-9 for a, b in zip(best, best[1:]))
    assert means[-1] >= means[0]


def test_diverse_beam_keeps_distinct_poses():
    # The top score is still retained under MMR, but the diverse beam should
    # cover more distinct backbones than a pure top-k beam.
    from abproto.beam_search import backbone_similarity
    spec = make_spec()
    gen = MockGenerator()
    pool = gen.generate(spec, 40, seed=1)
    from abproto.beam_search import diverse_topk, _default_score
    div = diverse_topk(pool, _default_score, k=8, diversity_weight=0.8)
    topk = sorted(pool, key=_default_score, reverse=True)[:8]

    def mean_pairwise_sim(ds):
        sims = [backbone_similarity(ds[i], ds[j]) for i in range(len(ds)) for j in range(i + 1, len(ds))]
        return sum(sims) / max(len(sims), 1)

    assert mean_pairwise_sim(div) <= mean_pairwise_sim(topk) + 1e-9


def test_smc_runs_and_ranks():
    from abproto.beam_search import SMCConfig, smc_search
    spec = make_spec()
    gen = MockGenerator()
    trace: list[dict] = []
    out = smc_search(spec, gen, SMCConfig(n_particles=24, steps=3, seed=0), trace=trace)
    assert len(out) == 24
    assert all("ess" in row for row in trace)
    # final best should be >= initial best (reweighting concentrates on good modes)
    assert trace[-1]["best_score"] >= trace[0]["best_score"] - 1e-9


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


def test_interface_confidence_preferred():
    from abproto.confidence import confidence_score
    from abproto.types import Confidence
    base = dict(plddt=80.0, pae=5.0, pde=4.0, iptm=0.7)
    good_iface = Confidence(**base, interface_pae=2.0, interface_plddt=85.0, contact_prob=0.9)
    bad_iface = Confidence(**base, interface_pae=12.0, interface_plddt=50.0, contact_prob=0.2)
    assert good_iface.has_interface()
    # identical global terms => ranking decided by the interface heads
    assert confidence_score(good_iface) > confidence_score(bad_iface)


def test_interface_rmsd_recovery_tracks_quality():
    spec = make_spec()
    gen, folder = MockGenerator(), MockFolder()
    # higher-quality designs should recover their pose better (lower RMSD)
    ds = sorted(gen.generate(spec, 20, seed=2),
                key=lambda d: d.confidence.iptm)
    low = pose_recovery(ds[0], folder, spec, seed=0)[0]      # worst quality
    high = pose_recovery(ds[-1], folder, spec, seed=0)[0]    # best quality
    assert high <= low
    assert high >= 0.0


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
