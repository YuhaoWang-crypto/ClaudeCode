"""End-to-end tests using the mock netMHCpan stand-in (no license needed).

Run with:  python -m pytest tests/  (or: python tests/test_pipeline.py)
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pipeline import aggregate, fasta                       # noqa: E402
from pipeline.models import REGISTRY, get_models            # noqa: E402
from pipeline.runner import run_models                      # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FASTA = os.path.join(ROOT, "examples", "example.fasta")
MOCK = os.path.join(HERE, "mock_netMHCpan")


def test_fasta_parse():
    recs = fasta.parse_fasta(FASTA)
    assert len(recs) == 1
    assert recs[0][0].startswith("sp|P0DTC2")
    assert recs[0][1].startswith("MFVFLVLL")


def test_missing_tool_is_skipped_not_fatal():
    specs = get_models(["netmhciipan"])          # not installed
    res = run_models(specs, FASTA, [], ["DRB1_0101"], None, (15,),
                     "/tmp/ep_test_skip", {}, max_workers=2)
    assert len(res) == 1
    assert res[0].status == "not_installed"
    assert res[0].records == []


def test_full_run_with_mock():
    specs = get_models(["netmhcpan"])
    res = run_models(specs, FASTA, ["HLA-A02:01", "HLA-B07:02"], [],
                     (9,), None, "/tmp/ep_test_full",
                     {"netMHCpan": MOCK}, max_workers=2)
    assert res[0].status == "ok"
    assert len(res[0].records) > 0
    # every record carries the essentials
    r0 = res[0].records[0]
    assert r0.peptide and r0.allele and r0.rank is not None
    # consensus produces binder-level hits
    cons = aggregate.consensus(res, top_n=25)
    assert cons and cons[0]["n_models"] >= 1
    assert cons[0]["best_rank_pct"] is not None


def test_biolib_backend_dispatch(monkeypatch=None):
    """BioLib-backed model: stub the cloud call, verify dispatch + CSV parse."""
    from pipeline import biolib_backend, runner
    from pipeline.models import REGISTRY

    def fake_run_app(app_uri, args, outdir, timeout=1800):
        os.makedirs(outdir, exist_ok=True)
        with open(os.path.join(outdir, "raw_output.csv"), "w") as f:
            f.write("Accession,Residue,BepiPred-3.0 linear epitope score,"
                    "BepiPred-3.0 score\n")
            for i, r in enumerate("MFVFLVLLPLVSSQ"):
                sc = 0.05 + (i % 7) * 0.03
                f.write(f"test,{r},{sc:.3f},{sc:.3f}\n")
        return os.path.abspath(outdir)

    orig = biolib_backend.run_app
    biolib_backend.run_app = fake_run_app
    try:
        spec = REGISTRY["bepipred-cloud"]
        assert spec.backend == "biolib" and spec.app_uri == "DTU/BepiPred-3"
        res = runner.run_models([spec], FASTA, [], [], None, None,
                                "/tmp/ep_test_biolib", {}, max_workers=2)
        assert res[0].status == "ok"
        assert len(res[0].records) == 14
        assert any(r.bind_level == "epitope" for r in res[0].records)
        assert "DTU/BepiPred-3" in res[0].command
    finally:
        biolib_backend.run_app = orig


def test_biolib_needs_structure_is_skipped():
    from pipeline import runner
    from pipeline.models import REGISTRY
    spec = REGISTRY["discotope-cloud"]      # needs a PDB
    res = runner.run_models([spec], FASTA, [], [], None, None,
                            "/tmp/ep_test_ds", {}, max_workers=2, structure=None)
    assert res[0].status == "needs_structure"


def test_registry_shape():
    for key, spec in REGISTRY.items():
        assert spec.category in ("mhc_i", "mhc_ii", "bcell_linear",
                                 "bcell_conf", "aux")
        assert callable(spec.build_command) and callable(spec.parse)


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for fn in fns:
        try:
            fn()
            print(f"PASS {fn.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"FAIL {fn.__name__}: {e}")
    print(f"\n{len(fns) - failed}/{len(fns)} passed")
    sys.exit(1 if failed else 0)
