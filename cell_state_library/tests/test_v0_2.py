"""Software tests for v0.2. Run: python -m unittest discover -s cell_state_library/tests -t .

DB-dependent tests are skipped if the library has not been built.
"""
import math
import random
import sqlite3
import unittest

import numpy as np

from cell_state_library import benchmark, build, identifiability as ident, panel_select, query
from cell_state_library.nnls import score_nnls
from cell_state_library.pair_code import STRENGTH_GRID, class_loglik, encode, hamming, score_pairs


def v01_score(a, g, b):
    """Core of v0.1 query.py (two-reference closed form), copied for regression."""
    norm = lambda x: math.sqrt(sum(t * t for t in x))
    norms = [norm(x) for x in (a, g, b)]
    a, g, b = [[t / n for t in x] for x, n in zip((a, g, b), norms)]
    dot = lambda x, y: sum(u * v for u, v in zip(x, y))
    ag, ab, gb = dot(a, g), dot(a, b), dot(g, b)
    det = 1 - ag * ag
    cands = [(0., 0.), (max(0., ab), 0.), (0., max(0., gb))]
    xa, xg = (ab - ag * gb) / det, (gb - ag * ab) / det
    if xa >= 0 and xg >= 0:
        cands.append((xa, xg))
    loss = lambda x: sum((x[0] * u + x[1] * v - w) ** 2 for u, v, w in zip(a, g, b))
    x = min(cands, key=loss)
    return x, math.sqrt(loss(x))


def rand_profile(rng, genes):
    return {gname: rng.gauss(0, 1) for gname in genes}


class TestNNLS(unittest.TestCase):
    def test_matches_v01_on_random_problems(self):
        rng = random.Random(1)
        genes = [f"G{i}" for i in range(60)]
        for _ in range(200):
            a, g, b = (rand_profile(rng, genes) for _ in range(3))
            r = score_nnls({"IFNA": a, "IFNG": g}, b)
            x, res = v01_score([a[k] for k in genes], [g[k] for k in genes], [b[k] for k in genes])
            self.assertAlmostEqual(r["scores"]["IFNA"], x[0], places=9)
            self.assertAlmostEqual(r["scores"]["IFNG"], x[1], places=9)
            self.assertAlmostEqual(r["relative_residual_norm"], res, places=9)

    def test_k3_satisfies_kkt(self):
        rng = np.random.default_rng(2)
        genes = [f"G{i}" for i in range(50)]
        for _ in range(50):
            refs = {f"R{k}": dict(zip(genes, rng.normal(size=50))) for k in range(3)}
            b = dict(zip(genes, rng.normal(size=50)))
            x = score_nnls(refs, b)["scores"]
            A = np.column_stack([np.array([refs[n][gn] for gn in genes]) for n in sorted(refs)])
            A /= np.linalg.norm(A, axis=0)
            bv = np.array([b[gn] for gn in genes])
            bv /= np.linalg.norm(bv)
            xv = np.array([x[n] for n in sorted(refs)])
            grad = A.T @ (A @ xv - bv)
            for xi, gi in zip(xv, grad):
                if xi > 1e-9:
                    self.assertAlmostEqual(gi, 0, places=7)
                else:
                    self.assertGreaterEqual(gi, -1e-7)

    def test_vectorised_twin_matches(self):
        rng = np.random.default_rng(3)
        a, g = rng.normal(size=40), rng.normal(size=40)
        B = rng.normal(size=(30, 40))
        B[rng.random(B.shape) < 0.2] = np.nan
        xa, xg, res = benchmark.nnls2(B, a, g)
        genes = [f"G{i}" for i in range(40)]
        for row in range(30):
            target = {genes[i]: B[row, i] for i in range(40) if not np.isnan(B[row, i])}
            r = score_nnls({"IFNA": dict(zip(genes, a)), "IFNG": dict(zip(genes, g))}, target, enforce_gate=False)
            self.assertAlmostEqual(r["scores"]["IFNA"], xa[row], places=8)
            self.assertAlmostEqual(r["scores"]["IFNG"], xg[row], places=8)
            self.assertAlmostEqual(r["relative_residual_norm"], res[row], places=6)

    def test_rejects_bad_input(self):
        refs = {"IFNA": {f"G{i}": 1.0 + i for i in range(40)}, "IFNG": {f"G{i}": (-1) ** i for i in range(40)}}
        with self.assertRaises(ValueError):
            score_nnls(refs, {"G1": 1.0})                      # coverage gate
        with self.assertRaises(ValueError):
            score_nnls(refs, {f"G{i}": 0.0 for i in range(40)})  # zero signal
        with self.assertRaises(ValueError):
            score_nnls(refs, {"G1": float("inf")})


class TestPairCode(unittest.TestCase):
    panel = [{"gene_i": "A", "gene_j": "B", "d_ifna": 2.0, "d_ifng": -1.5},
             {"gene_i": "C", "gene_j": "D", "d_ifna": 1.5, "d_ifng": 1.0},
             {"gene_i": "E", "gene_j": "F", "d_ifna": 0.2, "d_ifng": -2.0}]

    def test_shift_invariance(self):
        rng = random.Random(4)
        for _ in range(100):
            b = {k: rng.gauss(0, 2) for k in "ABCDEF"}
            c = rng.uniform(-5, 5)
            shifted = {k: v + c for k, v in b.items()}
            self.assertEqual(encode(self.panel, b, 0.7)[:3], encode(self.panel, shifted, 0.7)[:3])

    def test_missing_genes_drop_pairs(self):
        r = score_pairs(self.panel, {"A": 3, "B": 0, "C": 1}, 0.5, 0.7, ["IFNA", "IFNG"])
        self.assertEqual(r["code"], "+..")
        self.assertEqual(r["pairs_measured"], 1)
        with self.assertRaises(ValueError):
            score_pairs(self.panel, {"Z": 1.0}, 0.5, 0.7, ["IFNA"])

    def test_hamming(self):
        s1, p1, n1, _ = encode(self.panel, {"A": 3, "B": 0, "C": 0, "D": 0, "E": 0, "F": 3}, 0.7)
        s2, p2, n2, _ = encode(self.panel, {"A": 0, "B": 3, "C": 0, "D": 0, "E": 0, "F": 3}, 0.7)
        self.assertEqual((s1, s2), ([1, 0, -1], [-1, 0, -1]))
        self.assertEqual(hamming((p1, n1), (p2, n2)), 1)

    def test_numpy_twin_matches_stdlib(self):
        rng = np.random.default_rng(5)
        gi = {k: i for i, k in enumerate("ABCDEF")}
        B = rng.normal(0, 1.5, (50, 6))
        B[rng.random(B.shape) < 0.15] = np.nan
        ll_np = panel_select.loglik_np(self.panel, B, gi, 0.5, 0.7, "d_ifng", STRENGTH_GRID)
        for row in range(50):
            target = {k: B[row, i] for k, i in gi.items() if not np.isnan(B[row, i])}
            symbols = encode(self.panel, target, 0.7)[0]
            ll = class_loglik(self.panel, symbols, "d_ifng", 0.5, 0.7, STRENGTH_GRID)
            self.assertAlmostEqual(ll, ll_np[row], places=4)

    def test_selection_uses_disjoint_genes(self):
        rng = np.random.default_rng(6)
        genes = [f"G{i}" for i in range(40)]
        panel = panel_select.select_panel(genes, rng.normal(size=40), rng.normal(size=40), ["IFNA", "IFNG"], 8, 0.5, 0.7)
        used = [p["gene_i"] for p in panel] + [p["gene_j"] for p in panel]
        self.assertEqual(len(used), len(set(used)))
        self.assertTrue(all(p["d_ifna"] >= 0 for p in panel))


class TestIdentifiability(unittest.TestCase):
    def test_collinear_weak_gamma_fails(self):
        rng = np.random.default_rng(7)
        a = rng.normal(1, 1, 200)
        g = 0.05 * a + rng.normal(0, 0.01, 200)
        self.assertEqual(ident.check(a, g)["status"], "ifna_only")

    def test_strong_distinct_passes(self):
        rng = np.random.default_rng(8)
        self.assertEqual(ident.check(rng.normal(0, 1, 200), rng.normal(0, 1, 200))["status"], "both")

    def test_too_few_genes(self):
        self.assertEqual(ident.check([5, 1, 2] * 5, [1, 5, 2] * 5)["status"], "none")


@unittest.skipUnless(build.DB.exists(), "library not built")
class TestBuiltLibrary(unittest.TestCase):
    def setUp(self):
        self.con = sqlite3.connect(build.DB)

    def tearDown(self):
        self.con.close()

    def test_thirteen_active_subsets(self):
        n = self.con.execute("SELECT COUNT(DISTINCT cell_type_id) FROM reference_profiles WHERE status LIKE 'active%'").fetchone()[0]
        self.assertEqual(n, 13)
        self.assertEqual(self.con.execute("SELECT COUNT(*) FROM identifiability").fetchone()[0], 13)
        self.assertEqual(self.con.execute("SELECT COUNT(*) FROM pair_panels").fetchone()[0], 13 * build.PANEL_K)

    def test_long_format_round_trip(self):
        wide = self.con.execute("SELECT ifna_log2fc, ifng_log2fc FROM v_reference_wide WHERE cell_type_id='Monocyte' AND gene='CXCL9'").fetchone()
        self.assertAlmostEqual(wide[0], -0.0, delta=0.05)
        self.assertGreater(wide[1], 4)
        self.assertEqual(self.con.execute("SELECT COUNT(*) FROM v_reference_wide WHERE cell_type_id='Monocyte'").fetchone()[0], 824)

    def test_query_withholds_unidentifiable_and_rejects_staged(self):
        rows = dict(self.con.execute(
            "SELECT v.gene, v.log2fc FROM reference_values v JOIN reference_profiles p USING(profile_id) "
            "WHERE p.cell_type_id='NK.CD56dim' AND p.condition_id='IFNG'"))
        r = query.score({"cell_type": "NK.CD56dim", "celltype_level": "L2", "contrast": "test", "log2fc": rows}, build.DB)
        self.assertIsNone(r["nnls"]["scores"]["IFNg"]["value"])
        self.assertNotIn("IFNG", r["pair_code"]["loglik"])
        with self.assertRaises(ValueError):
            query.score({"cell_type": "Tcell", "celltype_level": "L1", "contrast": "test", "log2fc": rows}, build.DB)

    def test_alias_resolution(self):
        rows = dict(self.con.execute(
            "SELECT v.gene, v.log2fc FROM reference_values v JOIN reference_profiles p USING(profile_id) "
            "WHERE p.cell_type_id='Monocyte' AND p.condition_id='IFNA'"))
        for level in ("L1", "L2"):
            r = query.score({"cell_type": "Monocyte", "celltype_level": level, "contrast": "test", "log2fc": rows}, build.DB)
            self.assertAlmostEqual(r["nnls"]["scores"]["IFNa"]["value"], 1.0, places=9)
            self.assertEqual(r["pair_code"]["best_class"], "IFNA")


if __name__ == "__main__":
    unittest.main()
