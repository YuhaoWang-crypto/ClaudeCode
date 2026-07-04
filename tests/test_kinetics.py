"""Regression tests for the mechanistic delivery-kinetics model."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lipidlib.kinetics import DeliveryParams, metrics, analytic_auc


def test_numeric_matches_analytic_auc():
    m = metrics(DeliveryParams())
    rel = abs(m["auc_protein_numeric"] - m["auc_protein_analytic"]) / m["auc_protein_analytic"]
    assert rel < 0.02, f"numeric vs analytic AUC differ by {rel:.1%}"


def test_escape_efficiency_bottleneck():
    # default rates give the known ~1-2% endosomal-escape bottleneck
    assert 0.01 < DeliveryParams().escape_efficiency() < 0.03


def test_auc_linear_in_escape():
    # escape is rate-limiting -> AUC ~ proportional to escape rate
    base = analytic_auc(DeliveryParams())
    p2 = DeliveryParams(); p2.k_escape *= 2
    assert 1.8 < analytic_auc(p2) / base < 2.05


def test_protein_peaks_then_decays():
    m = metrics(DeliveryParams())
    assert 2 < m["t_peak_h"] < 48 and m["peak_protein"] > 0


if __name__ == "__main__":
    test_numeric_matches_analytic_auc()
    test_escape_efficiency_bottleneck()
    test_auc_linear_in_escape()
    test_protein_peaks_then_decays()
    print("all kinetics tests passed")
