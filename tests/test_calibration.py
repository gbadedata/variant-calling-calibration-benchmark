"""Tests for the calibration layer.

The key tests build inputs with KNOWN calibration properties and assert the
ECE behaves correctly: near zero for a well-calibrated caller, large for an
overconfident one. This proves the metric measures what it claims.
"""

from __future__ import annotations

import numpy as np

from src.benchmark.calibration import compute_calibration
from src.data_loader import load_synthetic, qual_to_confidence
from src.data_models import MatchStatus, VariantCall


def _call(qual: float, status: MatchStatus, i: int) -> VariantCall:
    return VariantCall.make("chr1", 100 + i, "C", "T", status, qual=qual)


def _well_calibrated_set(n: int = 2000, seed: int = 0) -> list[VariantCall]:
    """Build calls where empirical precision matches stated confidence.

    For each call we pick a QUAL, derive its confidence p, then draw TP with
    probability exactly p. By construction the caller is honest and ECE -> 0.
    """
    rng = np.random.default_rng(seed)
    calls = []
    for i in range(n):
        qual = float(rng.uniform(5, 60))
        p = qual_to_confidence(qual)
        status = MatchStatus.TP if rng.random() < p else MatchStatus.FP
        calls.append(_call(qual, status, i))
    return calls


def _overconfident_set(n: int = 2000, seed: int = 0) -> list[VariantCall]:
    """Build calls where true precision is far below stated confidence."""
    rng = np.random.default_rng(seed)
    calls = []
    for i in range(n):
        qual = float(rng.uniform(30, 60))  # high stated confidence
        true_p = 0.5                       # but only half are correct
        status = MatchStatus.TP if rng.random() < true_p else MatchStatus.FP
        calls.append(_call(qual, status, i))
    return calls


class TestECEBehaviour:
    def test_well_calibrated_low_ece(self) -> None:
        report = compute_calibration(_well_calibrated_set(), n_bins=10)
        assert report.ece < 0.05  # honest caller: small ECE

    def test_overconfident_high_ece(self) -> None:
        report = compute_calibration(_overconfident_set(), n_bins=10)
        assert report.ece > 0.2   # overconfident caller: large ECE

    def test_overconfident_positive_gap(self) -> None:
        report = compute_calibration(_overconfident_set(), n_bins=10)
        # stated confidence exceeds empirical precision overall
        assert report.overall_mean_confidence > report.overall_empirical_precision


class TestStructure:
    def test_excludes_fn(self) -> None:
        calls = [
            _call(50, MatchStatus.TP, 0),
            VariantCall.make("chr1", 999, "C", "T", MatchStatus.FN, qual=None),
        ]
        report = compute_calibration(calls)
        assert report.n_scored == 1

    def test_empty_input_safe(self) -> None:
        report = compute_calibration([])
        assert report.ece == 0.0
        assert report.n_scored == 0

    def test_bins_within_zero_one(self) -> None:
        report = compute_calibration(_well_calibrated_set(), n_bins=10)
        for b in report.bins:
            assert 0.0 <= b.lo < b.hi <= 1.0
            assert 0.0 <= b.empirical_precision <= 1.0

    def test_to_dict_serialisable(self) -> None:
        import json
        report = compute_calibration(load_synthetic())
        json.dumps(report.to_dict())


class TestFixtureCalibration:
    def test_fixture_is_overconfident(self) -> None:
        """The synthetic fixture encodes overconfidence in difficult regions,
        so its overall stated confidence should exceed empirical precision."""
        report = compute_calibration(load_synthetic())
        assert report.overall_mean_confidence > report.overall_empirical_precision
        assert report.ece > 0.0
