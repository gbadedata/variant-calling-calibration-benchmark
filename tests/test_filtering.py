"""Tests for the filtering/abstention layer."""

from __future__ import annotations

from src.benchmark.filtering import sweep_filter_thresholds
from src.data_loader import load_synthetic
from src.data_models import MatchStatus, VariantCall


def _call(qual, status, i):
    return VariantCall.make("chr1", 100 + i, "C", "T", status, qual=qual)


class TestSweepMechanics:
    def test_raising_threshold_filters_more(self) -> None:
        report = sweep_filter_thresholds(load_synthetic(), steps=10)
        fracs = [p.filtered_fraction for p in report.points]
        # filtered fraction is non-decreasing as threshold rises
        assert all(fracs[i] <= fracs[i + 1] + 1e-9 for i in range(len(fracs) - 1))

    def test_lowest_threshold_keeps_everything(self) -> None:
        report = sweep_filter_thresholds(load_synthetic(), steps=10)
        assert abs(report.points[0].filtered_fraction) < 1e-9

    def test_precision_improves_with_threshold(self) -> None:
        """Filtering low-QUAL calls should raise retained precision when QUAL
        carries real signal."""
        # Build a set where low QUAL = mostly FP, high QUAL = mostly TP
        calls = []
        for i in range(100):
            calls.append(_call(10.0, MatchStatus.FP, i))       # low qual, false
        for i in range(100, 200):
            calls.append(_call(60.0, MatchStatus.TP, i))       # high qual, true
        report = sweep_filter_thresholds(calls, steps=10)
        # at the top threshold, precision should be near 1.0
        assert report.points[-1].retained_precision > 0.95
        # at the bottom (keep all), precision is ~0.5
        assert report.points[0].retained_precision < 0.6


class TestReferenceThresholds:
    def test_safe_threshold_reaches_target(self) -> None:
        calls = []
        for i in range(100):
            calls.append(_call(10.0, MatchStatus.FP, i))
        for i in range(100, 200):
            calls.append(_call(60.0, MatchStatus.TP, i))
        report = sweep_filter_thresholds(calls, steps=20, target_precision=0.99)
        assert report.safe_threshold is not None
        # the safe threshold must actually achieve the target precision
        pt = next(p for p in report.points if abs(p.threshold - report.safe_threshold) < 1e-6)
        assert pt.retained_precision >= 0.99

    def test_unreachable_target_returns_none(self) -> None:
        # all calls false: no threshold reaches precision 0.99
        calls = [_call(50.0, MatchStatus.FP, i) for i in range(50)]
        report = sweep_filter_thresholds(calls, steps=10, target_precision=0.99)
        assert report.safe_threshold is None

    def test_max_f1_identified(self) -> None:
        report = sweep_filter_thresholds(load_synthetic(), steps=20)
        assert 0 <= report.max_f1 <= 1
        # max_f1 must be the actual maximum over points
        assert abs(report.max_f1 - max(p.retained_f1 for p in report.points)) < 1e-9


class TestStructure:
    def test_to_dict_serialisable(self) -> None:
        import json
        report = sweep_filter_thresholds(load_synthetic())
        json.dumps(report.to_dict())

    def test_empty_safe(self) -> None:
        report = sweep_filter_thresholds([])
        assert report.points == []
        assert report.safe_threshold is None
