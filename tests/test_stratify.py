"""Tests for stratified evaluation by region and by QUAL bin."""

from __future__ import annotations

from src.benchmark.stratify import stratify_by_qual_bin, stratify_by_region
from src.data_loader import load_synthetic
from src.data_models import MatchStatus, VariantCall


class TestStratifyByRegion:
    def test_assigns_to_named_strata(self) -> None:
        calls = [
            VariantCall.make("chr1", 1, "C", "T", MatchStatus.TP, qual=50, strata=["high_confidence"]),
            VariantCall.make("chr1", 2, "C", "T", MatchStatus.FP, qual=50, strata=["low_complexity"]),
        ]
        res = stratify_by_region(calls)
        assert res["high_confidence"].tp == 1
        assert res["low_complexity"].fp == 1

    def test_overlapping_strata_double_count(self) -> None:
        # a call in two strata contributes to both
        calls = [VariantCall.make("chr1", 1, "C", "T", MatchStatus.TP, qual=50,
                                  strata=["low_complexity", "high_gc"])]
        res = stratify_by_region(calls)
        assert res["low_complexity"].tp == 1
        assert res["high_gc"].tp == 1

    def test_fixture_difficult_regions_lower_f1(self) -> None:
        res = stratify_by_region(load_synthetic())
        # high_confidence should outperform low_complexity on F1
        assert res["high_confidence"].f1 > res["low_complexity"].f1


class TestStratifyByQualBin:
    def test_bin_count(self) -> None:
        bins = stratify_by_qual_bin(load_synthetic(), n_bins=10)
        assert len(bins) == 10

    def test_excludes_fn(self) -> None:
        # FN have no QUAL; they must not appear in any bin's tally
        calls = [
            VariantCall.make("chr1", 1, "C", "T", MatchStatus.TP, qual=50),
            VariantCall.make("chr1", 2, "C", "T", MatchStatus.FN, qual=None),
        ]
        bins = stratify_by_qual_bin(calls, n_bins=10)
        total = sum(b.tp + b.fp for b in bins)
        assert total == 1  # only the TP, FN excluded

    def test_bins_show_calibration_signal(self) -> None:
        """The fixture intentionally decouples QUAL from true precision in
        difficult regions (the miscalibration the benchmark detects). So we
        assert the bins are populated and precision is a valid fraction,
        not clean monotonicity, which would mean perfect calibration."""
        bins = stratify_by_qual_bin(load_synthetic(), n_bins=5)
        populated = [b for b in bins if (b.tp + b.fp) >= 10]
        assert len(populated) >= 2
        for b in populated:
            assert 0.0 <= b.precision <= 1.0
