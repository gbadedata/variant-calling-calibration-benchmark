"""Tests for data models, QUAL conversion, and the synthetic loader."""

from __future__ import annotations


from src.data_loader import load_synthetic, match_against_truth, qual_to_confidence
from src.data_models import (
    CalibrationBin,
    MatchStatus,
    StratumResult,
    VariantType,
)


class TestVariantType:
    def test_snv(self) -> None:
        assert VariantType.from_alleles("C", "T") == VariantType.SNV

    def test_insertion_is_indel(self) -> None:
        assert VariantType.from_alleles("A", "AT") == VariantType.INDEL

    def test_deletion_is_indel(self) -> None:
        assert VariantType.from_alleles("AT", "A") == VariantType.INDEL


class TestQualToConfidence:
    def test_qual_30_is_999(self) -> None:
        # QUAL 30 -> P_wrong 0.001 -> confidence 0.999
        assert abs(qual_to_confidence(30) - 0.999) < 1e-6

    def test_qual_10_is_90(self) -> None:
        assert abs(qual_to_confidence(10) - 0.9) < 1e-6

    def test_qual_20_is_99(self) -> None:
        assert abs(qual_to_confidence(20) - 0.99) < 1e-6

    def test_zero_qual_is_zero(self) -> None:
        assert qual_to_confidence(0) == 0.0

    def test_monotonic(self) -> None:
        assert qual_to_confidence(50) > qual_to_confidence(40) > qual_to_confidence(30)


class TestStratumResult:
    def test_precision_recall_f1(self) -> None:
        s = StratumResult(name="x", tp=80, fp=20, fn=20)
        assert abs(s.precision - 0.8) < 1e-9
        assert abs(s.recall - 0.8) < 1e-9
        assert abs(s.f1 - 0.8) < 1e-9

    def test_zero_denominators_safe(self) -> None:
        s = StratumResult(name="empty", tp=0, fp=0, fn=0)
        assert s.precision == 0.0
        assert s.recall == 0.0
        assert s.f1 == 0.0


class TestCalibrationBin:
    def test_empirical_precision(self) -> None:
        b = CalibrationBin(lo=0.9, hi=1.0, n=100, tp=95, fp=5, mean_confidence=0.99)
        assert abs(b.empirical_precision - 0.95) < 1e-9

    def test_gap_is_signed(self) -> None:
        # Overconfident: stated 0.99, empirical 0.95 -> positive gap
        b = CalibrationBin(lo=0.9, hi=1.0, n=100, tp=95, fp=5, mean_confidence=0.99)
        assert b.gap > 0
        assert abs(b.gap - 0.04) < 1e-9


class TestSyntheticLoader:
    def test_deterministic(self) -> None:
        a = load_synthetic(seed=7)
        b = load_synthetic(seed=7)
        assert [c.key for c in a] == [c.key for c in b]

    def test_has_all_statuses(self) -> None:
        calls = load_synthetic()
        statuses = {c.status for c in calls}
        assert MatchStatus.TP in statuses
        assert MatchStatus.FP in statuses
        assert MatchStatus.FN in statuses

    def test_fn_have_no_qual(self) -> None:
        calls = load_synthetic()
        for c in calls:
            if c.status == MatchStatus.FN:
                assert c.qual is None

    def test_tp_fp_have_qual(self) -> None:
        calls = load_synthetic()
        for c in calls:
            if c.status in (MatchStatus.TP, MatchStatus.FP):
                assert c.qual is not None

    def test_difficult_regions_lower_precision(self) -> None:
        """The fixture must encode lower precision in difficult regions."""
        calls = load_synthetic()
        def prec(stratum):
            sub = [c for c in calls if stratum in c.strata and c.qual is not None]
            tp = sum(1 for c in sub if c.status == MatchStatus.TP)
            return tp / len(sub) if sub else 0
        # high_confidence should be cleaner than low_complexity
        assert prec("high_confidence") > prec("low_complexity")


class TestMatchAgainstTruth:
    def test_tp_fp_fn_assignment(self) -> None:
        query = [("chr1", 100, "A", "T", 50.0), ("chr1", 200, "C", "G", 40.0)]
        truth = [("chr1", 100, "A", "T", None), ("chr1", 300, "G", "A", None)]
        calls = match_against_truth(query, truth)
        by_status = {c.key: c.status for c in calls}
        assert by_status["chr1:100:A:T"] == MatchStatus.TP   # in both
        assert by_status["chr1:200:C:G"] == MatchStatus.FP   # query only
        assert by_status["chr1:300:G:A"] == MatchStatus.FN   # truth only
