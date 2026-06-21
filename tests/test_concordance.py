"""Tests for the concordance engine, verified against hand-computed values."""

from __future__ import annotations

from src.benchmark.concordance import concordance_by_type, concordance_overall
from src.data_loader import load_synthetic
from src.data_models import MatchStatus, VariantCall


def _calls(spec):
    """spec: list of (ref, alt, status) -> VariantCalls at distinct positions."""
    out = []
    for i, (ref, alt, status) in enumerate(spec):
        out.append(VariantCall.make("chr1", 100 + i, ref, alt, status,
                                    qual=50.0 if status != MatchStatus.FN else None))
    return out


class TestOverall:
    def test_known_precision_recall(self) -> None:
        # 3 TP, 1 FP, 1 FN -> precision 3/4, recall 3/4
        calls = _calls([
            ("C", "T", MatchStatus.TP),
            ("C", "T", MatchStatus.TP),
            ("C", "T", MatchStatus.TP),
            ("C", "T", MatchStatus.FP),
            ("C", "T", MatchStatus.FN),
        ])
        r = concordance_overall(calls)
        assert r.tp == 3 and r.fp == 1 and r.fn == 1
        assert abs(r.precision - 0.75) < 1e-9
        assert abs(r.recall - 0.75) < 1e-9
        assert abs(r.f1 - 0.75) < 1e-9

    def test_perfect(self) -> None:
        calls = _calls([("C", "T", MatchStatus.TP)] * 5)
        r = concordance_overall(calls)
        assert r.precision == 1.0 and r.recall == 1.0 and r.f1 == 1.0


class TestByType:
    def test_snv_indel_separated(self) -> None:
        calls = [
            VariantCall.make("chr1", 1, "C", "T", MatchStatus.TP, qual=50),   # SNV TP
            VariantCall.make("chr1", 2, "C", "T", MatchStatus.FP, qual=50),   # SNV FP
            VariantCall.make("chr1", 3, "AT", "A", MatchStatus.TP, qual=50),  # indel TP
        ]
        by_type = concordance_by_type(calls)
        assert by_type["snv"].tp == 1 and by_type["snv"].fp == 1
        assert by_type["indel"].tp == 1 and by_type["indel"].fp == 0

    def test_indel_harder_in_fixture(self) -> None:
        """The synthetic fixture should show indels with lower F1 than SNVs."""
        calls = load_synthetic()
        by_type = concordance_by_type(calls)
        assert by_type["indel"].f1 < by_type["snv"].f1


class TestFixtureSanity:
    def test_overall_runs_on_fixture(self) -> None:
        r = concordance_overall(load_synthetic())
        assert 0 < r.precision < 1
        assert 0 < r.recall < 1
