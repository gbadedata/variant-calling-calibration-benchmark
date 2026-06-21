"""Tests for the hap.py output adapter.

Uses a small synthetic hap.py-style annotated VCF written to a temp file, so
no network or real hap.py run is needed. Verifies that BD decisions in the
TRUTH and QUERY sample columns map to the correct MatchStatus.
"""

from __future__ import annotations

from src.happy_adapter import parse_happy_vcf
from src.data_models import MatchStatus, VariantType


_HAPPY_VCF = """\
##fileformat=VCFv4.2
##FORMAT=<ID=BD,Number=1,Type=String,Description="Decision">
##FORMAT=<ID=BVT,Number=1,Type=String,Description="Variant type">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY
chr1\t100\t.\tC\tT\t55.0\tPASS\t.\tBD:BVT\tTP:SNP\tTP:SNP
chr1\t200\t.\tA\tG\t30.0\tPASS\t.\tBD:BVT\tN:NOCALL\tFP:SNP
chr1\t300\t.\tG\tA\t.\tPASS\t.\tBD:BVT\tFN:SNP\tN:NOCALL
chr1\t400\t.\tAT\tA\t42.0\tPASS\t.\tBD:BVT\tTP:INDEL\tTP:INDEL
"""


def _write(tmp_path):
    p = tmp_path / "happy.vcf"
    p.write_text(_HAPPY_VCF)
    return p


class TestParseHappyVcf:
    def test_counts_each_status(self, tmp_path) -> None:
        calls = parse_happy_vcf(_write(tmp_path))
        by_status = {}
        for c in calls:
            by_status[c.status] = by_status.get(c.status, 0) + 1
        assert by_status[MatchStatus.TP] == 2
        assert by_status[MatchStatus.FP] == 1
        assert by_status[MatchStatus.FN] == 1

    def test_tp_has_qual_fn_does_not(self, tmp_path) -> None:
        calls = parse_happy_vcf(_write(tmp_path))
        tp = next(c for c in calls if c.status == MatchStatus.TP and c.pos == 100)
        fn = next(c for c in calls if c.status == MatchStatus.FN)
        assert tp.qual == 55.0
        assert fn.qual is None

    def test_indel_typed(self, tmp_path) -> None:
        calls = parse_happy_vcf(_write(tmp_path))
        indel = next(c for c in calls if c.pos == 400)
        assert indel.variant_type == VariantType.INDEL

    def test_feeds_calibration(self, tmp_path) -> None:
        """The adapter output must work with the existing calibration layer."""
        from src.benchmark.calibration import compute_calibration
        calls = parse_happy_vcf(_write(tmp_path))
        report = compute_calibration(calls)
        assert report.n_scored == 3  # 2 TP + 1 FP (FN excluded)

    def test_feeds_concordance(self, tmp_path) -> None:
        from src.benchmark.concordance import concordance_overall
        calls = parse_happy_vcf(_write(tmp_path))
        r = concordance_overall(calls)
        assert r.tp == 2 and r.fp == 1 and r.fn == 1
