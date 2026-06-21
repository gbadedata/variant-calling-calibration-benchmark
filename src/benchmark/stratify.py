"""Stratified evaluation: concordance broken down by genomic context.

A single genome-wide F1 hides where a caller actually struggles. GIAB's
stratification BEDs exist precisely so that true/false positives and false
negatives can be split by region difficulty (low-complexity, segmental
duplications, GC extremes, low mappability). This module computes per-stratum
concordance, which is where the real story of a caller's behaviour lives:
strong in high-confidence regions, weaker in difficult ones.
"""

from __future__ import annotations

import logging

from src.data_models import MatchStatus, StratumResult, VariantCall

logger = logging.getLogger(__name__)


def stratify_by_region(calls: list[VariantCall]) -> dict[str, StratumResult]:
    """Compute concordance within each stratification region.

    A call contributes to every stratum it belongs to (strata can overlap, as
    GIAB's do). FN entries carry their stratum membership too, so recall is
    correctly attributed per region.
    """
    # Collect all stratum names present
    names: set[str] = set()
    for c in calls:
        names.update(c.strata)

    out: dict[str, StratumResult] = {}
    for name in sorted(names):
        subset = [c for c in calls if name in c.strata]
        tp = sum(1 for c in subset if c.status == MatchStatus.TP)
        fp = sum(1 for c in subset if c.status == MatchStatus.FP)
        fn = sum(1 for c in subset if c.status == MatchStatus.FN)
        out[name] = StratumResult(name=name, tp=tp, fp=fp, fn=fn)

    logger.info("stratified_by_region: %d strata", len(out))
    return out


def stratify_by_qual_bin(
    calls: list[VariantCall], n_bins: int = 10,
    qual_min: float = 0.0, qual_max: float = 100.0,
) -> list[StratumResult]:
    """Split calls into equal-width QUAL bins and score concordance per bin.

    Only TP/FP carry QUAL, so FN are excluded here (recall is not a function
    of the caller's emitted quality). This is the precursor to the calibration
    analysis: it shows whether higher-QUAL bins are in fact more precise.
    """
    scored = [c for c in calls if c.qual is not None]
    width = (qual_max - qual_min) / n_bins
    bins: list[StratumResult] = []

    for i in range(n_bins):
        lo = qual_min + i * width
        hi = lo + width
        # Last bin is inclusive of the upper edge
        if i == n_bins - 1:
            members = [c for c in scored if lo <= c.qual <= hi]
        else:
            members = [c for c in scored if lo <= c.qual < hi]
        tp = sum(1 for c in members if c.status == MatchStatus.TP)
        fp = sum(1 for c in members if c.status == MatchStatus.FP)
        bins.append(StratumResult(name=f"QUAL[{lo:.0f},{hi:.0f})", tp=tp, fp=fp, fn=0))

    return bins
