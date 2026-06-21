"""Filtering as abstention: the safe-vs-decisive tradeoff for QUAL thresholds.

This carries the abstention philosophy from clinical interpretation into
variant calling. Applying a QUAL filter is a deferral decision: every call
below the threshold is one the caller chooses NOT to commit to. Raising the
threshold removes false positives (good) but also discards true positives
(costly). The question is identical to the one a clinical model faces when it
returns "uncertain": where is the threshold that keeps you safe without
making you uselessly indecisive?

This module sweeps the QUAL threshold and, at each level, reports:
  - retained precision  (of the calls kept, what fraction are true)
  - retained recall     (of all true variants, what fraction are kept)
  - retained F1
  - filtered fraction   (how much the caller "abstained" on)

It then identifies two reference points:
  - max_f1_threshold: the threshold maximising retained F1 (best overall)
  - safe_threshold:   the lowest threshold achieving a target precision
                      (e.g. 0.99), i.e. the point where retained calls are
                      trustworthy enough to act on without review.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

from src.data_loader import qual_to_confidence
from src.data_models import MatchStatus, VariantCall

logger = logging.getLogger(__name__)


@dataclass
class FilterPoint:
    threshold: float
    retained_precision: float
    retained_recall: float
    retained_f1: float
    filtered_fraction: float
    n_retained: int

    def to_dict(self) -> dict:
        return {
            "threshold": round(self.threshold, 2),
            "retained_precision": round(self.retained_precision, 4),
            "retained_recall": round(self.retained_recall, 4),
            "retained_f1": round(self.retained_f1, 4),
            "filtered_fraction": round(self.filtered_fraction, 4),
            "n_retained": self.n_retained,
        }


@dataclass
class FilterReport:
    points: list[FilterPoint]
    max_f1_threshold: float
    max_f1: float
    safe_threshold: float | None
    target_precision: float

    def to_dict(self) -> dict:
        return {
            "max_f1_threshold": round(self.max_f1_threshold, 2),
            "max_f1": round(self.max_f1, 4),
            "safe_threshold": (
                round(self.safe_threshold, 2) if self.safe_threshold is not None else None
            ),
            "target_precision": self.target_precision,
            "points": [p.to_dict() for p in self.points],
        }


def sweep_filter_thresholds(
    calls: list[VariantCall],
    steps: int = 20,
    target_precision: float = 0.99,
) -> FilterReport:
    """Sweep QUAL thresholds and compute the precision/recall tradeoff.

    At each threshold t, calls with QUAL >= t are "retained" (committed to)
    and calls below t are "filtered" (deferred). Recall is computed against
    the full truth set: total true variants = retained TP + filtered TP + FN.
    Filtering a true positive therefore costs recall, exactly as a clinical
    abstention on a true-positive variant does.

    Args:
        calls: matched VariantCall list.
        steps: number of thresholds to sweep across the observed QUAL range.
        target_precision: the precision defining the "safe" threshold.

    Returns:
        A FilterReport with the full sweep and the two reference thresholds.
    """
    scored = [c for c in calls if c.qual is not None]
    fn_total = sum(1 for c in calls if c.status == MatchStatus.FN)
    total_true = sum(1 for c in scored if c.status == MatchStatus.TP) + fn_total

    if not scored:
        return FilterReport([], 0.0, 0.0, None, target_precision)

    quals = [c.qual for c in scored]
    lo, hi = min(quals), max(quals)
    if hi == lo:
        hi = lo + 1.0

    points: list[FilterPoint] = []
    for i in range(steps + 1):
        t = lo + (hi - lo) * i / steps
        retained = [c for c in scored if c.qual >= t]
        rtp = sum(1 for c in retained if c.status == MatchStatus.TP)
        rfp = sum(1 for c in retained if c.status == MatchStatus.FP)

        precision = rtp / (rtp + rfp) if (rtp + rfp) else 0.0
        recall = rtp / total_true if total_true else 0.0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
        filtered_fraction = 1.0 - len(retained) / len(scored)

        points.append(FilterPoint(
            threshold=t, retained_precision=precision, retained_recall=recall,
            retained_f1=f1, filtered_fraction=filtered_fraction, n_retained=len(retained),
        ))

    # Best overall F1
    best = max(points, key=lambda p: p.retained_f1)
    # Lowest threshold reaching target precision (keeps the most calls while safe)
    safe = next((p.threshold for p in points if p.retained_precision >= target_precision), None)

    report = FilterReport(
        points=points, max_f1_threshold=best.threshold, max_f1=best.retained_f1,
        safe_threshold=safe, target_precision=target_precision,
    )
    logger.info(
        "filter_sweep: max_f1=%.4f @ qual>=%.1f; safe(p>=%.2f) @ %s",
        best.retained_f1, best.threshold, target_precision,
        f"qual>={safe:.1f}" if safe is not None else "unreachable",
    )
    return report


def qual_threshold_as_confidence(threshold: float) -> float:
    """Helper: express a QUAL threshold as its implied confidence floor."""
    return qual_to_confidence(threshold)
