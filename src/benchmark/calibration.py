"""Calibration layer: is the caller's stated confidence honest?

Every variant caller emits a QUAL score that encodes a claimed probability
the call is correct (phred scale). Calibration asks whether that claim holds
empirically: among all calls the caller stamped at ~99% confidence, are 99%
actually true positives? A caller can have decent overall precision yet be
badly miscalibrated, systematically overconfident in difficult regions, for
instance, which is dangerous because downstream filters trust QUAL.

This module:
  1. Bins TP/FP calls by their QUAL-implied confidence.
  2. For each bin, compares mean stated confidence to empirical precision.
  3. Computes Expected Calibration Error (ECE): the support-weighted mean
     absolute gap between stated confidence and empirical precision. ECE = 0
     means perfectly honest confidence; larger ECE means more miscalibration.

ECE is the standard scalar summary of a calibration curve, used widely in
the ML literature for exactly this purpose. Reporting it puts variant-caller
confidence on the same footing as any probabilistic classifier.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

from src.data_loader import qual_to_confidence
from src.data_models import CalibrationBin, MatchStatus, VariantCall

logger = logging.getLogger(__name__)


@dataclass
class CalibrationReport:
    bins: list[CalibrationBin]
    ece: float
    n_scored: int
    overall_mean_confidence: float
    overall_empirical_precision: float

    def to_dict(self) -> dict:
        return {
            "ece": round(self.ece, 4),
            "n_scored": self.n_scored,
            "overall_mean_confidence": round(self.overall_mean_confidence, 4),
            "overall_empirical_precision": round(self.overall_empirical_precision, 4),
            "overall_gap": round(
                self.overall_mean_confidence - self.overall_empirical_precision, 4
            ),
            "bins": [b.to_dict() for b in self.bins],
        }


def compute_calibration(
    calls: list[VariantCall], n_bins: int = 10
) -> CalibrationReport:
    """Build the calibration curve and ECE from TP/FP calls.

    Calls are binned by their QUAL-implied confidence into n_bins equal-width
    bins over [0, 1]. FN are excluded (calibration concerns emitted
    confidence, which FN do not have).

    Args:
        calls: matched VariantCall list.
        n_bins: number of confidence bins.

    Returns:
        A CalibrationReport with per-bin detail and the scalar ECE.
    """
    scored = [c for c in calls if c.qual is not None and c.status in (MatchStatus.TP, MatchStatus.FP)]
    n = len(scored)
    if n == 0:
        return CalibrationReport([], 0.0, 0, 0.0, 0.0)

    # Precompute each call's stated confidence
    conf = [(c, qual_to_confidence(c.qual)) for c in scored]

    width = 1.0 / n_bins
    bins: list[CalibrationBin] = []
    ece = 0.0

    for i in range(n_bins):
        lo = i * width
        hi = lo + width
        if i == n_bins - 1:
            members = [(c, p) for (c, p) in conf if lo <= p <= hi]
        else:
            members = [(c, p) for (c, p) in conf if lo <= p < hi]

        if not members:
            continue

        tp = sum(1 for (c, _p) in members if c.status == MatchStatus.TP)
        fp = sum(1 for (c, _p) in members if c.status == MatchStatus.FP)
        mean_conf = sum(p for (_c, p) in members) / len(members)

        cbin = CalibrationBin(
            lo=lo, hi=hi, n=len(members), tp=tp, fp=fp, mean_confidence=mean_conf
        )
        bins.append(cbin)
        # ECE contribution: support-weighted absolute gap
        ece += (len(members) / n) * abs(cbin.gap)

    overall_conf = sum(p for (_c, p) in conf) / n
    tp_total = sum(1 for (c, _p) in conf if c.status == MatchStatus.TP)
    overall_prec = tp_total / n

    report = CalibrationReport(
        bins=bins, ece=ece, n_scored=n,
        overall_mean_confidence=overall_conf,
        overall_empirical_precision=overall_prec,
    )
    logger.info(
        "calibration: ece=%.4f mean_conf=%.4f empirical_prec=%.4f n=%d",
        ece, overall_conf, overall_prec, n,
    )
    return report
