"""Concordance engine: precision, recall, F1 against GIAB truth.

This is the base scoring layer. It consumes the matched VariantCall list
(each already tagged TP / FP / FN) and produces concordance metrics overall
and split by variant type (SNV vs indel).

It is deliberately small and pure: the interesting analysis (stratification,
calibration, abstention) is layered on top in later modules. Keeping the
base layer simple makes every higher layer testable in isolation.
"""

from __future__ import annotations

import logging

from src.data_models import MatchStatus, StratumResult, VariantCall, VariantType

logger = logging.getLogger(__name__)


def _tally(calls: list[VariantCall], name: str) -> StratumResult:
    """Count TP/FP/FN over a list of calls into a StratumResult."""
    tp = sum(1 for c in calls if c.status == MatchStatus.TP)
    fp = sum(1 for c in calls if c.status == MatchStatus.FP)
    fn = sum(1 for c in calls if c.status == MatchStatus.FN)
    return StratumResult(name=name, tp=tp, fp=fp, fn=fn)


def concordance_overall(calls: list[VariantCall]) -> StratumResult:
    """Overall precision/recall/F1 across all calls."""
    result = _tally(calls, "overall")
    logger.info(
        "concordance_overall: tp=%d fp=%d fn=%d p=%.4f r=%.4f f1=%.4f",
        result.tp, result.fp, result.fn,
        result.precision, result.recall, result.f1,
    )
    return result


def concordance_by_type(calls: list[VariantCall]) -> dict[str, StratumResult]:
    """Concordance split by variant type (SNV vs indel)."""
    out: dict[str, StratumResult] = {}
    for vtype in VariantType:
        subset = [c for c in calls if c.variant_type == vtype]
        out[vtype.value] = _tally(subset, vtype.value)
    return out
