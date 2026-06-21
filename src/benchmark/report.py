"""Unified benchmark report: combines all four evaluation layers.

Brings together concordance, stratification, calibration, and filtering into
one result object and serialises to a single JSON artefact plus a human
readable summary. This is the one place that answers, end to end: how
accurate is the caller, where does it struggle, is its confidence honest, and
where should it stop trusting itself?
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path

from src.benchmark.calibration import CalibrationReport, compute_calibration
from src.benchmark.concordance import concordance_by_type, concordance_overall
from src.benchmark.filtering import FilterReport, sweep_filter_thresholds
from src.benchmark.stratify import stratify_by_qual_bin, stratify_by_region
from src.data_models import StratumResult, VariantCall

logger = logging.getLogger(__name__)


@dataclass
class BenchmarkReport:
    caller_name: str
    n_calls: int
    overall: StratumResult
    by_type: dict[str, StratumResult]
    by_region: dict[str, StratumResult]
    by_qual_bin: list[StratumResult]
    calibration: CalibrationReport
    filtering: FilterReport

    def to_dict(self) -> dict:
        return {
            "caller_name": self.caller_name,
            "n_calls": self.n_calls,
            "concordance_overall": self.overall.to_dict(),
            "concordance_by_type": {k: v.to_dict() for k, v in self.by_type.items()},
            "stratified_by_region": {k: v.to_dict() for k, v in self.by_region.items()},
            "stratified_by_qual_bin": [b.to_dict() for b in self.by_qual_bin],
            "calibration": self.calibration.to_dict(),
            "filtering": self.filtering.to_dict(),
        }

    def write(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)
        logger.info("report_written: %s", path)
        return path


def build_report(
    caller_name: str, calls: list[VariantCall], n_bins: int = 10,
) -> BenchmarkReport:
    """Run all four evaluation layers and assemble the report."""
    return BenchmarkReport(
        caller_name=caller_name,
        n_calls=len(calls),
        overall=concordance_overall(calls),
        by_type=concordance_by_type(calls),
        by_region=stratify_by_region(calls),
        by_qual_bin=stratify_by_qual_bin(calls, n_bins=n_bins),
        calibration=compute_calibration(calls, n_bins=n_bins),
        filtering=sweep_filter_thresholds(calls),
    )


def print_summary(report: BenchmarkReport) -> None:
    """Human-readable summary to stdout."""
    o = report.overall
    print("\n" + "=" * 68)
    print("VARIANT-CALLING CALIBRATION BENCHMARK")
    print("=" * 68)
    print(f"  Caller:    {report.caller_name}")
    print(f"  Calls:     {report.n_calls}")
    print()
    print("  Concordance vs GIAB truth")
    print(f"    Overall:   P={o.precision:.4f}  R={o.recall:.4f}  F1={o.f1:.4f}")
    for vtype, r in report.by_type.items():
        print(f"    {vtype:8s}   P={r.precision:.4f}  R={r.recall:.4f}  F1={r.f1:.4f}  (tp={r.tp})")
    print()
    print("  By region (precision / recall / F1)")
    for name, r in report.by_region.items():
        print(f"    {name:18s} {r.precision:.3f} / {r.recall:.3f} / {r.f1:.3f}")
    print()
    c = report.calibration
    print("  Calibration (is the caller's confidence honest?)")
    print(f"    Expected Calibration Error (ECE): {c.ece:.4f}")
    print(f"    Mean stated confidence:           {c.overall_mean_confidence:.4f}")
    print(f"    Empirical precision:              {c.overall_empirical_precision:.4f}")
    gap = c.overall_mean_confidence - c.overall_empirical_precision
    verdict = "OVERCONFIDENT" if gap > 0.02 else ("UNDERCONFIDENT" if gap < -0.02 else "WELL CALIBRATED")
    print(f"    Verdict:                          {verdict} (gap {gap:+.4f})")
    print()
    f = report.filtering
    print("  Filtering as abstention (safe vs decisive)")
    print(f"    Best F1:        {f.max_f1:.4f} at QUAL >= {f.max_f1_threshold:.1f}")
    if f.safe_threshold is not None:
        print(f"    Safe threshold: QUAL >= {f.safe_threshold:.1f} reaches precision {f.target_precision}")
    else:
        print(f"    Safe threshold: precision {f.target_precision} not reachable at any QUAL")
    print("=" * 68)
