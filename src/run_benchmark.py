"""Benchmark runner: load calls, run all layers, write the report.

Usage:
    python3 -m src.run_benchmark                      # synthetic, offline
    python3 -m src.run_benchmark --query q.vcf.gz --truth giab.vcf.gz
                                                      # live, on real VCFs

The synthetic path needs no network and is what CI runs. The live path parses
a query VCF and the GIAB truth VCF, matches them, and (optionally) intersects
with stratification BEDs supplied via --strata-bed name=path pairs.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import structlog

from config.settings import settings
from src.benchmark.report import build_report, print_summary
from src.data_loader import load_synthetic, match_against_truth, parse_vcf

structlog.configure(
    processors=[
        structlog.processors.add_log_level,
        structlog.processors.TimeStamper(fmt="iso"),
        structlog.processors.JSONRenderer(),
    ],
    wrapper_class=structlog.make_filtering_bound_logger(logging.INFO),
)
log = structlog.get_logger()


def main() -> None:
    parser = argparse.ArgumentParser(description="Variant-calling calibration benchmark")
    parser.add_argument("--query", type=str, default=None, help="query VCF (live mode)")
    parser.add_argument("--truth", type=str, default=None, help="GIAB truth VCF (live mode)")
    parser.add_argument("--happy-vcf", type=str, default=None,
                        help="hap.py annotated output VCF (recommended for real data)")
    parser.add_argument("--caller-name", type=str, default=None, help="label for the report")
    parser.add_argument("--bins", type=int, default=settings.n_calibration_bins)
    args = parser.parse_args()

    if args.happy_vcf:
        from src.happy_adapter import parse_happy_vcf
        log.info("happy_mode", happy_vcf=args.happy_vcf)
        calls = parse_happy_vcf(Path(args.happy_vcf))
        caller_name = args.caller_name or Path(args.happy_vcf).stem
    elif args.query and args.truth:
        log.info("live_mode", query=args.query, truth=args.truth)
        query = parse_vcf(Path(args.query))
        truth = parse_vcf(Path(args.truth))
        calls = match_against_truth(query, truth)
        caller_name = args.caller_name or Path(args.query).stem
    else:
        log.info("synthetic_mode")
        calls = load_synthetic()
        caller_name = args.caller_name or "synthetic-caller"

    report = build_report(caller_name, calls, n_bins=args.bins)
    out = settings.reports_dir / f"benchmark_{caller_name.replace('/', '_')}.json"
    report.write(out)
    print_summary(report)
    print(f"\nReport written to {out}")


if __name__ == "__main__":
    main()
