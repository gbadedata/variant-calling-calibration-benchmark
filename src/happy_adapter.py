"""Adapter: ingest hap.py output into the framework's VariantCall model.

The framework's own exact matcher (data_loader.match_against_truth) is fine
for synthetic fixtures and normalised VCFs. For real benchmarking against
GIAB, the field standard is hap.py (Krusche et al. 2019), which performs
sophisticated variant normalisation and haplotype-aware matching, and handles
the messy realities (multi-allelic representation, complex variants, the MHC
region) that a naive matcher cannot.

The right division of labour is therefore:
  - hap.py does the matching (the hard, solved problem).
  - this framework adds the calibration and abstention analysis hap.py lacks.

This adapter reads hap.py's annotated output VCF and converts each query call
into a VariantCall with its truth status (from hap.py's BD/BVT FORMAT fields)
and its QUAL, so compute_calibration() and sweep_filter_thresholds() run on
real, properly matched data unchanged.

hap.py annotated VCF FORMAT fields used:
  BD  : decision  -> TP / FP / FN / N  (per sample: TRUTH and QUERY)
  BVT : variant type -> SNP / INDEL / NOCALL
The QUERY sample column carries the query call's decision; the TRUTH sample
column carries FN (truth present, query missed). QUAL is taken from the
query VCF's QUAL column, preserved by hap.py in the annotated output.
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path

from src.data_models import MatchStatus, VariantCall, VariantType

logger = logging.getLogger(__name__)


def _open(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _decision_to_status(bd: str) -> MatchStatus | None:
    bd = bd.upper()
    if bd == "TP":
        return MatchStatus.TP
    if bd == "FP":
        return MatchStatus.FP
    if bd == "FN":
        return MatchStatus.FN
    return None  # "N" / "." / NOCALL: not a scored call


def parse_happy_vcf(path: Path) -> list[VariantCall]:
    """Parse a hap.py annotated VCF (`*.vcf.gz`) into VariantCalls.

    hap.py writes two sample columns, TRUTH and QUERY, each with a BD
    (benchmarking decision) subfield. A query false positive is BD=FP in the
    QUERY sample; a true positive is BD=TP in both; a false negative is BD=FN
    in the TRUTH sample with the QUERY sample uncalled.
    """
    calls: list[VariantCall] = []
    sample_names: list[str] = []

    with _open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                sample_names = cols[9:]  # TRUTH, QUERY
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 10:
                continue
            chrom, pos, _id, ref, alt_field, qual_s, _filt, _info, fmt = cols[:9]
            samples = cols[9:]
            fmt_keys = fmt.split(":")
            try:
                bd_idx = fmt_keys.index("BD")
            except ValueError:
                continue

            # Map sample name -> its fields
            sample_fields = dict(zip(sample_names, samples))

            try:
                qual = None if qual_s in (".", "") else float(qual_s)
            except ValueError:
                qual = None

            for alt in alt_field.split(","):
                vtype = VariantType.from_alleles(ref, alt)
                # Determine status: prefer QUERY decision; fall back to TRUTH for FN
                status = None
                query = sample_fields.get("QUERY", "")
                truth = sample_fields.get("TRUTH", "")
                for sample in (query, truth):
                    parts = sample.split(":")
                    if len(parts) > bd_idx:
                        st = _decision_to_status(parts[bd_idx])
                        if st is not None:
                            status = st
                            break
                if status is None:
                    continue
                # FN entries have no meaningful query QUAL
                call_qual = qual if status in (MatchStatus.TP, MatchStatus.FP) else None
                calls.append(VariantCall(
                    chrom=chrom, pos=int(pos), ref=ref, alt=alt, qual=call_qual,
                    variant_type=vtype, status=status, strata=[],
                ))

    logger.info("happy_vcf_parsed: %d scored calls from %s", len(calls), path)
    return calls
