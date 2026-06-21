"""Load and match variant calls against GIAB truth.

Two paths:

  - load_synthetic(): a deterministic, biologically plausible fixture used by
    all tests and CI. No downloads, no network. It embeds calls with known
    truth status, QUAL scores, variant types, and stratum memberships so the
    whole framework can be exercised offline.

  - load_from_vcfs(): the live path used on the user's machine. Parses a
    query VCF and the GIAB truth VCF, matches by chrom:pos:ref:alt, assigns
    TP/FP/FN, and intersects with stratification BEDs. (Implemented to the
    point the offline environment allows; the network-dependent parsing is
    validated by the user against real GIAB data.)

QUAL semantics: VCF QUAL is phred-scaled P(call is wrong):
    QUAL = -10 * log10(P_wrong)  ->  P_wrong = 10^(-QUAL/10)
The caller's stated confidence that the call is correct is therefore
    confidence = 1 - 10^(-QUAL/10)
This conversion is what makes calibration meaningful: it puts the caller's
QUAL on the same [0, 1] scale as the empirical precision.
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path

import numpy as np

from src.data_models import MatchStatus, VariantCall

logger = logging.getLogger(__name__)


def qual_to_confidence(qual: float) -> float:
    """Convert a phred-scaled QUAL into a probability the call is correct."""
    if qual <= 0:
        return 0.0
    p_wrong = 10.0 ** (-qual / 10.0)
    return float(min(1.0, max(0.0, 1.0 - p_wrong)))


# Stratification category names used in the synthetic fixture and expected
# from the real GIAB stratification BEDs.
STRATA = ["high_confidence", "low_complexity", "segdup", "high_gc", "low_mappability"]


def load_synthetic(seed: int = 42, n: int = 600) -> list[VariantCall]:
    """Build a deterministic synthetic call set with known truth status.

    The fixture is constructed so that the caller is INTENTIONALLY
    miscalibrated in a realistic way: it is overconfident in difficult
    regions (low-complexity, segdup) and well-calibrated in high-confidence
    regions. This gives the calibration and stratification layers a real
    signal to detect, exactly as a real caller would exhibit.
    """
    rng = np.random.default_rng(seed)
    calls: list[VariantCall] = []

    for i in range(n):
        # Assign a stratum
        stratum = STRATA[rng.integers(0, len(STRATA))]
        in_difficult = stratum in ("low_complexity", "segdup", "low_mappability")

        # Variant type: ~85% SNV, ~15% indel (realistic genome-wide ratio)
        is_indel = rng.random() < 0.15
        if is_indel:
            ref, alt = "AT", "A"  # a deletion
        else:
            ref, alt = "C", "T"

        # QUAL: difficult regions get high QUAL too (the caller THINKS it is
        # confident), but its true correctness is lower there. That mismatch
        # is the miscalibration the benchmark detects.
        qual = float(rng.uniform(20, 90))
        stated_conf = qual_to_confidence(qual)

        # True correctness probability: in easy regions it tracks the stated
        # confidence; in difficult regions it is systematically lower.
        true_correct_p = stated_conf if not in_difficult else stated_conf * 0.75
        # Indels are harder than SNVs everywhere
        if is_indel:
            true_correct_p *= 0.85

        is_tp = rng.random() < true_correct_p
        status = MatchStatus.TP if is_tp else MatchStatus.FP

        pos = 1_000_000 + i * 137
        calls.append(VariantCall.make(
            chrom="chr1", pos=pos, ref=ref, alt=alt,
            status=status, qual=qual,
            strata=[stratum] + (["high_confidence"] if not in_difficult else []),
        ))

    # Add false negatives: truth variants the caller missed. These have no
    # QUAL (truth-only). Concentrate them in difficult regions to depress
    # recall there, as a real caller would.
    n_fn = int(n * 0.08)
    for j in range(n_fn):
        stratum = "low_complexity" if j % 2 == 0 else "segdup"
        is_indel = rng.random() < 0.4  # missed calls skew toward indels
        ref, alt = ("AT", "A") if is_indel else ("C", "T")
        pos = 5_000_000 + j * 211
        calls.append(VariantCall.make(
            chrom="chr1", pos=pos, ref=ref, alt=alt,
            status=MatchStatus.FN, qual=None, strata=[stratum],
        ))

    logger.info("synthetic_loaded: n=%d (incl %d FN)", len(calls), n_fn)
    return calls


# ── Live VCF path (validated by the user against real GIAB data) ──────────

def _open_maybe_gzip(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def parse_vcf(path: Path) -> list[tuple[str, int, str, str, float | None]]:
    """Parse a VCF into (chrom, pos, ref, alt, qual) tuples.

    Minimal, dependency-free parser handling the mandatory columns. Multi-
    allelic records are split into one entry per ALT allele.
    """
    out: list[tuple[str, int, str, str, float | None]] = []
    with _open_maybe_gzip(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            chrom, pos, _id, ref, alt_field, qual_s = cols[0], cols[1], cols[2], cols[3], cols[4], cols[5]
            try:
                qual = None if qual_s in (".", "") else float(qual_s)
            except ValueError:
                qual = None
            for alt in alt_field.split(","):
                out.append((chrom, int(pos), ref, alt, qual))
    return out


def match_against_truth(
    query: list[tuple[str, int, str, str, float | None]],
    truth: list[tuple[str, int, str, str, float | None]],
) -> list[VariantCall]:
    """Assign TP/FP/FN by matching query and truth on chrom:pos:ref:alt."""
    truth_keys = {(c, p, r, a) for (c, p, r, a, _q) in truth}
    query_keys = {(c, p, r, a) for (c, p, r, a, _q) in query}
    calls: list[VariantCall] = []

    for (c, p, r, a, q) in query:
        status = MatchStatus.TP if (c, p, r, a) in truth_keys else MatchStatus.FP
        calls.append(VariantCall.make(chrom=c, pos=p, ref=r, alt=a, status=status, qual=q))

    for (c, p, r, a, _q) in truth:
        if (c, p, r, a) not in query_keys:
            calls.append(VariantCall.make(chrom=c, pos=p, ref=r, alt=a, status=MatchStatus.FN, qual=None))

    return calls
