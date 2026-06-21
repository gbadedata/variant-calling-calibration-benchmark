"""Core data types for the variant-calling calibration benchmark.

A variant CALL (from the query VCF) is matched against GIAB truth and
assigned a MatchStatus. Each call also carries a QUAL score (the caller's
stated confidence) and a set of stratum memberships (which difficult-region
categories it falls in). These three things, truth status, confidence, and
genomic context, are everything the benchmark needs.
"""

from __future__ import annotations

from enum import Enum

from pydantic import BaseModel, Field


class VariantType(str, Enum):
    """SNV vs indel, the primary stratification axis."""

    SNV = "snv"
    INDEL = "indel"

    @staticmethod
    def from_alleles(ref: str, alt: str) -> "VariantType":
        """A call is an SNV iff ref and alt are both single bases."""
        if len(ref) == 1 and len(alt) == 1:
            return VariantType.SNV
        return VariantType.INDEL


class MatchStatus(str, Enum):
    """Concordance status of a call against the GIAB truth set.

    TP: call is present in truth (correct).
    FP: call is absent from truth (false positive).
    FN: truth variant the caller missed (false negative).

    FN entries originate from the truth set, not the query, so they carry no
    QUAL score. They matter for recall but are excluded from calibration,
    which is a property of the caller's emitted confidence.
    """

    TP = "tp"
    FP = "fp"
    FN = "fn"


class VariantCall(BaseModel):
    """A single variant, from the query VCF or (for FN) the truth set."""

    chrom: str
    pos: int
    ref: str
    alt: str
    qual: float | None = None  # None for FN (truth-only) entries
    variant_type: VariantType
    status: MatchStatus
    strata: list[str] = Field(default_factory=list)

    @property
    def key(self) -> str:
        return f"{self.chrom}:{self.pos}:{self.ref}:{self.alt}"

    @staticmethod
    def make(
        chrom: str, pos: int, ref: str, alt: str,
        status: MatchStatus, qual: float | None = None,
        strata: list[str] | None = None,
    ) -> "VariantCall":
        return VariantCall(
            chrom=chrom, pos=pos, ref=ref, alt=alt, qual=qual,
            variant_type=VariantType.from_alleles(ref, alt),
            status=status, strata=strata or [],
        )


class StratumResult(BaseModel):
    """Concordance metrics within one stratum (or overall)."""

    name: str
    tp: int
    fp: int
    fn: int

    @property
    def precision(self) -> float:
        d = self.tp + self.fp
        return self.tp / d if d else 0.0

    @property
    def recall(self) -> float:
        d = self.tp + self.fn
        return self.tp / d if d else 0.0

    @property
    def f1(self) -> float:
        p, r = self.precision, self.recall
        return 2 * p * r / (p + r) if (p + r) else 0.0

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "tp": self.tp, "fp": self.fp, "fn": self.fn,
            "precision": round(self.precision, 4),
            "recall": round(self.recall, 4),
            "f1": round(self.f1, 4),
        }


class CalibrationBin(BaseModel):
    """One QUAL-score bin for the calibration analysis.

    For calls whose QUAL falls in [lo, hi), we compare the mean stated
    confidence (derived from QUAL) against the empirical precision (the
    fraction that are true positives). A well-calibrated caller has these
    two close together in every bin.
    """

    lo: float
    hi: float
    n: int
    tp: int
    fp: int
    mean_confidence: float  # implied by QUAL, in [0, 1]

    @property
    def empirical_precision(self) -> float:
        d = self.tp + self.fp
        return self.tp / d if d else 0.0

    @property
    def gap(self) -> float:
        """Signed calibration gap: stated confidence minus empirical."""
        return self.mean_confidence - self.empirical_precision

    def to_dict(self) -> dict:
        return {
            "lo": round(self.lo, 2), "hi": round(self.hi, 2),
            "n": self.n, "tp": self.tp, "fp": self.fp,
            "mean_confidence": round(self.mean_confidence, 4),
            "empirical_precision": round(self.empirical_precision, 4),
            "gap": round(self.gap, 4),
        }
