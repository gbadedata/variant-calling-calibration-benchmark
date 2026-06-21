"""Configuration for the variant-calling calibration benchmark.

All parameters are overridable via environment variables prefixed VCB_ or a
.env file. Defaults let the framework run offline against the bundled
synthetic fixture with no downloads.
"""

from __future__ import annotations

from pathlib import Path

from pydantic_settings import BaseSettings, SettingsConfigDict

_ROOT = Path(__file__).resolve().parent.parent


class Settings(BaseSettings):
    model_config = SettingsConfigDict(env_prefix="VCB_", env_file=".env", extra="ignore")

    data_dir: Path = _ROOT / "data"
    reports_dir: Path = _ROOT / "evidence" / "reports"
    figures_dir: Path = _ROOT / "evidence" / "figures"

    # GIAB truth source (used by the live loader on the user's machine)
    giab_truth_vcf: str = "HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
    query_vcf: str = "query.vcf.gz"

    # Calibration
    n_calibration_bins: int = 10
    qual_min: float = 0.0
    qual_max: float = 100.0

    # Abstention / filtering sweep
    filter_threshold_steps: int = 20

    random_seed: int = 42


settings = Settings()
settings.reports_dir.mkdir(parents=True, exist_ok=True)
settings.figures_dir.mkdir(parents=True, exist_ok=True)
