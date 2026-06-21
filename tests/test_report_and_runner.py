"""Tests for the unified report."""

from __future__ import annotations

import json

from src.benchmark.report import build_report
from src.data_loader import load_synthetic


class TestBuildReport:
    def test_has_all_layers(self) -> None:
        d = build_report("c", load_synthetic()).to_dict()
        assert "concordance_overall" in d
        assert "stratified_by_region" in d
        assert "calibration" in d
        assert "filtering" in d

    def test_json_serialisable(self) -> None:
        json.dumps(build_report("c", load_synthetic()).to_dict())

    def test_writes_file(self, tmp_path) -> None:
        r = build_report("c", load_synthetic())
        out = tmp_path / "r.json"
        r.write(out)
        assert out.exists()
        loaded = json.loads(out.read_text())
        assert loaded["caller_name"] == "c"

    def test_fixture_overconfident_in_report(self) -> None:
        r = build_report("c", load_synthetic())
        # the fixture is overconfident, so ECE > 0 and gap positive
        assert r.calibration.ece > 0
        assert r.calibration.overall_mean_confidence > r.calibration.overall_empirical_precision
