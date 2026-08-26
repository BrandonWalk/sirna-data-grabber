"""Tests for sirna_data.rank_confidence_plot_cli (`sirna-rank-confidence-plot`)
-- skipped entirely if matplotlib (the optional `plot` extra) isn't
installed, same as test_rank_confidence_plot.py."""
from __future__ import annotations

import argparse
from pathlib import Path

import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")  # headless backend -- no display needed to run tests

import sirna_data.rank_confidence_plot_cli as cli_module  # noqa: E402
from sirna_data.rank_confidence_plot_cli import main  # noqa: E402


def _spy_on_plot_call(monkeypatch: pytest.MonkeyPatch) -> dict:
    """Wraps the real plot_probability_vs_num_tests so tests can inspect
    exactly what kwargs main() resolved and passed through to it, while
    still exercising the real matplotlib call underneath."""
    real_plot = cli_module.plot_probability_vs_num_tests
    captured: dict = {}

    def spy(*args, **kwargs):
        captured["args"] = args
        captured["kwargs"] = kwargs
        ax = real_plot(*args, **kwargs)
        captured["ax"] = ax
        return ax

    monkeypatch.setattr(cli_module, "plot_probability_vs_num_tests", spy)
    return captured


def test_writes_a_png_file(monkeypatch: pytest.MonkeyPatch, tmp_path: Path, capsys):
    out_path = tmp_path / "curves.png"
    monkeypatch.setattr(
        "sys.argv",
        [
            "sirna-rank-confidence-plot",
            "--pcc", "0.3686",
            "--n-items", "1000",
            "--save-path", str(out_path),
        ],
    )
    main()
    assert out_path.exists()
    assert out_path.stat().st_size > 0
    assert str(out_path) in capsys.readouterr().out


def test_default_marker_is_a_dot(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    captured = _spy_on_plot_call(monkeypatch)
    monkeypatch.setattr(
        "sys.argv",
        [
            "sirna-rank-confidence-plot",
            "--pcc", "0.3",
            "--n-items", "200",
            "--save-path", str(tmp_path / "out.png"),
        ],
    )
    main()
    assert captured["kwargs"]["marker"] == "o"
    assert captured["ax"].get_lines()[0].get_marker() == "o"


def test_marker_none_disables_per_point_markers(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    captured = _spy_on_plot_call(monkeypatch)
    monkeypatch.setattr(
        "sys.argv",
        [
            "sirna-rank-confidence-plot",
            "--pcc", "0.3",
            "--n-items", "200",
            "--marker", "none",
            "--save-path", str(tmp_path / "out.png"),
        ],
    )
    main()
    assert captured["kwargs"]["marker"] is None
    assert captured["ax"].get_lines()[0].get_marker() == "None"


def test_marker_none_is_case_insensitive(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    captured = _spy_on_plot_call(monkeypatch)
    monkeypatch.setattr(
        "sys.argv",
        [
            "sirna-rank-confidence-plot",
            "--pcc", "0.3",
            "--n-items", "200",
            "--marker", "None",
            "--save-path", str(tmp_path / "out.png"),
        ],
    )
    main()
    assert captured["kwargs"]["marker"] is None


def test_custom_marker_passed_through(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    captured = _spy_on_plot_call(monkeypatch)
    monkeypatch.setattr(
        "sys.argv",
        [
            "sirna-rank-confidence-plot",
            "--pcc", "0.3",
            "--n-items", "200",
            "--marker", "s",
            "--save-path", str(tmp_path / "out.png"),
        ],
    )
    main()
    assert captured["kwargs"]["marker"] == "s"


def test_custom_marker_size_passed_through(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    captured = _spy_on_plot_call(monkeypatch)
    monkeypatch.setattr(
        "sys.argv",
        [
            "sirna-rank-confidence-plot",
            "--pcc", "0.3",
            "--n-items", "200",
            "--marker-size", "9",
            "--save-path", str(tmp_path / "out.png"),
        ],
    )
    main()
    assert captured["kwargs"]["markersize"] == 9


def test_default_linestyle_is_solid(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    captured = _spy_on_plot_call(monkeypatch)
    monkeypatch.setattr(
        "sys.argv",
        [
            "sirna-rank-confidence-plot",
            "--pcc", "0.3",
            "--n-items", "200",
            "--save-path", str(tmp_path / "out.png"),
        ],
    )
    main()
    assert captured["kwargs"]["linestyle"] == "-"


def test_linestyle_none_disables_connecting_line(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    captured = _spy_on_plot_call(monkeypatch)
    monkeypatch.setattr(
        "sys.argv",
        [
            "sirna-rank-confidence-plot",
            "--pcc", "0.3",
            "--n-items", "200",
            "--linestyle", "none",
            "--save-path", str(tmp_path / "out.png"),
        ],
    )
    main()
    assert captured["kwargs"]["linestyle"] == "None"


def test_custom_linestyle_passed_through(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    captured = _spy_on_plot_call(monkeypatch)
    monkeypatch.setattr(
        "sys.argv",
        [
            "sirna-rank-confidence-plot",
            "--pcc", "0.3",
            "--n-items", "200",
            # "--" alone would be parsed by argparse as the end-of-options
            # marker, not an option value -- ":" exercises the same
            # passthrough behavior without that collision.
            "--linestyle", ":",
            "--save-path", str(tmp_path / "out.png"),
        ],
    )
    main()
    assert captured["kwargs"]["linestyle"] == ":"


def test_multiple_pccs_produce_multiple_curves(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    captured = _spy_on_plot_call(monkeypatch)
    monkeypatch.setattr(
        "sys.argv",
        [
            "sirna-rank-confidence-plot",
            "--pcc", "0.2", "0.5", "0.8",
            "--n-items", "500",
            "--save-path", str(tmp_path / "out.png"),
        ],
    )
    main()
    assert len(captured["ax"].get_lines()) == 3


def test_labels_applied_to_legend(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    captured = _spy_on_plot_call(monkeypatch)
    monkeypatch.setattr(
        "sys.argv",
        [
            "sirna-rank-confidence-plot",
            "--pcc", "0.2", "0.5",
            "--n-items", "500",
            "--labels", "Model A", "Model B",
            "--save-path", str(tmp_path / "out.png"),
        ],
    )
    main()
    legend_texts = [t.get_text() for t in captured["ax"].get_legend().get_texts()]
    assert legend_texts == ["Model A", "Model B"]


def test_labels_count_mismatch_errors(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    monkeypatch.setattr(
        "sys.argv",
        [
            "sirna-rank-confidence-plot",
            "--pcc", "0.2", "0.5",
            "--n-items", "500",
            "--labels", "Only One Label",
            "--save-path", str(tmp_path / "out.png"),
        ],
    )
    with pytest.raises(SystemExit):
        main()


def test_requires_pcc(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    monkeypatch.setattr(
        "sys.argv",
        [
            "sirna-rank-confidence-plot",
            "--n-items", "500",
            "--save-path", str(tmp_path / "out.png"),
        ],
    )
    with pytest.raises(SystemExit):
        main()


def test_requires_n_items(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    monkeypatch.setattr(
        "sys.argv",
        [
            "sirna-rank-confidence-plot",
            "--pcc", "0.3",
            "--save-path", str(tmp_path / "out.png"),
        ],
    )
    with pytest.raises(SystemExit):
        main()


def test_requires_save_path(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(
        "sys.argv", ["sirna-rank-confidence-plot", "--pcc", "0.3", "--n-items", "500"]
    )
    with pytest.raises(SystemExit):
        main()


def test_k_max_cannot_exceed_n_items(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    monkeypatch.setattr(
        "sys.argv",
        [
            "sirna-rank-confidence-plot",
            "--pcc", "0.3",
            "--n-items", "100",
            "--k-max", "200",
            "--save-path", str(tmp_path / "out.png"),
        ],
    )
    with pytest.raises(argparse.ArgumentTypeError):
        main()


def test_explicit_k_values_override_k_max(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    captured = _spy_on_plot_call(monkeypatch)
    monkeypatch.setattr(
        "sys.argv",
        [
            "sirna-rank-confidence-plot",
            "--pcc", "0.3",
            "--n-items", "500",
            "--k-values", "1", "50", "100",
            "--save-path", str(tmp_path / "out.png"),
        ],
    )
    main()
    assert list(captured["kwargs"]["k_values"]) == [1, 50, 100]
