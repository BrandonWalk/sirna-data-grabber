"""Tests for sirna_data.rank_confidence_plot -- skipped entirely if
matplotlib (the optional `plot` extra) isn't installed."""
from __future__ import annotations

from pathlib import Path

import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")  # headless backend -- no display needed to run tests

from sirna_data.rank_confidence import _default_k_values, probability_curves_for_pccs  # noqa: E402
from sirna_data.rank_confidence_plot import plot_probability_vs_num_tests  # noqa: E402


def test_plot_returns_axes_with_one_line_per_pcc():
    pccs = [0.1, 0.4, 0.8]
    ax = plot_probability_vs_num_tests(pccs, 300)
    assert len(ax.get_lines()) == len(pccs)


def test_plot_line_data_matches_probability_curves_for_pccs():
    pccs = [0.2, 0.6]
    n_items = 250
    k_values = [1, 50, 100, 200, 250]
    ax = plot_probability_vs_num_tests(pccs, n_items, k_values=k_values)
    expected = probability_curves_for_pccs(pccs, n_items, k_values)

    for line, p in zip(ax.get_lines(), pccs, strict=True):
        xdata = list(line.get_xdata())
        ydata = list(line.get_ydata())
        assert xdata == list(k_values)
        assert ydata == pytest.approx(expected[p])


def test_plot_uses_default_k_values_when_not_given():
    n_items = 800  # > _default_k_values' num_points=60, so the "large n_items" branch fires
    ax = plot_probability_vs_num_tests([0.3], n_items)
    xdata = list(ax.get_lines()[0].get_xdata())
    assert xdata == _default_k_values(n_items)


def test_plot_axis_labels_mention_top_n_when_greater_than_one():
    ax_default = plot_probability_vs_num_tests([0.3], 200, k_values=[10])
    ax_top5 = plot_probability_vs_num_tests([0.3], 200, k_values=[10], top_n=5)
    assert "the true best item" in ax_default.get_ylabel()
    assert "top 5" in ax_top5.get_ylabel()


def test_plot_legend_labels_include_pcc_values():
    ax = plot_probability_vs_num_tests([0.1, 0.5], 200, k_values=[10])
    legend_texts = [t.get_text() for t in ax.get_legend().get_texts()]
    assert any("0.1" in t for t in legend_texts)
    assert any("0.5" in t for t in legend_texts)


def test_plot_draws_on_existing_axes_when_given():
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    returned = plot_probability_vs_num_tests([0.3], 200, k_values=[10], ax=ax)
    assert returned is ax
    plt.close(fig)


def test_plot_save_path_writes_a_file(tmp_path: Path):
    out_path = tmp_path / "curves.png"
    plot_probability_vs_num_tests([0.2, 0.7], 300, save_path=str(out_path))
    assert out_path.exists()
    assert out_path.stat().st_size > 0
