"""Plot probability-of-capture vs. number-of-tests curves for one or more
model PCCs, using `sirna_data.rank_confidence`'s model.

Requires the optional `plot` extra for matplotlib, which is NOT installed
by the core package (`pip install sirna-data-grabber` alone):

    pip install sirna-data-grabber[plot]

This is a separate module (rather than living in `rank_confidence.py`
itself, or being imported by `sirna_data/__init__.py`) specifically so that
`import sirna_data` and everything in `rank_confidence.py` keep working
without matplotlib installed -- only importing this module requires it.
"""
from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING

from .rank_confidence import _default_k_values, probability_curves_for_pccs

if TYPE_CHECKING:
    from matplotlib.axes import Axes

try:
    import matplotlib.pyplot as plt
except ImportError as e:  # pragma: no cover -- exercised via the [plot] extra, or its absence
    raise ImportError(
        "sirna_data.rank_confidence_plot requires matplotlib, which is not installed. "
        "Install it with: pip install sirna-data-grabber[plot]"
    ) from e


def plot_probability_vs_num_tests(
    pccs: Sequence[float],
    n_items: int,
    *,
    top_n: int = 1,
    k_values: Sequence[int] | None = None,
    ax: Axes | None = None,
    save_path: str | None = None,
) -> Axes:
    """Line plot, one curve per PCC in `pccs`, of P(at least one of the true
    top `top_n` items is in the predicted top K) against K (number of
    top-predicted items checked) -- x-axis is K, y-axis is that
    probability, for `n_items` total ranked items.

    Parameters
    ----------
    pccs : one Pearson correlation per model to compare, each in [-1, 1].
        Must be non-empty.
    n_items : total number of ranked items (>= 1), same for every curve.
    top_n : passed through to every underlying probability calculation
        (default 1, the single true best item) -- see
        `sirna_data.rank_confidence`'s module docstring for what this
        generalizes and its extra caveat.
    k_values : which K's to plot each curve at. Defaults to
        `_default_k_values(n_items)` (a modest, evenly-spaced set -- see
        there for why); pass an explicit sequence (e.g.
        `range(1, n_items + 1)`) for the literal full curve or a denser
        sweep.
    ax : an existing `matplotlib.axes.Axes` to draw on. If not given, a new
        figure and axes are created.
    save_path : if given, the figure is also written to this path (e.g.
        "curves.png") via `Axes.figure.savefig` before returning.

    Returns
    -------
    The `matplotlib.axes.Axes` the curves were drawn on.
    """
    resolved_k_values = list(k_values) if k_values is not None else _default_k_values(n_items)
    curves = probability_curves_for_pccs(pccs, n_items, resolved_k_values, top_n=top_n)

    if ax is None:
        _, ax = plt.subplots()

    for p in pccs:
        ax.plot(resolved_k_values, curves[p], marker="o", markersize=3, label=f"PCC = {p:g}")

    captured_label = "the true best item" if top_n == 1 else f"one of the true top {top_n} items"
    ax.set_xlabel("Number of top-predicted items checked (K)")
    ax.set_ylabel(f"P({captured_label} is in the predicted top K)")
    title = f"n_items = {n_items}"
    if top_n != 1:
        title += f", top_n = {top_n}"
    ax.set_title(title)
    ax.set_ylim(0.0, 1.02)
    ax.legend()
    ax.grid(True, alpha=0.3)

    if save_path is not None:
        figure = ax.figure
        # ax.figure is typed as Figure | SubFigure (only relevant to the
        # nested-subfigures API, which this module never uses) -- narrow it
        # so mypy knows savefig is available.
        assert isinstance(figure, plt.Figure)
        figure.savefig(save_path, bbox_inches="tight")

    return ax
