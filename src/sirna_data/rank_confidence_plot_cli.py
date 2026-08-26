"""Command-line entry point: `sirna-rank-confidence-plot` (registered via
[project.scripts] in pyproject.toml). Given one or more Pearson correlations
between a predicted and true ranking, plots P(the true top `--top-n` item(s)
are captured within the top K predicted) against K, one curve per `--pcc` --
see `sirna_data.rank_confidence_plot.plot_probability_vs_num_tests` for the
underlying model, and `sirna_data.rank_confidence_cli`'s `sirna-rank-
confidence` for the single-value, table-only (no matplotlib needed)
counterpart to this command.

Requires the optional `plot` extra for matplotlib, same as the module this
wraps:

    pip install sirna-data-grabber[plot]

Examples:

    # Compare three models' PCCs, top-5 capture, over a PPP2CA-length
    # transcript's ~4596 candidate sites, but only plot K up to 800:
    sirna-rank-confidence-plot --pcc 0.3521 0.3686 0.3927 \\
        --n-items 4596 --top-n 5 --k-max 800 \\
        --labels "B-Graph" "OligoFormer" "Best blend (40% B-Graph)" \\
        --save-path bgraph_vs_oligoformer_ppp2ca.png

    # Single curve, full K range, default 60-point spacing:
    sirna-rank-confidence-plot --pcc 0.3686 --n-items 4596 --save-path oligoformer.png

    # Plain lines, no per-point dots -- handy when --num-points is dense:
    sirna-rank-confidence-plot --pcc 0.3686 --n-items 4596 --marker none \\
        --save-path oligoformer_lines.png
"""
from __future__ import annotations

import argparse

import matplotlib.pyplot as plt

from .rank_confidence_plot import plot_probability_vs_num_tests


def _k_values_for_args(args: argparse.Namespace) -> list[int]:
    """Resolves --k-values/--k-max/--num-points into the actual list of K's
    to evaluate/plot. --k-values (if given) wins outright. Otherwise builds
    an evenly-spaced --num-points set across [1, k_max] (k_max defaults to
    --n-items, i.e. the full range, unless --k-max narrows it) -- same
    spacing scheme as rank_confidence._default_k_values, just bounded to
    k_max instead of always n_items, so "--k-max 800" genuinely limits the
    plotted x-axis rather than just its label.
    """
    if args.k_values is not None:
        return args.k_values

    k_max = args.k_max if args.k_max is not None else args.n_items
    if k_max < 1:
        raise argparse.ArgumentTypeError(f"--k-max must be >= 1, got {k_max}")
    if k_max > args.n_items:
        raise argparse.ArgumentTypeError(
            f"--k-max ({k_max}) cannot exceed --n-items ({args.n_items})"
        )

    if k_max <= args.num_points:
        return list(range(1, k_max + 1))
    step = (k_max - 1) / (args.num_points - 1)
    return sorted({max(1, min(k_max, round(1 + i * step))) for i in range(args.num_points)})


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="sirna-rank-confidence-plot",
        description=(
            "Plot P(true top --top-n item(s) captured in predicted top K) vs. K, one "
            "curve per --pcc, for a ranking of --n-items total items."
        ),
    )
    parser.add_argument(
        "--pcc",
        type=float,
        nargs="+",
        required=True,
        metavar="R",
        help="One or more Pearson correlations, each in [-1, 1] -- one curve per value.",
    )
    parser.add_argument(
        "--n-items",
        type=int,
        required=True,
        metavar="N",
        help="Total number of ranked items (e.g. candidate siRNA sites scanned across a "
        "transcript).",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=1,
        metavar="N",
        help="How many of the true top items count as 'captured' (default: 1, the single true "
        "best item).",
    )
    parser.add_argument(
        "--k-max",
        type=int,
        default=None,
        metavar="K",
        help="Limit the plotted K range to [1, K_MAX] instead of the full [1, N_ITEMS] -- e.g. "
        "'--k-max 800' to zoom in on the first 800 candidate sites. Default: no limit, the "
        "full [1, N_ITEMS] range.",
    )
    parser.add_argument(
        "--num-points",
        type=int,
        default=60,
        metavar="N",
        help="How many K values to evaluate/plot across the K range (default: 60, evenly "
        "spaced, always including both endpoints). Ignored if --k-values is given.",
    )
    parser.add_argument(
        "--k-values",
        type=int,
        nargs="+",
        default=None,
        metavar="K",
        help="Explicit K values to evaluate/plot, overriding --k-max/--num-points entirely.",
    )
    parser.add_argument(
        "--labels",
        type=str,
        nargs="+",
        default=None,
        metavar="LABEL",
        help="Legend label per --pcc value, same order/count (default: 'PCC = <value>' per curve, "
        "from the underlying plotting function).",
    )
    parser.add_argument(
        "--marker",
        type=str,
        default="o",
        metavar="MARKER",
        help="Per-point marker style (matplotlib marker character, e.g. 'o', 's', '^'). "
        "Pass 'none' for plain lines with no per-point markers (default: 'o').",
    )
    parser.add_argument(
        "--marker-size",
        type=float,
        default=3,
        metavar="SIZE",
        help="Marker size (default: 3). Ignored when --marker is 'none'.",
    )
    parser.add_argument(
        "--linestyle",
        type=str,
        default="-",
        metavar="STYLE",
        help="Line style (matplotlib linestyle string, e.g. '-', '--', ':', '-.', or 'none' for "
        "markers with no connecting line; default: '-', a solid line).",
    )
    parser.add_argument(
        "--save-path",
        type=str,
        required=True,
        metavar="PATH",
        help="Where to write the chart (e.g. chart.png).",
    )
    args = parser.parse_args()

    if args.labels is not None and len(args.labels) != len(args.pcc):
        parser.error(
            f"--labels has {len(args.labels)} value(s) but --pcc has {len(args.pcc)} -- "
            "must match 1:1"
        )

    k_values = _k_values_for_args(args)

    # "none" (any case) means "no marker"/"no connecting line" -- matplotlib
    # itself wants that as an actual None for marker, but as the literal
    # string "None" for linestyle, so these aren't handled identically.
    marker = None if args.marker.strip().lower() == "none" else args.marker
    linestyle = "None" if args.linestyle.strip().lower() == "none" else args.linestyle

    ax = plot_probability_vs_num_tests(
        args.pcc,
        args.n_items,
        top_n=args.top_n,
        k_values=k_values,
        marker=marker,
        markersize=args.marker_size,
        linestyle=linestyle,
    )

    if args.labels is not None:
        for line, label in zip(ax.get_lines(), args.labels, strict=True):
            line.set_label(label)
        ax.legend()

    figure = ax.figure
    # ax.figure is typed as Figure | SubFigure (only relevant to the nested-
    # subfigures API, which this module never uses) -- narrow it so mypy
    # knows savefig is available, same pattern as rank_confidence_plot.py.
    assert isinstance(figure, plt.Figure)
    figure.savefig(args.save_path, bbox_inches="tight")
    print(
        f"Saved chart to {args.save_path}  "
        f"(K range: [{k_values[0]}, {k_values[-1]}], {len(k_values)} points)"
    )


if __name__ == "__main__":
    main()
