"""Command-line entry point: `sirna-rank-confidence` (registered via
[project.scripts] in pyproject.toml). Given a correlation between a
predicted and true ranking, prints how many top-predicted items need to be
checked to hit each requested confidence level -- see `sirna_data.rank_confidence`
for the underlying model and its caveats.

    sirna-rank-confidence --pcc 0.3686 --n-items 4561 \\
        --confidence 0.99 0.95 0.9 0.8 0.7 0.6 0.5

    sirna-rank-confidence --spcc 0.5 --n-items 1000
"""
from __future__ import annotations

import argparse

from .rank_confidence import min_top_k_for_confidence

DEFAULT_CONFIDENCE_LEVELS = [0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5]


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="sirna-rank-confidence",
        description=(
            "Given a correlation (SPCC or PCC) between a predicted and true ranking of "
            "N_ITEMS total items, print the smallest number of top-predicted items that "
            "need to be checked to be X%% confident the true best item is among them, "
            "for each requested confidence level."
        ),
    )
    correlation_group = parser.add_mutually_exclusive_group(required=True)
    correlation_group.add_argument(
        "--spcc",
        type=float,
        metavar="RHO",
        help="Spearman's Rank Correlation Coefficient, in [-1, 1]. Mutually exclusive with --pcc.",
    )
    correlation_group.add_argument(
        "--pcc",
        type=float,
        metavar="R",
        help=(
            "Pearson's Correlation Coefficient, in [-1, 1], used directly (no "
            "Spearman-to-Pearson conversion). Mutually exclusive with --spcc."
        ),
    )
    parser.add_argument(
        "--n-items",
        type=int,
        required=True,
        metavar="N",
        help="Total number of ranked items.",
    )
    parser.add_argument(
        "--confidence",
        type=float,
        nargs="+",
        default=DEFAULT_CONFIDENCE_LEVELS,
        metavar="LEVEL",
        help=(
            "One or more target confidence levels in (0, 1) (default: "
            f"{' '.join(str(c) for c in DEFAULT_CONFIDENCE_LEVELS)})."
        ),
    )
    args = parser.parse_args()

    if args.spcc is not None:
        correlation_kind, correlation_value = "spcc", args.spcc
    else:
        correlation_kind, correlation_value = "pcc", args.pcc

    print(f"{correlation_kind} = {correlation_value}, n_items = {args.n_items}\n")
    print(f"{'confidence':>10}  {'min top-K to check':>19}")
    for confidence in args.confidence:
        kwargs = {correlation_kind: correlation_value}
        k = min_top_k_for_confidence(n_items=args.n_items, confidence=confidence, **kwargs)
        print(f"{confidence:>10.4g}  {k:>19}")


if __name__ == "__main__":
    main()
