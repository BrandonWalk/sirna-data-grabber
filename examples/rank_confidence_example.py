"""Example: how many top-predicted items need checking to be confident the
true best item is among them, given a Pearson correlation (PCC) between a
predicted and true ranking.

    python examples/rank_confidence_example.py
"""
from __future__ import annotations

from sirna_data import min_top_k_for_confidence

PCC = 0.3686
N_ITEMS = 4561
CONFIDENCE_LEVELS = [0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5]


def main() -> None:
    print(f"PCC = {PCC}, n_items = {N_ITEMS}\n")
    print(f"{'confidence':>10}  {'min top-K to check':>19}")
    for confidence in CONFIDENCE_LEVELS:
        k = min_top_k_for_confidence(n_items=N_ITEMS, confidence=confidence, pcc=PCC)
        print(f"{confidence:>10.2f}  {k:>19}")


if __name__ == "__main__":
    main()
