"""Estimate how many top-ranked predictions must be checked to be confident
at least one true top item is among them, given only a rank/linear
correlation between a predicted and a true ranking.

Motivating question: a predicted ranking (e.g. of candidate siRNAs/genes by
predicted potency) correlates with the true ranking. If you could only
afford to check the single top *predicted* item, how confident could you be
that it's also one of the true best items? Usually not very. This module
answers the practical follow-up: how many of the top *predicted* items (K)
do you need to check to be, say, 95% confident that at least one of the true
top `top_n` items is among them?

`top_n` defaults to 1 -- "is the single true best item in the predicted top
K" -- the original, narrower question this module started out answering.
Set `top_n` higher to ask the easier, often more practically useful
question "is at least one of the true top 10 (or whatever) items in the
predicted top K" -- naturally requiring a smaller K for the same
confidence, since there are more chances to be captured.

Every public function here takes the correlation as *either* of two mutually
exclusive keyword arguments:

- `spcc`: Spearman's Rank Correlation Coefficient, in [-1, 1]. Converted to
  the model's underlying Pearson `r` via Greiner's relation (step 1 below).
- `pcc`: Pearson's Correlation Coefficient, in [-1, 1], used directly as the
  model's `r` with no conversion. Use this if you already computed Pearson's
  r on your own (e.g. on the raw predicted/true scores rather than their
  ranks) -- it feeds the model more directly than round-tripping through
  Spearman's rho.

Exactly one of `spcc`/`pcc` must be given; passing both or neither raises
`ValueError`. `probability_curves_for_pccs` / `min_top_k_for_confidence_multi`
run either function once per PCC in a list, e.g. to compare several models'
Pearson correlations side by side -- see their docstrings, and
`sirna_data.rank_confidence_plot` for a function that turns
`probability_curves_for_pccs`' output into a plot (one curve per PCC).

Model (a Gaussian-copula approximation -- see the Caveats section below):

1. Treat the true and predicted "scores" underlying the two rankings as a
   bivariate normal pair (T, P) with standard normal marginals and Pearson
   correlation r. Spearman's rho relates to r via Greiner's relation for a
   bivariate normal copula:

       rho_s = (6 / pi) * arcsin(r / 2)   =>   r = 2 * sin(pi * rho_s / 6)

   This is exact for that copula; the approximation in the overall method
   is treating a real predicted-vs-true score relationship as if it came
   from this copula in the first place. See `spearman_to_pearson`.

2. For each of the true top `top_n` items (true rank i = 1, ..., top_n), its
   true score is approximated by the expected value of the i-th largest of
   n iid standard normal draws, via the same simple Blom-style plug-in
   already used for i=1 (the maximum): t_hat_i = Phi^-1((n-i+1) / (n+1)).
   This point-estimate simplification ignores the (typically minor, for
   moderate/large n) extra spread from that order statistic's own sampling
   variability around its expectation.

3. Conditional on that true score, item i's predicted score P_i is
   Normal(r * t_hat_i, 1 - r**2) (standard bivariate-normal conditioning).

4. The OTHER n-1 items' predicted scores are approximated as iid
   Normal(0, 1) (their marginal distribution -- this ignores the mild
   dependence induced by conditioning that their true scores are all below
   the max, which is this model's main source of bias; see Caveats). Given
   P_i = p, the number of them exceeding p is Binomial(n-1, 1 - Phi(p)), so
   item i's PREDICTED rank is 1 + that count.

5. P(item i's predicted rank <= K) is obtained by integrating step 4's
   binomial CDF over the distribution of P_i from step 3 (numerically, via
   Simpson's rule).

6. For `top_n > 1`, step 5 is repeated for each of the `top_n` true top
   items, and their "predicted rank <= K" events are treated as
   INDEPENDENT to combine into "at least one captured":
   P(at least one) = 1 - product_i(1 - P(item i's predicted rank <= K)).
   This independence assumption is an additional approximation on top of
   step 4's (see Caveats) -- in reality these events are positively
   correlated (all `top_n` items share the same background of "other"
   items), so this tends to slightly OVERstate the combined probability,
   partially offsetting step 4's conservative bias. `min_top_k_for_confidence`
   finds the smallest K clearing the requested confidence via binary search
   (the probability is monotone non-decreasing in K).

Caveats -- this is a planning heuristic, not a certified statistical bound:

- Step 4's marginal (rather than order-statistic-conditioned) treatment of
  the other n-1 items means this model tends to be CONSERVATIVE -- it
  generally recommends a K at least as large as, and often larger than,
  strictly necessary -- especially at high correlation and/or small n. For
  the "how many do I need to check to be safe" use case this is the right
  direction to be wrong in, but don't treat the returned K as tight.
- Step 6's independence assumption (only relevant for `top_n > 1`) pulls in
  the opposite direction (slightly overstates the combined probability), so
  the two biases partially cancel rather than stack -- but neither is
  designed to exactly offset the other, so still treat the result as
  approximate, more so as `top_n` grows.
- Accuracy improves as n_items grows; at very small n_items (say, under
  ~20) or confidence extremely close to 1, treat the output as a rough
  guide rather than precise.
- The two rankings are assumed to have no ties.
- Cost grows roughly linearly with `top_n` (each additional true top item
  is another numerical integration) -- keep `top_n` modest (tens, not
  thousands) for interactive use, especially inside `min_top_k_for_confidence`'s
  binary search or a wide `k_values` sweep in `probability_curves_for_pccs`.
"""
from __future__ import annotations

import math
from collections.abc import Sequence
from statistics import NormalDist

_NORMAL = NormalDist()


def spearman_to_pearson(spcc: float) -> float:
    """Convert a Spearman rank correlation to the Pearson correlation `r` of
    the underlying bivariate-normal copula, via Greiner's relation
    `rho_s = (6/pi) * arcsin(r/2)`. Exact for that copula; an approximation
    for real rank data (see module docstring). Result is clamped to
    [-1, 1] to absorb floating-point drift at the +-1 endpoints.
    """
    if not -1.0 <= spcc <= 1.0:
        raise ValueError(f"spcc must be in [-1, 1], got {spcc!r}")
    r = 2.0 * math.sin(math.pi * spcc / 6.0)
    return max(-1.0, min(1.0, r))


def _resolve_pearson_r(spcc: float | None, pcc: float | None) -> float:
    """Validate that exactly one of `spcc`/`pcc` was given and return the
    corresponding Pearson `r` for the model (converting `spcc` via
    `spearman_to_pearson`, or validating `pcc`'s range directly)."""
    if (spcc is None) == (pcc is None):
        raise ValueError(
            "Provide exactly one of `spcc` or `pcc`, not both or neither "
            f"(got spcc={spcc!r}, pcc={pcc!r})."
        )
    if pcc is not None:
        if not -1.0 <= pcc <= 1.0:
            raise ValueError(f"pcc must be in [-1, 1], got {pcc!r}")
        return pcc
    assert spcc is not None  # for type-checkers; guaranteed by the check above
    return spearman_to_pearson(spcc)


def _expected_order_statistic_of_n_standard_normal(rank: int, n: int) -> float:
    """Blom-style plug-in approximation for E[the `rank`-th largest of n iid
    N(0, 1)] -- `rank=1` is the maximum (E[max] = Phi^-1(n/(n+1)))."""
    return _NORMAL.inv_cdf((n - rank + 1) / (n + 1))


def _log_binomial_cdf(k: int, n_trials: int, p: float) -> float:
    """log P(X <= k) for X ~ Binomial(n_trials, p).

    Computed via the standard term-by-term recurrence
    P(X=i+1) = P(X=i) * (n_trials-i)/(i+1) * p/(1-p), but accumulated in
    log-space with a running log-sum-exp. Doing this in log-space (rather
    than exact integer math.comb, or a plain floating-point recurrence)
    avoids two failure modes that matter for realistically large n_items:
    math.comb's exact-integer cost blowing up, and the plain recurrence's
    leading term silently underflowing to a *stuck* 0.0 for large n_trials
    (which would then incorrectly report 0 probability mass even in the
    distribution's real support region).
    """
    if k < 0:
        return -math.inf
    if k >= n_trials:
        return 0.0  # log(1)
    if p <= 0.0:
        return 0.0  # log(1): X is always 0 <= k
    if p >= 1.0:
        return -math.inf  # log(0): X is always n_trials > k

    log_p, log_1mp = math.log(p), math.log1p(-p)
    log_term = n_trials * log_1mp  # log P(X = 0)
    log_cdf = log_term
    for i in range(1, k + 1):
        log_term += math.log(n_trials - i + 1) - math.log(i) + log_p - log_1mp
        hi, lo = (log_cdf, log_term) if log_cdf >= log_term else (log_term, log_cdf)
        log_cdf = hi + math.log1p(math.exp(lo - hi))
    return log_cdf


def _probability_single_item_rank_leq_k(
    k: int, background_size: int, mean_p: float, sd_p: float
) -> float:
    """P(a single item's predicted rank is <= k), where the item's predicted
    score P is Normal(mean_p, sd_p**2) and its predicted rank among
    `background_size` OTHER items (each with predicted score approximated
    as iid marginal N(0, 1), independent of this item -- see module
    docstring step 4) is 1 + Binomial(background_size, 1 - Phi(P)).

    Internal helper shared by every top_n=1..N order statistic in
    `_probability_any_of_top_n_in_top_k_given_r`, and by
    `min_top_k_for_confidence`'s binary search (indirectly, via that
    function), so `spcc` isn't reconverted to `r` on every iteration.
    """
    if k > background_size:
        return 1.0
    if k < 1:
        return 0.0

    if sd_p == 0.0:
        # r = +-1: P is deterministic given mean_p, no integral needed.
        # Clamp like the general branch below -- log-sum-exp accumulation
        # in _log_binomial_cdf can overshoot 1.0 by a floating-point hair
        # (e.g. 1.0000000000000004) since it should mathematically sum to
        # exactly <= 1 but isn't guaranteed to in finite precision.
        q = 1.0 - _NORMAL.cdf(mean_p)
        return min(1.0, max(0.0, math.exp(_log_binomial_cdf(k - 1, background_size, q))))

    # Integrate P(rank <= k | P = p) * density(p) over p, reparameterized as
    # p = mean_p + sd_p * z with z ~ N(0, 1), via Simpson's rule. phi(z) is
    # negligible outside [-8, 8] to far beyond float precision.
    z_lo, z_hi, num_intervals = -8.0, 8.0, 400  # even, required for Simpson's rule
    h = (z_hi - z_lo) / num_intervals

    def integrand(z: float) -> float:
        p = mean_p + sd_p * z
        q = 1.0 - _NORMAL.cdf(p)
        return _NORMAL.pdf(z) * math.exp(_log_binomial_cdf(k - 1, background_size, q))

    total = integrand(z_lo) + integrand(z_hi)
    for i in range(1, num_intervals):
        total += (4.0 if i % 2 else 2.0) * integrand(z_lo + i * h)
    return min(1.0, max(0.0, total * h / 3.0))


def _probability_any_of_top_n_in_top_k_given_r(
    k: int, n_items: int, top_n: int, r: float
) -> float:
    """P(at least one of the true top `top_n` items' predicted rank is <=
    k), given the already-converted Pearson correlation `r` (see module
    docstring for the full model, including the extra independence
    assumption this introduces for `top_n > 1`). Internal helper shared by
    the two public functions below.
    """
    if n_items < 1:
        raise ValueError(f"n_items must be >= 1, got {n_items!r}")
    if not 1 <= top_n <= n_items:
        raise ValueError(f"top_n must be in [1, n_items] (n_items={n_items!r}), got {top_n!r}")
    if n_items == 1:
        return 1.0
    if top_n == n_items:
        # every item is one of the "true top n_items" by definition -- a
        # nonempty predicted top K trivially contains at least one of them.
        return 1.0 if k >= 1 else 0.0
    if k >= n_items:
        return 1.0
    if k < 1:
        return 0.0

    sd_p = math.sqrt(max(0.0, 1.0 - r * r))
    background_size = n_items - 1

    prob_none_captured = 1.0
    for rank in range(1, top_n + 1):
        t_hat = _expected_order_statistic_of_n_standard_normal(rank, n_items)
        mean_p = r * t_hat
        p_rank = _probability_single_item_rank_leq_k(k, background_size, mean_p, sd_p)
        prob_none_captured *= 1.0 - p_rank
    return min(1.0, max(0.0, 1.0 - prob_none_captured))


def probability_true_top_in_predicted_top_k(
    k: int,
    n_items: int,
    *,
    top_n: int = 1,
    spcc: float | None = None,
    pcc: float | None = None,
) -> float:
    """P(at least one of the true top `top_n` items' predicted rank is <=
    k), under this module's model (see module docstring) for n_items total
    items correlated at either Spearman's rho (`spcc`) or Pearson's r
    (`pcc`) -- exactly one of the two must be given, each validated to be
    in [-1, 1] (via `_resolve_pearson_r`, which raises `ValueError`
    otherwise).

    `top_n` (default 1, the single true best item) must be in
    [1, n_items] -- see the module docstring for how it generalizes the
    model and the extra caveat it introduces.

    The returned probability is itself checked to be a valid probability
    (in [0, 1], not NaN) before returning -- `_probability_any_of_top_n_in_top_k_given_r`
    already clamps its numerical-integration result into [0, 1], so this
    should never actually fire, but it's cheap insurance against a future
    change to the model quietly returning something nonsensical.
    """
    r = _resolve_pearson_r(spcc, pcc)
    prob = _probability_any_of_top_n_in_top_k_given_r(k, n_items, top_n, r)
    if math.isnan(prob) or not 0.0 <= prob <= 1.0:
        raise RuntimeError(
            f"Computed probability {prob!r} is not a valid probability (must be in [0, 1]) "
            "-- this indicates a bug in the underlying model, not bad input; please report it."
        )
    return prob


def min_top_k_for_confidence(
    n_items: int,
    confidence: float,
    *,
    top_n: int = 1,
    spcc: float | None = None,
    pcc: float | None = None,
) -> int:
    """How many of the top *predicted* items need to be checked to be at
    least `confidence` confident that at least one of the true top `top_n`
    items is among them, given only a correlation between the predicted and
    true rankings of `n_items` items total.

    Parameters
    ----------
    n_items : total number of ranked items (>= 1).
    confidence : target probability, in (0, 1) -- e.g. 0.95 for 95%.
    top_n : how many of the true best items count as "captured" if any one
        of them shows up in the predicted top K (default 1, i.e. only the
        single true best item). Must be in [1, n_items]. See the module
        docstring for how this generalizes the model and its extra
        caveat, and note cost grows roughly linearly with `top_n`.
    spcc : Spearman's Rank Correlation Coefficient between the predicted and
        true rankings, in [-1, 1]. Mutually exclusive with `pcc` -- give
        exactly one.
    pcc : Pearson's Correlation Coefficient, in [-1, 1], used directly as
        the model's underlying `r` (no Spearman-to-Pearson conversion).
        Mutually exclusive with `spcc` -- give exactly one.

    Returns
    -------
    The smallest K in [1, n_items] such that, under this module's model (see
    module docstring), P(at least one true top-`top_n` item is within the
    predicted top K) is at least `confidence`. Returns 1 if even the single
    top predicted item already clears the threshold, and n_items if nothing
    short of checking every item does (e.g. a non-positive correlation, or
    confidence very close to 1 with a small n_items).

    This is a modeled estimate, not an exact combinatorial guarantee -- see
    the module docstring's Caveats section (it's generally conservative, so
    treat the result as a "check at least this many" planning figure).
    """
    if n_items < 1:
        raise ValueError(f"n_items must be >= 1, got {n_items!r}")
    if not 0.0 < confidence < 1.0:
        raise ValueError(f"confidence must be in (0, 1), got {confidence!r}")
    if not 1 <= top_n <= n_items:
        raise ValueError(f"top_n must be in [1, n_items] (n_items={n_items!r}), got {top_n!r}")

    r = _resolve_pearson_r(spcc, pcc)  # validate spcc/pcc even if n_items == 1

    if n_items == 1:
        return 1

    lo, hi = 1, n_items
    while lo < hi:
        mid = (lo + hi) // 2
        prob = _probability_any_of_top_n_in_top_k_given_r(mid, n_items, top_n, r)
        if prob >= confidence:
            hi = mid
        else:
            lo = mid + 1
    return lo


def _default_k_values(n_items: int, num_points: int = 60) -> list[int]:
    """A modest, roughly evenly-spaced set of K values across [1, n_items]
    (deduplicated, always including both endpoints), used as the default
    `k_values` for `probability_curves_for_pccs` / the plotting module.
    Evaluating every single K from 1 to n_items is usually unnecessary for
    a probability curve and, since each evaluation's cost scales with K
    (see `_log_binomial_cdf`), can be slow for large `n_items` -- pass an
    explicit `k_values` (e.g. `range(1, n_items + 1)`) for the literal full
    curve or any other specific set of K's instead.
    """
    if n_items < 1:
        raise ValueError(f"n_items must be >= 1, got {n_items!r}")
    if n_items <= num_points:
        return list(range(1, n_items + 1))
    step = (n_items - 1) / (num_points - 1)
    values = sorted({max(1, min(n_items, round(1 + i * step))) for i in range(num_points)})
    return values


def probability_curves_for_pccs(
    pccs: Sequence[float],
    n_items: int,
    k_values: Sequence[int] | None = None,
    *,
    top_n: int = 1,
) -> dict[float, list[float]]:
    """`probability_true_top_in_predicted_top_k`, evaluated once per PCC in
    `pccs` at every K in `k_values` -- lets you compare several models'
    Pearson correlations side by side (e.g. to see how the number of tests
    needed relates to each model's PCC), and is what
    `sirna_data.rank_confidence_plot.plot_probability_vs_num_tests` plots.

    Parameters
    ----------
    pccs : one Pearson correlation per model to compare, each in [-1, 1]
        (validated by the underlying probability calls; a `ValueError` from
        any one of them propagates immediately). Must be non-empty.
    n_items : total number of ranked items (>= 1), same for every PCC.
    k_values : which K's to evaluate the probability at, same for every
        PCC. Defaults to `_default_k_values(n_items)` -- see there for why
        (and how to override for a denser/custom sweep). Must be
        non-empty if given explicitly.
    top_n : passed through to `probability_true_top_in_predicted_top_k` for
        every evaluation (default 1, the single true best item).

    Returns
    -------
    `{pcc: [probability at each K in k_values, in order]}`, one entry per
    input PCC (in insertion order; duplicate PCCs in `pccs` simply overwrite
    -- pass a distinguishable set if you need per-model results kept
    separate but share the same PCC value).
    """
    if not pccs:
        raise ValueError("pccs must be non-empty")
    resolved_k_values = list(k_values) if k_values is not None else _default_k_values(n_items)
    if not resolved_k_values:
        raise ValueError("k_values must be non-empty if given explicitly")

    return {
        p: [
            probability_true_top_in_predicted_top_k(k, n_items, top_n=top_n, pcc=p)
            for k in resolved_k_values
        ]
        for p in pccs
    }


def min_top_k_for_confidence_multi(
    pccs: Sequence[float],
    n_items: int,
    confidence: float,
    *,
    top_n: int = 1,
) -> dict[float, int]:
    """`min_top_k_for_confidence`, evaluated once per PCC in `pccs` -- e.g.
    to see how many top-predicted items each of several models needs
    checked to hit the same confidence target.

    Returns `{pcc: min_top_k}`, one entry per input PCC (in insertion
    order; duplicate PCCs in `pccs` simply overwrite). `pccs` must be
    non-empty; every other parameter is passed straight through to
    `min_top_k_for_confidence` for each PCC.
    """
    if not pccs:
        raise ValueError("pccs must be non-empty")
    return {
        p: min_top_k_for_confidence(n_items, confidence, top_n=top_n, pcc=p) for p in pccs
    }
