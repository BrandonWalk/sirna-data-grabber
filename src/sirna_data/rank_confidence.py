"""Estimate how many top-ranked predictions must be checked to be confident
the true best item is among them, given only a rank/linear correlation
between a predicted and a true ranking.

Motivating question: a predicted ranking (e.g. of candidate siRNAs/genes by
predicted potency) correlates with the true ranking. If you could only
afford to check the single top *predicted* item, how confident could you be
that it's also the true rank-1 item? Usually not very. This module answers
the practical follow-up: how many of the top *predicted* items (K) do you
need to check to be, say, 95% confident the true best item is among them?

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
`ValueError`.

Model (a Gaussian-copula approximation -- see the Caveats section below):

1. Treat the true and predicted "scores" underlying the two rankings as a
   bivariate normal pair (T, P) with standard normal marginals and Pearson
   correlation r. Spearman's rho relates to r via Greiner's relation for a
   bivariate normal copula:

       rho_s = (6 / pi) * arcsin(r / 2)   =>   r = 2 * sin(pi * rho_s / 6)

   This is exact for that copula; the approximation in the overall method
   is treating a real predicted-vs-true score relationship as if it came
   from this copula in the first place. See `spearman_to_pearson`.

2. The true rank-1 item has the largest of n iid true scores. Its true
   score is approximated by the expected value of that maximum, via Blom's
   plug-in: t_hat = Phi^-1(1 - 1/(n+1)) -- a standard, simple approximation
   for the mean of a normal maximum. This point-estimate simplification
   ignores the (typically minor, for moderate/large n) extra spread from
   the true maximum's own sampling variability around that expectation.

3. Conditional on that true score, the item's predicted score P is
   Normal(r * t_hat, 1 - r**2) (standard bivariate-normal conditioning).

4. The OTHER n-1 items' predicted scores are approximated as iid
   Normal(0, 1) (their marginal distribution -- this ignores the mild
   dependence induced by conditioning that their true scores are all below
   the max, which is this model's main source of bias; see Caveats). Given
   P = p, the number of them exceeding p is Binomial(n-1, 1 - Phi(p)), so
   the true rank-1 item's PREDICTED rank is 1 + that count.

5. P(predicted rank <= K) is obtained by integrating step 4's binomial CDF
   over the distribution of P from step 3 (numerically, via Simpson's
   rule), and the smallest K clearing the requested confidence is found by
   binary search (the probability is monotone non-decreasing in K).

Caveats -- this is a planning heuristic, not a certified statistical bound:

- Step 4's marginal (rather than order-statistic-conditioned) treatment of
  the other n-1 items means this model tends to be CONSERVATIVE -- it
  generally recommends a K at least as large as, and often larger than,
  strictly necessary -- especially at high correlation and/or small n. For
  the "how many do I need to check to be safe" use case this is the right
  direction to be wrong in, but don't treat the returned K as tight.
- Accuracy improves as n_items grows; at very small n_items (say, under
  ~20) or confidence extremely close to 1, treat the output as a rough
  guide rather than precise.
- The two rankings are assumed to have no ties.
"""
from __future__ import annotations

import math
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


def _expected_max_of_n_standard_normal(n: int) -> float:
    """Blom's plug-in approximation for E[max of n iid N(0, 1)]."""
    return _NORMAL.inv_cdf(n / (n + 1))


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


def _probability_top_in_top_k_given_r(k: int, n_items: int, r: float) -> float:
    """P(the true rank-1 item's predicted rank is <= k), given the
    already-converted Pearson correlation `r` (see module docstring for the
    model). Internal helper shared by the two public functions below, so
    `min_top_k_for_confidence`'s binary search doesn't reconvert `spcc` to
    `r` on every iteration.
    """
    if n_items < 1:
        raise ValueError(f"n_items must be >= 1, got {n_items!r}")
    if n_items == 1:
        return 1.0
    if k >= n_items:
        return 1.0
    if k < 1:
        return 0.0

    t_hat = _expected_max_of_n_standard_normal(n_items)
    mean_p = r * t_hat
    sd_p = math.sqrt(max(0.0, 1.0 - r * r))

    if sd_p == 0.0:
        # r = +-1: P is deterministic given t_hat, no integral needed.
        q = 1.0 - _NORMAL.cdf(mean_p)
        return math.exp(_log_binomial_cdf(k - 1, n_items - 1, q))

    # Integrate P(rank <= k | P = p) * density(p) over p, reparameterized as
    # p = mean_p + sd_p * z with z ~ N(0, 1), via Simpson's rule. phi(z) is
    # negligible outside [-8, 8] to far beyond float precision.
    z_lo, z_hi, num_intervals = -8.0, 8.0, 400  # even, required for Simpson's rule
    h = (z_hi - z_lo) / num_intervals

    def integrand(z: float) -> float:
        p = mean_p + sd_p * z
        q = 1.0 - _NORMAL.cdf(p)
        return _NORMAL.pdf(z) * math.exp(_log_binomial_cdf(k - 1, n_items - 1, q))

    total = integrand(z_lo) + integrand(z_hi)
    for i in range(1, num_intervals):
        total += (4.0 if i % 2 else 2.0) * integrand(z_lo + i * h)
    return min(1.0, max(0.0, total * h / 3.0))


def probability_true_top_in_predicted_top_k(
    k: int,
    n_items: int,
    *,
    spcc: float | None = None,
    pcc: float | None = None,
) -> float:
    """P(the true rank-1 item's predicted rank is <= k), under this module's
    model (see module docstring) for n_items total items correlated at
    either Spearman's rho (`spcc`) or Pearson's r (`pcc`) -- exactly one of
    the two must be given.
    """
    r = _resolve_pearson_r(spcc, pcc)
    return _probability_top_in_top_k_given_r(k, n_items, r)


def min_top_k_for_confidence(
    n_items: int,
    confidence: float,
    *,
    spcc: float | None = None,
    pcc: float | None = None,
) -> int:
    """How many of the top *predicted* items need to be checked to be at
    least `confidence` confident the true best (true rank-1) item is among
    them, given only a correlation between the predicted and true rankings
    of `n_items` items total.

    Parameters
    ----------
    n_items : total number of ranked items (>= 1).
    confidence : target probability, in (0, 1) -- e.g. 0.95 for 95%.
    spcc : Spearman's Rank Correlation Coefficient between the predicted and
        true rankings, in [-1, 1]. Mutually exclusive with `pcc` -- give
        exactly one.
    pcc : Pearson's Correlation Coefficient, in [-1, 1], used directly as
        the model's underlying `r` (no Spearman-to-Pearson conversion).
        Mutually exclusive with `spcc` -- give exactly one.

    Returns
    -------
    The smallest K in [1, n_items] such that, under this module's model (see
    module docstring), P(true rank-1 item is within the predicted top K) is
    at least `confidence`. Returns 1 if even the single top predicted item
    already clears the threshold, and n_items if nothing short of checking
    every item does (e.g. a non-positive correlation, or confidence very
    close to 1 with a small n_items).

    This is a modeled estimate, not an exact combinatorial guarantee -- see
    the module docstring's Caveats section (it's generally conservative, so
    treat the result as a "check at least this many" planning figure).
    """
    if n_items < 1:
        raise ValueError(f"n_items must be >= 1, got {n_items!r}")
    if not 0.0 < confidence < 1.0:
        raise ValueError(f"confidence must be in (0, 1), got {confidence!r}")

    r = _resolve_pearson_r(spcc, pcc)  # validate spcc/pcc even if n_items == 1

    if n_items == 1:
        return 1

    lo, hi = 1, n_items
    while lo < hi:
        mid = (lo + hi) // 2
        prob = _probability_top_in_top_k_given_r(mid, n_items, r)
        if prob >= confidence:
            hi = mid
        else:
            lo = mid + 1
    return lo
