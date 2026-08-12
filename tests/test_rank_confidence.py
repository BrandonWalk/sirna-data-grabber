from __future__ import annotations

import math

import pytest

from sirna_data.rank_confidence import (
    _log_binomial_cdf,
    min_top_k_for_confidence,
    probability_true_top_in_predicted_top_k,
    spearman_to_pearson,
)

# --------------------------------------------------------------------------
# spearman_to_pearson
# --------------------------------------------------------------------------


def test_spearman_to_pearson_zero_is_zero():
    assert spearman_to_pearson(0.0) == pytest.approx(0.0)


def test_spearman_to_pearson_plus_one_is_plus_one():
    assert spearman_to_pearson(1.0) == pytest.approx(1.0)


def test_spearman_to_pearson_minus_one_is_minus_one():
    assert spearman_to_pearson(-1.0) == pytest.approx(-1.0)


def test_spearman_to_pearson_monotonically_increasing():
    values = [spearman_to_pearson(s) for s in (-1.0, -0.5, 0.0, 0.5, 1.0)]
    assert values == sorted(values)


@pytest.mark.parametrize("bad_spcc", [-1.5, 1.5, float("nan") + 2])
def test_spearman_to_pearson_rejects_out_of_range(bad_spcc):
    with pytest.raises(ValueError):
        spearman_to_pearson(bad_spcc)


# --------------------------------------------------------------------------
# _log_binomial_cdf: cross-check against a direct math.comb computation
# --------------------------------------------------------------------------


def _direct_binomial_cdf(k: int, n_trials: int, p: float) -> float:
    return sum(math.comb(n_trials, i) * p**i * (1 - p) ** (n_trials - i) for i in range(k + 1))


@pytest.mark.parametrize("n_trials,p", [(20, 0.3), (50, 0.05), (100, 0.5), (10, 0.9)])
def test_log_binomial_cdf_matches_direct_computation(n_trials, p):
    for k in (0, 1, n_trials // 4, n_trials // 2, n_trials - 1):
        expected = _direct_binomial_cdf(k, n_trials, p)
        actual = math.exp(_log_binomial_cdf(k, n_trials, p))
        assert actual == pytest.approx(expected, abs=1e-9)


def test_log_binomial_cdf_edge_cases():
    assert _log_binomial_cdf(-1, 10, 0.5) == -math.inf
    assert math.exp(_log_binomial_cdf(10, 10, 0.5)) == pytest.approx(1.0)
    assert math.exp(_log_binomial_cdf(5, 10, 0.0)) == pytest.approx(1.0)
    assert math.exp(_log_binomial_cdf(5, 10, 1.0)) == pytest.approx(0.0)


def test_log_binomial_cdf_handles_large_n_without_stalling_at_zero():
    # A plain (non-log-space) floating recurrence starting from P(X=0) can
    # underflow to a *stuck* zero here; this must not happen.
    n_trials, p = 50_000, 0.1
    mean = n_trials * p
    around_mean = math.exp(_log_binomial_cdf(int(mean), n_trials, p))
    assert 0.3 < around_mean < 0.7  # CDF at the mean should be near 0.5


# --------------------------------------------------------------------------
# probability_true_top_in_predicted_top_k
# --------------------------------------------------------------------------


def test_probability_zero_correlation_is_uniform_rank():
    # With spcc=0, the true top item's predicted score is exchangeable with
    # every other item's, so its predicted rank is exactly uniform on
    # 1..n_items -- P(rank <= k) == k / n_items, not just approximately.
    n_items = 200
    for k in (1, 10, 50, 100, 199):
        prob = probability_true_top_in_predicted_top_k(k, n_items, spcc=0.0)
        assert prob == pytest.approx(k / n_items, abs=1e-3)


def test_probability_is_monotone_nondecreasing_in_k():
    n_items = 100
    probs = [
        probability_true_top_in_predicted_top_k(k, n_items, spcc=0.4) for k in range(1, n_items)
    ]
    assert probs == sorted(probs)


def test_probability_at_k_equals_n_items_is_one():
    assert probability_true_top_in_predicted_top_k(50, 50, spcc=0.2) == pytest.approx(1.0)


def test_probability_increases_with_correlation():
    # Fixed K, n_items: higher SPCC should never make it less likely the
    # true top item lands in the predicted top K.
    n_items, k = 100, 10
    probs = [
        probability_true_top_in_predicted_top_k(k, n_items, spcc=s)
        for s in (0.0, 0.2, 0.4, 0.6, 0.8, 0.99)
    ]
    assert probs == sorted(probs)


# --------------------------------------------------------------------------
# min_top_k_for_confidence
# --------------------------------------------------------------------------


def test_min_top_k_single_item_is_always_one():
    assert min_top_k_for_confidence(spcc=-0.9, n_items=1, confidence=0.999) == 1
    assert min_top_k_for_confidence(spcc=0.9, n_items=1, confidence=0.01) == 1


def test_min_top_k_within_bounds():
    for spcc in (-0.5, 0.0, 0.3, 0.7, 0.99):
        k = min_top_k_for_confidence(spcc=spcc, n_items=500, confidence=0.9)
        assert 1 <= k <= 500


def test_min_top_k_result_actually_clears_the_requested_confidence():
    n_items, confidence = 300, 0.9
    for spcc in (0.1, 0.5, 0.9):
        k = min_top_k_for_confidence(spcc=spcc, n_items=n_items, confidence=confidence)
        assert probability_true_top_in_predicted_top_k(k, n_items, spcc=spcc) >= confidence
        # and one fewer should not (except at the k=1 boundary)
        if k > 1:
            assert (
                probability_true_top_in_predicted_top_k(k - 1, n_items, spcc=spcc) < confidence
            )


def test_min_top_k_decreases_as_correlation_improves():
    n_items, confidence = 400, 0.95
    ks = [
        min_top_k_for_confidence(spcc=s, n_items=n_items, confidence=confidence)
        for s in (0.0, 0.3, 0.6, 0.9, 0.99)
    ]
    assert ks == sorted(ks, reverse=True)


def test_min_top_k_increases_as_confidence_target_rises():
    n_items, spcc = 400, 0.5
    ks = [
        min_top_k_for_confidence(spcc=spcc, n_items=n_items, confidence=c)
        for c in (0.5, 0.7, 0.9, 0.95, 0.99)
    ]
    assert ks == sorted(ks)


def test_min_top_k_zero_correlation_matches_uniform_expectation():
    # spcc=0 -> predicted rank is exactly uniform on 1..n_items, so the
    # smallest K clearing `confidence` is ceil(confidence * n_items) --
    # allow +-1 since probability_true_top_in_predicted_top_k is a numerical
    # integration, not bit-exact, and confidence*n_items lands exactly on
    # the boundary here.
    n_items, confidence = 1000, 0.95
    k = min_top_k_for_confidence(spcc=0.0, n_items=n_items, confidence=confidence)
    assert abs(k - math.ceil(confidence * n_items)) <= 1


@pytest.mark.parametrize("bad_spcc", [-1.1, 1.1])
def test_min_top_k_rejects_bad_spcc(bad_spcc):
    with pytest.raises(ValueError):
        min_top_k_for_confidence(spcc=bad_spcc, n_items=100, confidence=0.9)


@pytest.mark.parametrize("bad_confidence", [0.0, 1.0, -0.1, 1.1])
def test_min_top_k_rejects_bad_confidence(bad_confidence):
    with pytest.raises(ValueError):
        min_top_k_for_confidence(spcc=0.5, n_items=100, confidence=bad_confidence)


def test_min_top_k_rejects_bad_n_items():
    with pytest.raises(ValueError):
        min_top_k_for_confidence(spcc=0.5, n_items=0, confidence=0.9)


# --------------------------------------------------------------------------
# spcc vs. pcc: mutual exclusivity, and pcc used directly (no conversion)
# --------------------------------------------------------------------------


def test_pcc_used_directly_matches_spcc_converted_to_same_r():
    # Passing spcc=s should give identical results to passing
    # pcc=spearman_to_pearson(s), since spcc is just converted to that same
    # r internally before anything else happens.
    n_items, confidence, spcc = 250, 0.9, 0.6
    r = spearman_to_pearson(spcc)

    k_via_spcc = min_top_k_for_confidence(n_items=n_items, confidence=confidence, spcc=spcc)
    k_via_pcc = min_top_k_for_confidence(n_items=n_items, confidence=confidence, pcc=r)
    assert k_via_spcc == k_via_pcc

    prob_via_spcc = probability_true_top_in_predicted_top_k(10, n_items, spcc=spcc)
    prob_via_pcc = probability_true_top_in_predicted_top_k(10, n_items, pcc=r)
    assert prob_via_spcc == pytest.approx(prob_via_pcc)


def test_pcc_is_not_converted_like_spcc():
    # pcc=0.5 should NOT behave like spcc=0.5 (which converts to a
    # different, larger r via the arcsin relation).
    n_items, k = 300, 20
    prob_pcc = probability_true_top_in_predicted_top_k(k, n_items, pcc=0.5)
    prob_spcc = probability_true_top_in_predicted_top_k(k, n_items, spcc=0.5)
    assert prob_pcc != pytest.approx(prob_spcc)


@pytest.mark.parametrize("bad_pcc", [-1.1, 1.1])
def test_min_top_k_rejects_bad_pcc(bad_pcc):
    with pytest.raises(ValueError):
        min_top_k_for_confidence(n_items=100, confidence=0.9, pcc=bad_pcc)


def test_min_top_k_requires_exactly_one_of_spcc_or_pcc():
    with pytest.raises(ValueError):
        min_top_k_for_confidence(n_items=100, confidence=0.9)  # neither given
    with pytest.raises(ValueError):
        min_top_k_for_confidence(n_items=100, confidence=0.9, spcc=0.5, pcc=0.5)  # both given


def test_probability_requires_exactly_one_of_spcc_or_pcc():
    with pytest.raises(ValueError):
        probability_true_top_in_predicted_top_k(10, 100)
    with pytest.raises(ValueError):
        probability_true_top_in_predicted_top_k(10, 100, spcc=0.5, pcc=0.5)


@pytest.mark.parametrize("bad_spcc", [-1.1, 1.1, float("inf"), float("-inf")])
def test_probability_rejects_bad_spcc(bad_spcc):
    with pytest.raises(ValueError):
        probability_true_top_in_predicted_top_k(10, 100, spcc=bad_spcc)


@pytest.mark.parametrize("bad_pcc", [-1.1, 1.1, float("inf"), float("-inf")])
def test_probability_rejects_bad_pcc(bad_pcc):
    with pytest.raises(ValueError):
        probability_true_top_in_predicted_top_k(10, 100, pcc=bad_pcc)


# --------------------------------------------------------------------------
# probability_true_top_in_predicted_top_k always returns a valid probability
# --------------------------------------------------------------------------


def test_probability_is_always_a_valid_probability():
    for n_items in (2, 5, 50, 1000):
        for k in (1, max(1, n_items // 2), n_items - 1, n_items):
            for corr_kwargs in (
                {"spcc": -0.99},
                {"spcc": 0.0},
                {"spcc": 0.5},
                {"spcc": 0.99},
                {"spcc": -1.0},
                {"spcc": 1.0},
                {"pcc": -0.99},
                {"pcc": 0.0},
                {"pcc": 0.5},
                {"pcc": 0.99},
                {"pcc": -1.0},
                {"pcc": 1.0},
            ):
                prob = probability_true_top_in_predicted_top_k(k, n_items, **corr_kwargs)
                assert not math.isnan(prob)
                assert 0.0 <= prob <= 1.0
