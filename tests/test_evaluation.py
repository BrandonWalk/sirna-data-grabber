from __future__ import annotations

import math

import pytest

from sirna_data.evaluation import evaluate_predictions

# --------------------------------------------------------------------------
# basic regression metrics
# --------------------------------------------------------------------------


def test_perfect_predictions():
    actual = [0.1, 0.5, 0.9, 0.3, 0.7]
    metrics = evaluate_predictions(actual, actual)
    assert metrics.n == 5
    assert metrics.pcc == pytest.approx(1.0)
    assert metrics.spcc == pytest.approx(1.0)
    assert metrics.mse == pytest.approx(0.0)
    assert metrics.rmse == pytest.approx(0.0)
    assert metrics.r2 == pytest.approx(1.0)


def test_mse_rmse_known_values():
    predicted = [0.0, 0.0, 0.0, 0.0]
    actual = [1.0, 1.0, 1.0, 1.0]
    metrics = evaluate_predictions(predicted, actual)
    assert metrics.mse == pytest.approx(1.0)
    assert metrics.rmse == pytest.approx(1.0)


def test_r2_undefined_when_actual_is_constant():
    predicted = [0.1, 0.2, 0.3]
    actual = [0.5, 0.5, 0.5]  # zero variance -> ss_tot == 0
    metrics = evaluate_predictions(predicted, actual)
    assert math.isnan(metrics.r2)


def test_pcc_spcc_perfect_negative_correlation():
    predicted = [0.1, 0.2, 0.3, 0.4]
    actual = [0.4, 0.3, 0.2, 0.1]
    metrics = evaluate_predictions(predicted, actual)
    assert metrics.pcc == pytest.approx(-1.0)
    assert metrics.spcc == pytest.approx(-1.0)


def test_pcc_undefined_with_single_point():
    metrics = evaluate_predictions([0.5], [0.6])
    assert math.isnan(metrics.pcc)
    assert math.isnan(metrics.spcc)


# --------------------------------------------------------------------------
# AUC
# --------------------------------------------------------------------------


def test_auc_perfect_separation_is_one():
    # every "hit" (actual >= 0.5) predicted-scores strictly above every non-hit
    predicted = [0.1, 0.2, 0.8, 0.9]
    actual = [0.1, 0.2, 0.8, 0.9]
    metrics = evaluate_predictions(predicted, actual, hit_threshold=0.5)
    assert metrics.auc == pytest.approx(1.0)


def test_auc_perfect_inversion_is_zero():
    # hits (actual >= 0.5, i.e. last two) get the LOWEST predicted scores
    predicted = [0.9, 0.8, 0.2, 0.1]
    actual = [0.1, 0.2, 0.8, 0.9]
    metrics = evaluate_predictions(predicted, actual, hit_threshold=0.5)
    assert metrics.auc == pytest.approx(0.0)


def test_auc_undefined_when_all_one_class():
    predicted = [0.1, 0.2, 0.3]
    actual = [0.1, 0.2, 0.3]  # all below hit_threshold -> no positives
    metrics = evaluate_predictions(predicted, actual, hit_threshold=0.9)
    assert math.isnan(metrics.auc)


def test_auc_ties_give_point_five_when_uninformative():
    # predicted scores identical across the board -> ranks tie -> AUC = 0.5
    predicted = [0.5, 0.5, 0.5, 0.5]
    actual = [0.1, 0.2, 0.8, 0.9]
    metrics = evaluate_predictions(predicted, actual, hit_threshold=0.5)
    assert metrics.auc == pytest.approx(0.5)


# --------------------------------------------------------------------------
# F1
# --------------------------------------------------------------------------


def test_f1_perfect_prediction_is_one():
    predicted = [0.1, 0.2, 0.8, 0.9]
    actual = [0.1, 0.2, 0.8, 0.9]
    metrics = evaluate_predictions(predicted, actual, hit_threshold=0.5)
    assert metrics.f1 == pytest.approx(1.0)


def test_f1_no_overlap_is_zero():
    # true hits are indices 2,3; predicted hits are indices 0,1 -- no overlap
    predicted = [0.9, 0.9, 0.1, 0.1]
    actual = [0.1, 0.1, 0.9, 0.9]
    metrics = evaluate_predictions(predicted, actual, hit_threshold=0.5)
    assert metrics.f1 == pytest.approx(0.0)


def test_f1_no_positives_at_all_is_zero_not_nan():
    predicted = [0.1, 0.2, 0.3]
    actual = [0.1, 0.2, 0.3]
    metrics = evaluate_predictions(predicted, actual, hit_threshold=0.9)
    assert metrics.f1 == pytest.approx(0.0)
    assert not math.isnan(metrics.f1)


def test_f1_known_precision_recall():
    # true hits: idx 0,1,2 ; predicted hits: idx 0,1,3
    # TP=2 (0,1), FP=1 (3), FN=1 (2) -> precision=2/3, recall=2/3, F1=2/3
    actual = [0.9, 0.9, 0.9, 0.1]
    predicted = [0.9, 0.9, 0.1, 0.9]
    metrics = evaluate_predictions(predicted, actual, hit_threshold=0.5)
    assert metrics.f1 == pytest.approx(2 / 3)


# --------------------------------------------------------------------------
# per-gene breakdown and weighted averages
# --------------------------------------------------------------------------


def test_no_gene_breakdown_when_genes_omitted():
    metrics = evaluate_predictions([0.1, 0.2], [0.1, 0.2])
    assert metrics.by_gene == {}
    assert metrics.weighted_pcc_by_gene is None
    assert metrics.weighted_spcc_by_gene is None


def test_gene_breakdown_groups_correctly():
    predicted = [0.1, 0.2, 0.5, 0.6, 0.9]
    actual = [0.1, 0.3, 0.4, 0.7, 0.8]
    genes = ["A", "A", "B", "B", "B"]
    metrics = evaluate_predictions(predicted, actual, genes)
    assert set(metrics.by_gene) == {"A", "B"}
    assert metrics.by_gene["A"].n == 2
    assert metrics.by_gene["B"].n == 3


def test_gene_with_single_sirna_has_nan_correlation_but_is_recorded():
    predicted = [0.1, 0.2, 0.3]
    actual = [0.1, 0.2, 0.3]
    genes = ["A", "A", "SOLO"]
    metrics = evaluate_predictions(predicted, actual, genes)
    assert metrics.by_gene["SOLO"].n == 1
    assert math.isnan(metrics.by_gene["SOLO"].pcc)
    assert math.isnan(metrics.by_gene["SOLO"].spcc)


def test_weighted_average_matches_manual_calculation():
    # Gene A: perfect correlation (pcc=spcc=1.0), 2 sirnas
    # Gene B: perfect anti-correlation (pcc=spcc=-1.0), 3 sirnas
    predicted = [0.1, 0.2, 0.9, 0.5, 0.1]
    actual = [0.1, 0.2, 0.1, 0.5, 0.9]
    genes = ["A", "A", "B", "B", "B"]
    metrics = evaluate_predictions(predicted, actual, genes)

    assert metrics.by_gene["A"].pcc == pytest.approx(1.0)
    assert metrics.by_gene["B"].pcc == pytest.approx(-1.0)
    expected = (2 * 1.0 + 3 * -1.0) / 5
    assert metrics.weighted_pcc_by_gene == pytest.approx(expected)
    assert metrics.weighted_spcc_by_gene == pytest.approx(expected)


def test_weighted_average_excludes_undefined_genes_from_denominator():
    # SOLO has n=1 -> NaN correlation, must not count toward the weighted
    # average's denominator (would otherwise silently drag it down).
    predicted = [0.1, 0.2, 0.3, 0.9]
    actual = [0.1, 0.2, 0.3, 0.1]
    genes = ["A", "A", "A", "SOLO"]
    metrics = evaluate_predictions(predicted, actual, genes)
    assert metrics.by_gene["A"].pcc == pytest.approx(1.0)
    assert math.isnan(metrics.by_gene["SOLO"].pcc)
    # weighted average should be exactly gene A's value (weight 3), SOLO
    # (weight 1, but undefined) excluded entirely
    assert metrics.weighted_pcc_by_gene == pytest.approx(1.0)


# --------------------------------------------------------------------------
# input validation
# --------------------------------------------------------------------------


def test_rejects_mismatched_lengths():
    with pytest.raises(ValueError):
        evaluate_predictions([0.1, 0.2], [0.1])


def test_rejects_mismatched_genes_length():
    with pytest.raises(ValueError):
        evaluate_predictions([0.1, 0.2], [0.1, 0.2], genes=["A"])


def test_rejects_empty_input():
    with pytest.raises(ValueError):
        evaluate_predictions([], [])
