"""Evaluate a set of predicted knockdown (KD) values against actual/measured
KD values: overall regression + hit-calling metrics, and (if genes are
given) a per-gene correlation breakdown plus gene-weighted averages.

    from sirna_data.evaluation import evaluate_predictions

    metrics = evaluate_predictions(predicted_kd, actual_kd, genes)
    metrics.pcc, metrics.spcc, metrics.mse, metrics.rmse, metrics.r2
    metrics.auc, metrics.f1
    metrics.by_gene["TP53"].pcc          # only set if `genes` was given
    metrics.weighted_pcc_by_gene         # only set if `genes` was given

AUC and F1 need binary hit/non-hit labels, not continuous KD values, so both
are derived from `hit_threshold` (default 0.7, override to match whatever
scale/convention your KD values use -- e.g. a 0-1 fraction vs. a 0-100
percent scale): a sample counts as a "hit" when its KD is >= hit_threshold.
- `actual_kd >= hit_threshold` gives the TRUE hit/non-hit labels, used by
  both AUC and F1.
- `predicted_kd >= hit_threshold` gives the PREDICTED hit/non-hit labels,
  used by F1 only -- AUC instead uses `predicted_kd` directly as a
  continuous ranking score against the true labels (the standard,
  threshold-free use of AUC), so the choice of hit_threshold only affects
  AUC through the true-label side.
"""
from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass, field

import pandas as pd


@dataclass
class GeneCorrelation:
    """Per-gene correlation breakdown entry (see `PredictionMetrics.by_gene`)."""

    n: int  # number of siRNAs tested against this gene
    pcc: float  # NaN if undefined (n < 2, or zero variance in either series)
    spcc: float  # NaN if undefined (n < 2, or zero variance in either series)


@dataclass
class PredictionMetrics:
    """Result of `evaluate_predictions()`. `by_gene`, `weighted_pcc_by_gene`,
    and `weighted_spcc_by_gene` are only populated when `genes` was given.
    """

    n: int
    pcc: float
    spcc: float
    mse: float
    rmse: float
    r2: float
    auc: float
    f1: float
    hit_threshold: float
    by_gene: dict[str, GeneCorrelation] = field(default_factory=dict)
    weighted_pcc_by_gene: float | None = None
    weighted_spcc_by_gene: float | None = None


def _pearson(x: Sequence[float], y: Sequence[float]) -> float:
    """Pearson correlation, NaN if undefined (n < 2 or either series is
    constant) -- pandas' own NaN behavior for these cases, used directly."""
    return float(pd.Series(x, dtype=float).corr(pd.Series(y, dtype=float), method="pearson"))


def _spearman(x: Sequence[float], y: Sequence[float]) -> float:
    """Spearman correlation, NaN if undefined (n < 2 or either series is
    constant), with ties handled via average ranks.

    Computed as the Pearson correlation of x's and y's ranks (average rank
    for ties) rather than via pandas' `.corr(method="spearman")`, which
    would pull in scipy as a transitive import -- this package intentionally
    has no dependencies beyond pandas and requests (see rank_confidence.py's
    similar scipy-avoidance). Ranking itself uses only pandas/numpy.
    """
    x_ranks = pd.Series(x, dtype=float).rank(method="average")
    y_ranks = pd.Series(y, dtype=float).rank(method="average")
    return _pearson(x_ranks.tolist(), y_ranks.tolist())


def _auc(scores: Sequence[float], positive: Sequence[bool]) -> float:
    """ROC AUC of `scores` (continuous) against `positive` (binary true
    labels), via the Mann-Whitney U / rank-sum identity
    `AUC = (R_pos - n_pos*(n_pos+1)/2) / (n_pos * n_neg)`, where R_pos is the
    sum of `scores`' ranks (average rank for ties) restricted to the
    positive-labelled samples. Equivalent to sklearn's `roc_auc_score` for
    the binary case, without adding sklearn as a dependency. NaN if there
    are no positives or no negatives (AUC is undefined without both).
    """
    n_pos = sum(positive)
    n_neg = len(positive) - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    ranks = pd.Series(scores, dtype=float).rank(method="average")
    r_pos = float(ranks[list(positive)].sum())
    u_stat = r_pos - n_pos * (n_pos + 1) / 2
    return u_stat / (n_pos * n_neg)


def _f1(true_positive_labels: Sequence[bool], pred_positive_labels: Sequence[bool]) -> float:
    """Binary F1 of predicted vs. true hit/non-hit labels. Follows the
    common zero-division convention (matching scikit-learn's default):
    precision and/or recall are treated as 0.0 when their denominator is 0
    (no predicted positives / no actual positives respectively), so F1 is
    0.0 whenever there's nothing correctly predicted as positive -- never
    NaN.
    """
    tp = sum(t and p for t, p in zip(true_positive_labels, pred_positive_labels, strict=True))
    fp = sum(
        (not t) and p for t, p in zip(true_positive_labels, pred_positive_labels, strict=True)
    )
    fn = sum(
        t and (not p) for t, p in zip(true_positive_labels, pred_positive_labels, strict=True)
    )
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    if precision + recall == 0.0:
        return 0.0
    return 2 * precision * recall / (precision + recall)


def evaluate_predictions(
    predicted_kd: Sequence[float],
    actual_kd: Sequence[float],
    genes: Sequence[str] | None = None,
    *,
    hit_threshold: float = 0.7,
) -> PredictionMetrics:
    """Evaluate `predicted_kd` against `actual_kd`.

    Parameters
    ----------
    predicted_kd : predicted knockdown value per siRNA.
    actual_kd : measured/actual knockdown value per siRNA, same length and
        order as `predicted_kd`.
    genes : optional gene each siRNA was tested on, same length/order as the
        two KD lists. If given, adds a per-gene PCC/SPCC breakdown
        (`by_gene`) and gene-weighted average PCC/SPCC (weighted by the
        number of siRNAs tested per gene).
    hit_threshold : KD value at/above which a sample counts as a "hit", used
        to derive the binary labels AUC and F1 need. Must be in the same
        units/scale as `predicted_kd`/`actual_kd` (default 0.7 assumes a
        0-1 fractional scale; pass e.g. 70.0 if your KD values are a 0-100
        percent scale instead). See the module docstring for exactly how
        it's applied to each metric.

    Returns
    -------
    A `PredictionMetrics` with the overall metrics always populated, and
    `by_gene`/`weighted_pcc_by_gene`/`weighted_spcc_by_gene` populated only
    when `genes` is given. Any individual metric that's undefined for the
    given data (e.g. PCC with fewer than 2 points, or AUC with no negative
    examples) is NaN rather than raising -- check with `math.isnan(...)`.
    """
    n = len(predicted_kd)
    if len(actual_kd) != n:
        raise ValueError(
            f"predicted_kd and actual_kd must be the same length, got {n} and {len(actual_kd)}"
        )
    if n == 0:
        raise ValueError("predicted_kd/actual_kd must not be empty")
    if genes is not None and len(genes) != n:
        raise ValueError(f"genes must be the same length as predicted_kd, got {len(genes)}")

    predicted_kd = [float(v) for v in predicted_kd]
    actual_kd = [float(v) for v in actual_kd]

    mse = sum((a - p) ** 2 for a, p in zip(actual_kd, predicted_kd, strict=True)) / n
    rmse = math.sqrt(mse)
    mean_actual = sum(actual_kd) / n
    ss_tot = sum((a - mean_actual) ** 2 for a in actual_kd)
    r2 = float("nan") if ss_tot == 0.0 else 1.0 - (mse * n) / ss_tot

    true_hits = [a >= hit_threshold for a in actual_kd]
    pred_hits = [p >= hit_threshold for p in predicted_kd]

    metrics = PredictionMetrics(
        n=n,
        pcc=_pearson(predicted_kd, actual_kd),
        spcc=_spearman(predicted_kd, actual_kd),
        mse=mse,
        rmse=rmse,
        r2=r2,
        auc=_auc(predicted_kd, true_hits),
        f1=_f1(true_hits, pred_hits),
        hit_threshold=hit_threshold,
    )

    if genes is None:
        return metrics

    indices_by_gene: dict[str, list[int]] = defaultdict(list)
    for i, gene in enumerate(genes):
        indices_by_gene[gene].append(i)

    by_gene: dict[str, GeneCorrelation] = {}
    weighted_pcc_num = weighted_pcc_den = 0.0
    weighted_spcc_num = weighted_spcc_den = 0.0
    for gene, idx in indices_by_gene.items():
        gene_predicted = [predicted_kd[i] for i in idx]
        gene_actual = [actual_kd[i] for i in idx]
        gene_pcc = _pearson(gene_predicted, gene_actual)
        gene_spcc = _spearman(gene_predicted, gene_actual)
        by_gene[gene] = GeneCorrelation(n=len(idx), pcc=gene_pcc, spcc=gene_spcc)

        # Genes with an undefined per-gene correlation (too few siRNAs, or
        # zero variance) contribute no term to the weighted average --
        # there's nothing to weight, and including their weight without a
        # value would just bias the average down.
        if not math.isnan(gene_pcc):
            weighted_pcc_num += len(idx) * gene_pcc
            weighted_pcc_den += len(idx)
        if not math.isnan(gene_spcc):
            weighted_spcc_num += len(idx) * gene_spcc
            weighted_spcc_den += len(idx)

    metrics.by_gene = by_gene
    metrics.weighted_pcc_by_gene = (
        weighted_pcc_num / weighted_pcc_den if weighted_pcc_den else float("nan")
    )
    metrics.weighted_spcc_by_gene = (
        weighted_spcc_num / weighted_spcc_den if weighted_spcc_den else float("nan")
    )
    return metrics
