"""Split a list of `SiRNARecord` into train/test sets, either grouped by
gene (no gene appears in both splits) or as a plain random split with no
regard for gene at all -- plus a leave-N-genes-out generator for gene-level
cross-validation.

Motivating problem: naive random splitting of siRNA/knockdown records can
put two siRNAs targeting the *same* gene into both the train and test sets.
A model can then partly "solve" a test example just by having seen a
different siRNA against the same gene/mRNA during training, inflating
apparent test performance relative to how the model will actually do on a
genuinely novel gene. Grouping by gene (the default here) avoids that;
splitting by individual record instead is occasionally still useful (e.g.
matching a baseline paper's own non-grouped split, or specifically measuring
within-gene generalization).

`train_test_split` deliberately mirrors `sklearn.model_selection
.train_test_split`'s name and parameter conventions (`test_size`,
`random_state`) so it's a familiar drop-in for anyone used to sklearn -- the
one deliberate difference is `by_gene`, which defaults to `True` and has no
sklearn equivalent (sklearn's version has no concept of grouping).
"""
from __future__ import annotations

import random
from collections.abc import Iterator, Sequence

from .raw_loader import SiRNARecord


def train_test_split(
    records: Sequence[SiRNARecord],
    test_size: float = 0.2,
    *,
    by_gene: bool = True,
    random_state: int | None = None,
) -> tuple[list[SiRNARecord], list[SiRNARecord]]:
    """Split `records` into (train, test) lists.

    Parameters
    ----------
    records : the records to split, e.g. straight from `load_records()`.
    test_size : target fraction in (0, 1) of whichever unit is being split
        -- genes if `by_gene=True`, individual records if `by_gene=False`.
        Rounded to the nearest whole gene/record, then clamped to at least
        1 and at most (count - 1) so neither split ends up empty.
    by_gene : if True (the default), split by gene -- every record for a
        given gene goes entirely into train or entirely into test, so no
        gene appears in both. This is almost always what you want when
        evaluating how well a model generalizes to a genuinely new gene,
        since siRNAs targeting the same gene/mRNA share sequence context
        that a naive per-record split would leak across train and test. If
        False, individual records are split at random with no regard for
        which gene they belong to -- the same gene can (and, with more than
        one siRNA per gene, usually will) appear in both splits. This is the
        one parameter with no sklearn equivalent -- everything else here
        follows `sklearn.model_selection.train_test_split`'s conventions.
    random_state : random seed for reproducibility, same meaning as
        sklearn's `random_state`. `None` (the default) uses unseeded
        randomness -- a different split each call.

    Returns
    -------
    (train_records, test_records) -- each a list of `SiRNARecord`, in the
    same relative order they appeared in `records` (not shuffled).

    Raises `ValueError` if `records` is empty, `test_size` is not in
    (0, 1), or there are too few genes/records to form two non-empty splits
    (e.g. `by_gene=True` with only a single gene present).
    """
    if not records:
        raise ValueError("records must be non-empty")
    if not 0.0 < test_size < 1.0:
        raise ValueError(f"test_size must be in (0, 1), got {test_size!r}")

    rng = random.Random(random_state)

    if by_gene:
        genes = sorted({r.gene for r in records})
        if len(genes) < 2:
            raise ValueError(
                f"by_gene=True needs at least 2 distinct genes to split, got {len(genes)!r}"
            )
        rng.shuffle(genes)
        n_test = min(len(genes) - 1, max(1, round(test_size * len(genes))))
        test_genes = set(genes[:n_test])
        train_records = [r for r in records if r.gene not in test_genes]
        test_records = [r for r in records if r.gene in test_genes]
        return train_records, test_records

    n = len(records)
    if n < 2:
        raise ValueError(f"by_gene=False needs at least 2 records to split, got {n!r}")
    indices = list(range(n))
    rng.shuffle(indices)
    n_test = min(n - 1, max(1, round(test_size * n)))
    test_idx = set(indices[:n_test])
    train_records = [r for i, r in enumerate(records) if i not in test_idx]
    test_records = [r for i, r in enumerate(records) if i in test_idx]
    return train_records, test_records


def leave_n_genes_out(
    records: Sequence[SiRNARecord],
    n: int,
    *,
    random_state: int | None = None,
) -> Iterator[tuple[list[SiRNARecord], list[SiRNARecord]]]:
    """Generator for gene-level cross-validation: partitions every distinct
    gene in `records` into disjoint groups of `n` genes, then yields one
    (train, test) fold per group -- test is every record for that group's
    genes, train is everything else. Every gene appears in exactly one
    test fold across the full iteration (never repeated, never skipped),
    generalizing leave-one-gene-out CV (`n=1`) to leave-`n`-genes-out.

    The gene order is shuffled once up front (seeded by `random_state`) before
    being cut into consecutive groups of `n`, so folds are reproducible given
    the same `random_state` but not tied to genes' original/alphabetical
    order. If the gene count isn't evenly divisible by `n`, the last fold
    holds out fewer than `n` genes.

    Parameters
    ----------
    records : the records to fold over, e.g. straight from `load_records()`.
    n : how many genes to hold out per fold. Must be a positive integer;
        can exceed the number of distinct genes, in which case a single
        fold is yielded with every gene (and therefore every record) in
        test and nothing in train.
    random_state : random seed controlling the shuffle order genes are
        grouped in, same meaning as `train_test_split`'s `random_state`.
        `None` (the default) uses unseeded randomness -- a different
        partition each call.

    Yields
    ------
    (train_records, test_records) tuples, `ceil(n_genes / n)` of them in
    total, each a list of `SiRNARecord` in the same relative order they
    appeared in `records` (not shuffled).

    Raises `ValueError` up front (before yielding anything) if `records` is
    empty or `n` is not a positive integer.
    """
    if not records:
        raise ValueError("records must be non-empty")
    if n < 1:
        raise ValueError(f"n must be a positive integer, got {n!r}")

    genes = sorted({r.gene for r in records})
    rng = random.Random(random_state)
    rng.shuffle(genes)

    for i in range(0, len(genes), n):
        test_genes = set(genes[i : i + n])
        train_records = [r for r in records if r.gene not in test_genes]
        test_records = [r for r in records if r.gene in test_genes]
        yield train_records, test_records
