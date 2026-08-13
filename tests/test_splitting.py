from __future__ import annotations

import math

import pytest

from sirna_data.raw_loader import SiRNARecord, load_records
from sirna_data.splitting import leave_n_genes_out, train_test_split


def _record(row_id: str, gene: str) -> SiRNARecord:
    return SiRNARecord(
        row_id=row_id,
        gene=gene,
        accession=f"ACC_{gene}",
        guide_seq="ACGUACGUACGUACGUACG",
        duplex_len=19,
        mrna_window="AAAAAACGUACGUACGUACGUACGCCCCC",
        site_start=5,
        site_len=19,
        has_flanking_context=True,
        label=50.0,
        technology="Luciferase reporter assay",
        source="TestSource",
    )


def _records_for_genes(genes_with_counts: dict[str, int]) -> list[SiRNARecord]:
    records = []
    for gene, count in genes_with_counts.items():
        for i in range(count):
            records.append(_record(f"{gene}_{i}", gene))
    return records


# --------------------------------------------------------------------------
# by_gene=True (default): no gene straddles both splits
# --------------------------------------------------------------------------


def test_by_gene_default_is_true():
    records = _records_for_genes({f"GENE{i}": 3 for i in range(10)})
    train, test = train_test_split(records, random_state=0)
    train_explicit, test_explicit = train_test_split(records, by_gene=True, random_state=0)
    assert [r.row_id for r in train] == [r.row_id for r in train_explicit]
    assert [r.row_id for r in test] == [r.row_id for r in test_explicit]


def test_by_gene_no_gene_appears_in_both_splits():
    records = _records_for_genes({f"GENE{i}": 4 for i in range(20)})
    train, test = train_test_split(records, test_size=0.3, random_state=1)
    train_genes = {r.gene for r in train}
    test_genes = {r.gene for r in test}
    assert train_genes.isdisjoint(test_genes)


def test_by_gene_covers_every_record_exactly_once():
    records = _records_for_genes({f"GENE{i}": 3 for i in range(15)})
    train, test = train_test_split(records, random_state=2)
    assert sorted(r.row_id for r in train + test) == sorted(r.row_id for r in records)
    assert len(train) + len(test) == len(records)


def test_by_gene_test_size_approximately_respected():
    records = _records_for_genes({f"GENE{i}": 2 for i in range(100)})
    train, test = train_test_split(records, test_size=0.2, random_state=3)
    test_genes = {r.gene for r in test}
    assert 15 <= len(test_genes) <= 25  # ~20 of 100 genes, some rounding slack


def test_by_gene_is_deterministic_with_same_random_state():
    records = _records_for_genes({f"GENE{i}": 3 for i in range(20)})
    train1, test1 = train_test_split(records, random_state=42)
    train2, test2 = train_test_split(records, random_state=42)
    assert [r.row_id for r in train1] == [r.row_id for r in train2]
    assert [r.row_id for r in test1] == [r.row_id for r in test2]


def test_by_gene_requires_at_least_two_genes():
    records = _records_for_genes({"ONLY_GENE": 5})
    with pytest.raises(ValueError):
        train_test_split(records, by_gene=True)


def test_by_gene_preserves_relative_order_within_each_split():
    records = _records_for_genes({f"GENE{i}": 1 for i in range(20)})
    train, test = train_test_split(records, random_state=5)
    original_order = [r.row_id for r in records]
    train_positions = [original_order.index(r.row_id) for r in train]
    test_positions = [original_order.index(r.row_id) for r in test]
    assert train_positions == sorted(train_positions)
    assert test_positions == sorted(test_positions)


# --------------------------------------------------------------------------
# by_gene=False: random split by individual record, ignoring gene
# --------------------------------------------------------------------------


def test_by_record_covers_every_record_exactly_once():
    records = _records_for_genes({f"GENE{i}": 5 for i in range(10)})
    train, test = train_test_split(records, by_gene=False, random_state=6)
    assert sorted(r.row_id for r in train + test) == sorted(r.row_id for r in records)
    assert len(train) + len(test) == len(records)


def test_by_record_test_size_approximately_respected():
    records = _records_for_genes({f"GENE{i}": 5 for i in range(20)})  # 100 records
    train, test = train_test_split(records, test_size=0.25, by_gene=False, random_state=7)
    assert len(test) == round(0.25 * 100)


def test_by_record_can_put_same_gene_in_both_splits():
    # With enough siRNAs per gene and a record-level (not gene-level) split,
    # at least one gene should end up represented on both sides.
    records = _records_for_genes({f"GENE{i}": 10 for i in range(5)})  # 50 records
    train, test = train_test_split(records, test_size=0.5, by_gene=False, random_state=8)
    train_genes = {r.gene for r in train}
    test_genes = {r.gene for r in test}
    assert train_genes & test_genes


def test_by_record_is_deterministic_with_same_random_state():
    records = _records_for_genes({f"GENE{i}": 5 for i in range(10)})
    train1, test1 = train_test_split(records, by_gene=False, random_state=9)
    train2, test2 = train_test_split(records, by_gene=False, random_state=9)
    assert [r.row_id for r in train1] == [r.row_id for r in train2]
    assert [r.row_id for r in test1] == [r.row_id for r in test2]


def test_by_record_requires_at_least_two_records():
    records = _records_for_genes({"GENE1": 1})
    with pytest.raises(ValueError):
        train_test_split(records, by_gene=False)


# --------------------------------------------------------------------------
# shared validation
# --------------------------------------------------------------------------


def test_rejects_empty_records():
    with pytest.raises(ValueError):
        train_test_split([], random_state=0)


@pytest.mark.parametrize("bad_test_size", [0.0, 1.0, -0.1, 1.1])
def test_rejects_bad_test_size(bad_test_size):
    records = _records_for_genes({f"GENE{i}": 2 for i in range(5)})
    with pytest.raises(ValueError):
        train_test_split(records, test_size=bad_test_size)


def test_unseeded_calls_still_return_valid_disjoint_splits():
    # No random_state given -- just check the invariants hold, not exact output.
    records = _records_for_genes({f"GENE{i}": 3 for i in range(10)})
    train, test = train_test_split(records)
    assert len(train) > 0
    assert len(test) > 0
    assert len(train) + len(test) == len(records)


# --------------------------------------------------------------------------
# integration: load_records() -> train_test_split()
# --------------------------------------------------------------------------


def test_train_test_split_on_loaded_records(patch_data_dir, fake_data_dir):
    # patch_data_dir/fake_data_dir (conftest.py) build a small, self-contained
    # fake data/raw/ directory so this exercises the real load_records() merge
    # path without touching the actual multi-source dataset.
    records = load_records(
        csv_path=fake_data_dir / "sirna_efficacy.csv",
        fasta_path=fake_data_dir / "mrna_transcripts.fasta",
    )
    assert len(records) > 0
    all_genes = {r.gene for r in records}
    assert len(all_genes) >= 2  # fixture ships 6 genes across its 5 sources

    train, test = train_test_split(records, test_size=0.3, random_state=0)

    # every record accounted for exactly once, no gene split across both
    assert sorted(r.row_id for r in train + test) == sorted(r.row_id for r in records)
    train_genes = {r.gene for r in train}
    test_genes = {r.gene for r in test}
    assert train_genes.isdisjoint(test_genes)
    assert train_genes | test_genes == all_genes


def test_train_test_split_on_loaded_records_by_record(patch_data_dir, fake_data_dir):
    records = load_records(
        csv_path=fake_data_dir / "sirna_efficacy.csv",
        fasta_path=fake_data_dir / "mrna_transcripts.fasta",
    )
    train, test = train_test_split(records, test_size=0.3, by_gene=False, random_state=0)
    assert sorted(r.row_id for r in train + test) == sorted(r.row_id for r in records)


# --------------------------------------------------------------------------
# leave_n_genes_out: gene-level cross-validation generator
# --------------------------------------------------------------------------


def test_leave_n_genes_out_is_a_generator():
    import types

    records = _records_for_genes({f"GENE{i}": 3 for i in range(6)})
    result = leave_n_genes_out(records, 2)
    assert isinstance(result, types.GeneratorType)


def test_leave_one_gene_out_yields_one_fold_per_gene():
    records = _records_for_genes({f"GENE{i}": 3 for i in range(6)})
    folds = list(leave_n_genes_out(records, 1, random_state=0))
    assert len(folds) == 6
    for _train, test in folds:
        test_genes = {r.gene for r in test}
        assert len(test_genes) == 1


def test_leave_n_genes_out_fold_count_matches_ceil_division():
    records = _records_for_genes({f"GENE{i}": 2 for i in range(10)})
    folds = list(leave_n_genes_out(records, 3, random_state=0))
    assert len(folds) == math.ceil(10 / 3)


def test_leave_n_genes_out_every_gene_appears_in_exactly_one_test_fold():
    records = _records_for_genes({f"GENE{i}": 3 for i in range(10)})
    all_genes = {r.gene for r in records}
    folds = list(leave_n_genes_out(records, 3, random_state=0))

    seen_test_genes: list[str] = []
    for train, test in folds:
        test_genes = {r.gene for r in test}
        train_genes = {r.gene for r in train}
        assert train_genes.isdisjoint(test_genes)  # never leaks within a fold
        seen_test_genes.extend(test_genes)

    # every gene shows up exactly once across all test folds combined
    assert sorted(seen_test_genes) == sorted(all_genes)


def test_leave_n_genes_out_every_record_covered_each_fold():
    records = _records_for_genes({f"GENE{i}": 4 for i in range(8)})
    for train, test in leave_n_genes_out(records, 2, random_state=0):
        assert sorted(r.row_id for r in train + test) == sorted(r.row_id for r in records)
        assert len(train) + len(test) == len(records)


def test_leave_n_genes_out_last_fold_can_be_smaller():
    # 7 genes, n=3 -> folds of size 3, 3, 1
    records = _records_for_genes({f"GENE{i}": 2 for i in range(7)})
    folds = list(leave_n_genes_out(records, 3, random_state=0))
    sizes = [len({r.gene for r in test}) for _, test in folds]
    assert sorted(sizes) == [1, 3, 3]


def test_leave_n_genes_out_is_deterministic_with_same_random_state():
    records = _records_for_genes({f"GENE{i}": 3 for i in range(10)})
    folds1 = list(leave_n_genes_out(records, 2, random_state=7))
    folds2 = list(leave_n_genes_out(records, 2, random_state=7))
    for (train1, test1), (train2, test2) in zip(folds1, folds2, strict=True):
        assert [r.row_id for r in train1] == [r.row_id for r in train2]
        assert [r.row_id for r in test1] == [r.row_id for r in test2]


def test_leave_n_genes_out_n_larger_than_gene_count_yields_single_fold():
    records = _records_for_genes({f"GENE{i}": 2 for i in range(3)})
    folds = list(leave_n_genes_out(records, 10, random_state=0))
    assert len(folds) == 1
    train, test = folds[0]
    assert train == []
    assert sorted(r.row_id for r in test) == sorted(r.row_id for r in records)


def test_leave_n_genes_out_n_equals_one_matches_leave_one_gene_out_semantics():
    records = _records_for_genes({f"GENE{i}": 2 for i in range(5)})
    folds = list(leave_n_genes_out(records, 1, random_state=3))
    assert len(folds) == 5
    tested_genes = [next(iter({r.gene for r in test})) for _, test in folds]
    assert sorted(tested_genes) == sorted({r.gene for r in records})


def test_leave_n_genes_out_rejects_empty_records():
    with pytest.raises(ValueError):
        list(leave_n_genes_out([], 1))


@pytest.mark.parametrize("bad_n", [0, -1])
def test_leave_n_genes_out_rejects_non_positive_n(bad_n):
    records = _records_for_genes({f"GENE{i}": 2 for i in range(3)})
    with pytest.raises(ValueError):
        list(leave_n_genes_out(records, bad_n))


def test_leave_n_genes_out_on_loaded_records(patch_data_dir, fake_data_dir):
    records = load_records(
        csv_path=fake_data_dir / "sirna_efficacy.csv",
        fasta_path=fake_data_dir / "mrna_transcripts.fasta",
    )
    all_genes = {r.gene for r in records}
    seen_test_genes: set[str] = set()
    fold_count = 0
    for train, test in leave_n_genes_out(records, 2, random_state=0):
        fold_count += 1
        train_genes = {r.gene for r in train}
        test_genes = {r.gene for r in test}
        assert train_genes.isdisjoint(test_genes)
        assert sorted(r.row_id for r in train + test) == sorted(r.row_id for r in records)
        seen_test_genes |= test_genes
    assert fold_count == math.ceil(len(all_genes) / 2)
    assert seen_test_genes == all_genes
