from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

from sirna_data.genes import describe_genes, list_genes
from sirna_data.raw_loader import SiRNARecord, load_records

FLANK = 5  # must match conftest.FLANK

# Every gene the fake fixture dataset contributes, in the order
# list_genes()/describe_genes() should return them (case-insensitive by
# gene name).
FIXTURE_GENES = [
    "GENEA",  # siRNAEfficacyDB
    "GENEB",  # siRNAEfficacyDB (unlocatable accession)
    "GENEC",  # Monopoli2023
    "GENED",  # Shabalina2006
    "GENEE",  # Davis2025 (2 rows)
    "GENEF",  # CMsiRNAdb_full
    "MARTGENE",  # Martinelli (2 rows)
    "PCSK9",  # CMsiRNAdb
    "REMOVED",  # REMOVED_REMOVED
]


def _record(*, gene: str, source: str) -> SiRNARecord:
    """A minimal SiRNARecord, for the cases that are about gene/source
    bookkeeping rather than about any real source's data."""
    return SiRNARecord(
        row_id="r",
        gene=gene,
        accession="ACC",
        guide_seq="ACGU",
        duplex_len=4,
        mrna_window="ACGU",
        site_start=0,
        site_len=4,
        has_flanking_context=False,
        label=50.0,
        technology="test",
        source=source,
    )


def test_list_genes_loads_and_lists_every_gene(patch_data_dir: Path):
    assert list_genes(flank_nt=FLANK) == FIXTURE_GENES


def test_list_genes_from_already_loaded_records(patch_data_dir: Path):
    records = load_records(flank_nt=FLANK)
    assert list_genes(records) == FIXTURE_GENES


def test_describe_genes_counts_records_per_gene(patch_data_dir: Path):
    by_gene = {info.gene: info for info in describe_genes(flank_nt=FLANK)}
    assert by_gene["GENEA"].n_records == 1
    # both multi-row fixture sources
    assert by_gene["MARTGENE"].n_records == 2
    assert by_gene["GENEE"].n_records == 2
    assert sum(info.n_records for info in by_gene.values()) == len(load_records(flank_nt=FLANK))


def test_describe_genes_reports_sources_and_licenses(patch_data_dir: Path):
    by_gene = {info.gene: info for info in describe_genes(flank_nt=FLANK)}
    assert by_gene["PCSK9"].sources == ("CMsiRNAdb",)
    assert by_gene["PCSK9"].licenses == ("CC BY-NC-ND 4.0",)
    assert by_gene["GENEA"].licenses == ("CC BY-NC",)
    # the REMOVED panel's license could not be established -- that must show
    # up here rather than being reported as permissive or dropped.
    assert by_gene["REMOVED"].licenses == ("unresolved",)


def test_describe_genes_reports_accessions(patch_data_dir: Path):
    by_gene = {info.gene: info for info in describe_genes(flank_nt=FLANK)}
    assert by_gene["GENEA"].accessions == ("ACC1",)
    assert by_gene["GENEB"].accessions == ("ACC_MISSING",)


def test_describe_genes_tolerates_missing_accessions(patch_data_dir: Path):
    # ~900 real CMsiRNAdb_full rows have no accession at all, which arrives
    # from pandas as NaN rather than a string -- it must not blow up the
    # sort, and must leave `accessions` empty rather than inventing one.
    records = load_records(flank_nt=FLANK)
    target = next(r for r in records if r.gene == "GENEF")
    records = [r for r in records if r.gene != "GENEF"] + [
        replace(target, accession=float("nan")),
        replace(target, accession="  "),
    ]
    by_gene = {info.gene: info for info in describe_genes(records)}
    assert by_gene["GENEF"].accessions == ()
    assert by_gene["GENEF"].n_records == 2


def test_describe_genes_and_list_genes_agree(patch_data_dir: Path):
    records = load_records(flank_nt=FLANK)
    assert list_genes(records) == [info.gene for info in describe_genes(records)]


def test_gene_strings_are_not_normalized():
    # Two textually-distinct raw gene strings (siRNAEfficacyDB's "EGFP "
    # has a trailing space) stay two entries, matching load_records().
    records = [
        _record(gene="EGFP", source="Martinelli_sirna_reproduction"),
        _record(gene="EGFP ", source="siRNAEfficacyDB"),
    ]
    infos = describe_genes(records)
    assert [info.gene for info in infos] == ["EGFP", "EGFP "]


def test_unknown_source_is_reported_as_unresolved():
    infos = describe_genes([_record(gene="GENEX", source="SomeSourceAddedLater")])
    assert infos[0].licenses == ("unresolved",)


# --------------------------------------------------------------------------
# argument handling
# --------------------------------------------------------------------------


def test_load_records_arguments_are_forwarded(patch_data_dir: Path):
    # Only the permissively-licensed fixture sources: Monopoli (GENEC),
    # Shabalina (GENED), Davis2025 (GENEE).
    assert list_genes(flank_nt=FLANK, licenses=["CC BY 4.0", "CC BY 2.0"]) == [
        "GENEC",
        "GENED",
        "GENEE",
    ]
    assert "PCSK9" not in list_genes(flank_nt=FLANK, include_cmsirnadb=False)


def test_records_and_load_arguments_are_mutually_exclusive(patch_data_dir: Path):
    records = load_records(flank_nt=FLANK)
    with pytest.raises(TypeError, match="one or the other"):
        list_genes(records, licenses=["CC BY 4.0"])
    with pytest.raises(TypeError, match="one or the other"):
        describe_genes(records, include_monopoli=False)


def test_empty_record_list_describes_nothing():
    assert list_genes([]) == []
    assert describe_genes([]) == []
