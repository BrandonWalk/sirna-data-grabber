"""Fetcher for the Martinelli 2023 / sirna-reproduction subset.

The upstream table has no gene column, so gene identity is recovered by
matching each row against the known transcripts. These tests pin that rule,
including the two cases it must refuse: a sequence in no transcript, and one
in more than one (the two firefly luciferase reporters overlap).
"""
from __future__ import annotations

import csv
from pathlib import Path

import pytest

from sirna_data.fetch import martinelli
from sirna_data.raw_loader import _revcomp

SITE_A = "GACGUAAACGGCCACAAGU"          # 19nt, in the EGFP transcript below
SITE_B = "UGGCCAAUGCCAAGGAGAU"          # 19nt, in the APOB transcript below
SHARED = "AAAACCCCGGGGUUUUAAA"          # in BOTH luciferase transcripts
ORPHAN = "CAUCAUCAUCAUCAUCAUC"          # in none

TRANSCRIPTS = {
    "U55763": "AAAAA" + SITE_A + "CCCCC",
    "NM_000384": "GGGGG" + SITE_B + "UUUUU",
    "U47296": "AAA" + SHARED + "CCC",
    "X65324": "GGG" + SHARED + "UUU",
}

COLUMNS = martinelli.SOURCE_COLUMNS


def _row(identifier: str, sense: str, antisense: str = "", pct: str = "84.9") -> dict[str, str]:
    return {
        COLUMNS["id"]: identifier, COLUMNS["pmid"]: "19282453",
        COLUMNS["sense"]: sense, COLUMNS["sense_modification"]: "hexitol nucleic acid",
        COLUMNS["antisense"]: antisense, COLUMNS["antisense_modification"]: "0",
        COLUMNS["pct"]: pct,
    }


def test_places_a_row_by_its_sense_strand():
    rows = martinelli.assign_genes([_row("SM1", SITE_A + "TT")], TRANSCRIPTS)
    assert len(rows) == 1
    assert rows[0]["Gene"] == "EGFP"
    assert rows[0]["Accession_number"] == "U55763"


def test_sequences_are_stored_in_the_rna_alphabet():
    """The source writes DNA overhangs as T; the corpus is RNA throughout,
    and the chemistry is recorded in the modification columns instead."""
    rows = martinelli.assign_genes(
        [_row("SM1", SITE_A + "TT", antisense=_revcomp(SITE_A) + "TT")], TRANSCRIPTS
    )
    assert rows[0]["Sequence"] == SITE_A + "UU"
    assert rows[0]["Sequence_antisense"] == _revcomp(SITE_A) + "UU"
    assert rows[0]["Modification_sense"] == "hexitol nucleic acid"


def test_falls_back_to_the_antisense_strand_when_the_sense_does_not_match():
    rows = martinelli.assign_genes(
        [_row("SM2", ORPHAN + "TT", antisense=_revcomp(SITE_B) + "TT")], TRANSCRIPTS
    )
    assert len(rows) == 1
    assert rows[0]["Gene"] == "APOB"
    assert rows[0]["Sequence"] == SITE_B   # the verified target site, not the orphan strand


def test_drops_a_row_that_matches_no_transcript():
    assert martinelli.assign_genes([_row("SM3", ORPHAN + "TT")], TRANSCRIPTS) == []


def test_drops_a_row_that_matches_more_than_one_transcript():
    """Ambiguous between the two firefly reporters -- dropped, not guessed."""
    assert martinelli.assign_genes([_row("SM4", SHARED)], TRANSCRIPTS) == []


def test_drops_a_row_too_short_to_place():
    assert martinelli.assign_genes([_row("SM5", "ACGUACGU")], TRANSCRIPTS) == []


def test_carries_both_strands_and_their_modifications_through():
    rows = martinelli.assign_genes(
        [_row("SM6", SITE_A + "TT", antisense=_revcomp(SITE_A) + "TT")], TRANSCRIPTS
    )
    assert rows[0]["Sequence_antisense"] == _revcomp(SITE_A) + "UU"
    assert rows[0]["Modification_sense"] == "hexitol nucleic acid"
    assert rows[0]["Modification_antisense"] == "0"
    assert float(rows[0]["PCT"]) == 84.9
    assert rows[0]["Experiment_ID"] == "SM6"


@pytest.fixture
def only_these_transcripts(monkeypatch: pytest.MonkeyPatch):
    """fetch() insists NCBI returns every accession it asked for, so the
    fixture's four-transcript world needs a matching gene map."""
    monkeypatch.setattr(
        martinelli, "GENE_MAP",
        {accession: martinelli.GENE_MAP[accession] for accession in TRANSCRIPTS},
    )


def test_fetch_writes_csv_and_fasta(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, only_these_transcripts):
    monkeypatch.setattr(martinelli, "fetch_table", lambda: [_row("SM1", SITE_A + "TT")])
    monkeypatch.setattr(martinelli, "fetch_transcripts", lambda _a: TRANSCRIPTS)
    martinelli.fetch(tmp_path)

    rows = list(csv.DictReader(open(tmp_path / "martinelli_extra.csv")))
    assert len(rows) == 1 and list(rows[0]) == martinelli.COLUMNS
    fasta = (tmp_path / "martinelli_transcripts.fasta").read_text()
    assert fasta.count(">") == len(TRANSCRIPTS)
    sequence_lines = [line for line in fasta.splitlines() if not line.startswith(">")]
    assert "U" not in "".join(sequence_lines)   # FASTA is written in the DNA alphabet


def test_fetch_raises_when_nothing_can_be_placed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, only_these_transcripts):
    monkeypatch.setattr(martinelli, "fetch_table", lambda: [_row("SM9", ORPHAN)])
    monkeypatch.setattr(martinelli, "fetch_transcripts", lambda _a: TRANSCRIPTS)
    with pytest.raises(ValueError, match="no rows could be placed"):
        martinelli.fetch(tmp_path)


def test_fetch_table_rejects_an_upstream_file_missing_a_column(
    monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(martinelli, "_download", lambda _url: b"id,sequence\n1,ACGU\n")
    with pytest.raises(ValueError, match="missing expected column"):
        martinelli.fetch_table()


# --------------------------------------------------------------------------
# canonical form: what `fetch()` writes, and what `normalize_file` rewrites a
# hand-built file into, must be the same thing
# --------------------------------------------------------------------------


MIXED_ROW = {
    "Experiment_ID": "SM10", "PMID": "123", "Gene": "EGFP",
    "Accession_number": "U55763", "Sequence": SITE_A + "TT",
    "Modification_sense": "0", "Sequence_antisense": _revcomp(SITE_A) + "tt",
    "Modification_antisense": "0", "PCT": "84",
}


def test_canonical_row_normalizes_alphabet_and_numbers():
    row = martinelli.canonical_row(MIXED_ROW)
    assert row["Sequence"] == SITE_A + "UU"
    assert row["Sequence_antisense"] == _revcomp(SITE_A) + "UU"
    assert row["PCT"] == "84.0"


def test_canonical_row_is_idempotent():
    once = martinelli.canonical_row(MIXED_ROW)
    assert martinelli.canonical_row(once) == once


def test_sort_key_orders_by_experiment_number_not_lexically():
    rows = [{"Experiment_ID": name} for name in ("SM100", "SM9", "SM1000")]
    assert [r["Experiment_ID"] for r in sorted(rows, key=martinelli.sort_key)] == [
        "SM9", "SM100", "SM1000"
    ]


def _write(path: Path, rows: list[dict[str, str]]) -> None:
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=martinelli.COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def test_normalize_file_rewrites_and_reports_changed_rows(tmp_path: Path):
    path = tmp_path / "martinelli_extra.csv"
    _write(path, [MIXED_ROW])
    assert martinelli.normalize_file(path) == 1
    row = next(iter(csv.DictReader(open(path))))
    assert row["Sequence"] == SITE_A + "UU"
    assert row["PCT"] == "84.0"
    # running it again changes nothing
    assert martinelli.normalize_file(path) == 0


def test_normalize_file_puts_rows_in_canonical_order(tmp_path: Path):
    path = tmp_path / "martinelli_extra.csv"
    _write(path, [dict(MIXED_ROW, Experiment_ID=name) for name in ("SM100", "SM9")])
    martinelli.normalize_file(path)
    assert [row["Experiment_ID"] for row in csv.DictReader(open(path))] == ["SM9", "SM100"]


def test_normalize_file_can_rederive_the_stored_target_site(tmp_path: Path):
    """A hand-picked short window is replaced by the canonical one."""
    path = tmp_path / "martinelli_extra.csv"
    _write(path, [dict(MIXED_ROW, Sequence=SITE_A[:19], Sequence_antisense="")])
    martinelli.normalize_file(path, TRANSCRIPTS)
    row = next(iter(csv.DictReader(open(path))))
    assert row["Sequence"] == SITE_A   # re-derived against the transcript


def test_normalize_file_refuses_to_reassign_a_gene(tmp_path: Path):
    path = tmp_path / "martinelli_extra.csv"
    _write(path, [dict(MIXED_ROW, Accession_number="NM_000384")])  # claims APOB
    with pytest.raises(ValueError, match="refusing to reassign a gene"):
        martinelli.normalize_file(path, TRANSCRIPTS)


def test_normalize_file_refuses_a_row_it_cannot_place(tmp_path: Path):
    path = tmp_path / "martinelli_extra.csv"
    _write(path, [dict(MIXED_ROW, Sequence=ORPHAN, Sequence_antisense="")])
    with pytest.raises(ValueError, match="no longer matches any known transcript"):
        martinelli.normalize_file(path, TRANSCRIPTS)


def test_normalize_file_rejects_a_file_missing_a_column(tmp_path: Path):
    path = tmp_path / "martinelli_extra.csv"
    path.write_text("Experiment_ID,Sequence\nSM1,ACGU\n")
    with pytest.raises(ValueError, match="missing expected column"):
        martinelli.normalize_file(path)
