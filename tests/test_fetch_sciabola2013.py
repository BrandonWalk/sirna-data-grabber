"""Fetcher for Sciabola et al. 2013's in-house Hep3B panel.

The Word decoding itself is covered in tests/test_fetch_ole.py; here the
document is replaced by a synthetic Word text stream so these tests stay
about the fetcher's own logic: table geometry, strand normalization,
transcript verification and the per-dose columns.
"""
from __future__ import annotations

import csv
import io
import zipfile
from pathlib import Path

import pytest

from sirna_data.fetch import sciabola2013
from sirna_data.fetch._ole import CELL_MARK

SENSE = "ACGUACGUACGUACGUACGUA"          # 21nt, listed sense in the fake table
ANTISENSE_TARGET = "GGCAUGCAUGCAUGCAUGCAU"  # 21nt, listed as the guide strand
UNKNOWN = "UUUUUUUUUUUUUUUUUUUUU"          # in no transcript -> dropped
HIF1A = "NM_001530.4"


def _row(cells: list[str]) -> str:
    """One Word table row: a cell mark after every cell, then the row mark."""
    return CELL_MARK.join(cells) + CELL_MARK + CELL_MARK


def _document() -> str:
    sequences = _row(["Gene", "19p2name", "19p2_seq", "Rname", "R-Dicer_seq",
                      "Lname", "L-Dicer_seq"])
    sequences += _row(["HIF1a", "1", SENSE, "1R", "x", "1L", "y"])
    sequences += _row(["HIF1a", "2", sciabola2013.revcomp(ANTISENSE_TARGET), "2R", "x", "2L", "y"])
    sequences += _row(["HIF1a", "3", UNKNOWN, "3R", "x", "3L", "y"])
    values = CELL_MARK.join(["19p2", "R-Dicer", "L-Dicer"]) + CELL_MARK + CELL_MARK
    values += _row(["Name", "0.08nM", "0.4nM", "2nM", "10nM",
                    "Name", "0.08nM", "0.4nM", "2nM", "Name", "0.08nM", "0.4nM", "2nM"])
    values += _row(["1", "10", "20", "30", "40", "1R", "", "", "", "1L", "", "", ""])
    values += _row(["2", "", "", "", "80", "2R", "", "", "", "2L", "", "", ""])
    values += _row(["3", "5", "15", "", "", "3R", "", "", "", "3L", "", "", ""])
    return f"Table S3\r{sequences}\rTable S4\r{values}\r"


@pytest.fixture
def fake_document(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(sciabola2013._ole, "word_text", lambda _blob: _document())
    monkeypatch.setattr(sciabola2013, "fetch_document", lambda: b"fake-doc")


@pytest.fixture
def fake_transcripts(monkeypatch: pytest.MonkeyPatch):
    transcripts = {
        accession: ("AAAAA" + SENSE + "CCCCC" + ANTISENSE_TARGET + "GGGGG")
        for accession in set(sciabola2013.ACCESSIONS.values())
        | set(sciabola2013.VARIANT_ACCESSIONS.values())
    }
    monkeypatch.setattr(sciabola2013, "fetch_transcripts", lambda _a: transcripts)
    return transcripts


def test_parse_tables_recovers_both_grids(fake_document):
    sequences, values = sciabola2013.parse_tables(b"fake-doc")
    assert [row[1] for row in sequences] == ["1", "2", "3"]
    assert sequences[0][2] == SENSE
    # empty cells must survive: row 2 has only its 10 nM value
    assert values[1][:5] == ["2", "", "", "", "80"]


def test_build_rows_normalizes_strand_and_verifies_against_transcripts(
    fake_document, fake_transcripts):
    sequences, values = sciabola2013.parse_tables(b"fake-doc")
    rows, dropped = sciabola2013.build_rows(sequences, values, fake_transcripts)

    assert [row["Compound_Name"] for row in rows] == ["19p2_1", "19p2_2"]
    assert dropped == ["19p2_3"]

    listed_sense, listed_guide = rows
    assert listed_sense["Table_S3_Strand"] == "sense"
    assert listed_sense["Sense_21mer"] == SENSE
    # the row listed as a guide is stored as its reverse complement
    assert listed_guide["Table_S3_Strand"] == "antisense"
    assert listed_guide["Sense_21mer"] == ANTISENSE_TARGET
    assert listed_guide["Table_S3_Sequence"] == sciabola2013.revcomp(ANTISENSE_TARGET)
    # gene labels are normalized, the source's own spelling is kept
    assert listed_sense["Gene"] == "HIF1A" and listed_sense["Source_Gene_Label"] == "HIF1a"
    assert listed_sense["Accession_number"] == HIF1A


def test_build_rows_keeps_every_dose_column_including_absent_ones(
    fake_document, fake_transcripts):
    sequences, values = sciabola2013.parse_tables(b"fake-doc")
    rows, _dropped = sciabola2013.build_rows(sequences, values, fake_transcripts)
    assert [rows[0][c] for c in sciabola2013.DOSE_COLUMNS] == [10.0, 20.0, 30.0, 40.0]
    # only the top dose was measured for compound 2 -- the rest stay empty,
    # not zero, so a caller can tell "not measured" from "no knockdown"
    assert [rows[1][c] for c in sciabola2013.DOSE_COLUMNS] == [None, None, None, 80.0]


def test_fetch_writes_csv_and_fasta(tmp_path: Path, fake_document, fake_transcripts):
    sciabola2013.fetch(tmp_path)

    rows = list(csv.DictReader(open(tmp_path / "sciabola2013_extra.csv")))
    assert len(rows) == 2
    assert rows[0]["Sense_21mer"] == SENSE
    assert rows[1]["Pct_Inhibition_10nM"] == "80.0"
    assert rows[1]["Pct_Inhibition_008nM"] == ""

    fasta = (tmp_path / "sciabola2013_transcripts.fasta").read_text()
    assert fasta.startswith(f">{HIF1A}\n")   # only accessions actually used
    assert fasta.count(">") == 1


def test_fetch_document_picks_the_right_member_out_of_the_zip(
    monkeypatch: pytest.MonkeyPatch):
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr("gks1191f1p.jpg", b"not it")
        archive.writestr("supp_gks1191_nar-01814-n-2012-File002.doc", b"the document")
    monkeypatch.setattr(sciabola2013, "_download", lambda _url: buffer.getvalue())
    assert sciabola2013.fetch_document() == b"the document"


def test_fetch_document_raises_when_the_member_is_absent(monkeypatch: pytest.MonkeyPatch):
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr("gks1191f1p.jpg", b"not it")
    monkeypatch.setattr(sciabola2013, "_download", lambda _url: buffer.getvalue())
    with pytest.raises(FileNotFoundError, match="File002.doc"):
        sciabola2013.fetch_document()


def test_parse_tables_raises_when_a_table_is_missing(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(sciabola2013._ole, "word_text", lambda _blob: "no tables here")
    with pytest.raises(ValueError, match="no sequence and/or value table"):
        sciabola2013.parse_tables(b"fake-doc")
