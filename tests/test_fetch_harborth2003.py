"""Fetcher for Harborth et al. 2003's lamin A/C panel (via Ichihara 2007).

The Excel decoding itself is covered in tests/test_fetch_ole.py; here the
workbook is replaced by a synthetic cell grid so these tests stay about the
fetcher's own logic: picking the file out of the zip, filtering to the right
study, verifying against the transcript, and writing the two output files.
"""
from __future__ import annotations

import csv
import io
import zipfile
from pathlib import Path

import pytest

from sirna_data.fetch import harborth2003

SENSE = "GGAGGACCUGCAGGAGCUC"        # 19nt, present in the fake transcript
SENSE_2 = "CUGGACUUCCAGAAGAACA"
TRANSCRIPT = "AAAAA" + SENSE + "CCCCC" + SENSE_2 + "GGGGG"

HEADER = ["Dataset and number", "Training", "Subset", "Spare", "Authors", "Gene",
          "Accession number", "Antisense, 21 mer", "Sense, 19 mer", "% Inhibition"]


def _grid(rows: list[list[object]]) -> dict[tuple[int, int], object]:
    return {(r, c): value for r, row in enumerate(rows) for c, value in enumerate(row)}


@pytest.fixture
def fake_workbook(monkeypatch: pytest.MonkeyPatch):
    cells = _grid([
        HEADER,
        ["A1", "x", "x", "", "Huesken", "TC10", "BD135193", "AAAA", "CCCC", 46.2],
        ["B1", "", "", "", "Harborth", "Lamin A", "AH001498", SENSE + "uu", SENSE, 83.0],
        ["B2", "", "", "", "Harborth", "Lamin A", "AH001498", SENSE_2 + "uu", SENSE_2, 98.0],
    ])
    monkeypatch.setattr(harborth2003._ole, "xls_cells", lambda _blob: cells)
    return cells


@pytest.fixture
def no_network(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(harborth2003, "fetch_workbook", lambda: b"fake-xls")
    monkeypatch.setattr(harborth2003, "fetch_transcript", lambda: TRANSCRIPT)


def test_parse_panel_keeps_only_harborth_rows(fake_workbook):
    panel = harborth2003.parse_panel(b"fake-xls")
    assert [row["sense"] for row in panel] == [SENSE, SENSE_2]
    assert {row["authors"] for row in panel} == {"Harborth"}
    assert panel[0]["source_accession"] == "AH001498"
    assert panel[0]["pct_inhibition"] == 83.0


def test_parse_panel_raises_when_a_column_is_missing(monkeypatch: pytest.MonkeyPatch):
    cells = _grid([["Authors", "Gene"], ["Harborth", "Lamin A"]])
    monkeypatch.setattr(harborth2003._ole, "xls_cells", lambda _blob: cells)
    with pytest.raises(ValueError, match="missing expected column"):
        harborth2003.parse_panel(b"fake-xls")


def test_parse_panel_raises_when_the_study_is_absent(monkeypatch: pytest.MonkeyPatch):
    cells = _grid([HEADER, ["A1", "", "", "", "Huesken", "TC10", "X", "A", "C", 1.0]])
    monkeypatch.setattr(harborth2003._ole, "xls_cells", lambda _blob: cells)
    with pytest.raises(ValueError, match="no rows with Authors"):
        harborth2003.parse_panel(b"fake-xls")


def test_fetch_writes_csv_and_fasta(tmp_path: Path, fake_workbook, no_network):
    harborth2003.fetch(tmp_path)

    rows = list(csv.DictReader(open(tmp_path / "harborth2003_extra.csv")))
    assert len(rows) == 2
    assert [row["Compound_Name"] for row in rows] == ["B1", "B2"]
    first = rows[0]
    assert first["Gene"] == "Lamin A"
    assert first["Accession_number"] == harborth2003.ACCESSION   # the RefSeq located in
    assert first["Source_Accession"] == "AH001498"               # the paper's own accession
    assert first["Antisense_21mer"] == SENSE + "uu"              # both strands from the table
    assert first["Sense_19mer"] == SENSE
    assert float(first["Pct_Inhibition"]) == 83.0
    assert first["Cell"] == "HeLa"

    fasta = (tmp_path / "harborth2003_transcripts.fasta").read_text()
    assert fasta.startswith(f">{harborth2003.ACCESSION}\n")
    assert "".join(fasta.splitlines()[1:]) == TRANSCRIPT.replace("U", "T")


def test_fetch_raises_when_a_sense_strand_is_not_in_the_transcript(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, fake_workbook):
    monkeypatch.setattr(harborth2003, "fetch_workbook", lambda: b"fake-xls")
    monkeypatch.setattr(harborth2003, "fetch_transcript", lambda: "ACGUACGUACGU")
    with pytest.raises(ValueError, match="not in NM_170707"):
        harborth2003.fetch(tmp_path)


def test_fetch_workbook_picks_the_right_member_out_of_the_zip(
    monkeypatch: pytest.MonkeyPatch):
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr("figure1.jpg", b"not it")
        archive.writestr("nar_gkm699_nar-01617-met-s-2007-File007.xls", b"the workbook")
    monkeypatch.setattr(harborth2003, "_download", lambda _url: buffer.getvalue())
    assert harborth2003.fetch_workbook() == b"the workbook"


def test_fetch_workbook_raises_when_the_member_is_absent(monkeypatch: pytest.MonkeyPatch):
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr("figure1.jpg", b"not it")
    monkeypatch.setattr(harborth2003, "_download", lambda _url: buffer.getvalue())
    with pytest.raises(FileNotFoundError, match="File007.xls"):
        harborth2003.fetch_workbook()
