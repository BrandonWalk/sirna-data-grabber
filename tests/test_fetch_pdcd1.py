"""Fetcher for the REMOVED REMOVED panel.

The upstream CSV's header labels are undocumented, so the parser identifies
columns by content; these tests pin that behaviour, including the shapes it
must refuse rather than guess at.
"""
from __future__ import annotations

import csv
from pathlib import Path

import pytest

from sirna_data.fetch import REMOVED

SEQUENCES = [
    "CGGAGAGCTTCGTGCTAAA", "GTGGTTGGTAGCGACAGAA", "GCCTGAGGGCTGGCTGTAA",
    "TCGTCTGGGCGGTGCTACA", "CCAGGATGGTTCTTAGACT", "GTCACCTGGGCCATGGCCA",
    "CAGTTCCAAACCCTGGTGG", "TGCAGCTTCTCCAACACAT",
]
TRANSCRIPT = "AAAA" + "".join(SEQUENCES) + "TTTT"

UPSTREAM = "\n".join(
    ["id,sequence,luc,qpcr"]
    + [f"163-{i},{seq},0.97{i},0.96{i}" for i, seq in enumerate(SEQUENCES, start=1)]
)


@pytest.fixture
def no_network(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(REMOVED, "_download", lambda _url: UPSTREAM.encode())
    monkeypatch.setattr(REMOVED, "fetch_transcript", lambda: TRANSCRIPT)


def test_parse_panel_finds_columns_by_content():
    panel = REMOVED.parse_panel(UPSTREAM)
    assert len(panel) == len(SEQUENCES)
    assert [row["Sequence"] for row in panel] == SEQUENCES
    assert panel[0]["Experiment_ID"] == "163-1"
    assert panel[0]["Gene"] == "REMOVED"
    assert panel[0]["Accession_number"] == "NM_005018"


def test_fractional_efficiencies_become_percentages():
    panel = REMOVED.parse_panel(UPSTREAM)
    assert panel[0]["Efficiency_LUC_Pct"] == pytest.approx(97.1)
    assert panel[0]["Efficiency_QPCR_Pct"] == pytest.approx(96.1)


def test_values_already_in_percent_are_left_alone():
    text = "id,sequence,luc,qpcr\n163-1," + SEQUENCES[0] + ",97.3,96.5"
    panel = REMOVED.parse_panel(text)
    assert panel[0]["Efficiency_LUC_Pct"] == pytest.approx(97.3)


def test_refuses_a_file_with_no_sequence_column():
    with pytest.raises(ValueError, match="looks like an siRNA sequence"):
        REMOVED.parse_panel("id,value\n1,0.5\n")


def test_refuses_a_file_without_two_efficiency_columns():
    text = "id,sequence,luc\n163-1," + SEQUENCES[0] + ",0.97"
    with pytest.raises(ValueError, match="expected two efficiency columns"):
        REMOVED.parse_panel(text)


def test_fetch_writes_csv_and_fasta(tmp_path: Path, no_network):
    REMOVED.fetch(tmp_path)

    rows = list(csv.DictReader(open(tmp_path / "REMOVED")))
    assert len(rows) == 8
    assert list(rows[0]) == REMOVED.COLUMNS
    assert rows[0]["Sequence"] == SEQUENCES[0]

    fasta = (tmp_path / "REMOVED").read_text()
    assert fasta.startswith(">NM_005018\n")
    assert "".join(fasta.splitlines()[1:]) == TRANSCRIPT


def test_fetch_refuses_a_panel_of_the_wrong_size(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    short = "\n".join(UPSTREAM.splitlines()[:4])
    monkeypatch.setattr(REMOVED, "_download", lambda _url: short.encode())
    monkeypatch.setattr(REMOVED, "fetch_transcript", lambda: TRANSCRIPT)
    with pytest.raises(ValueError, match="expected 8 siRNAs"):
        REMOVED.fetch(tmp_path)


def test_fetch_refuses_sequences_that_are_not_in_the_transcript(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(REMOVED, "_download", lambda _url: UPSTREAM.encode())
    monkeypatch.setattr(REMOVED, "fetch_transcript", lambda: "ACGTACGTACGTACGTACGT")
    with pytest.raises(ValueError, match="not in NM_005018"):
        REMOVED.fetch(tmp_path)
