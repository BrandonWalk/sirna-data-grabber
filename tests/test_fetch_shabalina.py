from __future__ import annotations

import csv
import io
import zipfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from sirna_data.fetch import shabalina


def _fasta_bytes(entries: dict[str, str]) -> bytes:
    lines = []
    for acc, seq in entries.items():
        lines.append(f">{acc} description")
        lines.append(seq)
    return ("\n".join(lines) + "\n").encode()


def _mock_urlopen(payload: bytes):
    cm = MagicMock()
    cm.__enter__.return_value.read.return_value = payload
    return cm


def _s1a_row(accession, seq, activity, reference=""):
    """Build one 26-column raw source row: [acc, start, end, seq, ...19
    filler cols..., activity (idx 23), filler (idx 24), reference (idx 25)]."""
    return [accession, "1", "19", seq] + [""] * 19 + [activity] + [""] + [reference]


def _make_zip_payload(data_rows: list[list[str]]) -> bytes:
    buf = io.StringIO()
    writer = csv.writer(buf)
    writer.writerow(["header"] * 26)
    for row in data_rows:
        writer.writerow(row)
    zip_buf = io.BytesIO()
    with zipfile.ZipFile(zip_buf, "w") as zf:
        zf.writestr(shabalina.TABLE_S1A_FILENAME, buf.getvalue())
    return zip_buf.getvalue()


def test_revcomp_reverses_and_complements_rna_bases():
    assert shabalina.revcomp("ACGU") == "ACGU"  # palindromic under revcomp
    assert shabalina.revcomp("AAUUCCGG") == "CCGGAAUU"
    assert shabalina.revcomp("GGGG") == "CCCC"


def test_fetch_table_s1a_parses_zip_and_forward_fills_reference():
    rows = [
        _s1a_row("NM_000321", "AAAAAAAAAAAAAAAAAAA", "42.0", "Author et al. 2001"),
        _s1a_row("M25346", "CCCCCCCCCCCCCCCCCCC", "10.5"),  # no reference -> forward filled
        ["", "", "", "", ""],  # blank accession row, must be skipped
        # short row (fewer than 26 raw columns) -- exercises the padding
        # branch that pads ragged rows out to 26 columns before indexing
        # activity_remaining (23) / reference (25).
        ["M84918", "1", "19", "GGGGGGGGGGGGGGGGGGG", "extra", "extra2"],
    ]
    payload = _make_zip_payload(rows)
    with patch.object(shabalina.urllib.request, "urlopen", return_value=_mock_urlopen(payload)):
        parsed = shabalina.fetch_table_s1a()

    assert len(parsed) == 3
    assert parsed[2]["accession"] == "M84918"
    assert parsed[2]["antisense_19mer"] == "GGGGGGGGGGGGGGGGGGG"
    assert parsed[2]["activity_remaining"] == ""  # padded-in column, not present in source row
    assert parsed[2]["reference"] == "Author et al. 2001"  # still forward-filled from row 1
    assert parsed[0]["accession"] == "NM_000321"
    assert parsed[0]["antisense_19mer"] == "AAAAAAAAAAAAAAAAAAA"
    assert parsed[0]["activity_remaining"] == "42.0"
    assert parsed[0]["reference"] == "Author et al. 2001"
    # second row had no References entry -> carries forward the prior one
    assert parsed[1]["reference"] == "Author et al. 2001"


def test_fetch_transcripts_parses_fasta_and_strips_version():
    payload = _fasta_bytes({"NM_000321.3": "ACGT", "M25346.1": "TTTT"})
    with patch.object(
        shabalina.urllib.request, "urlopen", return_value=_mock_urlopen(payload)
    ) as mock_urlopen:
        result = shabalina.fetch_transcripts(["NM_000321", "M25346"])

    mock_urlopen.assert_called_once()
    assert result == {"NM_000321": "ACGT", "M25346": "TTTT"}


def test_fetch_filters_already_present_genes_and_applies_revcomp(tmp_path: Path, capsys):
    all_rows = [
        {
            "accession": "NM_000321",
            "start": "1",
            "end": "19",
            "antisense_19mer": "AAAAAAAAAAAAAAAAAAA",
            "activity_remaining": "42.0",
            "reference": "ref1",
        },
        {
            "accession": "J03132",  # maps to None (ICAM-1 duplicate) -> dropped
            "start": "1",
            "end": "19",
            "antisense_19mer": "CCCCCCCCCCCCCCCCCCC",
            "activity_remaining": "5.0",
            "reference": "ref2",
        },
    ]
    with patch.object(shabalina, "fetch_table_s1a", return_value=all_rows), patch.object(
        shabalina, "fetch_transcripts", return_value={"NM_000321": "ACGT"}
    ) as mock_fetch_transcripts:
        shabalina.fetch(tmp_path)

    mock_fetch_transcripts.assert_called_once_with(["NM_000321"])

    csv_path = tmp_path / "shabalina_extra.csv"
    fasta_path = tmp_path / "shabalina_transcripts.fasta"
    assert csv_path.exists() and fasta_path.exists()

    rows = list(csv.reader(csv_path.read_text().splitlines()))
    assert rows[0] == ["Sequence", "Gene", "Accession_number", "Activity_Remaining_Pct"]
    assert len(rows) == 2  # header + the one kept row (J03132 dropped)
    assert rows[1][0] == shabalina.revcomp("AAAAAAAAAAAAAAAAAAA")
    assert rows[1][1] == "RB1"
    assert rows[1][2] == "NM_000321"

    out = capsys.readouterr().out
    assert "Keeping 1 rows for 1 new genes" in out


def test_fetch_raises_on_unmapped_accession(tmp_path: Path):
    all_rows = [
        {
            "accession": "TOTALLY_UNKNOWN",
            "start": "1",
            "end": "19",
            "antisense_19mer": "AAAAAAAAAAAAAAAAAAA",
            "activity_remaining": "1.0",
            "reference": "ref",
        }
    ]
    with patch.object(shabalina, "fetch_table_s1a", return_value=all_rows):
        with pytest.raises(RuntimeError, match="Unmapped accessions"):
            shabalina.fetch(tmp_path)
