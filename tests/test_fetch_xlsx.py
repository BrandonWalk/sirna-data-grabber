"""The stdlib .xlsx reader behind the Davis 2025 fetcher."""
from __future__ import annotations

import io
import zipfile

import pytest

from sirna_data.fetch import _xlsx

RELS = """<?xml version="1.0"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rId1" Target="worksheets/sheet1.xml"/></Relationships>"""
WORKBOOK = """<?xml version="1.0"?><workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"
 xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
<sheets><sheet name="Table 1" sheetId="1" r:id="rId1"/></sheets></workbook>"""
SHARED = """<?xml version="1.0"?><sst xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
<si><t>Compound</t></si><si><t>APP</t></si></sst>"""
SHEET = """<?xml version="1.0"?><worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
<sheetData>
<row r="1"><c r="A1" t="s"><v>0</v></c><c r="C1" t="inlineStr"><is><t>Inline</t></is></c></row>
<row r="2"><c r="A2" t="s"><v>1</v></c><c r="B2"><v>37.5</v></c></row>
</sheetData></worksheet>"""


def _workbook(**overrides: str) -> bytes:
    members = {
        "xl/_rels/workbook.xml.rels": RELS,
        "xl/workbook.xml": WORKBOOK,
        "xl/sharedStrings.xml": SHARED,
        "xl/worksheets/sheet1.xml": SHEET,
    }
    members.update(overrides)
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        for name, content in members.items():
            archive.writestr(name, content)
    return buffer.getvalue()


def test_reads_shared_inline_and_numeric_cells():
    rows = _xlsx.sheet_rows(_workbook(), "Table 1")
    assert rows[0] == ["Compound", "", "Inline"]   # B1 is absent -> padded, not shifted
    assert rows[1] == ["APP", "37.5"]


def test_column_letters_place_values_at_their_own_index():
    assert _xlsx._column_index("A1") == 0
    assert _xlsx._column_index("B2") == 1
    assert _xlsx._column_index("AA10") == 26
    assert _xlsx._column_index("AB7") == 27


def test_missing_sheet_names_what_is_available():
    with pytest.raises(_xlsx.UnsupportedWorkbookError, match="no sheet named 'Nope'"):
        _xlsx.sheet_rows(_workbook(), "Nope")


def test_rejects_something_that_is_not_a_workbook():
    with pytest.raises(_xlsx.UnsupportedWorkbookError, match="not a zip"):
        _xlsx.sheet_rows(b"\xd0\xcf\x11\xe0 an old .xls", "Table 1")

    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr("hello.txt", "not a workbook")
    with pytest.raises(_xlsx.UnsupportedWorkbookError, match="no xl/workbook.xml"):
        _xlsx.sheet_rows(buffer.getvalue(), "Table 1")


def test_works_without_a_shared_string_table():
    sheet = SHEET.replace('t="s"><v>0</v>', '><v>1</v>').replace('t="s"><v>1</v>', '><v>2</v>')
    blob = _workbook(**{"xl/worksheets/sheet1.xml": sheet})
    rows = _xlsx.sheet_rows(blob, "Table 1")
    assert rows[0][0] == "1"
