"""A reader for the modern Excel format (`.xlsx`), which Davis et al. 2025
ships its supplementary tables in.

Like `_ole.py`, this exists so `sirna-data-fetch` needs nothing beyond
pandas and the standard library -- openpyxl would do the job but cannot be a
dependency of a package that promises to reconstruct the dataset from a bare
`pip install`. An .xlsx is a zip of XML, so the standard library is enough.

Only cell values are read. Blank cells matter here (a row's columns must
line up with its header), so cells are placed by their column letter rather
than by the order they happen to appear in the file.
"""
from __future__ import annotations

import io
import zipfile
from xml.etree import ElementTree as ET

_MAIN = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_RELS = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"


class UnsupportedWorkbookError(ValueError):
    """The bytes handed in are not a readable .xlsx workbook."""


def _column_index(reference: str) -> int:
    """'A' -> 0, 'B' -> 1, ... 'AA' -> 26, from a cell reference like 'AB12'."""
    index = 0
    for char in reference:
        if not char.isalpha():
            break
        index = index * 26 + (ord(char.upper()) - ord("A") + 1)
    return index - 1


def sheet_rows(blob: bytes, sheet_name: str) -> list[list[str]]:
    """Every row of one named worksheet, as lists of strings.

    Rows are padded so that a value always sits at its own column index --
    an empty cell reads as "" rather than shifting its neighbours left.
    """
    try:
        archive = zipfile.ZipFile(io.BytesIO(blob))
    except zipfile.BadZipFile as error:
        raise UnsupportedWorkbookError("not a .xlsx workbook (not a zip)") from error

    with archive:
        names = archive.namelist()
        if "xl/workbook.xml" not in names:
            raise UnsupportedWorkbookError("zip has no xl/workbook.xml -- not a workbook")

        targets = {
            rel.get("Id"): rel.get("Target", "").lstrip("/").replace("xl/", "")
            for rel in ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        }
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        sheets = {
            str(sheet.get("name")): targets.get(str(sheet.get(f"{_RELS}id")), "")
            for sheet in workbook.find(f"{_MAIN}sheets") or []
        }
        if sheet_name not in sheets:
            raise UnsupportedWorkbookError(
                f"no sheet named {sheet_name!r} (found: {sorted(sheets)})"
            )

        shared: list[str] = []
        if "xl/sharedStrings.xml" in names:
            shared = [
                "".join(node.text or "" for node in item.iter(f"{_MAIN}t"))
                for item in ET.fromstring(archive.read("xl/sharedStrings.xml"))
            ]

        sheet = ET.fromstring(archive.read(f"xl/{sheets[sheet_name]}"))

    def value(cell: ET.Element) -> str:
        kind = cell.get("t")
        if kind == "inlineStr":
            return "".join(node.text or "" for node in cell.iter(f"{_MAIN}t"))
        node = cell.find(f"{_MAIN}v")
        if node is None or node.text is None:
            return ""
        if kind == "s":
            index = int(node.text)
            return shared[index] if index < len(shared) else ""
        return node.text

    rows: list[list[str]] = []
    for row in sheet.find(f"{_MAIN}sheetData") or []:
        cells: list[str] = []
        for cell in row:
            index = _column_index(cell.get("r", ""))
            if index < 0:
                index = len(cells)
            cells.extend([""] * (index - len(cells) + 1))
            cells[index] = value(cell)
        rows.append(cells)
    return rows
