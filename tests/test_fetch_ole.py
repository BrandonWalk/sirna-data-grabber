"""The legacy-binary readers behind the Sciabola 2013 and Harborth 2003
fetchers (`sirna_data.fetch._ole`).

Building a whole synthetic .doc/.xls in a fixture would test little beyond
the fixture, so these cover the parts that actually encode business rules:
Word's table geometry (an empty cell is indistinguishable from a row end by
delimiter alone), Excel's packed numbers, and its shared-string table --
including the CONTINUE spill that is the format's classic trap.
"""
from __future__ import annotations

import struct

import pytest

from sirna_data.fetch import _ole
from sirna_data.fetch._ole import CELL_MARK


def _row(cells: list[str]) -> str:
    return CELL_MARK.join(cells) + CELL_MARK + CELL_MARK


class TestWordTableGrid:
    def test_recovers_rows_from_the_column_count(self):
        block = _row(["a", "b", "c"]) + _row(["d", "e", "f"])
        assert _ole.word_table_grid(block, 3) == [["a", "b", "c"], ["d", "e", "f"]]

    def test_empty_cells_do_not_end_a_row(self):
        # the whole reason rows are counted rather than split: "b" missing
        # leaves two adjacent marks, exactly like a row end
        block = _row(["a", "", "c"]) + _row(["", "e", ""])
        assert _ole.word_table_grid(block, 3) == [["a", "", "c"], ["", "e", ""]]

    def test_skip_cells_drops_a_merged_header_row(self):
        block = CELL_MARK.join(["merged"]) + CELL_MARK + CELL_MARK + _row(["a", "b"])
        assert _ole.word_table_grid(block, 2, skip_cells=2) == [["a", "b"]]

    def test_cells_are_stripped(self):
        assert _ole.word_table_grid(_row([" a ", "b\t"]), 2) == [["a", "b"]]


class TestExcelNumbers:
    @pytest.mark.parametrize(
        "rk, expected",
        [
            (0b1000_1110, 35.0),            # integer, no cents flag
            (0b1000_1111, 0.35),            # integer, cents flag -> /100
            ((struct.unpack("<I", struct.pack("<d", 1.5)[4:])[0] & 0xFFFFFFFC), 1.5),
        ],
    )
    def test_rk_decoding(self, rk, expected):
        assert _ole._rk_value(rk) == pytest.approx(expected)

    def test_negative_integers_round_trip(self):
        assert _ole._rk_value(((-7 & 0x3FFFFFFF) << 2) | 0x02) == -7.0


class TestSharedStrings:
    @staticmethod
    def _string(text: str, wide: bool = False) -> bytes:
        payload = text.encode("utf-16-le" if wide else "cp1252")
        return struct.pack("<HB", len(text), 1 if wide else 0) + payload

    def test_reads_plain_and_wide_strings(self):
        body = self._string("Harborth") + self._string("Lamin A", wide=True)
        payload = struct.pack("<II", 2, 2) + body
        assert _ole._shared_strings(payload, []) == ["Harborth", "Lamin A"]

    def test_reads_a_string_split_across_a_continue_record(self):
        # "GAGCUCC" broken mid-string: the CONTINUE record restates the width
        payload = struct.pack("<II", 1, 1) + struct.pack("<HB", 7, 0) + b"GAG"
        continuation = b"\x00" + b"CUCC"
        assert _ole._shared_strings(payload, [continuation]) == ["GAGCUCC"]

    def test_stops_cleanly_when_data_runs_out(self):
        payload = struct.pack("<II", 5, 5) + self._string("only one")
        assert _ole._shared_strings(payload, []) == ["only one"]


class TestContainerErrors:
    def test_rejects_non_ole_bytes(self):
        with pytest.raises(_ole.UnsupportedFileError, match="not an OLE2"):
            _ole.streams(b"PK\x03\x04 this is a zip")

    def test_word_text_rejects_a_document_without_the_word_stream(
        self, monkeypatch: pytest.MonkeyPatch):
        monkeypatch.setattr(_ole, "streams", lambda _data: {"Workbook": b""})
        with pytest.raises(_ole.UnsupportedFileError, match="not a Word 97"):
            _ole.word_text(b"whatever")

    def test_xls_cells_rejects_a_file_without_a_workbook_stream(
        self, monkeypatch: pytest.MonkeyPatch):
        monkeypatch.setattr(_ole, "streams", lambda _data: {"WordDocument": b""})
        with pytest.raises(_ole.UnsupportedFileError, match="not an Excel 97"):
            _ole.xls_cells(b"whatever")
