"""Readers for the two legacy Microsoft binary formats this package's
supplementary sources ship in: Word 97 (`.doc`) and Excel 97 (`.xls`).

Both sit inside an OLE2 compound file, so the container reader is shared.
This module exists so `sirna-data-fetch` can rebuild every source from its
primary files with **no dependency beyond pandas and the standard library**
(see pyproject.toml) -- LibreOffice, antiword, xlrd and openpyxl would each
do the job, but none of them can be a dependency of a pip-installable
package that promises to work from a bare `pip install`.

Only what these two sources actually need is implemented: whole-stream reads
from the container, Word's text stream with table cell/row marks intact, and
the value-bearing Excel cell records. Anything else (formatting, formulas,
encryption, revision marks) is ignored by design.
"""
from __future__ import annotations

import struct

_OLE_SIGNATURE = b"\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1"
_END_OF_CHAIN = 0xFFFFFFFE
_FREE_SECTOR = 0xFFFFFFFF

# Word cell/row marks, kept in the text stream so table structure survives.
CELL_MARK = "\x07"


class UnsupportedFileError(ValueError):
    """The bytes handed in are not the legacy format they were meant to be."""


def streams(data: bytes) -> dict[str, bytes]:
    """Every named stream in an OLE2 compound file, as {name: bytes}."""
    if data[:8] != _OLE_SIGNATURE:
        raise UnsupportedFileError("not an OLE2 compound file")
    sector_size = 1 << struct.unpack_from("<H", data, 30)[0]
    mini_size = 1 << struct.unpack_from("<H", data, 32)[0]
    (n_fat, dir_start, _cutoff, mini_start, _n_mini,
     difat_start, n_difat) = struct.unpack_from("<II4xIIIII", data, 44)

    def sector(i: int) -> bytes:
        return data[(i + 1) * sector_size: (i + 2) * sector_size]

    difat = list(struct.unpack_from("<109I", data, 76))
    nxt = difat_start
    for _ in range(n_difat):
        if nxt >= _END_OF_CHAIN:
            break
        block = sector(nxt)
        difat += list(struct.unpack_from(f"<{sector_size // 4 - 1}I", block, 0))
        nxt = struct.unpack_from("<I", block, sector_size - 4)[0]

    fat: list[int] = []
    for s in difat[:n_fat]:
        if s >= _END_OF_CHAIN:
            continue
        fat += list(struct.unpack_from(f"<{sector_size // 4}I", sector(s), 0))

    def chain(start: int, table: list[int]) -> list[int]:
        out: list[int] = []
        seen: set[int] = set()
        s = start
        while s < _END_OF_CHAIN and s not in seen:
            seen.add(s)
            out.append(s)
            s = table[s] if s < len(table) else _END_OF_CHAIN
        return out

    directory = b"".join(sector(s) for s in chain(dir_start, fat))
    entries = []
    for off in range(0, len(directory) - 127, 128):
        raw = directory[off: off + 128]
        name_len = struct.unpack_from("<H", raw, 64)[0]
        name = raw[: max(0, name_len - 2)].decode("utf-16-le", "ignore")
        entry_type = raw[66]
        start, size = struct.unpack_from("<I", raw, 116)[0], struct.unpack_from("<I", raw, 120)[0]
        entries.append((name, entry_type, start, size))

    root = next((e for e in entries if e[1] == 5), None)
    if root is None:
        raise UnsupportedFileError("OLE2 file has no root directory entry")
    mini_stream = b"".join(sector(s) for s in chain(root[2], fat))[: root[3]]

    mini_fat: list[int] = []
    if mini_start < _END_OF_CHAIN:
        for s in chain(mini_start, fat):
            mini_fat += list(struct.unpack_from(f"<{sector_size // 4}I", sector(s), 0))

    out: dict[str, bytes] = {}
    for name, entry_type, start, size in entries:
        if entry_type != 2:  # streams only, not storages
            continue
        if size < _cutoff:
            blob = b"".join(
                mini_stream[s * mini_size: (s + 1) * mini_size] for s in chain(start, mini_fat)
            )
        else:
            blob = b"".join(sector(s) for s in chain(start, fat))
        out[name] = blob[:size]
    return out


def word_text(data: bytes) -> str:
    """The text of a Word 97 `.doc`, with table cell marks (\\x07) intact.

    Walks the FIB to the piece table, so it handles both the 8-bit (CP1252)
    and UTF-16 pieces Word mixes within one document.
    """
    found = streams(data)
    if "WordDocument" not in found:
        raise UnsupportedFileError("no WordDocument stream -- not a Word 97 .doc")
    document = found["WordDocument"]
    flags = struct.unpack_from("<H", document, 10)[0]
    table_name = "1Table" if flags & 0x0200 else "0Table"
    if table_name not in found:
        raise UnsupportedFileError(f"missing {table_name} stream")
    table = found[table_name]

    clx_start, clx_len = struct.unpack_from("<II", document, 0x01A2)
    clx = table[clx_start: clx_start + clx_len]
    i = 0
    while i < len(clx) and clx[i] == 1:  # skip the Prc (formatting) blocks
        i += 3 + struct.unpack_from("<h", clx, i + 1)[0]
    if i >= len(clx) or clx[i] != 2:
        raise UnsupportedFileError("no piece table in the Word document")
    piece_len = struct.unpack_from("<I", clx, i + 1)[0]
    pieces = clx[i + 5: i + 5 + piece_len]

    n_pieces = (piece_len - 4) // 12
    positions = list(struct.unpack_from(f"<{n_pieces + 1}I", pieces, 0))
    chunks = []
    for p in range(n_pieces):
        fc = struct.unpack_from("<I", pieces, 4 * (n_pieces + 1) + 8 * p + 2)[0]
        n_chars = positions[p + 1] - positions[p]
        if fc & 0x40000000:  # 8-bit (CP1252) piece
            start = (fc & ~0x40000000) // 2
            chunks.append(document[start: start + n_chars].decode("cp1252", "replace"))
        else:
            chunks.append(document[fc: fc + 2 * n_chars].decode("utf-16-le", "replace"))
    return "".join(chunks)


def word_table_grid(block: str, n_columns: int, skip_cells: int = 0) -> list[list[str]]:
    """Rows of a Word table, from one `\\x07`-delimited text block.

    Word writes a cell mark after every cell AND a row mark after the last
    one, and an empty cell is just two adjacent marks -- so rows cannot be
    found by splitting on a delimiter. They are recovered instead from the
    table's known column count: every row is `n_columns + 1` marks wide.
    `skip_cells` drops a leading merged header row.
    """
    cells = [c.strip() for c in block.split(CELL_MARK)][skip_cells:]
    stride = n_columns + 1
    return [cells[i: i + n_columns] for i in range(0, len(cells) - n_columns + 1, stride)]


def _rk_value(rk: int) -> float:
    """Decode Excel's packed RK number."""
    value = (
        float(rk >> 2 if (rk >> 2) < 2 ** 29 else (rk >> 2) - 2 ** 30)
        if rk & 0x02
        else struct.unpack("<d", struct.pack("<Q", (rk & 0xFFFFFFFC) << 32))[0]
    )
    return value / 100 if rk & 0x01 else value


def _shared_strings(payload: bytes, continuations: list[bytes]) -> list[str]:
    """Excel's shared string table, which spills across CONTINUE records."""
    _total, unique = struct.unpack_from("<II", payload, 0)
    buffer, index = payload, 8
    pending = list(continuations)
    out: list[str] = []

    def next_buffer() -> bool:
        nonlocal buffer, index
        if not pending:
            return False
        buffer, index = pending.pop(0), 0
        return True

    while len(out) < unique:
        while index + 3 > len(buffer):
            if not next_buffer():
                return out
        length, grbit = struct.unpack_from("<HB", buffer, index)
        index += 3
        n_runs = 0
        ext_size = 0
        if grbit & 0x08:  # rich text
            n_runs = struct.unpack_from("<H", buffer, index)[0]
            index += 2
        if grbit & 0x04:  # far-east extension
            ext_size = struct.unpack_from("<I", buffer, index)[0]
            index += 4
        wide = grbit & 0x01
        parts, remaining = [], length
        while remaining:
            width = 2 if wide else 1
            take = min((len(buffer) - index) // width, remaining)
            raw = buffer[index: index + take * width]
            parts.append(raw.decode("utf-16-le" if wide else "cp1252", "replace"))
            index += take * width
            remaining -= take
            if remaining:
                if not next_buffer():
                    break
                wide = buffer[index] & 0x01  # each CONTINUE restates the width
                index += 1
        index += n_runs * 4 + ext_size
        out.append("".join(parts))
    return out


def xls_cells(data: bytes) -> dict[tuple[int, int], str | float]:
    """{(row, column): value} for an Excel 97 `.xls` workbook's cells.

    Only the value-bearing records are read (LABELSST, NUMBER, RK, MULRK);
    formatting, formulas and everything else are skipped.
    """
    found = streams(data)
    if "Workbook" not in found and "Book" not in found:
        raise UnsupportedFileError("no Workbook stream -- not an Excel 97 .xls")
    book = found.get("Workbook") or found["Book"]

    records: list[tuple[int, bytes]] = []
    pos = 0
    while pos + 4 <= len(book):
        number, size = struct.unpack_from("<HH", book, pos)
        records.append((number, book[pos + 4: pos + 4 + size]))
        pos += 4 + size

    strings: list[str] = []
    for i, (number, payload) in enumerate(records):
        if number != 0x00FC:  # SST
            continue
        continuations = []
        for next_number, next_payload in records[i + 1:]:
            if next_number != 0x003C:  # CONTINUE
                break
            continuations.append(next_payload)
        strings = _shared_strings(payload, continuations)
        break

    cells: dict[tuple[int, int], str | float] = {}
    for number, payload in records:
        if number == 0x00FD:  # LABELSST
            row, column, _xf, index = struct.unpack_from("<HHHI", payload, 0)
            cells[(row, column)] = strings[index] if index < len(strings) else ""
        elif number == 0x0203:  # NUMBER
            row, column, _xf, value = struct.unpack_from("<HHHd", payload, 0)
            cells[(row, column)] = value
        elif number == 0x027E:  # RK
            row, column, _xf, rk = struct.unpack_from("<HHHI", payload, 0)
            cells[(row, column)] = _rk_value(rk)
        elif number == 0x00BD:  # MULRK
            row, first_column = struct.unpack_from("<HH", payload, 0)
            for k in range((len(payload) - 6) // 6):
                rk = struct.unpack_from("<I", payload, 4 + k * 6 + 2)[0]
                cells[(row, first_column + k)] = _rk_value(rk)
    return cells
