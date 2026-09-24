"""Fetch Davis et al. 2025's fully chemically modified siRNA panel -- APP,
MAPT, BACE1 and SNCA, measured by QuantiGene 2.0 in SH-SY5Y cells.

Davis, Hildebrand, MacMillan, Monopoli, ... Pai & Khvorova 2025, Nucleic
Acids Research 53(12):gkaf479 (doi:10.1093/nar/gkaf479, PMC12205987),
**CC BY 4.0**. Supplemental Table 1, filtered to the paper's own "Included
in Filtered Dataset" subset.

Oxford's own supplementary CDN link for this article sits behind a
Cloudflare Turnstile challenge, which is why this data originally had to be
downloaded by hand. Europe PMC serves the same zip through its public
supplementaryFiles API with no challenge, so no bot-check is being worked
around here -- see ../../../data/DATA_SOURCES.md.
"""
from __future__ import annotations

import csv
import io
import urllib.request
import zipfile
from pathlib import Path

from . import _xlsx

SUPPL_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12205987/supplementaryFiles"
INNER_ZIP_SUFFIX = "gkaf479_supplemental_files.zip"
WORKBOOK_SUFFIX = "Supplemental Tables.xlsx"
SHEET = "Supplemental Table 1"

# Header labels in the sheet -> the columns this package stores. The two
# "Average (% Untreated Control)" columns are the native assay and the
# reporter assay, in that order, and are told apart by position.
HEADERS = {
    "Compound Name": "Compound_Name",
    "Target": "Gene",
    "NCBI Accession Used for siRNA Design": "Accession_number",
    "20mer Target Sequence": "Sense_20mer",
    "50mer Consensus Sequence Across mRNA Variants Expressed in SH-SY5Y cells": "MRNA_50mer_Window",
    "Dataset": "Dataset",
}
NATIVE_COLUMNS = ("Native_Avg_Pct_Untreated", "Native_STDEV")
REPORTER_COLUMNS = ("Reporter_Avg_Pct_Untreated", "Reporter_STDEV")
AVERAGE_LABEL = "Average (% Untreated Control)"
FILTER_COLUMN = "Included in Filtered Dataset"
FILTER_VALUE = "Yes"

COLUMNS = [
    "Compound_Name", "Gene", "Accession_number", "Sense_20mer", "MRNA_50mer_Window",
    *NATIVE_COLUMNS, *REPORTER_COLUMNS, "Dataset",
]


def _as_float(cell: str) -> str:
    """Excel keeps more digits than a float can hold; store the double."""
    return "" if not cell else repr(float(cell))


def _download(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=180) as response:
        return response.read()


def _member(archive: zipfile.ZipFile, suffix: str, where: str) -> bytes:
    names = [n for n in archive.namelist() if n.endswith(suffix)]
    if not names:
        raise FileNotFoundError(f"{suffix} not in {where} (found: {sorted(archive.namelist())})")
    return archive.read(names[0])


def fetch_workbook() -> bytes:
    """The supplementary workbook, inside a zip inside Europe PMC's zip."""
    with zipfile.ZipFile(io.BytesIO(_download(SUPPL_URL))) as outer:
        inner_blob = _member(outer, INNER_ZIP_SUFFIX, "Europe PMC's supplementary files")
    with zipfile.ZipFile(io.BytesIO(inner_blob)) as inner:
        return _member(inner, WORKBOOK_SUFFIX, INNER_ZIP_SUFFIX)


def parse_table(workbook: bytes) -> list[dict[str, str]]:
    """The filtered rows of Supplemental Table 1, keyed by this module's column names.

    The sheet carries several title/legend rows above the real header, so the
    header is found by content rather than by a fixed row number.
    """
    rows = _xlsx.sheet_rows(workbook, SHEET)
    header_index = next(
        (i for i, row in enumerate(rows) if "Compound Name" in row and "Target" in row),
        None,
    )
    if header_index is None:
        raise ValueError(f"no header row in {SHEET!r} -- the workbook's layout has changed")

    header = rows[header_index]
    position = {label: i for i, label in enumerate(header)}
    missing = [label for label in (*HEADERS, FILTER_COLUMN) if label not in position]
    if missing:
        raise ValueError(f"{SHEET!r} is missing expected column(s): {missing}")

    averages = [i for i, label in enumerate(header) if label == AVERAGE_LABEL]
    if len(averages) != 2:
        raise ValueError(
            f"expected the native and reporter {AVERAGE_LABEL!r} columns, found {len(averages)}"
        )
    native, reporter = averages  # the STDEV of each sits immediately to its right

    out: list[dict[str, str]] = []
    def cell(row: list[str], index: int) -> str:
        return row[index].strip() if index < len(row) else ""

    for row in rows[header_index + 1:]:
        if cell(row, position[FILTER_COLUMN]) != FILTER_VALUE:
            continue
        record = {name: cell(row, position[label]) for label, name in HEADERS.items()}
        # The sheet cites versioned accessions; this package keys on the
        # bare accession, as the other sources do.
        record["Accession_number"] = record["Accession_number"].split(".")[0]
        for columns, first in ((NATIVE_COLUMNS, native), (REPORTER_COLUMNS, reporter)):
            record[columns[0]] = _as_float(cell(row, first))
            record[columns[1]] = _as_float(cell(row, first + 1))
        out.append(record)
    if not out:
        raise ValueError(f"no rows with {FILTER_COLUMN} == {FILTER_VALUE!r}")
    return out


def fetch(dest: Path) -> None:
    """Write davis2025_extra.csv. This source ships its own mRNA window per
    row, so unlike every other fetcher here there is no transcript FASTA."""
    dest.mkdir(parents=True, exist_ok=True)
    table = parse_table(fetch_workbook())

    misplaced = [
        row["Compound_Name"] for row in table
        if row["Sense_20mer"] and row["Sense_20mer"] not in row["MRNA_50mer_Window"]
    ]
    if misplaced:
        raise ValueError(
            f"{len(misplaced)} target site(s) are not inside their own 50mer window "
            f"({misplaced[:3]}) -- the source table has changed"
        )

    csv_path = dest / "davis2025_extra.csv"
    with open(csv_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(table)
    print(f"  davis2025: {len(table)} rows -> {csv_path.name}")
