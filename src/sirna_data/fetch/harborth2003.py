"""Fetch Harborth et al. 2003's lamin A/C panel -- 44 standard 21-nt duplexes
tiling lamin A/C, measured as protein knockdown in HeLa cells.

Harborth et al. 2003 (Antisense Nucleic Acid Drug Dev. 13:83-105) is NOT open
access and carries no reuse license, so the values are taken from its
republication in Ichihara et al. 2007 (Nucleic Acids Research 35(18):e123,
doi:10.1093/nar/gkm699, PMC2094068, CC BY-NC 2.0 UK), whose supplementary
workbook tabulates every dataset behind its i-Score model with a per-row
`Authors` column. Retrieved through Europe PMC's public supplementaryFiles
API, the same documented route `monopoli.py` and `shabalina.py` use.

Cite both papers when using this subset -- Harborth for the experiments,
Ichihara for the compilation these numbers come from. See
../../../data/DATA_SOURCES.md.
"""
from __future__ import annotations

import csv
import io
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

from . import _ole

SUPPL_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC2094068/supplementaryFiles"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
WORKBOOK_SUFFIX = "File007.xls"

AUTHORS = "Harborth"  # the workbook's own per-row study label
GENE = "Lamin A"  # kept as the source spells it, to match the corpus's gene group
CELL_LINE = "HeLa"  # from Harborth et al. 2003's Methods; not a column in the table

# The table cites AH001498, a segmented GenBank gene record in which only 42 of
# the 44 sense strands can be located; all 44 are in this RefSeq mRNA, so that
# is what records carry (the cited accession is preserved per row).
ACCESSION = "NM_170707"

# Header labels in the workbook's first sheet -> our column names.
COLUMNS = {
    "Authors": "authors",
    "Gene": "gene",
    "Accession number": "source_accession",
    "Antisense, 21 mer": "antisense",
    "Sense, 19 mer": "sense",
    "% Inhibition": "pct_inhibition",
}


def _download(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=120) as response:
        return response.read()


def fetch_workbook() -> bytes:
    """The supplementary .xls, out of Europe PMC's zip of supplementary files."""
    with zipfile.ZipFile(io.BytesIO(_download(SUPPL_URL))) as archive:
        names = [n for n in archive.namelist() if n.endswith(WORKBOOK_SUFFIX)]
        if not names:
            raise FileNotFoundError(
                f"{WORKBOOK_SUFFIX} not in Europe PMC's supplementary files for "
                f"PMC2094068 (found: {sorted(archive.namelist())})"
            )
        return archive.read(names[0])


def parse_panel(workbook: bytes) -> list[dict[str, str | float | None]]:
    """The 44 Harborth rows of the workbook, keyed by our column names."""
    cells = _ole.xls_cells(workbook)
    grid: dict[int, dict[int, str | float]] = {}
    for (row, column), value in cells.items():
        grid.setdefault(row, {})[column] = value

    header_row = min(grid)
    header = {
        str(value).strip(): column
        for column, value in grid[header_row].items()
        if str(value).strip() in COLUMNS
    }
    missing = set(COLUMNS) - set(header)
    if missing:
        raise ValueError(f"workbook is missing expected column(s): {sorted(missing)}")

    panel: list[dict[str, str | float | None]] = []
    for row in sorted(grid):
        if row == header_row:
            continue
        values = grid[row]
        if str(values.get(header["Authors"], "")).strip() != AUTHORS:
            continue
        panel.append({
            name: values.get(header[label])
            for label, name in COLUMNS.items()
        })
    if not panel:
        raise ValueError(f"no rows with Authors == {AUTHORS!r} in the workbook")
    return panel


def fetch_transcript() -> str:
    """The RefSeq mRNA the panel's target sites are located in."""
    query = urllib.parse.urlencode(
        {"db": "nuccore", "id": ACCESSION, "rettype": "fasta", "retmode": "text"}
    )
    text = _download(f"{EFETCH_URL}?{query}").decode()
    return "".join(line.strip() for line in text.splitlines() if not line.startswith(">"))


def fetch(dest: Path) -> None:
    """Write harborth2003_extra.csv and harborth2003_transcripts.fasta."""
    dest.mkdir(parents=True, exist_ok=True)
    panel = parse_panel(fetch_workbook())
    transcript = fetch_transcript().upper().replace("T", "U")

    unlocated = [
        str(row["sense"]).upper() for row in panel
        if str(row["sense"]).upper() not in transcript
    ]
    if unlocated:
        raise ValueError(
            f"{len(unlocated)} of {len(panel)} sense strands are not in {ACCESSION}: "
            f"{unlocated[:3]} -- the source table or the transcript has changed"
        )

    csv_path = dest / "harborth2003_extra.csv"
    with open(csv_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "Compound_Name", "Gene", "Accession_number", "Antisense_21mer",
            "Sense_19mer", "Pct_Inhibition", "Source_Accession", "Cell",
        ])
        for i, row in enumerate(panel, start=1):
            inhibition = row["pct_inhibition"]
            if inhibition is None:
                raise ValueError(f"row {i} has no % inhibition value")
            writer.writerow([
                f"B{i}", GENE, ACCESSION, row["antisense"], row["sense"],
                float(inhibition), row["source_accession"], CELL_LINE,
            ])

    fasta_path = dest / "harborth2003_transcripts.fasta"
    with open(fasta_path, "w") as handle:
        handle.write(f">{ACCESSION}\n")
        for i in range(0, len(transcript), 70):
            handle.write(transcript[i: i + 70].replace("U", "T") + "\n")

    print(f"  harborth2003: {len(panel)} rows -> {csv_path.name}, {fasta_path.name}")
