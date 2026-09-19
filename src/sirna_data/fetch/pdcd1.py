"""Fetch the 8-siRNA REMOVED panel from the REMOVED repo (Xu, Zhao et al.
2024, GENE, doi:10.1016/j.gene.2024.148330).

**License unresolved**: the repo ships no LICENSE file, so it is
all-rights-reserved by default under GitHub's terms, and the associated
Elsevier paper is not confirmed open access. That is exactly why this data
is fetched rather than redistributed -- `REMOVED` is gitignored, and
running this fetcher is the user's own decision about their own copy. See
../../../NOTICE.md before relying on it.

The file was deleted from `main` (commit d2ad931) but is still present in
the parent commit, which is what this pins.

The upstream CSV's header labels are not documented anywhere, so columns are
identified by their content -- the sequence column by its alphabet, the two
efficiency columns by being numeric -- and the result is verified against
the real REMOVED transcript before anything is written.
"""
from __future__ import annotations

import csv
import io
import re
import urllib.parse
import urllib.request
from pathlib import Path

TABLE_URL = (
    "https://raw.githubusercontent.com/REMOVED/REMOVED/f3254ee/"
    "data/ExperimentalData/REMOVED_8.csv"
)
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

GENE = "REMOVED"
ACCESSION = "NM_005018"
EXPECTED_ROWS = 8
SEQUENCE_PATTERN = re.compile(r"^[ACGTUacgtu]{18,25}$")
COLUMNS = [
    "Experiment_ID", "Sequence", "Gene", "Accession_number",
    "Efficiency_LUC_Pct", "Efficiency_QPCR_Pct",
]


def _download(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=120) as response:
        return response.read()


def _as_percent(value: str) -> float:
    """The source reports fractional knockdown (0-1); this package uses %."""
    number = float(value)
    return number * 100 if abs(number) <= 1 else number


def parse_panel(text: str) -> list[dict[str, str | float]]:
    """The panel's rows, with columns identified by their content."""
    rows = [row for row in csv.reader(io.StringIO(text)) if any(cell.strip() for cell in row)]
    if not rows:
        raise ValueError(f"no rows in {TABLE_URL}")

    def is_number(cell: str) -> bool:
        try:
            float(cell)
        except ValueError:
            return False
        return True

    body = [row for row in rows if any(SEQUENCE_PATTERN.match(cell.strip()) for cell in row)]
    if not body:
        raise ValueError("no row contains anything that looks like an siRNA sequence")

    first = body[0]
    sequence_column = next(
        i for i, cell in enumerate(first) if SEQUENCE_PATTERN.match(cell.strip())
    )
    numeric_columns = [
        i for i, cell in enumerate(first) if i != sequence_column and is_number(cell.strip())
    ]
    if len(numeric_columns) < 2:
        raise ValueError(
            "expected two efficiency columns (luciferase and qPCR), "
            f"found {len(numeric_columns)} numeric column(s) in {first!r}"
        )
    luciferase, qpcr = numeric_columns[:2]
    identifier = next(
        (i for i in range(len(first)) if i not in (sequence_column, luciferase, qpcr)),
        None,
    )

    panel: list[dict[str, str | float]] = []
    for i, row in enumerate(body, start=1):
        panel.append({
            "Experiment_ID": row[identifier].strip() if identifier is not None else str(i),
            "Sequence": row[sequence_column].strip().upper(),
            "Gene": GENE,
            "Accession_number": ACCESSION,
            "Efficiency_LUC_Pct": _as_percent(row[luciferase]),
            "Efficiency_QPCR_Pct": _as_percent(row[qpcr]),
        })
    return panel


def fetch_transcript() -> str:
    query = urllib.parse.urlencode(
        {"db": "nuccore", "id": ACCESSION, "rettype": "fasta", "retmode": "text"}
    )
    text = _download(f"{EFETCH_URL}?{query}").decode()
    return "".join(line.strip() for line in text.splitlines() if not line.startswith(">"))


def fetch(dest: Path) -> None:
    """Write REMOVED and REMOVED."""
    dest.mkdir(parents=True, exist_ok=True)
    panel = parse_panel(_download(TABLE_URL).decode("utf-8-sig"))
    if len(panel) != EXPECTED_ROWS:
        raise ValueError(
            f"expected {EXPECTED_ROWS} siRNAs in the REMOVED panel, parsed {len(panel)} -- "
            "the upstream file has changed"
        )

    transcript = fetch_transcript().upper()
    rna = transcript.replace("T", "U")
    unlocated = [
        str(row["Sequence"]) for row in panel
        if str(row["Sequence"]).replace("T", "U") not in rna
    ]
    if unlocated:
        raise ValueError(
            f"{len(unlocated)} sequence(s) are not in {ACCESSION}: {unlocated[:3]}"
        )

    csv_path = dest / "REMOVED"
    with open(csv_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(panel)

    fasta_path = dest / "REMOVED"
    with open(fasta_path, "w") as handle:
        handle.write(f">{ACCESSION}\n")
        for i in range(0, len(transcript), 70):
            handle.write(transcript[i: i + 70] + "\n")

    print(f"  REMOVED: {len(panel)} rows -> {csv_path.name}, {fasta_path.name}")
