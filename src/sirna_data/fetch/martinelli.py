"""Fetch the gene-identified subset of the Martinelli 2023 /
`sirna-reproduction` siRNAmod-derived table -- 577 of its 907 rows.

Martinelli 2023 (bioRxiv 10.1101/2023.08.16.553554), CC BY-NC 4.0. The
upstream file is a 907-row table of chemically modified siRNAs with a
sequence, both strands' modifications and a % inhibition -- but **no gene or
accession column**, which is why only part of it is loadable here.

Gene identity is recovered mechanically: every row's 19nt sense core is
searched against the transcripts of the 13 accessions below, and a row is
kept only when it lands in exactly one of them. Those accessions are the
output of the one-off source-document tracing described in
../../../data/DATA_SOURCES.md; with them known, the assignment itself is
reproducible, which is what this fetcher does.
"""
from __future__ import annotations

import csv
import io
import urllib.parse
import urllib.request
from pathlib import Path

from ..raw_loader import _revcomp

TABLE_URL = (
    "https://raw.githubusercontent.com/mmartinelli-bio/sirna-reproduction/"
    "main/data/sirna_reproduction_martinelli_907.csv"
)
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
CORE_LENGTH = 19  # the duplex core expected to match the transcript

# accession -> gene/reporter label, recovered by tracing the source documents.
GENE_MAP = {
    "AF025846": "Luciferase_renilla",
    "NM_000384": "APOB",
    "NM_001009022": "SOD2_chimp",
    "NM_001025366": "VEGFA",
    "NM_001111035": "ACP5",
    "NM_004064": "CDKN1B",
    "NM_012614": "NPY",
    "NM_201628": "KAZRIN",
    "NR_001458": "MIR155HG",
    "U47296": "Luciferase_firefly",
    "U55763": "EGFP",
    "X65324": "Luciferase_firefly",
    "XM_001111972": "CASR_rhesus",
}

SOURCE_COLUMNS = {
    "id": "siRNAmodDB ID",
    "pmid": "PMID",
    "sense": "Sequence of sense strand",
    "sense_modification": "modifications (sense strand)",
    "antisense": "Sequence of antisense strand",
    "antisense_modification": "modifications (antisense strand)",
    "pct": "PCT",
}
COLUMNS = [
    "Experiment_ID", "PMID", "Gene", "Accession_number", "Sequence",
    "Modification_sense", "Sequence_antisense", "Modification_antisense", "PCT",
]


def _as_float(value: str) -> str:
    """Percentages as numbers, so "0" and "0.0" cannot both appear."""
    text = (value or "").strip()
    return repr(float(text)) if text else ""


SEQUENCE_FIELDS = ("Sequence", "Sequence_antisense")


def sort_key(row: dict[str, str]) -> tuple[int, str]:
    """Canonical row order: by the source's own experiment number.

    Deterministic and independent of the upstream file's ordering, so a
    re-fetch and a normalized hand-built file agree line for line.
    """
    identifier = str(row.get("Experiment_ID", "")).strip()
    digits = "".join(character for character in identifier if character.isdigit())
    return (int(digits) if digits else 0, identifier)


def canonical_row(row: dict[str, str]) -> dict[str, str]:
    """One row in this package's canonical form.

    The upstream table mixes conventions -- DNA overhangs written T on some
    rows and U on others, percentages as both "0" and "0.0" -- so every row
    that leaves this module goes through here: sequences in the RNA
    alphabet (the chemistry lives in the modification columns, not in the
    letters), percentages as numbers. Idempotent, so it can be applied to
    an already-canonical file without changing it.
    """
    out = {column: str(row.get(column, "")).strip() for column in COLUMNS}
    for field in SEQUENCE_FIELDS:
        out[field] = out[field].upper().replace("T", "U")
    out["PCT"] = _as_float(out["PCT"])
    return out


def normalize_file(path: Path, transcripts: dict[str, str] | None = None) -> int:
    """Rewrite an existing martinelli_extra.csv in canonical form.

    Returns the number of rows whose contents changed. This is what brings a
    hand-built file into line with what `fetch()` produces, so re-fetching
    is a no-op rather than a diff.

    Without `transcripts` this only canonicalizes formatting (alphabet and
    number formatting), which needs no inputs. Given them, each row's stored
    target site is also re-derived by the same longest-uniquely-matching-
    window rule `fetch()` uses, so a hand-picked shorter window is replaced
    by the canonical one. A row that re-places to a DIFFERENT gene is an
    error, not a silent rewrite.
    """
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"{path} has no rows")
    missing = [column for column in COLUMNS if column not in rows[0]]
    if missing:
        raise ValueError(f"{path} is missing expected column(s): {missing}")

    canonical = sorted((canonical_row(row) for row in rows), key=sort_key)
    if transcripts is not None:
        for row in canonical:
            placed = _place(row["Sequence"], row["Sequence_antisense"], transcripts)
            if placed is None:
                raise ValueError(
                    f"{row['Experiment_ID']} no longer matches any known transcript"
                )
            accession, stored = placed
            if accession != row["Accession_number"]:
                raise ValueError(
                    f"{row['Experiment_ID']} re-places to {accession}, not "
                    f"{row['Accession_number']} -- refusing to reassign a gene"
                )
            row["Sequence"] = stored
    before_by_id = {str(row["Experiment_ID"]).strip(): row for row in rows}
    changed = sum(
        1 for after in canonical
        if any(
            str(before_by_id.get(after["Experiment_ID"], {}).get(column, "")).strip()
            != after[column]
            for column in COLUMNS
        )
    )
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(canonical)
    return changed


def _download(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=120) as response:
        return response.read()


def fetch_table() -> list[dict[str, str]]:
    """The upstream 907-row table."""
    text = _download(TABLE_URL).decode("utf-8-sig")
    rows = list(csv.DictReader(io.StringIO(text)))
    if not rows:
        raise ValueError(f"no rows in {TABLE_URL}")
    missing = [label for label in SOURCE_COLUMNS.values() if label not in rows[0]]
    if missing:
        raise ValueError(f"upstream table is missing expected column(s): {missing}")
    return rows


def fetch_transcripts(accessions: list[str]) -> dict[str, str]:
    """{accession: RNA sequence} from NCBI, version stripped from the key."""
    query = urllib.parse.urlencode(
        {"db": "nuccore", "id": ",".join(accessions), "rettype": "fasta", "retmode": "text"}
    )
    text = _download(f"{EFETCH_URL}?{query}").decode()
    out: dict[str, str] = {}
    name = None
    for line in text.splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0].split(".")[0]
            out[name] = ""
        elif name:
            out[name] += line.strip().upper().replace("T", "U")
    return out


def _place(sense: str, antisense: str, transcripts: dict[str, str]) -> tuple[str, str] | None:
    """(accession, stored sense sequence) for a row, or None if unplaceable.

    The longest window of either strand that occurs in exactly ONE of the
    known transcripts decides the row's gene; a window matching two
    transcripts (the two firefly luciferase reporters overlap) is ambiguous
    and the row is dropped.

    What gets stored depends on which strand placed the row. Placed by its
    sense strand, the verified window is spliced back into that strand, so
    the row keeps whatever the source gives beyond the duplex core (a 2nt
    overhang, usually). Placed by its antisense strand -- the source's sense
    column is wrong or missing for those rows -- there is no sense strand to
    keep, so the target site is the reverse complement's last CORE_LENGTH
    bases: the 2nt that reverse-complement the guide's 3' overhang are not
    part of the target.

    Sequences are stored in the RNA alphabet, as everything else in this
    corpus is: the source writes DNA overhangs as T, but the chemistry lives
    in the modification columns, not in the letters.
    """
    literal = sense.strip().upper().replace("T", "U")
    candidates = (("sense", literal.replace("T", "U")),
                  ("antisense", _revcomp(antisense.strip().upper().replace("T", "U"))))
    for origin, candidate in candidates:
        for length in range(len(candidate), CORE_LENGTH - 1, -1):
            for start in range(0, len(candidate) - length + 1):
                window = candidate[start: start + length]
                hits = [a for a, seq in transcripts.items() if window in seq]
                if len(hits) != 1:
                    continue
                if origin == "sense":
                    return hits[0], literal[:start] + window + literal[start + length:]
                return hits[0], candidate[-CORE_LENGTH:]
    return None


def assign_genes(
    table: list[dict[str, str]], transcripts: dict[str, str]
) -> list[dict[str, str]]:
    """Rows whose sequence lands in exactly one of the known transcripts."""
    out: list[dict[str, str]] = []
    for row in table:
        sense = row[SOURCE_COLUMNS["sense"]] or ""
        antisense = row[SOURCE_COLUMNS["antisense"]] or ""
        if len(sense.strip()) < CORE_LENGTH:
            continue
        placed = _place(sense, antisense, transcripts)
        if placed is None:
            continue
        accession, stored = placed
        out.append(canonical_row({
            "Experiment_ID": row[SOURCE_COLUMNS["id"]],
            "PMID": row[SOURCE_COLUMNS["pmid"]],
            "Gene": GENE_MAP[accession],
            "Accession_number": accession,
            "Sequence": stored,
            "Modification_sense": row[SOURCE_COLUMNS["sense_modification"]],
            "Sequence_antisense": antisense,
            "Modification_antisense": row[SOURCE_COLUMNS["antisense_modification"]],
            "PCT": row[SOURCE_COLUMNS["pct"]],
        }))
    return out


def fetch(dest: Path) -> None:
    """Write martinelli_extra.csv and martinelli_transcripts.fasta."""
    dest.mkdir(parents=True, exist_ok=True)
    table = fetch_table()
    transcripts = fetch_transcripts(sorted(GENE_MAP))
    missing = [a for a in GENE_MAP if a not in transcripts]
    if missing:
        raise ValueError(f"NCBI returned no sequence for: {missing}")

    rows = sorted(assign_genes(table, transcripts), key=sort_key)
    if not rows:
        raise ValueError("no rows could be placed in any of the known transcripts")

    csv_path = dest / "martinelli_extra.csv"
    with open(csv_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    fasta_path = dest / "martinelli_transcripts.fasta"
    with open(fasta_path, "w") as handle:
        for accession in sorted(transcripts):
            sequence = transcripts[accession].replace("U", "T")
            handle.write(f">{accession}\n")
            for i in range(0, len(sequence), 70):
                handle.write(sequence[i: i + 70] + "\n")

    print(
        f"  martinelli: {len(rows)} of {len(table)} rows placed -> "
        f"{csv_path.name}, {fasta_path.name}"
    )


def main() -> None:
    """`python -m sirna_data.fetch.martinelli --normalize <csv>` rewrites an
    existing martinelli_extra.csv in canonical form (see `canonical_row`)."""
    import argparse

    parser = argparse.ArgumentParser(
        prog="python -m sirna_data.fetch.martinelli",
        description="Normalize an existing martinelli_extra.csv in place.",
    )
    parser.add_argument("--normalize", metavar="CSV", type=Path, required=True)
    parser.add_argument(
        "--transcripts", metavar="FASTA", type=Path, default=None,
        help="also re-derive each row's stored target site from these transcripts",
    )
    arguments = parser.parse_args()
    transcripts = None
    if arguments.transcripts is not None:
        from ..raw_loader import _dna_to_rna, read_fasta

        transcripts = {
            accession.split(".")[0]: _dna_to_rna(sequence)
            for accession, sequence in read_fasta(arguments.transcripts).items()
        }
    changed = normalize_file(arguments.normalize, transcripts)
    print(f"{arguments.normalize}: {changed} row(s) rewritten")


if __name__ == "__main__":
    main()
