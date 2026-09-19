"""Fetch Sciabola et al. 2013's in-house siRNA panel -- 21-nt duplexes against
ten hepatocellular-carcinoma-relevant genes, measured as mRNA knockdown in
Hep3B cells by QuantiGene 2.0 at up to four doses.

Sciabola, Cao, Orozco, Faustino & Stanton 2013, Nucleic Acids Research
41(3):1383-1394 (doi:10.1093/nar/gks1191, PMC3561943), CC BY-NC 3.0. The
sequences and per-dose values are Supplementary Tables S3 and S4, inside a
Word document retrieved through Europe PMC's public supplementaryFiles API
(Oxford's own supplementary CDN returns HTTP 403).

The table's sequence column mixes strand orientations, so every row is tested
in both orientations against its gene's RefSeq transcript and stored
sense-normalized; rows that match neither are dropped. See
../../../data/DATA_SOURCES.md.
"""
from __future__ import annotations

import csv
import io
import re
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

from . import _ole

SUPPL_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3561943/supplementaryFiles"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
DOCUMENT_SUFFIX = "File002.doc"

SEQUENCE_COLUMNS = 7  # Gene, 19p2name, 19p2_seq, Rname, R-Dicer_seq, Lname, L-Dicer_seq
VALUE_COLUMNS = 13  # Name + 4 doses (21-mer), then Name + 3 doses per Dicer design
VALUE_HEADER_CELLS = 4  # the merged "19p2 | R-Dicer | L-Dicer" row above the real header
CORE_LENGTH = 19  # duplex core; the extra 2 bases are the 3' overhang

# The table's own gene spellings -> current official symbols, so these rows
# group with the same genes contributed by other sources.
GENE_SYMBOLS = {
    "HIF1a": "HIF1A", "HK2": "HK2", "HPSE": "HPSE", "SURVIVIN": "BIRC5",
    "c-Myc": "MYC", "EZH2": "EZH2", "FRAP1": "MTOR", "bRAF": "BRAF",
    "CTNNB1": "CTNNB1", "PIK3CA": "PIK3CA",
}
ACCESSIONS = {
    "HIF1A": "NM_001530.4", "HK2": "NM_000189.5", "HPSE": "NM_006665.6",
    "BIRC5": "NM_001168.3", "EZH2": "NM_004456.5", "MTOR": "NM_004958.4",
    "BRAF": "NM_004333.6", "CTNNB1": "NM_001904.4", "MYC": "NM_002467.6",
    "PIK3CA": "NM_006218.4",
}
# One HPSE row locates only in this transcript variant, not the one above.
VARIANT_ACCESSIONS = {"HPSE": "NM_001098540.3"}

DOSE_COLUMNS = ["Pct_Inhibition_008nM", "Pct_Inhibition_04nM",
                "Pct_Inhibition_2nM", "Pct_Inhibition_10nM"]
NOT_MEASURED = {"", "n.d.", "nd", "-", "na", "n/a"}


def _download(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=120) as response:
        return response.read()


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGU", "UGCA"))[::-1]


def fetch_document() -> bytes:
    """The supplementary .doc, out of Europe PMC's zip of supplementary files."""
    with zipfile.ZipFile(io.BytesIO(_download(SUPPL_URL))) as archive:
        names = [n for n in archive.namelist() if n.endswith(DOCUMENT_SUFFIX)]
        if not names:
            raise FileNotFoundError(
                f"{DOCUMENT_SUFFIX} not in Europe PMC's supplementary files for "
                f"PMC3561943 (found: {sorted(archive.namelist())})"
            )
        return archive.read(names[0])


def parse_tables(document: bytes) -> tuple[list[list[str]], list[list[str]]]:
    """(sequence rows, per-dose value rows) from the supplementary document."""
    text = _ole.word_text(document)
    blocks = [b for b in re.split(r"(?:\r|\x0c)+", text) if _ole.CELL_MARK in b]
    sequences: list[list[str]] = []
    values: list[list[str]] = []
    for block in blocks:
        if "19p2_seq" in block:
            sequences = _ole.word_table_grid(block, SEQUENCE_COLUMNS)
        elif "0.08nM" in block:
            values = _ole.word_table_grid(block, VALUE_COLUMNS, skip_cells=VALUE_HEADER_CELLS)
    if not sequences or not values:
        raise ValueError("supplementary document has no sequence and/or value table")
    return sequences[1:], values[1:]  # drop each table's header row


def fetch_transcripts(accessions: list[str]) -> dict[str, str]:
    """{accession: RNA sequence} from NCBI, one request for the whole set."""
    query = urllib.parse.urlencode(
        {"db": "nuccore", "id": ",".join(accessions), "rettype": "fasta", "retmode": "text"}
    )
    text = _download(f"{EFETCH_URL}?{query}").decode()
    out: dict[str, str] = {}
    name = None
    for line in text.splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = ""
        elif name:
            out[name] += line.strip().upper().replace("T", "U")
    return out


def _number(cell: str) -> float | None:
    cell = (cell or "").strip()
    return None if cell.lower() in NOT_MEASURED else float(cell)


def build_rows(
    sequences: list[list[str]], values: list[list[str]], transcripts: dict[str, str]
) -> tuple[list[dict[str, str | float | None]], list[str]]:
    """Sense-normalized, transcript-verified rows, plus the compounds dropped."""
    doses = {row[0].strip(): row[1:5] for row in values if row and row[0].strip()}
    rows: list[dict[str, str | float | None]] = []
    dropped: list[str] = []
    for entry in sequences:
        label, number, listed = entry[0].strip(), entry[1].strip(), entry[2].strip().upper()
        gene = GENE_SYMBOLS.get(label)
        if gene is None or not re.fullmatch(r"[ACGU]+", listed):
            continue
        candidates = [ACCESSIONS[gene]]
        if gene in VARIANT_ACCESSIONS:
            candidates.append(VARIANT_ACCESSIONS[gene])
        placed = None
        for accession in candidates:
            transcript = transcripts.get(accession, "")
            for sense, strand in ((listed, "sense"), (revcomp(listed), "antisense")):
                if sense in transcript or sense[-CORE_LENGTH:] in transcript:
                    placed = (accession, sense, strand)
                    break
            if placed:
                break
        if placed is None:
            dropped.append(f"19p2_{number}")
            continue
        accession, sense, strand = placed
        measured = [_number(cell) for cell in doses.get(number, [""] * 4)]
        if not any(value is not None for value in measured):
            dropped.append(f"19p2_{number}")
            continue
        rows.append({
            "Compound_Name": f"19p2_{number}", "Gene": gene, "Source_Gene_Label": label,
            "Accession_number": accession, "Sense_21mer": sense,
            "Table_S3_Sequence": listed, "Table_S3_Strand": strand,
            **dict(zip(DOSE_COLUMNS, measured, strict=True)),
        })
    return rows, dropped


def fetch(dest: Path) -> None:
    """Write sciabola2013_extra.csv and sciabola2013_transcripts.fasta."""
    dest.mkdir(parents=True, exist_ok=True)
    sequences, values = parse_tables(fetch_document())
    accessions = sorted(set(ACCESSIONS.values()) | set(VARIANT_ACCESSIONS.values()))
    transcripts = fetch_transcripts(accessions)
    missing = [a for a in accessions if a not in transcripts]
    if missing:
        raise ValueError(f"NCBI returned no sequence for: {missing}")

    rows, dropped = build_rows(sequences, values, transcripts)
    if not rows:
        raise ValueError("no rows could be located in their transcripts")

    csv_path = dest / "sciabola2013_extra.csv"
    with open(csv_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    used = {str(row["Accession_number"]) for row in rows}
    fasta_path = dest / "sciabola2013_transcripts.fasta"
    with open(fasta_path, "w") as handle:
        for accession in sorted(used):
            sequence = transcripts[accession].replace("U", "T")
            handle.write(f">{accession}\n")
            for i in range(0, len(sequence), 70):
                handle.write(sequence[i: i + 70] + "\n")

    note = f" ({len(dropped)} dropped: not locatable)" if dropped else ""
    print(f"  sciabola2013: {len(rows)} rows{note} -> {csv_path.name}, {fasta_path.name}")
