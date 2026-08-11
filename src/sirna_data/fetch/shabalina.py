"""Fetch the siRNA training database from Shabalina, Spiridonov & Ogurtsov
2006 (BMC Bioinformatics 7:65, CC BY 2.0) -- Additional File 4 ("TableS1A"),
a 653-siRNA / 52-gene heterogeneous compilation used to train their
"ThermoComposition" method.

Retrieved via Europe PMC's public supplementaryFiles API (same legitimate,
documented bulk-access endpoint used for the Monopoli et al. 2023 data --
see ../../../data/DATA_SOURCES.md), not scraping.

We only keep rows targeting genes *not already present* in siRNAEfficacyDB
(our primary dataset): roughly half of the 653 rows (mostly from Khvorova
et al. 2003 and a few other sources) turned out to be exact-sequence
duplicates of genes we already have, confirmed by direct antisense-sequence
matching against data/raw/sirna_efficacy.csv before this module was written.
The gene->accession mapping and the handful of exclusions/corrections below
were derived by hand from NCBI esummary/efetch lookups on each of the 52
accessions in the source table -- see ../../../data/DATA_SOURCES.md for the
full reasoning (in particular: NM_000314/PTEN is excluded as a duplicate of
the already-present MMAC1 gene under its old name; NM_004351 is CBLB, a
distinct paralog from NM_005188/CBL, not the same gene; U47298 is the pGL3
luciferase reporter vector, i.e. the same "Firefly luciferase" gene already
present; M25346 (a puromycin-resistance marker, "PAC") and the two
tissue-factor orthologs are kept as legitimate distinct non-endogenous
targets, the same way "Firefly luciferase"/"SEAP"/"EGFP" already are in the
primary dataset).
"""
from __future__ import annotations

import csv
import io
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

SUPPL_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC1431570/supplementaryFiles"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
TABLE_S1A_FILENAME = "1471-2105-7-65-S4.csv"

# accession -> gene symbol, for every accession in the source table.
# `None` means "same gene as one already in siRNAEfficacyDB" -- excluded.
GENE_MAP: dict[str, str | None] = {
    "AF493916": "HRAS",
    "AJ272212": "PSKH1",
    "AK122643": "FLJ16071",
    "J03132": None,  # ICAM-1, already present
    "M16553": "F3_human",
    "M25346": "PAC",
    "M26071": "F3_mouse",
    "M33197": None,  # GAPDH, already present
    "M60857": None,  # Cyclophilin B, already present
    "M84918": "MyoD",
    "NM_000314": None,  # PTEN == MMAC1 (old name), already present via U92436
    "NM_000321": "RB1",
    "NM_000368": "TSC1",
    "NM_000389": "CDKN1A",
    "NM_000548": "TSC2",
    "NM_000875": "IGF1R",
    "NM_001010": "RPS6",
    "NM_001315": "MAPK14",
    "NM_001344": "DAD1",
    "NM_001626": "AKT2",
    "NM_002015": "FOXO1",
    "NM_002037": "FYN",
    "NM_002046": None,  # GAPDH, already present
    "NM_002093": "GSK3B",
    "NM_002211": "ITGB1",
    "NM_002467": "MYC",
    "NM_002613": "PDPK1",
    "NM_002870": "RAB13",
    "NM_002953": "RPS6KA1",
    "NM_004095": "EIF4EBP1",
    "NM_004351": "CBLB",  # distinct paralog from NM_005188/CBL -- not a duplicate
    "NM_004383": "CSK",
    "NM_004404": "SEPTIN2",
    "NM_004517": "ILK",
    "NM_004586": "RPS6KA3",
    "NM_005027": "PIK3R2",
    "NM_005163": "AKT1",
    "NM_005188": "CBL",
    "NM_005544": "IRS1",
    "NM_005572": None,  # Lamin A, already present
    "NM_005938": "FOXO4",
    "NM_006218": "PIK3CA",
    "NM_006930": "SKP1",
    "NM_019884": "GSK3A",
    "NM_020548": None,  # DBI, already present
    "NM_031313": "ALPG",
    "NM_144586": "LYPD1",
    "U47298": None,  # pGL3 vector == Firefly luciferase, already present
    "U92436": None,  # MMAC1, already present
    "X75932": None,  # PLK, already present
    "XM_043865": "PIK3R1",
    "LaminA": None,  # Lamin A, already present
}


def revcomp(seq: str) -> str:
    comp = {"A": "U", "U": "A", "G": "C", "C": "G"}
    return "".join(comp[b] for b in reversed(seq))


def fetch_table_s1a() -> list[dict]:
    """Downloads the supplementary-files zip and parses Additional File 4.

    The CSV has no header names for its first four columns (accession,
    start, end, antisense 19-mer sequence) and a sparse "References" column
    that's only populated on the first row of each citation block (forward
    filled here to attribute every row to its source publication).
    """
    with urllib.request.urlopen(SUPPL_URL, timeout=60) as resp:
        blob = resp.read()
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:
        with zf.open(TABLE_S1A_FILENAME) as fh:
            text = fh.read().decode("utf-8", errors="replace")

    reader = csv.reader(io.StringIO(text))
    next(reader)  # header
    rows = []
    last_ref = None
    for r in reader:
        if len(r) < 6 or not r[0].strip():
            continue
        while len(r) < 26:
            r.append("")
        if r[25].strip():
            last_ref = r[25].strip()
        rows.append(
            {
                "accession": r[0].strip(),
                "start": r[1].strip(),
                "end": r[2].strip(),
                "antisense_19mer": r[3].strip().upper(),
                "activity_remaining": r[23].strip(),
                "reference": last_ref,
            }
        )
    return rows


def fetch_transcripts(accessions: list[str]) -> dict[str, str]:
    params = urllib.parse.urlencode(
        {"db": "nuccore", "id": ",".join(accessions), "rettype": "fasta", "retmode": "text"}
    )
    with urllib.request.urlopen(f"{EFETCH_URL}?{params}", timeout=60) as resp:
        text = resp.read().decode()
    sequences: dict[str, str] = {}
    header: str | None = None
    chunks: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                sequences[header] = "".join(chunks)
            header = line[1:].split()[0].split(".")[0]  # drop version suffix
            chunks = []
        elif line.strip():
            chunks.append(line.strip())
    if header is not None:
        sequences[header] = "".join(chunks)
    return sequences


def fetch(dest: Path) -> None:
    """Writes shabalina_extra.csv and shabalina_transcripts.fasta into `dest`."""
    dest.mkdir(parents=True, exist_ok=True)

    all_rows = fetch_table_s1a()
    print(f"Fetched {len(all_rows)} rows from Additional File 4")

    new_rows = []
    unknown_accessions = set()
    for r in all_rows:
        acc = r["accession"]
        if acc not in GENE_MAP:
            unknown_accessions.add(acc)
            continue
        gene = GENE_MAP[acc]
        if gene is None:
            continue
        new_rows.append({**r, "gene": gene})
    if unknown_accessions:
        raise RuntimeError(f"Unmapped accessions found in source table: {unknown_accessions}")

    print(f"Keeping {len(new_rows)} rows for {len(set(r['gene'] for r in new_rows))} new genes")

    csv_path = dest / "shabalina_extra.csv"
    with open(csv_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Sequence", "Gene", "Accession_number", "Activity_Remaining_Pct"])
        for r in new_rows:
            sense = revcomp(r["antisense_19mer"])
            w.writerow([sense, r["gene"], r["accession"], r["activity_remaining"]])
    print(f"Wrote {csv_path} ({len(new_rows)} rows)")

    accessions = sorted({r["accession"] for r in new_rows})
    sequences = fetch_transcripts(accessions)
    fasta_path = dest / "shabalina_transcripts.fasta"
    with open(fasta_path, "w") as fh:
        for acc, seq in sequences.items():
            fh.write(f">{acc}\n{seq}\n")
    print(f"Wrote {fasta_path} ({len(sequences)} transcripts)")
