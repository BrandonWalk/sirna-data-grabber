#!/usr/bin/env python3
"""Fetch the small supplementary dataset from Monopoli et al. 2023 (Molecular
Therapy Nucleic Acids, CC BY 4.0) -- Table S3, 20 chemically-modified siRNAs
against 4 genes (APP, MAPT, BACE1, SNCA) not present in siRNAEfficacyDB.

Retrieved via Europe PMC's public supplementaryFiles API (no bot-detection
evasion involved -- see data/DATA_SOURCES.md for how this differs from the
ThermoFisher catalog we declined to scrape).

This is a distinct, smaller augmentation to the primary dataset: different
assay chemistry (cholesterol-conjugated, heavily 2'-F/2'-OMe/phosphorothioate
modified "sdRNA", not a standard unmodified siRNA duplex) -- see
data/DATA_SOURCES.md for the caveats before trusting this data the same way
as the primary set.
"""
from __future__ import annotations

import urllib.parse
import urllib.request
from pathlib import Path

RAW_DIR = Path(__file__).resolve().parents[1] / "data" / "raw"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# Table S3 of Monopoli et al. 2023 (doi:10.1016/j.omtn.2023.06.010), transcribed
# by hand from the extracted supplementary PDF text. "sequence" is the mRNA
# target-site sequence (verified by exact match against the transcripts below,
# NOT the antisense guide -- the guide is derived as its reverse complement in
# raw_loader.py). "reporter_remaining_pct" is % luciferase reporter expression
# remaining vs. untreated control (lower = more potent knockdown).
TABLE_S3 = [
    # sequence, gene, reporter_remaining_pct, sd_pct
    ("UUCAAUAUGCUAAAGAAGUA", "APP", 19.2, 5.7),
    ("GUCCAAGUGUGGCUCAAAGG", "MAPT", 61.9, 8.5),
    ("GGUCCUAAGCCCACAAUCAU", "MAPT", 6.6, 2.7),
    ("UGAUCGGGCCCGAAAACGAA", "BACE1", 26.9, 8.1),
    ("UUUUGAAAGGCUUUCCUCAG", "MAPT", 12.5, 2.0),
    ("CUUUGUGAUUCCCUACCGCU", "APP", 51.6, 6.5),
    ("CAUUGAGACUUCAAGCUUUU", "APP", 10.2, 2.6),
    ("UAGUGCAUGAAUAGAUUCUC", "APP", 4.8, 0.6),
    ("GUGGGAGUUCAGCUGCUUCU", "APP", 10.7, 1.7),
    ("GUCACCUUAAAGGAGAUCAA", "SNCA", 11.2, 3.2),
    ("UGCUGCCAUGAUUUUGGCCA", "MAPT", 50.3, 6.6),
    ("AGCCUCUGAAGUUGGACAGC", "APP", 45.9, 9.9),
    ("AUGGUUUCUGGCUAGGAGAG", "BACE1", 99.2, 1.3),
    ("AUGAUCGCUUUCUACACUGU", "APP", 12.2, 1.9),
    ("ACUUUCAGAACUGCUACCAU", "BACE1", 19.7, 2.5),
    ("AUGGGUGCUGAAAAUAAACU", "SNCA", 29.2, 3.2),
    ("AAGCAGCAUAUUUUAAAAAU", "SNCA", 72.2, 6.8),
    ("CAAGUGACAAAUGUUGGAGG", "SNCA", 47.1, 5.1),
    ("CAAAGUCCAGGCACAAGAGU", "MAPT", 90.6, 5.3),
    ("AUUCUCCAAAACAAUUUUCU", "APP", 32.2, 7.3),
]

# Canonical RefSeq accession used per gene (verified: each Table S3 sequence
# above is an exact substring of the corresponding transcript).
GENE_ACCESSION = {
    "APP": "NM_000484",
    "MAPT": "NM_005910",
    "BACE1": "NM_012104",
    "SNCA": "NM_000345",
}


def fetch_transcripts() -> dict[str, str]:
    accs = list(GENE_ACCESSION.values())
    params = urllib.parse.urlencode(
        {"db": "nuccore", "id": ",".join(accs), "rettype": "fasta", "retmode": "text"}
    )
    with urllib.request.urlopen(f"{EFETCH_URL}?{params}", timeout=60) as resp:
        text = resp.read().decode()
    sequences: dict[str, str] = {}
    header, chunks = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                sequences[header] = "".join(chunks)
            header = line[1:].split()[0].split(".")[0]
            chunks = []
        elif line.strip():
            chunks.append(line.strip())
    if header is not None:
        sequences[header] = "".join(chunks)
    return sequences


def main() -> None:
    RAW_DIR.mkdir(parents=True, exist_ok=True)

    csv_path = RAW_DIR / "monopoli_extra.csv"
    with open(csv_path, "w") as fh:
        fh.write("Sequence,Gene,Accession_number,Reporter_Remaining_Pct,SD_Pct\n")
        for seq, gene, remaining, sd in TABLE_S3:
            fh.write(f"{seq},{gene},{GENE_ACCESSION[gene]},{remaining},{sd}\n")
    print(f"Wrote {csv_path} ({len(TABLE_S3)} rows)")

    sequences = fetch_transcripts()
    fasta_path = RAW_DIR / "monopoli_transcripts.fasta"
    with open(fasta_path, "w") as fh:
        for acc, seq in sequences.items():
            fh.write(f">{acc}\n{seq}\n")
    print(f"Wrote {fasta_path} ({len(sequences)} transcripts)")


if __name__ == "__main__":
    main()
