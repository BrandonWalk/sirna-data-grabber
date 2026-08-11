"""Fetch siRNA knockdown-efficacy records and their matching full-length
mRNA transcripts.

Sources (see ../../../data/DATA_SOURCES.md for full attribution/license notes):
  - siRNAEfficacyDB (Zhang et al. 2024, IET Systems Biology), CC BY-NC.
  - NCBI Nucleotide (efetch), public domain sequence records.
"""
from __future__ import annotations

import tempfile
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

EFFICACY_DB_URL = "https://cellknowledge.com.cn/siRNAEfficacy/download/Gene_all.txt"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
BATCH_SIZE = 40  # keep each efetch request small and polite
REQUEST_DELAY_S = 0.4  # NCBI allows 3 req/s without an API key


def download_efficacy_table() -> pd.DataFrame:
    """Downloads to a temp file, not `dest` -- only the cleaned
    `sirna_efficacy.csv` fetch() writes at the end is meant to persist;
    keeping this pre-cleaning snapshot around too was just a redundant
    duplicate."""
    with tempfile.NamedTemporaryFile(suffix=".tsv") as tmp:
        urllib.request.urlretrieve(EFFICACY_DB_URL, tmp.name)
        df = pd.read_csv(tmp.name, sep="\t")
    print(f"Downloaded {len(df)} siRNA records covering {df['Gene'].nunique()} genes")
    return df


def fetch_mrna_fasta(accessions: list[str]) -> dict[str, str]:
    """Batch-fetch full mRNA sequences from NCBI by RefSeq/GenBank accession."""
    sequences: dict[str, str] = {}
    accessions = sorted(set(accessions))
    for i in range(0, len(accessions), BATCH_SIZE):
        batch = accessions[i : i + BATCH_SIZE]
        params = urllib.parse.urlencode(
            {"db": "nuccore", "id": ",".join(batch), "rettype": "fasta", "retmode": "text"}
        )
        with urllib.request.urlopen(f"{EFETCH_URL}?{params}", timeout=60) as resp:
            text = resp.read().decode()
        header: str | None = None
        chunks: list[str] = []
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
        print(f"  fetched {min(i + BATCH_SIZE, len(accessions))}/{len(accessions)} accessions")
        time.sleep(REQUEST_DELAY_S)
    return sequences


def fetch(dest: Path) -> None:
    """Writes sirna_efficacy.csv and mrna_transcripts.fasta into `dest`."""
    dest.mkdir(parents=True, exist_ok=True)
    df = download_efficacy_table()

    # Drop the Renilla luciferase assay control rows (not an endogenous gene target).
    before = len(df)
    df = df[df["Accession_number"] != "-"].reset_index(drop=True)
    print(f"Dropped {before - len(df)} control rows without a gene accession")

    accessions = df["Accession_number"].unique().tolist()
    print(f"Fetching {len(accessions)} unique mRNA transcripts from NCBI...")
    sequences = fetch_mrna_fasta(accessions)
    missing = set(accessions) - set(sequences)
    if missing:
        print(f"WARNING: {len(missing)} accessions did not resolve: {sorted(missing)}")

    fasta_path = dest / "mrna_transcripts.fasta"
    with open(fasta_path, "w") as fh:
        for acc, seq in sequences.items():
            fh.write(f">{acc}\n{seq}\n")

    csv_path = dest / "sirna_efficacy.csv"
    df.to_csv(csv_path, index=False)
    print(f"Wrote {csv_path} ({len(df)} rows) and {fasta_path} ({len(sequences)} sequences)")
