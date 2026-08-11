"""Fetch CMsiRNAdb's original, unmodified bulk TSV release, plus the NCBI
RefSeq transcripts raw_loader.py's CMsiRNAdb loaders need.

Sources (see ../../../data/DATA_SOURCES.md and
../../../data/CMSIRNADB_FULL_RETRIEVAL.md for full attribution/license notes):
  - CMsiRNAdb (He et al. 2026, BMC Bioinformatics), CC BY-NC-ND 4.0. Its
    "No Derivatives" term means only the ORIGINAL, unmodified download may be
    redistributed -- this module writes exactly that file
    (cmsirnadb_full_raw.tsv) and nothing derived from it. All filtering/
    collapsing happens later, in code, at load time
    (_load_cmsirnadb_records/_load_cmsirnadb_full_records in raw_loader.py).
  - NCBI Nucleotide (efetch), public domain sequence records.
"""
from __future__ import annotations

import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

CMSIRNADB_URL = "https://www.cellknowledge.com.cn/CMsiRNAdb/download/CMsiRNA_data_update.tsv"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
BATCH_SIZE = 40  # keep each efetch request small and polite
REQUEST_DELAY_S = 0.4  # NCBI allows 3 req/s without an API key

# Must match CMSIRNADB_PCSK9_ACCESSION in raw_loader.py: the one canonical
# human PCSK9 coding transcript every PCSK9 row is located against,
# regardless of what accession the row itself cites (many cite a non-coding
# variant or no accession at all -- see ../../../data/DATA_SOURCES.md).
PCSK9_ACCESSION = "NM_174936.4"


def download_cmsirnadb_table(dest: Path) -> pd.DataFrame:
    tsv_path = dest / "cmsirnadb_full_raw.tsv"
    urllib.request.urlretrieve(CMSIRNADB_URL, tsv_path)
    # dtype=str: several downstream columns (accessions, sequences) must
    # round-trip as exact strings, matching raw_loader.py's own read.
    df = pd.read_csv(tsv_path, sep="\t", dtype=str)
    print(f"Downloaded {len(df)} CMsiRNAdb rows covering {df['Target_Gene'].nunique()} genes")
    return df


def fetch_transcript_fasta(accessions: list[str]) -> dict[str, str]:
    """Batch-fetch full mRNA sequences from NCBI by RefSeq accession.

    Keys are kept as the FULL versioned accession (e.g. "NM_174936.4"), not
    stripped -- raw_loader.py looks transcripts up by the exact
    Accession_number string CMsiRNAdb's own table uses, which includes the
    version suffix.
    """
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
                header = line[1:].split()[0]
                chunks = []
            elif line.strip():
                chunks.append(line.strip())
        if header is not None:
            sequences[header] = "".join(chunks)
        print(f"  fetched {min(i + BATCH_SIZE, len(accessions))}/{len(accessions)} accessions")
        time.sleep(REQUEST_DELAY_S)
    return sequences


def write_fasta(path: Path, sequences: dict[str, str]) -> None:
    with open(path, "w") as fh:
        for acc, seq in sequences.items():
            fh.write(f">{acc}\n{seq}\n")


def fetch(dest: Path) -> None:
    """Writes cmsirnadb_full_raw.tsv, cmsirnadb_transcripts.fasta, and
    cmsirnadb_full_transcripts.fasta into `dest`."""
    dest.mkdir(parents=True, exist_ok=True)
    df = download_cmsirnadb_table(dest)

    # PCSK9: fetch only the one canonical accession every PCSK9 row is
    # located against (see _load_cmsirnadb_records) -- not derived from the
    # table's own (inconsistent) Accession_number values for PCSK9 rows.
    print(f"Fetching PCSK9 canonical transcript ({PCSK9_ACCESSION})...")
    pcsk9_seq = fetch_transcript_fasta([PCSK9_ACCESSION])
    if PCSK9_ACCESSION not in pcsk9_seq:
        print(f"WARNING: {PCSK9_ACCESSION} did not resolve")
    pcsk9_fasta_path = dest / "cmsirnadb_transcripts.fasta"
    write_fasta(pcsk9_fasta_path, pcsk9_seq)
    print(f"Wrote {pcsk9_fasta_path} ({len(pcsk9_seq)} sequence)")

    # The other 12 genes: fetch every distinct accession the table cites for
    # them. A handful may not resolve (e.g. superseded/predicted RefSeq
    # records) -- that's fine, _load_cmsirnadb_full_records falls back to
    # duplex-only context for any row whose accession doesn't resolve, the
    # same graceful degradation every other source in this project uses.
    others = df[df["Target_Gene"] != "PCSK9"]
    accessions = others["Accession_number"].dropna().unique().tolist()
    print(f"Fetching {len(accessions)} unique transcripts for the other 12 genes from NCBI...")
    sequences = fetch_transcript_fasta(accessions)
    missing = set(accessions) - set(sequences)
    if missing:
        print(f"WARNING: {len(missing)} accessions did not resolve: {sorted(missing)}")
    full_fasta_path = dest / "cmsirnadb_full_transcripts.fasta"
    write_fasta(full_fasta_path, sequences)
    print(f"Wrote {full_fasta_path} ({len(sequences)} sequences)")

    print(f"Wrote {dest / 'cmsirnadb_full_raw.tsv'} ({len(df)} rows)")
