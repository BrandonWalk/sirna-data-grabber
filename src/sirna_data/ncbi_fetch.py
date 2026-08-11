"""Look up a gene's RefSeq mRNA transcript from NCBI, for callers who want to
pass a gene symbol instead of a raw sequence.

Uses NCBI's public E-utilities REST API (no key required for this call
volume) -- this is the same programmatic access documented at
https://www.ncbi.nlm.nih.gov/books/NBK25501/, not scraping: esearch finds the
best-matching RefSeq mRNA record in the nuccore database, efetch retrieves
its FASTA, and the DNA-alphabet sequence is converted to RNA (T -> U).
"""
from __future__ import annotations

from dataclasses import dataclass

import requests

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL_PARAMS = {"tool": "sirna-data-grabber", "email": "sirna-data-grabber@example.com"}
REQUEST_TIMEOUT = 15


class GeneNotFoundError(Exception):
    pass


@dataclass
class FetchedTranscript:
    accession: str
    description: str
    sequence: str  # RNA alphabet (T -> U already applied)


def fetch_mrna_by_gene(gene: str, organism: str = "Homo sapiens") -> FetchedTranscript:
    gene = gene.strip()
    if not gene:
        raise GeneNotFoundError("gene symbol is empty")

    term = f'{gene}[Gene Name] AND "{organism}"[Organism] AND biomol_mrna[PROP] AND refseq[filter]'
    search_resp = requests.get(
        f"{EUTILS_BASE}/esearch.fcgi",
        params={**TOOL_PARAMS, "db": "nuccore", "term": term, "retmode": "json", "retmax": "5"},
        timeout=REQUEST_TIMEOUT,
    )
    search_resp.raise_for_status()
    ids = search_resp.json().get("esearchresult", {}).get("idlist", [])
    if not ids:
        raise GeneNotFoundError(
            f"no RefSeq mRNA found for gene {gene!r} in {organism} -- check the symbol "
            "(e.g. official HGNC symbol) and organism"
        )

    fetch_resp = requests.get(
        f"{EUTILS_BASE}/efetch.fcgi",
        params={
            **TOOL_PARAMS,
            "db": "nuccore",
            "id": ids[0],
            "rettype": "fasta",
            "retmode": "text",
        },
        timeout=REQUEST_TIMEOUT,
    )
    fetch_resp.raise_for_status()
    fasta = fetch_resp.text

    lines = fasta.splitlines()
    if not lines or not lines[0].startswith(">"):
        raise GeneNotFoundError(f"unexpected response from NCBI for gene {gene!r}")

    header = lines[0][1:]
    accession = header.split()[0]
    description = header[len(accession) :].strip()
    dna_seq = "".join(lines[1:])
    if not dna_seq:
        raise GeneNotFoundError(f"empty sequence returned for gene {gene!r}")

    return FetchedTranscript(
        accession=accession,
        description=description,
        sequence=dna_seq.upper().replace("T", "U"),
    )
