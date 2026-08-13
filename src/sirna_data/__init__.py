"""sirna_data: reusable loader for the siRNA knockdown-efficacy dataset.

    from sirna_data import load_records, SiRNARecord, fetch_mrna_by_gene

`load_records()` returns the full merged, provenance-documented dataset
(siRNAEfficacyDB + supplementary sources) as a list of `SiRNARecord`, each
pairing an siRNA guide with its real local mRNA target context and
experimentally measured knockdown label. See ../../data/DATA_SOURCES.md for
where every row comes from.

`fetch_mrna_by_gene()` looks up a gene's RefSeq mRNA transcript live from
NCBI, for callers that only have a gene symbol.
"""
from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

from .evaluation import GeneCorrelation, PredictionMetrics, evaluate_predictions
from .ncbi_fetch import FetchedTranscript, GeneNotFoundError, fetch_mrna_by_gene
from .rank_confidence import (
    min_top_k_for_confidence,
    probability_true_top_in_predicted_top_k,
    spearman_to_pearson,
)
from .raw_loader import DATA_DIR, SiRNARecord, load_records, read_fasta
from .splitting import leave_n_genes_out, train_test_split

__all__ = [
    "load_records",
    "SiRNARecord",
    "read_fasta",
    "DATA_DIR",
    "fetch_mrna_by_gene",
    "FetchedTranscript",
    "GeneNotFoundError",
    "min_top_k_for_confidence",
    "probability_true_top_in_predicted_top_k",
    "spearman_to_pearson",
    "evaluate_predictions",
    "PredictionMetrics",
    "GeneCorrelation",
    "train_test_split",
    "leave_n_genes_out",
]

try:
    __version__ = version("sirna-data-grabber")
except PackageNotFoundError:
    # Not installed as a package (e.g. running from a raw checkout without
    # `pip install -e .`) -- avoid hardcoding a string here that would just
    # go stale at the next release like the old __version__ = "0.1.0" did.
    __version__ = "0.0.0+unknown"
