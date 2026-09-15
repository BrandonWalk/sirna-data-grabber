"""List the genes available to load, and what each one comes with.

`load_records()` answers "give me every record"; this module answers the
question you usually have first -- "what genes are even in here, and which
of them can I actually use?" -- without making you group the records
yourself:

    from sirna_data import list_genes, describe_genes

    list_genes()                              # ['ACP5', 'AGT', 'AKT1', ...]
    list_genes(licenses=["CC BY 4.0"])        # just the permissive subset
    describe_genes()[0].sources               # provenance per gene

Both functions derive their answer from the actual fetched data files (via
`load_records`), not from a hardcoded table -- so they can't go stale as
sources are added or their filtering changes, and they honor every
`load_records()` argument (`include_*`, `licenses=`, `data_dir=`). The cost
is that they do a full load: seconds, and they need the data present. Pass
records you've already loaded to skip that:

    records = load_records()
    genes = list_genes(records)
"""
from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

from .licenses import LICENSE_UNRESOLVED, license_for_source
from .raw_loader import SiRNARecord, load_records


@dataclass(frozen=True)
class GeneInfo:
    """What one gene contributes to a loaded dataset.

    gene : the gene string exactly as the raw sources spell it. Not
        normalized: siRNAEfficacyDB's `Gene` column has a trailing space on
        `"EGFP "` where Martinelli's says `"EGFP"`, and those stay two
        distinct entries here because that's how the raw data actually
        compares (see the README's gene table footnote).
    n_records : how many siRNA records target it.
    sources : the distinct `SiRNARecord.source` values contributing to it,
        sorted -- a gene can come from more than one source (e.g. `APP`).
    licenses : the distinct license ids of those sources, sorted. More than
        one means this gene's records are NOT uniformly usable: check
        per-record provenance, or re-load with `licenses=` to get only the
        subset you can use. A source with no registry entry is reported as
        `"unresolved"` rather than silently dropped.
    accessions : the distinct transcript accessions its target sites were
        located in, sorted -- more than one where sources disagree on the
        canonical transcript. Records whose source left the accession blank
        (a raw-data gap on ~900 CMsiRNAdb_full rows, which then fall back to
        duplex-only context) contribute nothing here, so this can be empty
        for a gene that still has records.
    """

    gene: str
    n_records: int
    sources: tuple[str, ...]
    licenses: tuple[str, ...]
    accessions: tuple[str, ...]


def _resolve_records(
    records: Sequence[SiRNARecord] | None, load_kwargs: dict[str, Any]
) -> Sequence[SiRNARecord]:
    if records is None:
        return load_records(**load_kwargs)
    if load_kwargs:
        unexpected = ", ".join(sorted(load_kwargs))
        raise TypeError(
            f"got both `records` and load_records() argument(s) ({unexpected}) -- pass "
            "one or the other: `records` describes an already-loaded list as-is, while "
            "load_records() arguments are only used when loading for you."
        )
    return records


def describe_genes(
    records: Sequence[SiRNARecord] | None = None, **load_kwargs: Any
) -> list[GeneInfo]:
    """Every gene in the dataset as a `GeneInfo` (record count, sources,
    licenses, accessions), sorted case-insensitively by gene name.

    With no arguments, loads the full default dataset first. `**load_kwargs`
    are passed straight through to `load_records()`, so any subset can be
    described the same way it's loaded -- `describe_genes(licenses=["CC BY
    4.0"])`, `describe_genes(include_cmsirnadb_full=False)`,
    `describe_genes(data_dir="./my_data")`. Pass an already-loaded record
    list as `records` to describe it without re-loading (in which case no
    `load_records()` arguments are accepted, since nothing is being loaded).
    """
    resolved = _resolve_records(records, load_kwargs)

    counts: dict[str, int] = {}
    sources: dict[str, set[str]] = {}
    licenses: dict[str, set[str]] = {}
    accessions: dict[str, set[str]] = {}
    for record in resolved:
        gene = record.gene
        counts[gene] = counts.get(gene, 0) + 1
        sources.setdefault(gene, set()).add(record.source)
        entry = license_for_source(record.source)
        licenses.setdefault(gene, set()).add(
            entry.license_id if entry is not None else LICENSE_UNRESOLVED
        )
        # Guarded rather than added blindly: a blank accession comes
        # through from pandas as NaN, not a string, on the raw rows that
        # don't state one.
        accession = record.accession
        accessions.setdefault(gene, set())
        if isinstance(accession, str) and accession.strip():
            accessions[gene].add(accession)

    return [
        GeneInfo(
            gene=gene,
            n_records=counts[gene],
            sources=tuple(sorted(sources[gene])),
            licenses=tuple(sorted(licenses[gene])),
            accessions=tuple(sorted(accessions[gene])),
        )
        for gene in sorted(counts, key=lambda g: (g.strip().lower(), g))
    ]


def list_genes(
    records: Sequence[SiRNARecord] | None = None, **load_kwargs: Any
) -> list[str]:
    """Just the gene names available to load, sorted case-insensitively.

    Takes the same arguments as `describe_genes` (which is what to call for
    counts, per-gene provenance, and licenses):

        list_genes()                          # every gene in the default load
        list_genes(licenses=["CC BY 4.0"])    # only the permissive sources' genes
        list_genes(records)                   # from records you already loaded
    """
    return [info.gene for info in describe_genes(records, **load_kwargs)]
