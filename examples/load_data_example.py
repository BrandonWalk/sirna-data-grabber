"""Example: load the siRNA knockdown-efficacy dataset and take a first look
at it -- see what genes and licenses are available, inspect a record,
filter by source/gene, and do a gene-grouped train/test split for a model.

    python examples/load_data_example.py

See ../README.md and ../data/DATA_SOURCES.md for what's actually in the
dataset and where each source comes from.
"""
from __future__ import annotations

from collections import Counter

from sirna_data import (
    SiRNARecord,
    describe_genes,
    list_licenses,
    list_sources,
    load_records,
    train_test_split,
)


def describe_record(r: SiRNARecord) -> None:
    print(f"  row_id           = {r.row_id}")
    print(f"  gene / accession = {r.gene} / {r.accession}")
    print(f"  guide_seq        = {r.guide_seq} ({r.duplex_len}nt duplex)")
    print(f"  label (%KD)      = {r.label:.1f}")
    print(f"  source / assay   = {r.source} / {r.technology}")
    print(f"  modified?        = {r.is_modified} ({r.modification_chemistry})")
    if r.has_flanking_context:
        site = r.mrna_window[r.site_start : r.site_start + r.site_len]
        print(f"  mrna_window      = {r.mrna_window}")
        print(f"  target site      = {site} (at offset {r.site_start})")
    else:
        print("  mrna_window      = <no transcript match -- duplex-only context>")


def main() -> None:
    # load_records() pulls every source in data/raw/ (siRNAEfficacyDB,
    # CMsiRNAdb, Shabalina, Martinelli, Davis2025, Monopoli,
    # Sciabola2013, Harborth2003) into one merged list. Each include_* flag defaults to
    # True; pass include_cmsirnadb_full=False etc. to drop a source, or
    # data_dir="/path/to/data" to point at a different copy of data/raw/.
    records = load_records()
    print(f"Loaded {len(records)} records total\n")

    # What's in here, before doing anything with it: every gene available
    # to load, with its record count, provenance and license(s).
    genes = describe_genes(records)  # or describe_genes() to load for you
    print(f"{len(genes)} genes available; the 5 with the most records:")
    for info in sorted(genes, key=lambda g: -g.n_records)[:5]:
        print(f"  {info.gene:<10} {info.n_records:>6} records  {'+'.join(info.licenses)}")

    # Every record is a SiRNARecord: siRNA guide + local mRNA target
    # context + experimentally measured knockdown label. Look at one.
    print("First record:")
    describe_record(records[0])

    # Basic breakdowns.
    by_source = Counter(r.source for r in records)
    print("\nRecords per source:")
    for source, n in by_source.most_common():
        print(f"  {source:<28} {n:>6}")

    n_genes = len({r.gene for r in records})
    n_with_context = sum(r.has_flanking_context for r in records)
    n_modified = sum(r.is_modified for r in records)
    print(f"\n{n_genes} distinct genes")
    print(f"{n_with_context}/{len(records)} records have real flanking mRNA context")
    print(f"{n_modified}/{len(records)} records are chemically modified constructs")

    # The data's licenses vary per source and most are non-commercial, so
    # "load only what I'm allowed to use" is a first-class filter rather
    # than something to reimplement per project.
    print(f"\nLicenses present: {', '.join(list_licenses())}")
    permissive = [entry.license_id for entry in list_sources() if entry.commercial_use]
    commercial_ok = load_records(licenses=permissive)
    print(
        f"{len(commercial_ok)}/{len(records)} records across "
        f"{len({r.gene for r in commercial_ok})} genes are from sources that permit "
        "commercial use (attribution still required -- see NOTICE.md)"
    )

    # Filter down to one gene, e.g. to inspect every siRNA tested against
    # PCSK9 across all sources that cover it.
    pcsk9_records = [r for r in records if r.gene == "PCSK9"]
    print(f"\n{len(pcsk9_records)} records target PCSK9")

    # Gene-grouped train/test split -- no gene appears in both splits, so a
    # model can't "cheat" on a held-out gene just by having seen a
    # different siRNA against the same mRNA during training. Use
    # by_gene=False for a plain per-record split instead, or
    # leave_n_genes_out() for full gene-level cross-validation.
    train, test = train_test_split(records, test_size=0.2, by_gene=True, random_state=0)
    print(f"\nGene-grouped split: {len(train)} train / {len(test)} test records")
    print(f"  {len({r.gene for r in train})} train genes, {len({r.gene for r in test})} test genes")


if __name__ == "__main__":
    main()
