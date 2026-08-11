"""Fetch the raw data `sirna_data.load_records()` reads, from each source's
original location, into a local directory -- so a plain `pip install
sirna-data-grabber` install can reconstruct the dataset without cloning the
git repo.

    from sirna_data.fetch import cmsirnadb, monopoli, shabalina, sirna_efficacy
    sirna_efficacy.fetch(Path("./my_data"))

Or from the command line (installed automatically, no extras required --
every fetcher here only needs pandas, already a core dependency, plus the
standard library):

    sirna-data-fetch --dest ./my_data
    export SIRNA_DATA_DIR=./my_data

This only covers the four sources `load_records()` actually reads
(siRNAEfficacyDB, Monopoli 2023, Shabalina 2006, CMsiRNAdb). See
../../../data/DATA_SOURCES.md for full provenance/license notes per source,
and NOTICE.md before using the fetched data commercially -- most of it is
non-commercial only.
"""
from __future__ import annotations
