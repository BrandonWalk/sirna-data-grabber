"""The guard that keeps the README honest: everything `load_records()` can
load must be reachable by `sirna-data-fetch`.

Without this, a source can be added with no fetcher and nothing fails --
which is exactly how this package came to advertise a dataset a
`pip install` could not reconstruct.
"""
from __future__ import annotations

import inspect

from sirna_data.fetch.cli import SOURCES
from sirna_data.raw_loader import load_records

# The one fetcher that feeds two loader flags: `cmsirnadb` downloads the
# single raw TSV that both the PCSK9 subset and the 12-gene addition are
# derived from at load time.
FETCHER_FOR_FLAG = {"cmsirnadb_full": "cmsirnadb"}


def _loader_flags() -> set[str]:
    return {
        name[len("include_"):]
        for name in inspect.signature(load_records).parameters
        if name.startswith("include_")
    }


def test_every_loadable_source_has_a_fetcher():
    uncovered = {
        flag for flag in _loader_flags()
        if FETCHER_FOR_FLAG.get(flag, flag) not in SOURCES
    }
    assert not uncovered, (
        f"{sorted(uncovered)} can be loaded but not fetched -- a pip install "
        "cannot reconstruct the dataset the README advertises. Add a fetcher, "
        "or document the omission explicitly."
    )


def test_every_fetcher_feeds_a_loadable_source():
    flags = _loader_flags()
    covered = {FETCHER_FOR_FLAG.get(flag, flag) for flag in flags}
    orphans = set(SOURCES) - covered
    assert not orphans, f"{sorted(orphans)} fetches data no loader reads"


def test_every_fetcher_is_callable_with_a_destination():
    for name, fetch in SOURCES.items():
        parameters = list(inspect.signature(fetch).parameters)
        assert parameters == ["dest"], f"{name}.fetch should take exactly (dest), got {parameters}"
