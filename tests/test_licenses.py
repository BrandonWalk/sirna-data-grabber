from __future__ import annotations

import inspect
from pathlib import Path

import pytest

from sirna_data.licenses import (
    LICENSE_UNRESOLVED,
    SOURCE_LICENSES,
    license_for_source,
    list_licenses,
    list_sources,
    normalize_license,
    source_keys_for_licenses,
)
from sirna_data.raw_loader import load_records

FLANK = 5  # must match conftest.FLANK


# --------------------------------------------------------------------------
# the registry itself: it's only useful if it stays in sync with the loaders
# --------------------------------------------------------------------------


def test_every_registry_key_is_a_real_load_records_flag():
    flags = set(inspect.signature(load_records).parameters)
    missing = {k for k in SOURCE_LICENSES if f"include_{k}" not in flags}
    assert not missing, f"registry keys with no include_* flag on load_records: {missing}"


def test_every_load_records_source_flag_is_in_the_registry():
    # The reverse direction: a source added to load_records without a
    # license entry would silently become unselectable by license.
    flags = {
        name[len("include_") :]
        for name in inspect.signature(load_records).parameters
        if name.startswith("include_")
    }
    assert flags == set(SOURCE_LICENSES)


def test_every_loaded_record_source_resolves_to_a_license(
    patch_data_dir: Path, fake_data_dir: Path
):
    records = load_records(flank_nt=FLANK)
    unresolvable = {r.source for r in records if license_for_source(r.source) is None}
    assert not unresolvable, f"record sources with no registry entry: {unresolvable}"


def test_registry_keys_match_their_own_key_field():
    assert all(key == entry.key for key, entry in SOURCE_LICENSES.items())


def test_unresolved_sources_report_unknown_rather_than_permissive():
    for entry in list_sources():
        if entry.license_id == LICENSE_UNRESOLVED:
            # None means "nobody established this", which must never be
            # confused with False ("established as forbidden") or True.
            assert entry.commercial_use is None
            assert entry.derivatives_redistributable is None
        else:
            assert entry.commercial_use is not None


def test_no_derivatives_terms_are_flagged_as_such():
    # CMsiRNAdb's "ND" is the reason its subsets are derived in code at load
    # time instead of shipped as CSVs -- the registry must say so.
    for entry in list_sources():
        if "ND" in entry.license_id.split():
            assert entry.derivatives_redistributable is False


# --------------------------------------------------------------------------
# list_licenses / normalize_license
# --------------------------------------------------------------------------


def test_list_licenses_is_sorted_deduped_and_includes_unresolved():
    licenses = list_licenses()
    assert licenses == sorted(set(licenses))
    assert LICENSE_UNRESOLVED in licenses
    # both CMsiRNAdb entries share one license -> one entry, not two
    assert licenses.count("CC BY-NC-ND 4.0") == 1


def test_list_sources_preserves_load_order():
    assert [e.key for e in list_sources()] == list(SOURCE_LICENSES)


@pytest.mark.parametrize(
    "given",
    [
        "CC BY-NC-ND 4.0",
        "cc by-nc-nd 4.0",
        "CC-BY-NC-ND-4.0",
        "cc_by_nc_nd_4.0",
        "  CC BY NC ND 4.0  ",
    ],
)
def test_normalize_license_accepts_punctuation_and_case_variants(given: str):
    assert normalize_license(given) == "CC BY-NC-ND 4.0"


def test_normalize_license_keeps_versions_distinct():
    # Two sources here state two genuinely different things; folding the
    # version away would silently widen a selection.
    assert normalize_license("CC BY-NC") == "CC BY-NC"
    assert normalize_license("CC BY-NC 4.0") == "CC BY-NC 4.0"


def test_normalize_license_rejects_unknown_and_lists_valid_ids():
    with pytest.raises(ValueError) as excinfo:
        normalize_license("MIT")
    message = str(excinfo.value)
    assert "MIT" in message
    for license_id in list_licenses():
        assert license_id in message


# --------------------------------------------------------------------------
# source_keys_for_licenses
# --------------------------------------------------------------------------


def test_source_keys_for_licenses_returns_load_order():
    keys = source_keys_for_licenses(["CC BY 4.0", "CC BY 2.0"])
    assert keys == ["monopoli", "shabalina", "davis2025"]


def test_source_keys_for_licenses_groups_sources_sharing_a_license():
    assert source_keys_for_licenses(["CC BY-NC-ND 4.0"]) == ["cmsirnadb", "cmsirnadb_full"]


def test_source_keys_for_licenses_never_includes_unresolved_implicitly():
    unresolved_keys = {e.key for e in list_sources() if e.license_id == LICENSE_UNRESOLVED}
    selected = set(
        source_keys_for_licenses([lic for lic in list_licenses() if lic != LICENSE_UNRESOLVED])
    )
    assert not (selected & unresolved_keys)
    # ...but is selectable when asked for explicitly.
    assert set(source_keys_for_licenses([LICENSE_UNRESOLVED])) == unresolved_keys


def test_source_keys_for_licenses_rejects_unknown_license():
    with pytest.raises(ValueError):
        source_keys_for_licenses(["CC BY 4.0", "WTFPL"])


def test_source_keys_for_licenses_accepts_empty_selection():
    assert source_keys_for_licenses([]) == []


# --------------------------------------------------------------------------
# license_for_source
# --------------------------------------------------------------------------


def test_license_for_source_exact_match():
    assert license_for_source("siRNAEfficacyDB").key == "sirna_efficacy"


def test_license_for_source_resolves_discriminated_sources_by_prefix():
    # A `<registered source>_<discriminator>` string resolves to its entry.
    assert license_for_source("CMsiRNAdb_some_future_subset").key == "cmsirnadb"


def test_license_for_source_prefers_exact_over_prefix():
    # "CMsiRNAdb_full" also prefix-matches "CMsiRNAdb" -- its own entry wins.
    assert license_for_source("CMsiRNAdb_full").key == "cmsirnadb_full"
    assert license_for_source("CMsiRNAdb").key == "cmsirnadb"


def test_license_for_source_returns_none_for_unknown_source():
    assert license_for_source("SomeSourceAddedLater") is None
