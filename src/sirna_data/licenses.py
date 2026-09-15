"""Machine-readable per-source data-license registry, plus the machinery
behind `load_records(licenses=[...])` -- selecting sources by license
instead of by name.

The CODE in this package is MIT licensed; the DATA it loads is not. Every
source in `data/raw/` carries its original publisher's own terms -- a few
permissive (CC BY), most non-commercial (CC BY-NC / CC BY-NC-ND), and two
genuinely unresolved (no LICENSE file in the source repo at all, so
all-rights-reserved by default). That split is documented in prose in
NOTICE.md and data/DATA_SOURCES.md; this module is the same table in code,
so a caller can express "load only what my project's license terms allow"
directly instead of hand-maintaining their own mapping of `include_*` flags:

    from sirna_data import list_licenses, list_sources, load_records

    list_licenses()                        # every license id in the dataset
    load_records(licenses=["CC BY 4.0"])   # only the CC BY 4.0 sources

This registry is documentation, not legal advice, and it is only as current
as NOTICE.md -- read each source's real terms before relying on it. In
particular `commercial_use=None` on the `unresolved` sources means *nobody
established what the terms are*, not that anything goes.
"""
from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import dataclass

#: `license_id` for a source whose terms could not be established at all
#: (the source repo ships no LICENSE file, so it is all-rights-reserved by
#: default under GitHub's terms). Deliberately a normal license id rather
#: than None, so it can be passed to `load_records(licenses=...)` like any
#: other -- but only ever explicitly, never as part of a real CC selection.
LICENSE_UNRESOLVED = "unresolved"


@dataclass(frozen=True)
class SourceLicense:
    """The license facts for one loadable source.

    key : the source's `load_records()` flag suffix -- `"cmsirnadb_full"`
        here is `include_cmsirnadb_full=` there. This is the identity a
        license selection ultimately resolves to.
    name : human-readable source name, as used in NOTICE.md's table.
    license_id : canonical license id, e.g. `"CC BY-NC-ND 4.0"`, or
        `LICENSE_UNRESOLVED`. Version-specific: `"CC BY-NC"` (what
        siRNAEfficacyDB states) and `"CC BY-NC 4.0"` (what Martinelli
        states) are deliberately distinct ids, because that is how the two
        sources actually state their terms.
    commercial_use : True if the license permits commercial use, False if
        it forbids it, None if unresolved (do NOT read None as "probably
        fine").
    derivatives_redistributable : True if adaptations of the data may be
        redistributed, False for the "ND" (No Derivatives) terms that are
        why CMsiRNAdb is filtered in code at load time rather than shipped
        as a derived CSV (see raw_loader.py's CMsiRNAdb note), None if
        unresolved.
    record_sources : the `SiRNARecord.source` values this source produces.
        A record's source may also be one of these with a `_`-suffixed
        discriminator (`REMOVED` -> `REMOVED_Sciabola2013`), which
        `license_for_source` resolves by prefix.
    url : the source's own landing page/paper, where there is one.
    notes : the short caveat worth seeing next to the license id.
    """

    key: str
    name: str
    license_id: str
    commercial_use: bool | None
    derivatives_redistributable: bool | None
    record_sources: tuple[str, ...]
    url: str | None
    notes: str


#: Every source `load_records()` can load, keyed by its `include_*` flag
#: suffix, in load order. Mirrors NOTICE.md's table -- keep the two in sync.
SOURCE_LICENSES: dict[str, SourceLicense] = {
    "sirna_efficacy": SourceLicense(
        key="sirna_efficacy",
        name="siRNAEfficacyDB (Zhang et al. 2024)",
        license_id="CC BY-NC",
        commercial_use=False,
        derivatives_redistributable=True,
        record_sources=("siRNAEfficacyDB",),
        url="https://cellknowledge.com.cn/siRNAEfficacy",
        notes="Non-commercial only. The primary source (3,532 of the records).",
    ),
    "monopoli": SourceLicense(
        key="monopoli",
        name="Monopoli, Korkin & Khvorova 2023",
        license_id="CC BY 4.0",
        commercial_use=True,
        derivatives_redistributable=True,
        record_sources=("Monopoli2023",),
        url="https://doi.org/10.1093/nar/gkad665",
        notes="Attribution required. 20 chemically modified sdRNAs.",
    ),
    "REMOVED": SourceLicense(
        key="REMOVED",
        name="REMOVED panel (Xu, Zhao et al. 2024 / REMOVED repo)",
        license_id=LICENSE_UNRESOLVED,
        commercial_use=None,
        derivatives_redistributable=None,
        record_sources=("REMOVED_REMOVED",),
        url="https://github.com/REMOVED/REMOVED",
        notes=(
            "No LICENSE file in the source repo (all-rights-reserved by "
            "default) and the data was recovered from a deleted file in its "
            "git history -- verify redistribution rights with the authors "
            "before relying on this subset. See NOTICE.md."
        ),
    ),
    "shabalina": SourceLicense(
        key="shabalina",
        name="Shabalina, Spiridonov & Ogurtsov 2006",
        license_id="CC BY 2.0",
        commercial_use=True,
        derivatives_redistributable=True,
        record_sources=("Shabalina2006",),
        url="https://doi.org/10.1186/1471-2105-7-65",
        notes="Attribution required. Heterogeneous published compilation.",
    ),
    "martinelli": SourceLicense(
        key="martinelli",
        name="Martinelli et al. / sirna-reproduction 2023",
        license_id="CC BY-NC 4.0",
        commercial_use=False,
        derivatives_redistributable=True,
        record_sources=("Martinelli_sirna_reproduction",),
        url="https://github.com/mmartinelli-bio/sirna-reproduction",
        notes=(
            "Non-commercial only. Gene identity for the included rows was "
            "recovered here, not given by the source -- see "
            "data/DATA_SOURCES.md."
        ),
    ),
    "REMOVED": SourceLicense(
        key="REMOVED",
        name="REMOVED repo (Sciabola et al. 2013 + Harborth et al. 2001)",
        license_id=LICENSE_UNRESOLVED,
        commercial_use=None,
        derivatives_redistributable=None,
        record_sources=("REMOVED",),
        url="https://github.com/drugparadigm/REMOVED",
        notes=(
            "No LICENSE file in the source repo (all-rights-reserved by "
            "default). Its `label` column's scale is also undocumented -- "
            "see _load_REMOVED_records's caveat before trusting the "
            "labels quantitatively."
        ),
    ),
    "cmsirnadb": SourceLicense(
        key="cmsirnadb",
        name="CMsiRNAdb, human PCSK9 subset (He et al. 2026)",
        license_id="CC BY-NC-ND 4.0",
        commercial_use=False,
        derivatives_redistributable=False,
        record_sources=("CMsiRNAdb",),
        url="https://cellknowledge.com.cn/CMsiRNAdb/",
        notes=(
            "Non-commercial only, and No Derivatives: only the original "
            "unmodified download may be redistributed, which is why this "
            "subset is derived in code at load time."
        ),
    ),
    "cmsirnadb_full": SourceLicense(
        key="cmsirnadb_full",
        name="CMsiRNAdb, the other 12 genes (He et al. 2026)",
        license_id="CC BY-NC-ND 4.0",
        commercial_use=False,
        derivatives_redistributable=False,
        record_sources=("CMsiRNAdb_full",),
        url="https://cellknowledge.com.cn/CMsiRNAdb/",
        notes=(
            "Same terms as the PCSK9 subset above -- same underlying "
            "database and the same load-time-derivation reason."
        ),
    ),
    "davis2025": SourceLicense(
        key="davis2025",
        name="Davis, Monopoli et al. 2025 (NAR gkaf479)",
        license_id="CC BY 4.0",
        commercial_use=True,
        derivatives_redistributable=True,
        record_sources=("Davis2025",),
        url="https://doi.org/10.1093/nar/gkaf479",
        notes="Attribution required. Largest single permissively-licensed source here.",
    ),
}


def _normalize(license_id: str) -> str:
    """Fold a license id to a comparison key, so callers don't have to
    match the registry's exact punctuation: `"cc-by-nc-nd-4.0"`,
    `"CC_BY_NC_ND_4.0"` and `"CC BY-NC-ND 4.0"` all compare equal. Version
    numbers are NOT folded away -- `"CC BY-NC"` and `"CC BY-NC 4.0"` stay
    distinct ids, because two different sources here state exactly those
    two different things."""
    return re.sub(r"[\s_-]+", " ", license_id).strip().upper()


_BY_NORMALIZED_LICENSE: dict[str, str] = {
    _normalize(entry.license_id): entry.license_id for entry in SOURCE_LICENSES.values()
}


def list_sources() -> list[SourceLicense]:
    """Every loadable source and its license facts, in `load_records()`'s
    own load order."""
    return list(SOURCE_LICENSES.values())


def list_licenses() -> list[str]:
    """Every distinct license id present in this dataset, sorted -- the
    exact strings `load_records(licenses=...)` accepts (case- and
    punctuation-insensitively; see `normalize_license`).

    Includes `LICENSE_UNRESOLVED` ("unresolved"), which is a real, selectable
    id covering the two sources whose terms could not be established.
    """
    return sorted(_BY_NORMALIZED_LICENSE.values())


def normalize_license(license_id: str) -> str:
    """Resolve `license_id` to the canonical id used in the registry,
    accepting any capitalization/punctuation variant (`"cc by 4.0"` ->
    `"CC BY 4.0"`).

    Raises ValueError -- listing every valid id -- for anything not present
    in this dataset, so a typo'd or wishfully-broad selection fails loudly
    instead of silently loading nothing.
    """
    canonical = _BY_NORMALIZED_LICENSE.get(_normalize(license_id))
    if canonical is None:
        valid = ", ".join(list_licenses())
        raise ValueError(
            f"unknown license {license_id!r} -- no source in this dataset carries it. "
            f"Valid license ids (see list_licenses()): {valid}"
        )
    return canonical


def source_keys_for_licenses(licenses: Iterable[str]) -> list[str]:
    """The `load_records()` flag suffixes of every source carrying one of
    `licenses`, in load order.

    Each id is resolved through `normalize_license`, so an unknown license
    raises rather than quietly contributing nothing.
    """
    wanted = {normalize_license(lic) for lic in licenses}
    return [key for key, entry in SOURCE_LICENSES.items() if entry.license_id in wanted]


def license_for_source(source: str) -> SourceLicense | None:
    """The license entry for a `SiRNARecord.source` value, or None if that
    source isn't in the registry (which, for anything `load_records()`
    produced, means a source was added without a registry entry -- the
    completeness test in tests/test_licenses.py exists to catch that).

    Handles the discriminated sources too: `"REMOVED_Sciabola2013"` and
    `"REMOVED_Harborth2001"` both resolve to the REMOVED entry, while
    exact ids always win over a prefix (so `"CMsiRNAdb_full"` resolves to
    its own entry rather than `"CMsiRNAdb"`'s).
    """
    for entry in SOURCE_LICENSES.values():
        if source in entry.record_sources:
            return entry
    for entry in SOURCE_LICENSES.values():
        if any(source.startswith(f"{rs}_") for rs in entry.record_sources):
            return entry
    return None
