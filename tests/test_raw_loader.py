from __future__ import annotations

from pathlib import Path

import pytest

from sirna_data.raw_loader import (
    CMSIRNADB_PCSK9_ACCESSION,
    SiRNARecord,
    _dna_to_rna,
    _load_cmsirnadb_full_records,
    _load_cmsirnadb_records,
    _load_monopoli_records,
    _load_shabalina_records,
    _load_sirnaefficacydb_records,
    _locate_window,
    _revcomp,
    load_records,
    read_fasta,
)

# FLANK must match conftest.FLANK (the value fake_data_dir was built with).
FLANK = 5


# --------------------------------------------------------------------------
# small pure-function helpers
# --------------------------------------------------------------------------


def test_dna_to_rna_converts_and_uppercases():
    assert _dna_to_rna("acgtACGT") == "ACGUACGU"


def test_dna_to_rna_leaves_existing_u_alone():
    assert _dna_to_rna("acgu") == "ACGU"


def test_revcomp():
    # complement(A A G G) = U U C C, reversed = C C U U
    assert _revcomp("AAGG") == "CCUU"


def test_revcomp_is_involutive():
    seq = "ACGUACGUACGU"
    assert _revcomp(_revcomp(seq)) == seq


def test_read_fasta_multi_record_and_multiline(tmp_path: Path):
    fasta = tmp_path / "test.fasta"
    fasta.write_text(
        ">ACC1 some description here\n"
        "ACGT\n"
        "ACGT\n"
        ">ACC2\n"
        "TTTT\n"
    )
    sequences = read_fasta(fasta)
    # header key is only the first whitespace-separated token (accession),
    # multi-line sequences are concatenated.
    assert sequences == {"ACC1": "ACGTACGT", "ACC2": "TTTT"}


def test_locate_window_found_with_full_flank():
    transcript = "AAAAA" + "ACGUACGUACGU" + "CCCCC"  # 5 + 12 + 5
    window, site_start, has_context = _locate_window("ACGUACGUACGU", transcript, flank_nt=5)
    assert has_context is True
    assert site_start == 5
    assert window == transcript  # both flanks fully available


def test_locate_window_clips_at_left_boundary():
    site = "ACGUACGUACGU"
    transcript = site + "CCCCC"  # nothing to the left of the site
    window, site_start, has_context = _locate_window(site, transcript, flank_nt=5)
    assert has_context is True
    assert site_start == 0
    assert window == transcript


def test_locate_window_clips_at_right_boundary():
    site = "ACGUACGUACGU"
    transcript = "AAAAA" + site  # nothing to the right of the site
    window, site_start, has_context = _locate_window(site, transcript, flank_nt=5)
    assert has_context is True
    assert site_start == 5
    assert window == transcript


def test_locate_window_not_found_falls_back_to_duplex_only():
    window, site_start, has_context = _locate_window("ACGUACGUACGU", "GGGGGGGGGGGG", flank_nt=5)
    assert has_context is False
    assert site_start == 0
    assert window == "ACGUACGUACGU"


def test_locate_window_handles_missing_transcript():
    window, site_start, has_context = _locate_window("ACGUACGUACGU", None, flank_nt=5)
    assert has_context is False
    assert site_start == 0
    assert window == "ACGUACGUACGU"


# --------------------------------------------------------------------------
# siRNAEfficacyDB (primary source) loader
# --------------------------------------------------------------------------


def test_load_sirnaefficacydb_records(fake_data_dir: Path, fixture_constants):
    records = _load_sirnaefficacydb_records(
        fake_data_dir / "sirna_efficacy.csv", fake_data_dir / "mrna_transcripts.fasta", FLANK
    )
    assert len(records) == 2
    assert all(isinstance(r, SiRNARecord) for r in records)

    located, missing = records
    assert located.gene == "GENEA"
    assert located.accession == "ACC1"
    assert located.has_flanking_context is True
    assert located.site_start == 5
    assert located.label == 55.5
    assert located.technology == "Luciferase reporter assay"
    assert located.source == "siRNAEfficacyDB"
    assert located.mrna_window == fixture_constants.transcript(fixture_constants.sites["primary"])

    # ACC_MISSING has no matching FASTA record -> falls back to duplex-only.
    assert missing.gene == "GENEB"
    assert missing.has_flanking_context is False
    assert missing.mrna_window == "UUUUUUUUUUUU"


# --------------------------------------------------------------------------
# supplementary-source loaders (each reads DATA_DIR directly)
# --------------------------------------------------------------------------


def test_load_monopoli_records(patch_data_dir: Path, fixture_constants):
    records = _load_monopoli_records(FLANK)
    assert len(records) == 1
    r = records[0]
    assert r.gene == "GENEC"
    assert r.source == "Monopoli2023"
    assert r.label == pytest.approx(100.0 - 30.0)
    assert r.has_flanking_context is True
    assert r.guide_seq == _revcomp(fixture_constants.sites["monopoli"])


def test_load_shabalina_records(patch_data_dir: Path):
    records = _load_shabalina_records(FLANK)
    assert len(records) == 1
    r = records[0]
    assert r.gene == "GENED"
    assert r.source == "Shabalina2006"
    assert r.label == pytest.approx(100.0 - 40.0)


def test_load_cmsirnadb_records(patch_data_dir: Path, fixture_constants):
    """PCSK9 subset: fixture has 4 raw rows -- 1 kept, 3 excluded (mouse
    accession, non-human cell type, inclisiran-matching antisense)."""
    records = _load_cmsirnadb_records(FLANK)
    assert len(records) == 1
    r = records[0]
    assert r.gene == "PCSK9"
    assert r.source == "CMsiRNAdb"
    assert r.label == pytest.approx(65.0)
    assert "HepG2" in r.technology
    # every surviving row is normalized to the canonical accession,
    # regardless of what the raw row itself stated.
    assert r.accession == CMSIRNADB_PCSK9_ACCESSION
    assert r.has_flanking_context is True
    assert r.guide_seq == _revcomp(fixture_constants.cmsirnadb_pcsk9_site)


def test_load_cmsirnadb_records_excludes_species_and_inclisiran(patch_data_dir: Path):
    records = _load_cmsirnadb_records(FLANK)
    # the 3 excluded rows (mouse accession 50.0, non-human cell type 55.0,
    # inclisiran-matching 72.0) must not appear.
    labels = {r.label for r in records}
    assert 50.0 not in labels
    assert 55.0 not in labels
    assert 72.0 not in labels


def test_load_cmsirnadb_full_records(patch_data_dir: Path, fixture_constants):
    """Non-PCSK9 genes: fixture has 4 raw GENEF rows -- 2 duplicate-duplex
    measurements (70.0, 90.0) collapsed into 1 record with the median label,
    1 outlier (500.0) dropped, 1 contaminated-sequence row dropped."""
    records = _load_cmsirnadb_full_records(FLANK)
    assert len(records) == 1
    r = records[0]
    assert r.gene == "GENEF"
    assert r.source == "CMsiRNAdb_full"
    assert r.label == pytest.approx(80.0)  # median(70.0, 90.0)
    assert "Hela" in r.technology
    assert r.has_flanking_context is True
    assert r.guide_seq == _revcomp(fixture_constants.cmsirnadb_other_site)


def test_load_cmsirnadb_full_records_dedup_against_existing(
    patch_data_dir: Path, fixture_constants):
    existing = frozenset({fixture_constants.cmsirnadb_other_site})
    records = _load_cmsirnadb_full_records(FLANK, existing_sequences=existing)
    assert records == []


@pytest.mark.parametrize(
    "loader",
    [
        _load_monopoli_records,
        _load_shabalina_records,
        _load_cmsirnadb_records,
        _load_cmsirnadb_full_records,
    ],
)
def test_supplementary_loaders_return_empty_when_files_missing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, loader):
    import sirna_data.raw_loader as raw_loader_module

    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()
    monkeypatch.setattr(raw_loader_module, "DATA_DIR", empty_dir)
    assert loader(FLANK) == []


# --------------------------------------------------------------------------
# load_records: the public, merged entry point
# --------------------------------------------------------------------------


def test_load_records_merges_every_source(patch_data_dir: Path, fake_data_dir: Path):
    records = load_records(
        csv_path=fake_data_dir / "sirna_efficacy.csv",
        fasta_path=fake_data_dir / "mrna_transcripts.fasta",
        flank_nt=FLANK,
    )
    # 2 primary + 1 each of 4 supplementary sources
    assert len(records) == 6
    sources = {r.source for r in records}
    assert sources == {
        "siRNAEfficacyDB",
        "Monopoli2023",
        "Shabalina2006",
        "CMsiRNAdb",
        "CMsiRNAdb_full",
    }


def test_load_records_respects_include_flags(patch_data_dir: Path, fake_data_dir: Path):
    records = load_records(
        csv_path=fake_data_dir / "sirna_efficacy.csv",
        fasta_path=fake_data_dir / "mrna_transcripts.fasta",
        flank_nt=FLANK,
        include_monopoli=False,
        include_shabalina=False,
        include_cmsirnadb=False,
        include_cmsirnadb_full=False,
    )
    # only the 2 primary rows
    assert len(records) == 2
    assert {r.source for r in records} == {"siRNAEfficacyDB"}


def test_load_records_can_exclude_primary_source(patch_data_dir: Path, fake_data_dir: Path):
    # the primary siRNAEfficacyDB set must be excludable too -- no source is
    # forced on unconditionally.
    records = load_records(
        csv_path=fake_data_dir / "sirna_efficacy.csv",
        fasta_path=fake_data_dir / "mrna_transcripts.fasta",
        flank_nt=FLANK,
        include_sirna_efficacy=False,
    )
    assert "siRNAEfficacyDB" not in {r.source for r in records}
    # the other 4 supplementary sources still load by default
    assert len(records) == 4


def test_load_records_all_flags_false_returns_nothing(patch_data_dir: Path, fake_data_dir: Path):
    records = load_records(
        csv_path=fake_data_dir / "sirna_efficacy.csv",
        fasta_path=fake_data_dir / "mrna_transcripts.fasta",
        flank_nt=FLANK,
        include_sirna_efficacy=False,
        include_monopoli=False,
        include_shabalina=False,
        include_cmsirnadb=False,
        include_cmsirnadb_full=False,
    )
    assert records == []


def test_load_records_defaults_to_data_dir(patch_data_dir: Path):
    # csv_path/fasta_path omitted -> should fall back to DATA_DIR/<default filenames>,
    # which patch_data_dir has already pointed at the fixture directory.
    records = load_records(flank_nt=FLANK)
    assert len(records) == 6


def test_load_records_data_dir_arg_without_env_var(
    fake_data_dir: Path, monkeypatch: pytest.MonkeyPatch):
    # No SIRNA_DATA_DIR env var, no monkeypatched module DATA_DIR -- just
    # pass the directory straight to load_records().
    monkeypatch.delenv("SIRNA_DATA_DIR", raising=False)
    records = load_records(flank_nt=FLANK, data_dir=fake_data_dir)
    assert len(records) == 6
    assert {r.source for r in records} == {
        "siRNAEfficacyDB",
        "Monopoli2023",
        "Shabalina2006",
        "CMsiRNAdb",
        "CMsiRNAdb_full",
    }


def test_load_records_data_dir_accepts_str(fake_data_dir: Path, monkeypatch: pytest.MonkeyPatch):
    monkeypatch.delenv("SIRNA_DATA_DIR", raising=False)
    records = load_records(flank_nt=FLANK, data_dir=str(fake_data_dir))
    assert len(records) == 6


# --------------------------------------------------------------------------
# SIRNA_DATA_DIR env var override (module-level constant, needs a reload)
# --------------------------------------------------------------------------


def test_data_dir_env_var_override(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, reload_raw_loader):
    custom_dir = tmp_path / "custom_raw"
    custom_dir.mkdir()
    monkeypatch.setenv("SIRNA_DATA_DIR", str(custom_dir))

    reloaded = reload_raw_loader()
    assert reloaded.DATA_DIR == custom_dir
