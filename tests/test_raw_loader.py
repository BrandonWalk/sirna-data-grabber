from __future__ import annotations

from pathlib import Path

import pytest

from sirna_data.raw_loader import (
    CMSIRNADB_PCSK9_ACCESSION,
    SiRNARecord,
    _cmsirnadb_align_modifications,
    _cmsirnadb_chemistry_summary,
    _cmsirnadb_parse_modification_types,
    _dna_to_rna,
    _load_cmsirnadb_full_records,
    _load_cmsirnadb_records,
    _load_martinelli_records,
    _load_monopoli_records,
    _load_pdcd1_records,
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
    # Dataset-level chemistry tag, not per-position (Monopoli's raw table
    # has no position-indexed modification columns the way CMsiRNAdb's does).
    assert r.is_modified is True
    assert r.modification_chemistry is not None
    assert "sdRNA" in r.modification_chemistry
    assert r.sense_modifications is None
    assert r.antisense_modifications is None


def test_load_pdcd1_records(patch_data_dir: Path, fixture_constants):
    records = _load_pdcd1_records(FLANK)
    assert len(records) == 1
    r = records[0]
    assert r.gene == "PDCD1"
    assert r.source == "siRNABERT_PDCD1"
    assert r.label == pytest.approx(96.5)
    assert r.has_flanking_context is True
    assert r.guide_seq == _revcomp(fixture_constants.sites["pdcd1"])
    # Unmodified source (standard synthetic siRNA) -- the schema's defaults.
    assert r.is_modified is False
    assert r.modification_chemistry is None
    assert r.sense_modifications is None
    assert r.antisense_modifications is None


def test_load_shabalina_records(patch_data_dir: Path):
    records = _load_shabalina_records(FLANK)
    assert len(records) == 1
    r = records[0]
    assert r.gene == "GENED"
    assert r.source == "Shabalina2006"
    assert r.label == pytest.approx(100.0 - 40.0)


def test_load_martinelli_records(patch_data_dir: Path, fixture_constants):
    records = _load_martinelli_records(FLANK)
    assert len(records) == 2
    by_id = {r.row_id: r for r in records}

    # Row 1: sense modified ("locked nucleic acid"), antisense the "0" =
    # unmodified sentinel -- still counts as is_modified overall.
    r1 = by_id["martinelli_SM1"]
    assert r1.gene == "MARTGENE"
    assert r1.accession == "MARTACC"
    assert r1.source == "Martinelli_sirna_repro"
    assert r1.label == pytest.approx(75.0)
    assert r1.has_flanking_context is True
    # both strands taken directly from the source, not revcomp-derived.
    assert r1.guide_seq == _revcomp(fixture_constants.martinelli_site) + "UU"
    assert r1.duplex_len == 19
    assert r1.site_len == 19  # 3' overhang excluded from the matched core
    assert r1.is_modified is True
    assert r1.modification_chemistry == "locked nucleic acid (per-molecule, Martinelli/sirna-repro)"
    # per-molecule, not per-position, annotation -- these stay unset.
    assert r1.sense_modifications is None
    assert r1.antisense_modifications is None

    # Row 2: "0"/"0" on both strands -- the fully-unmodified sentinel path.
    r2 = by_id["martinelli_SM2"]
    assert r2.label == pytest.approx(20.0)
    assert r2.is_modified is False
    assert r2.modification_chemistry is None


# --------------------------------------------------------------------------
# CMsiRNAdb modification-annotation parsing (pure functions)
# --------------------------------------------------------------------------


def test_cmsirnadb_parse_modification_types_mixed():
    # position 1 modified, 2-3 bare-letter (unmodified sentinel).
    field = "1*2'-O-Methylcytidine || 2*A || 3*G"
    assert _cmsirnadb_parse_modification_types(field, 3) == (
        "2'-O-Methylcytidine", None, None,
    )


def test_cmsirnadb_parse_modification_types_all_unmodified():
    field = "1*A || 2*C || 3*G"
    assert _cmsirnadb_parse_modification_types(field, 3) == (None, None, None)


def test_cmsirnadb_parse_modification_types_empty_or_missing():
    assert _cmsirnadb_parse_modification_types("", 3) is None
    assert _cmsirnadb_parse_modification_types("nan", 3) is None


def test_cmsirnadb_parse_modification_types_length_mismatch():
    # only 2 positions given but seq_len says 3 -- don't guess, return None.
    assert _cmsirnadb_parse_modification_types("1*A || 2*C", 3) is None


def test_cmsirnadb_align_modifications_slices_to_window():
    # full raw strand is 5nt, but we only stored a 3nt core starting at
    # offset 1 -- the returned tuple should be sliced/aligned to that core.
    seq_full = "AACGU"
    field = "1*A || 2*2'-O-Methyladenosine || 3*C || 4*G || 5*U"
    result = _cmsirnadb_align_modifications(seq_full, "ACG", field)
    assert result == ("2'-O-Methyladenosine", None, None)


def test_cmsirnadb_align_modifications_target_not_found():
    result = _cmsirnadb_align_modifications("AACGU", "UUUU", "1*A || 2*C || 3*G || 4*U || 5*A")
    assert result is None


def test_cmsirnadb_chemistry_summary_none_when_no_modifications():
    assert _cmsirnadb_chemistry_summary(None, None) == (False, None)
    assert _cmsirnadb_chemistry_summary((None, None), None) == (False, None)


def test_cmsirnadb_chemistry_summary_classifies_families():
    sense = ("2'-O-Methyladenosine", None)
    antisense = ("2'-Fluorocytidine", "2'-O-Methyl-3'-Phosphorothioate uridine")
    is_modified, summary = _cmsirnadb_chemistry_summary(sense, antisense)
    assert is_modified is True
    assert summary is not None
    assert "2'-OMe" in summary
    assert "2'-F" in summary
    assert "PS-backbone" in summary
    assert "CMsiRNAdb" in summary


def test_load_cmsirnadb_full_records_missing_modification_columns_is_safe(
    tmp_path: Path, fixture_constants):
    """Backward compat: an old-format raw TSV with no modification columns
    at all must still load cleanly, just with is_modified=False and no
    per-position data -- not an error."""
    data_dir = tmp_path / "raw"
    data_dir.mkdir()
    site = fixture_constants.cmsirnadb_other_site
    revcomp = _revcomp
    (data_dir / "cmsirnadb_full_raw.tsv").write_text(
        "Accession_number\tTarget_Gene\tAntisense_seqence\tSense_seqence\tInhibition\tCell_Type\n"
        f"NM_000001.1\tGENEF\t{revcomp(site)}\t{site}\t70.0\tHela\n"
    )
    with open(data_dir / "cmsirnadb_full_transcripts.fasta", "w") as fh:
        fh.write(f">NM_000001.1\n{('AAAAA' + site + 'CCCCC')}\n")
    records = _load_cmsirnadb_full_records(FLANK, data_dir=data_dir)
    assert len(records) == 1
    r = records[0]
    assert r.is_modified is False
    assert r.modification_chemistry is None
    assert r.sense_modifications is None
    assert r.antisense_modifications is None


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
    # Modification: sense positions 1-2 are 2'-OMe, rest unmodified -- a
    # real tuple of mostly-None, not a bare None ("checked, unmodified
    # here" differs from "no annotation at all"). Antisense in the fixture
    # is bare-letter throughout too, which is CMsiRNAdb's own sentinel for
    # "confirmed unmodified" (not the same as a missing/empty field, which
    # is what actually produces a bare None -- see
    # test_load_cmsirnadb_full_records_missing_modification_columns_is_safe).
    assert r.is_modified is True
    assert r.modification_chemistry == "2'-OMe (per-position, CMsiRNAdb)"
    assert r.sense_modifications is not None
    assert r.sense_modifications[0] == "2'-O-Methyladenosine"
    assert r.sense_modifications[1] == "2'-O-Methylcytidine"
    assert all(m is None for m in r.sense_modifications[2:])
    assert r.antisense_modifications is not None
    assert all(m is None for m in r.antisense_modifications)


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
    # Modification: fully 2'-OMe on sense, fully 2'-F on antisense in the
    # fixture (both duplicate rows share the same annotation -- the
    # representative row's data is what should end up on the collapsed
    # record).
    assert r.is_modified is True
    assert r.modification_chemistry == "2'-OMe/2'-F (per-position, CMsiRNAdb)"
    assert r.sense_modifications is not None
    assert all(m is not None for m in r.sense_modifications)
    assert r.antisense_modifications is not None
    assert all(m is not None for m in r.antisense_modifications)


def test_load_cmsirnadb_full_records_dedup_against_existing(
    patch_data_dir: Path, fixture_constants):
    existing = frozenset({fixture_constants.cmsirnadb_other_site})
    records = _load_cmsirnadb_full_records(FLANK, existing_sequences=existing)
    assert records == []


@pytest.mark.parametrize(
    "loader",
    [
        _load_monopoli_records,
        _load_pdcd1_records,
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
    # 2 primary + 1 each of monopoli/pdcd1/shabalina/cmsirnadb/cmsirnadb_full
    # + 2 martinelli
    assert len(records) == 9
    sources = {r.source for r in records}
    assert sources == {
        "siRNAEfficacyDB",
        "Monopoli2023",
        "siRNABERT_PDCD1",
        "Shabalina2006",
        "Martinelli_sirna_repro",
        "CMsiRNAdb",
        "CMsiRNAdb_full",
    }


def test_load_records_respects_include_flags(patch_data_dir: Path, fake_data_dir: Path):
    records = load_records(
        csv_path=fake_data_dir / "sirna_efficacy.csv",
        fasta_path=fake_data_dir / "mrna_transcripts.fasta",
        flank_nt=FLANK,
        include_monopoli=False,
        include_pdcd1=False,
        include_shabalina=False,
        include_martinelli=False,
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
    # the other supplementary sources still load by default (1 each of
    # monopoli/pdcd1/shabalina/cmsirnadb/cmsirnadb_full + 2 martinelli)
    assert len(records) == 7


def test_load_records_all_flags_false_returns_nothing(patch_data_dir: Path, fake_data_dir: Path):
    records = load_records(
        csv_path=fake_data_dir / "sirna_efficacy.csv",
        fasta_path=fake_data_dir / "mrna_transcripts.fasta",
        flank_nt=FLANK,
        include_sirna_efficacy=False,
        include_monopoli=False,
        include_pdcd1=False,
        include_shabalina=False,
        include_martinelli=False,
        include_cmsirnadb=False,
        include_cmsirnadb_full=False,
    )
    assert records == []


def test_load_records_defaults_to_data_dir(patch_data_dir: Path):
    # csv_path/fasta_path omitted -> should fall back to DATA_DIR/<default filenames>,
    # which patch_data_dir has already pointed at the fixture directory.
    records = load_records(flank_nt=FLANK)
    assert len(records) == 9


def test_load_records_data_dir_arg_without_env_var(
    fake_data_dir: Path, monkeypatch: pytest.MonkeyPatch):
    # No SIRNA_DATA_DIR env var, no monkeypatched module DATA_DIR -- just
    # pass the directory straight to load_records().
    monkeypatch.delenv("SIRNA_DATA_DIR", raising=False)
    records = load_records(flank_nt=FLANK, data_dir=fake_data_dir)
    assert len(records) == 9
    assert {r.source for r in records} == {
        "siRNAEfficacyDB",
        "Monopoli2023",
        "siRNABERT_PDCD1",
        "Shabalina2006",
        "Martinelli_sirna_repro",
        "CMsiRNAdb",
        "CMsiRNAdb_full",
    }


def test_load_records_data_dir_accepts_str(fake_data_dir: Path, monkeypatch: pytest.MonkeyPatch):
    monkeypatch.delenv("SIRNA_DATA_DIR", raising=False)
    records = load_records(flank_nt=FLANK, data_dir=str(fake_data_dir))
    assert len(records) == 9


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
