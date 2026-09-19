"""Sciabola et al. 2013 Supplementary Tables S3/S4 loader.

Every fixture path is the synthetic one built in conftest.fake_data_dir --
these tests never touch the real data/raw/ files.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from sirna_data.raw_loader import (
    SCIABOLA2013_CORE_LEN,
    _load_sciabola2013_records,
    _revcomp,
    load_records,
)

# Must match conftest.FLANK (the value fake_data_dir was built with).
FLANK = 5


def test_load_sciabola2013_records(patch_data_dir: Path, fixture_constants):
    records = _load_sciabola2013_records(FLANK)
    assert len(records) == 3
    by_id = {r.row_id: r for r in records}

    r = by_id["sciabola2013_19p2_1"]
    assert r.gene == "SCIAGENE"
    assert r.accession == "SCIACC"
    assert r.source == "Sciabola2013"
    # sense stored in the CSV, guide derived by revcomp like every other loader
    assert r.guide_seq == _revcomp(fixture_constants.sciabola2013_site)
    assert r.duplex_len == 21
    assert r.site_len == 21
    assert r.has_flanking_context is True
    site = r.mrna_window[r.site_start : r.site_start + r.site_len]
    assert site == fixture_constants.sciabola2013_site
    assert "Hep3B" in r.technology
    assert "QuantiGene" in r.technology


def test_label_is_the_mean_of_the_doses_present(patch_data_dir: Path):
    by_id = {r.row_id: r for r in _load_sciabola2013_records(FLANK)}
    # 4 doses: mean(10, 20, 30, 40)
    assert by_id["sciabola2013_19p2_1"].label == pytest.approx(25.0)
    assert "mean of 4 dose(s)" in by_id["sciabola2013_19p2_1"].technology
    # 10 nM only -- the follow-up-gene pattern
    assert by_id["sciabola2013_19p2_2"].label == pytest.approx(80.0)
    assert "mean of 1 dose(s)" in by_id["sciabola2013_19p2_2"].technology
    # 2 of 4 doses present: mean(5, 15), NOT a mean over four columns
    assert by_id["sciabola2013_19p2_3"].label == pytest.approx(10.0)


def test_falls_back_to_the_19nt_core_when_the_overhang_is_not_in_the_transcript(
    patch_data_dir: Path, fixture_constants):
    r = {x.row_id: x for x in _load_sciabola2013_records(FLANK)}["sciabola2013_19p2_2"]
    assert r.has_flanking_context is True
    assert r.site_len == SCIABOLA2013_CORE_LEN
    site = r.mrna_window[r.site_start : r.site_start + r.site_len]
    assert site == fixture_constants.sciabola2013_core_only_site[2:]
    # the record still carries the full 21nt duplex
    assert r.duplex_len == 21


def test_unlocatable_row_falls_back_to_duplex_only_context(patch_data_dir: Path):
    r = {x.row_id: x for x in _load_sciabola2013_records(FLANK)}["sciabola2013_19p2_3"]
    assert r.has_flanking_context is False
    assert r.mrna_window == _revcomp(r.guide_seq)
    assert r.site_start == 0


def test_dedup_against_existing_sequences(patch_data_dir: Path, fixture_constants):
    existing = frozenset({fixture_constants.sciabola2013_site})
    records = _load_sciabola2013_records(FLANK, existing_sequences=existing)
    assert "sciabola2013_19p2_1" not in {r.row_id for r in records}
    assert len(records) == 2


def test_returns_empty_when_files_missing(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    import sirna_data.raw_loader as raw_loader_module

    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()
    monkeypatch.setattr(raw_loader_module, "DATA_DIR", empty_dir)
    assert _load_sciabola2013_records(FLANK) == []


def test_license_selection_includes_sciabola2013(patch_data_dir: Path):
    records = load_records(flank_nt=FLANK, licenses=["CC BY-NC 3.0"])
    assert {r.source for r in records} == {"Sciabola2013"}
    assert len(records) == 3
