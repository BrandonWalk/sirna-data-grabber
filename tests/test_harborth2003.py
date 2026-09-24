"""Harborth et al. 2003 lamin A/C panel (values via Ichihara et al. 2007),
plus the supersession it triggers: siRNAEfficacyDB's corrupted
`Lamin A` block.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from sirna_data.raw_loader import (
    _load_harborth2003_records,
    _load_sirnaefficacydb_records,
    _revcomp,
    load_records,
)

# Must match conftest.FLANK (the value fake_data_dir was built with).
FLANK = 5


def test_load_harborth2003_records(patch_data_dir: Path, fixture_constants):
    records = _load_harborth2003_records(FLANK)
    assert len(records) == 2
    r = {x.row_id: x for x in records}["harborth2003_B1"]
    assert r.gene == "LAMGENE"
    assert r.accession == "LAMACC"
    assert r.source == "Harborth2003"
    assert r.label == pytest.approx(83.0)
    # both strands come from the raw table; the guide is NOT derived by revcomp
    assert r.guide_seq == (_revcomp(fixture_constants.harborth2003_site) + "UU")
    assert r.duplex_len == 19  # 19nt core, the trailing 2nt are overhang
    assert r.has_flanking_context is True
    site = r.mrna_window[r.site_start : r.site_start + r.site_len]
    assert site == fixture_constants.harborth2003_site
    assert "protein" in r.technology
    assert "HeLa" in r.technology


def test_row_without_a_transcript_falls_back_to_duplex_only(patch_data_dir: Path):
    r = {x.row_id: x for x in _load_harborth2003_records(FLANK)}["harborth2003_B2"]
    assert r.has_flanking_context is False
    assert r.site_start == 0


def test_returns_empty_when_files_missing(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    import sirna_data.raw_loader as raw_loader_module

    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()
    monkeypatch.setattr(raw_loader_module, "DATA_DIR", empty_dir)
    assert _load_harborth2003_records(FLANK) == []


# --------------------------------------------------------------------------
# supersessions
# --------------------------------------------------------------------------


def test_sirna_efficacy_loader_skips_superseded_genes(fake_data_dir: Path):
    csv_path = fake_data_dir / "sirna_efficacy.csv"
    fasta_path = fake_data_dir / "mrna_transcripts.fasta"
    everything = _load_sirnaefficacydb_records(csv_path, fasta_path, FLANK)
    assert {r.gene for r in everything} == {"GENEA", "GENEB"}

    kept = _load_sirnaefficacydb_records(
        csv_path, fasta_path, FLANK, superseded_genes=frozenset({"GENEA"})
    )
    assert {r.gene for r in kept} == {"GENEB"}


def test_load_records_supersedes_the_sirna_efficacy_block(
    fake_data_dir: Path, monkeypatch: pytest.MonkeyPatch):
    """With Harborth2003 covering a gene, siRNAEfficacyDB's rows for that
    gene are dropped -- the real-world case being its 44 byte-identical
    `Lamin A` rows."""
    import sirna_data.raw_loader as raw_loader_module

    text = (fake_data_dir / "harborth2003_extra.csv").read_text()
    (fake_data_dir / "harborth2003_extra.csv").write_text(text.replace("LAMGENE", "GENEA"))
    monkeypatch.setattr(raw_loader_module, "DATA_DIR", fake_data_dir)

    records = load_records(flank_nt=FLANK)
    assert {r.source for r in records if r.gene == "GENEA"} == {"Harborth2003"}
    assert len([r for r in records if r.gene == "GENEA"]) == 2


def test_disabling_harborth2003_restores_the_old_sources(patch_data_dir: Path):
    records = load_records(flank_nt=FLANK, include_harborth2003=False)
    assert not [r for r in records if r.source == "Harborth2003"]
    assert {r.source for r in records if r.gene == "GENEA"} == {"siRNAEfficacyDB"}


def test_license_selection_includes_harborth2003(patch_data_dir: Path):
    records = load_records(flank_nt=FLANK, licenses=["CC BY-NC 2.0 UK"])
    assert {r.source for r in records} == {"Harborth2003"}
    assert len(records) == 2
