"""Shared fixtures: build a small, self-contained fake data/raw/ directory so
tests never touch the real (100MB+) dataset. Every target-site sequence uses
only A/C/G/U so DNA->RNA conversion (T->U) is a no-op and fixtures can write
the exact same string into both the CSV and the FASTA transcript.
"""
from __future__ import annotations

import importlib
from pathlib import Path

import pytest

import sirna_data.raw_loader as raw_loader_module

FLANK = 5
LEFT_FLANK = "AAAAA"
RIGHT_FLANK = "CCCCC"

# (source key) -> 12nt "site" sequence, distinct per source so fixtures don't
# collide with each other if ever merged.
SITES = {
    "primary": "ACGUACGUACGU",
    "monopoli": "GGGCCCGGGCCC",
    "shabalina": "AAACCCGGGUUU",
}

# CMsiRNAdb rows need a full 19nt core (CMSIRNADB_CORE_LEN), unlike the
# 12nt sites above.
CMSIRNADB_PCSK9_SITE = ("ACGU" * 5)[:19]
CMSIRNADB_OTHER_SITE = ("GCUA" * 5)[:19]


def _transcript(site: str) -> str:
    return LEFT_FLANK + site + RIGHT_FLANK


def _write_fasta(path: Path, records: dict[str, str]) -> None:
    with open(path, "w") as fh:
        for header, seq in records.items():
            fh.write(f">{header}\n{seq}\n")


@pytest.fixture
def fake_data_dir(tmp_path: Path) -> Path:
    """Build a full fake data/raw/ directory: primary + every supplementary
    source, one row each, with sites embedded in a short synthetic
    transcript ("AAAAA" + site + "CCCCC")."""
    data_dir = tmp_path / "raw"
    data_dir.mkdir()

    # Primary siRNAEfficacyDB set: row0 locates cleanly, row1's accession has
    # no matching transcript (exercises the has_flanking_context=False path).
    (data_dir / "sirna_efficacy.csv").write_text(
        "Gene,Accession_number,Antisense_21mer,Sense_19mer,%Inhibition,Technology\n"
        f"GENEA,ACC1,CGUACGUACGUACGUACGUA,{SITES['primary']},55.5,Luciferase reporter assay\n"
        "GENEB,ACC_MISSING,AAAAAAAAAAAAAAAAAAAAA,UUUUUUUUUUUU,10.0,Western blotting\n"
    )
    _write_fasta(data_dir / "mrna_transcripts.fasta", {"ACC1": _transcript(SITES["primary"])})

    (data_dir / "monopoli_extra.csv").write_text(
        "Sequence,Gene,Accession_number,Reporter_Remaining_Pct\n"
        f"{SITES['monopoli']},GENEC,ACC3,30.0\n"
    )
    _write_fasta(
        data_dir / "monopoli_transcripts.fasta", {"ACC3": _transcript(SITES["monopoli"])}
    )

    (data_dir / "shabalina_extra.csv").write_text(
        "Sequence,Gene,Accession_number,Activity_Remaining_Pct\n"
        f"{SITES['shabalina']},GENED,ACC4,40.0\n"
    )
    _write_fasta(
        data_dir / "shabalina_transcripts.fasta", {"ACC4": _transcript(SITES["shabalina"])}
    )

    # CMsiRNAdb: a single raw master TSV feeds both _load_cmsirnadb_records
    # (PCSK9) and _load_cmsirnadb_full_records (everything else) -- see
    # raw_loader.py's module-level CMsiRNAdb note for why there's no
    # pre-filtered CSV here (CC BY-NC-ND forbids redistributing a derivative).
    revcomp = raw_loader_module._revcomp
    incl_core = raw_loader_module.CMSIRNADB_INCLISIRAN_CORE
    tsv_columns = "Accession_number\tTarget_Gene\tAntisense_seqence\tSense_seqence\tInhibition\tCell_Type\n"
    tsv_rows = [
        # kept: normal PCSK9 row, accession deliberately "wrong" (not the
        # canonical one) to check every surviving row gets normalized to it.
        f"NR_110451.3\tPCSK9\t{revcomp(CMSIRNADB_PCSK9_SITE)}\t{CMSIRNADB_PCSK9_SITE}\t65.0\tHepG2\n",
        # excluded: mouse accession (species filter).
        f"NM_153565.2\tPCSK9\t{revcomp(CMSIRNADB_PCSK9_SITE)}\t{CMSIRNADB_PCSK9_SITE}\t50.0\tMus musculus\n",
        # excluded: non-human hepatocytes cell type (species filter).
        f"NR_110451.3\tPCSK9\t{revcomp(CMSIRNADB_PCSK9_SITE)}\t{CMSIRNADB_PCSK9_SITE}\t55.0\tNon-human hepatocytes\n",
        # excluded: antisense contains inclisiran's real target core.
        f"NR_110451.3\tPCSK9\t{incl_core}\t{CMSIRNADB_PCSK9_SITE}\t72.0\tHepG2\n",
        # kept: non-PCSK9 gene, first of two duplicate-duplex measurements
        # (collapsed to one row, label = median).
        f"NM_000001.1\tGENEF\t{revcomp(CMSIRNADB_OTHER_SITE)}\t{CMSIRNADB_OTHER_SITE}\t70.0\tHela\n",
        f"NM_000001.1\tGENEF\t{revcomp(CMSIRNADB_OTHER_SITE)}\t{CMSIRNADB_OTHER_SITE}\t90.0\tHela\n",
        # excluded: data-entry outlier (outside [-50, 100]).
        f"NM_000001.1\tGENEF\t{revcomp(CMSIRNADB_OTHER_SITE)}\t{CMSIRNADB_OTHER_SITE}\t500.0\tHela\n",
        # excluded: contaminated sequence (non-ACGU character).
        f"NM_000001.1\tGENEF\t{revcomp(CMSIRNADB_OTHER_SITE)}\t{CMSIRNADB_OTHER_SITE}N\t80.0\tHela\n",
    ]
    (data_dir / "cmsirnadb_full_raw.tsv").write_text(tsv_columns + "".join(tsv_rows))

    _write_fasta(
        data_dir / "cmsirnadb_transcripts.fasta",
        {"NM_174936.4": _transcript(CMSIRNADB_PCSK9_SITE)},
    )
    _write_fasta(
        data_dir / "cmsirnadb_full_transcripts.fasta",
        {"NM_000001.1": _transcript(CMSIRNADB_OTHER_SITE)},
    )

    return data_dir


@pytest.fixture
def patch_data_dir(monkeypatch: pytest.MonkeyPatch, fake_data_dir: Path) -> Path:
    """Point raw_loader.DATA_DIR at the fake fixture directory for the
    duration of a test (module-level constant, so patched directly rather
    than via env var)."""
    monkeypatch.setattr(raw_loader_module, "DATA_DIR", fake_data_dir)
    return fake_data_dir


class FixtureConstants:
    """Exposes the constants/helpers used to build fake_data_dir, so test
    modules can compute expected values without duplicating them."""

    flank = FLANK
    left_flank = LEFT_FLANK
    right_flank = RIGHT_FLANK
    sites = SITES
    cmsirnadb_pcsk9_site = CMSIRNADB_PCSK9_SITE
    cmsirnadb_other_site = CMSIRNADB_OTHER_SITE
    transcript = staticmethod(_transcript)


@pytest.fixture
def fixture_constants() -> type[FixtureConstants]:
    return FixtureConstants


@pytest.fixture
def reload_raw_loader(monkeypatch: pytest.MonkeyPatch):
    """Yields a function that reloads sirna_data.raw_loader (to pick up a
    freshly-set SIRNA_DATA_DIR env var), and reloads it back to its original
    state afterwards so later tests aren't affected by import-time caching."""
    yield lambda: importlib.reload(raw_loader_module)
    monkeypatch.delenv("SIRNA_DATA_DIR", raising=False)
    importlib.reload(raw_loader_module)
