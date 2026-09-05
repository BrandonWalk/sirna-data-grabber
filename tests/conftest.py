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
    "pdcd1": "CCGGAACCGGAA",
    "shabalina": "AAACCCGGGUUU",
}

# Davis2025 rows need a 20nt target site (unlike the 12nt sites above) plus
# its own 50nt local window (15nt flank + 20nt site + 15nt flank), since
# that source ships its own window instead of a full-length transcript.
DAVIS2025_SITE = ("ACGU" * 5)[:20]
DAVIS2025_LEFT_FLANK = "GGGGGGGGGGGGGGG"  # 15nt
DAVIS2025_RIGHT_FLANK = "UUUUUUUUUUUUUUU"  # 15nt

# CMsiRNAdb rows need a full 19nt core (CMSIRNADB_CORE_LEN), unlike the
# 12nt sites above.
CMSIRNADB_PCSK9_SITE = ("ACGU" * 5)[:19]
CMSIRNADB_OTHER_SITE = ("GCUA" * 5)[:19]

# Martinelli rows also need a full 19nt core -- only the first 19nt of the
# stored (21nt, overhang-included) "Sequence" field is expected to match the
# transcript; the trailing 2nt overhang is each row's own synthetic tail,
# not genomic sequence, so it must NOT appear in the transcript fixture.
MARTINELLI_SITE = ("GCAU" * 5)[:19]

_2OME_NAMES = {
    "A": "2'-O-Methyladenosine", "C": "2'-O-Methylcytidine",
    "G": "2'-O-Methylguanosine", "U": "2'-O-Methyluridine",
}
_2F_NAMES = {
    "A": "2'-Fluoroadenosine", "C": "2'-Fluorocytidine",
    "G": "2'-Fluoroguanosine", "U": "2'-Fluorouridine",
}


def _mod_types_field(seq: str, modified_positions: set[int], names: dict[str, str]) -> str:
    """Build a `Modification_Types_*_strand`-style cell: bare base letter
    (CMsiRNAdb's own "no annotation" sentinel) at every position not in
    `modified_positions` (1-indexed), the real modified-nucleoside name at
    every position that is."""
    parts = []
    for i, base in enumerate(seq, start=1):
        value = names[base] if i in modified_positions else base
        parts.append(f"{i}*{value}")
    return " || ".join(parts)


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

    (data_dir / "pdcd1_extra.csv").write_text(
        "Experiment_ID,Sequence,Gene,Accession_number,Efficiency_LUC_Pct,Efficiency_QPCR_Pct\n"
        f"163-1,{SITES['pdcd1']},PDCD1,NM_005018,97.0,96.5\n"
    )
    _write_fasta(
        data_dir / "pdcd1_transcripts.fasta", {"NM_005018": _transcript(SITES["pdcd1"])}
    )

    (data_dir / "shabalina_extra.csv").write_text(
        "Sequence,Gene,Accession_number,Activity_Remaining_Pct\n"
        f"{SITES['shabalina']},GENED,ACC4,40.0\n"
    )
    _write_fasta(
        data_dir / "shabalina_transcripts.fasta", {"ACC4": _transcript(SITES["shabalina"])}
    )

    # Martinelli/sirna-reproduction: both strands given directly (not derived by
    # revcomp), each with its own 2nt 3' overhang appended past the 19nt
    # core that's expected to match the transcript. Row 1 exercises the
    # modified path (sense modified, antisense the "0" = unmodified
    # sentinel); row 2 exercises the fully-unmodified ("0"/"0") path.
    martinelli_antisense = raw_loader_module._revcomp(MARTINELLI_SITE)
    (data_dir / "martinelli_extra.csv").write_text(
        "Experiment_ID,PMID,Gene,Accession_number,Sequence,Modification_sense,"
        "Sequence_antisense,Modification_antisense,PCT\n"
        f"SM1,12345678,MARTGENE,MARTACC,{MARTINELLI_SITE}UU,locked nucleic acid,"
        f"{martinelli_antisense}UU,0,75.0\n"
        f"SM2,12345678,MARTGENE,MARTACC,{MARTINELLI_SITE}UU,0,"
        f"{martinelli_antisense}UU,0,20.0\n"
    )
    _write_fasta(
        data_dir / "martinelli_transcripts.fasta", {"MARTACC": _transcript(MARTINELLI_SITE)}
    )

    # CMsiRNAdb: a single raw master TSV feeds both _load_cmsirnadb_records
    # (PCSK9) and _load_cmsirnadb_full_records (everything else) -- see
    # raw_loader.py's module-level CMsiRNAdb note for why there's no
    # pre-filtered CSV here (CC BY-NC-ND forbids redistributing a derivative).
    revcomp = raw_loader_module._revcomp
    incl_core = raw_loader_module.CMSIRNADB_INCLISIRAN_CORE

    # Modification annotation for the surviving rows -- see
    # test_load_cmsirnadb_records / test_load_cmsirnadb_full_records for
    # what these are expected to parse to. Excluded rows get an empty
    # annotation since it's never read for them.
    pcsk9_antisense = revcomp(CMSIRNADB_PCSK9_SITE)
    pcsk9_sense_mods = _mod_types_field(CMSIRNADB_PCSK9_SITE, {1, 2}, _2OME_NAMES)
    pcsk9_antisense_mods = _mod_types_field(pcsk9_antisense, set(), _2OME_NAMES)  # all bare = unmodified
    genef_antisense = revcomp(CMSIRNADB_OTHER_SITE)
    genef_sense_mods = _mod_types_field(
        CMSIRNADB_OTHER_SITE, set(range(1, len(CMSIRNADB_OTHER_SITE) + 1)), _2OME_NAMES
    )
    genef_antisense_mods = _mod_types_field(
        genef_antisense, set(range(1, len(genef_antisense) + 1)), _2F_NAMES
    )

    tsv_columns = (
        "Accession_number\tTarget_Gene\tAntisense_seqence\tSense_seqence\tInhibition\t"
        "Cell_Type\tModification_Types_Sense_strand\tModification_Types_Antisense_strand\n"
    )
    tsv_rows = [
        # kept: normal PCSK9 row, accession deliberately "wrong" (not the
        # canonical one) to check every surviving row gets normalized to it.
        # Sense strand: positions 1-2 modified (2'-OMe), rest unmodified.
        # Antisense strand: fully unannotated (bare-letter sentinel), same
        # as CMsiRNAdb's own real-data pattern for rows with no antisense
        # chemistry info.
        f"NR_110451.3\tPCSK9\t{pcsk9_antisense}\t{CMSIRNADB_PCSK9_SITE}\t65.0\tHepG2\t"
        f"{pcsk9_sense_mods}\t{pcsk9_antisense_mods}\n",
        # excluded: mouse accession (species filter).
        f"NM_153565.2\tPCSK9\t{pcsk9_antisense}\t{CMSIRNADB_PCSK9_SITE}\t50.0\tMus musculus\t\t\n",
        # excluded: non-human hepatocytes cell type (species filter).
        f"NR_110451.3\tPCSK9\t{pcsk9_antisense}\t{CMSIRNADB_PCSK9_SITE}\t55.0\tNon-human hepatocytes\t\t\n",
        # excluded: antisense contains inclisiran's real target core.
        f"NR_110451.3\tPCSK9\t{incl_core}\t{CMSIRNADB_PCSK9_SITE}\t72.0\tHepG2\t\t\n",
        # kept: non-PCSK9 gene, first of two duplicate-duplex measurements
        # (collapsed to one row, label = median). Fully modified on both
        # strands (different chemistry family per strand) so the collapsed
        # record's modification data is unambiguous to assert on.
        f"NM_000001.1\tGENEF\t{genef_antisense}\t{CMSIRNADB_OTHER_SITE}\t70.0\tHela\t"
        f"{genef_sense_mods}\t{genef_antisense_mods}\n",
        f"NM_000001.1\tGENEF\t{genef_antisense}\t{CMSIRNADB_OTHER_SITE}\t90.0\tHela\t"
        f"{genef_sense_mods}\t{genef_antisense_mods}\n",
        # excluded: data-entry outlier (outside [-50, 100]).
        f"NM_000001.1\tGENEF\t{genef_antisense}\t{CMSIRNADB_OTHER_SITE}\t500.0\tHela\t\t\n",
        # excluded: contaminated sequence (non-ACGU character).
        f"NM_000001.1\tGENEF\t{genef_antisense}\t{CMSIRNADB_OTHER_SITE}N\t80.0\tHela\t\t\n",
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

    # Davis2025: ships its own 50nt local window per row instead of a
    # full-length transcript (see _load_davis2025_records's docstring) --
    # no davis2025_transcripts.fasta. Row 1 locates cleanly; row 2's window
    # has a '?' consensus-ambiguity placeholder (exercises the
    # has_flanking_context=False fallback path).
    davis2025_site_2 = ("UGCA" * 5)[:20]
    davis2025_window_2 = DAVIS2025_LEFT_FLANK + davis2025_site_2 + DAVIS2025_RIGHT_FLANK
    davis2025_window_2 = "?" + davis2025_window_2[1:]  # inject one ambiguous base
    (data_dir / "davis2025_extra.csv").write_text(
        "Compound_Name,Gene,Accession_number,Sense_20mer,MRNA_50mer_Window,"
        "Native_Avg_Pct_Untreated,Native_STDEV,Reporter_Avg_Pct_Untreated,Reporter_STDEV,Dataset\n"
        f"GENEE_100_Blunt_2'-OMe/-F,GENEE,ACC5,{DAVIS2025_SITE},"
        f"{DAVIS2025_LEFT_FLANK}{DAVIS2025_SITE}{DAVIS2025_RIGHT_FLANK},25.0,3.0,,,Original\n"
        f"GENEE_200_Asymmetric_2'-OMe Rich,GENEE,ACC5,{davis2025_site_2},"
        f"{davis2025_window_2},40.0,4.0,60.0,5.0,Walk\n"
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
    martinelli_site = MARTINELLI_SITE
    davis2025_site = DAVIS2025_SITE
    davis2025_left_flank = DAVIS2025_LEFT_FLANK
    davis2025_right_flank = DAVIS2025_RIGHT_FLANK
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
