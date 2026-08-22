"""Load siRNAEfficacyDB (+ supplementary source) records and locate each
target site in its full-length mRNA transcript, producing the (siRNA, local
mRNA window, label) triples that downstream feature-engineering / model code
can turn into whatever representation it needs (e.g. graphs for a GNN,
feature vectors for a classical model, etc.).

This is the reusable "load the data" layer: it only depends on pandas, knows
nothing about torch/PyG/ViennaRNA/RNA-FM, and can be installed standalone by
any project that wants this dataset.
"""
from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from pathlib import Path
from statistics import median

import pandas as pd

# Resolution order for where the raw CSV/FASTA files live:
#   1. SIRNA_DATA_DIR env var, if set (point this at any copy of data/raw/).
#   2. data/raw/ next to this repo checkout (works out of the box when this
#      package is installed editable from the sirna-data-grabber repo, since
#      data/ and src/ are siblings at the repo root).
DATA_DIR = (
    Path(os.environ["SIRNA_DATA_DIR"])
    if os.environ.get("SIRNA_DATA_DIR")
    else Path(__file__).resolve().parents[2] / "data" / "raw"
)
FLANK_NT = 30  # nucleotides of mRNA context kept on each side of the target site

_LICENSE_NOTICE = (
    "sirna_data: this code is MIT licensed, but the DATA it loads is not -- "
    "most sources (siRNAEfficacyDB, CMsiRNAdb) are Creative Commons "
    "Non-Commercial and restrict commercial use. See NOTICE.md in the "
    "sirna-data-grabber repo for the per-source license summary before using "
    "this data commercially. Set SIRNA_DATA_QUIET=1 to silence this message."
)
_license_notice_shown = False


def _maybe_show_license_notice() -> None:
    """Print a one-time-per-process reminder that the data (unlike the code)
    carries mixed, mostly non-commercial licenses. Best-effort documentation
    aid, not a substitute for reading NOTICE.md / data/DATA_SOURCES.md."""
    global _license_notice_shown
    if _license_notice_shown or os.environ.get("SIRNA_DATA_QUIET"):
        return
    print(_LICENSE_NOTICE, file=sys.stderr)
    _license_notice_shown = True


def _dna_to_rna(seq: str) -> str:
    return seq.upper().replace("T", "U")


def read_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    header: str | None = None
    chunks: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    sequences[header] = "".join(chunks)
                header = line[1:].split()[0]
                chunks = []
            elif line:
                chunks.append(line)
        if header is not None:
            sequences[header] = "".join(chunks)
    return sequences


@dataclass
class SiRNARecord:
    row_id: str
    gene: str
    accession: str
    guide_seq: str  # antisense strand (may include a 2nt 3' overhang)
    duplex_len: int  # number of 5' guide positions that actually pair with the target
    mrna_window: str  # local mRNA context, RNA alphabet
    site_start: int  # index of the target site's first nt within mrna_window
    site_len: int
    has_flanking_context: bool
    label: float
    technology: str  # assay type (e.g. "Luciferase reporter assay") -- see data/DATA_SOURCES.md
    source: str  # provenance, e.g. "siRNAEfficacyDB" or "Monopoli2023" -- see data/DATA_SOURCES.md
    # Chemical modification -- unset (False/None) for every source unless a
    # loader below explicitly says otherwise. Most sources in this dataset
    # are standard/unmodified synthetic siRNA; a minority (CMsiRNAdb,
    # Monopoli2023) are chemically modified therapeutic-style constructs.
    # See "Chemical modification data" in data/DATA_SOURCES.md.
    is_modified: bool = False
    # Short human-readable summary of the chemistry class, e.g.
    # "2'-OMe/2'-F/PS (per-position)" or a dataset-level architecture note
    # when no per-position detail is available. None when is_modified=False.
    modification_chemistry: str | None = None
    # Per-position modification name for each nt of `guide_seq`'s
    # corresponding sense strand / of `guide_seq` itself, aligned 1:1 by
    # index (entry i describes position i of the sequence); None at a
    # position means that nt is an unmodified/natural ribonucleotide. Only
    # populated for sources with real per-position annotation in the raw
    # data (currently CMsiRNAdb); None (the whole field, not just entries)
    # when no per-position detail exists for this record.
    sense_modifications: tuple[str | None, ...] | None = None
    antisense_modifications: tuple[str | None, ...] | None = None


def _locate_window(
    sense: str, transcript: str | None, flank_nt: int) -> tuple[str, int, bool]:
    """Find `sense` in `transcript` and return (mrna_window, site_start_in_window,
    has_flanking_context), falling back to duplex-only context if not found."""
    if transcript is None:
        return sense, 0, False
    site_start = transcript.find(sense)
    if site_start == -1:
        return sense, 0, False
    left = max(0, site_start - flank_nt)
    right = min(len(transcript), site_start + len(sense) + flank_nt)
    return transcript[left:right], site_start - left, True


def _revcomp(seq: str) -> str:
    comp = {"A": "U", "U": "A", "G": "C", "C": "G"}
    return "".join(comp[b] for b in reversed(seq))


def _load_sirnaefficacydb_records(
    csv_path: Path, fasta_path: Path, flank_nt: int) -> list[SiRNARecord]:
    """siRNAEfficacyDB (Zhang et al. 2024) -- the primary source. Itself a
    compilation of classic published assays (Huesken et al. 2005 and
    others; see data/DATA_SOURCES.md), but the raw file has no per-row
    author/study column, so it can only be loaded as this one merged
    siRNAEfficacyDB source, not split back out by original study.
    """
    df = pd.read_csv(csv_path)
    transcripts = {acc: _dna_to_rna(seq) for acc, seq in read_fasta(fasta_path).items()}

    records: list[SiRNARecord] = []
    for i, row in df.iterrows():
        guide_seq = row["Antisense_21mer"].upper()
        sense = _dna_to_rna(row["Sense_19mer"])
        transcript = transcripts.get(row["Accession_number"])
        mrna_window, window_site_start, has_flanking_context = _locate_window(
            sense, transcript, flank_nt
        )
        records.append(
            SiRNARecord(
                row_id=f"row{i}",
                gene=row["Gene"],
                accession=row["Accession_number"],
                guide_seq=guide_seq,
                duplex_len=min(19, len(guide_seq)),
                mrna_window=mrna_window,
                site_start=window_site_start,
                site_len=len(sense),
                has_flanking_context=has_flanking_context,
                label=float(row["%Inhibition"]),
                technology=row["Technology"],
                source="siRNAEfficacyDB",
            )
        )
    return records


def _load_monopoli_records(flank_nt: int, data_dir: Path | None = None) -> list[SiRNARecord]:
    """Monopoli et al. 2023 Table S3: 20 chemically-modified siRNAs (a
    cholesterol-conjugated, heavily 2'-F/2'-OMe/phosphorothioate-modified
    "sdRNA" architecture, not a standard unmodified duplex) against 4 genes
    absent from siRNAEfficacyDB. See data/DATA_SOURCES.md for the caveats
    before trusting this the same way as the primary dataset.
    """
    data_dir = data_dir or DATA_DIR
    csv_path = data_dir / "monopoli_extra.csv"
    fasta_path = data_dir / "monopoli_transcripts.fasta"
    if not csv_path.exists() or not fasta_path.exists():
        return []

    df = pd.read_csv(csv_path)
    transcripts = {acc: _dna_to_rna(seq) for acc, seq in read_fasta(fasta_path).items()}

    records: list[SiRNARecord] = []
    for i, row in df.iterrows():
        sense = row["Sequence"].upper()  # verified to match the transcript directly
        guide_seq = _revcomp(sense)
        transcript = transcripts.get(row["Accession_number"])
        mrna_window, window_site_start, has_flanking_context = _locate_window(
            sense, transcript, flank_nt
        )
        records.append(
            SiRNARecord(
                row_id=f"monopoli_row{i}",
                gene=row["Gene"],
                accession=row["Accession_number"],
                guide_seq=guide_seq,
                duplex_len=len(guide_seq),  # fully paired, no 3' overhang in this design
                mrna_window=mrna_window,
                site_start=window_site_start,
                site_len=len(sense),
                has_flanking_context=has_flanking_context,
                label=100.0 - float(row["Reporter_Remaining_Pct"]),
                technology="Dual-luciferase reporter assay (modified sdRNA)",
                source="Monopoli2023",
                is_modified=True,
                # Dataset-level chemistry note, not per-position: the paper
                # describes a single fixed sdRNA architecture applied
                # uniformly to all 20 rows, but doesn't give a per-position
                # modification map the way CMsiRNAdb's raw table does.
                modification_chemistry=(
                    "sdRNA: heavy 2'-F/2'-OMe/phosphorothioate, "
                    "cholesterol-conjugated (dataset-level; no per-position "
                    "map in source -- Monopoli et al. 2023)"
                ),
            )
        )
    return records


def _load_pdcd1_records(flank_nt: int, data_dir: Path | None = None) -> list[SiRNARecord]:
    """PDCD1 (PD-1) 8-siRNA panel, recovered from a deleted
    `data/ExperimentalData/PDCD1_8.csv` in the git history of
    github.com/ChengkuiZhao/siRNABERT (commit f3254ee, before it was
    removed in d2ad931). Dual-readout luciferase + qPCR knockdown assay;
    every sense sequence verified by exact substring match against the
    real NM_005018.3 (PDCD1) mRNA transcript before being trusted -- see
    data/DATA_SOURCES.md. +1 gene (PDCD1), not present in any other
    integrated source.
    """
    data_dir = data_dir or DATA_DIR
    csv_path = data_dir / "pdcd1_extra.csv"
    fasta_path = data_dir / "pdcd1_transcripts.fasta"
    if not csv_path.exists() or not fasta_path.exists():
        return []

    df = pd.read_csv(csv_path)
    transcripts = {acc: _dna_to_rna(seq) for acc, seq in read_fasta(fasta_path).items()}

    records: list[SiRNARecord] = []
    for i, row in df.iterrows():
        sense = _dna_to_rna(row["Sequence"])  # verified to match the transcript directly
        guide_seq = _revcomp(sense)
        transcript = transcripts.get(row["Accession_number"])
        mrna_window, window_site_start, has_flanking_context = _locate_window(
            sense, transcript, flank_nt
        )
        records.append(
            SiRNARecord(
                row_id=f"pdcd1_row{i}",
                gene=row["Gene"],
                accession=row["Accession_number"],
                guide_seq=guide_seq,
                duplex_len=len(guide_seq),  # 19nt, fully paired, no 3' overhang given
                mrna_window=mrna_window,
                site_start=window_site_start,
                site_len=len(sense),
                has_flanking_context=has_flanking_context,
                label=float(row["Efficiency_QPCR_Pct"]),
                technology="Dual-readout (luciferase reporter + qPCR) knockdown assay",
                source="siRNABERT_PDCD1",
            )
        )
    return records


def _load_shabalina_records(flank_nt: int, data_dir: Path | None = None) -> list[SiRNARecord]:
    """Shabalina, Spiridonov & Ogurtsov 2006 (BMC Bioinformatics 7:65)
    Additional File 4: a 653-siRNA heterogeneous compilation, filtered down
    to the 269 rows (41 genes) targeting genes absent from siRNAEfficacyDB --
    see sirna_data.fetch.shabalina for how the kept/excluded genes
    were determined (roughly half of the 653 rows turned out to be
    exact-sequence duplicates of genes we already have) and
    data/DATA_SOURCES.md for the full provenance and caveats.
    """
    data_dir = data_dir or DATA_DIR
    csv_path = data_dir / "shabalina_extra.csv"
    fasta_path = data_dir / "shabalina_transcripts.fasta"
    if not csv_path.exists() or not fasta_path.exists():
        return []

    df = pd.read_csv(csv_path)
    transcripts = {acc: _dna_to_rna(seq) for acc, seq in read_fasta(fasta_path).items()}

    records: list[SiRNARecord] = []
    for i, row in df.iterrows():
        sense = row["Sequence"].upper()
        guide_seq = _revcomp(sense)
        transcript = transcripts.get(row["Accession_number"])
        mrna_window, window_site_start, has_flanking_context = _locate_window(
            sense, transcript, flank_nt
        )
        records.append(
            SiRNARecord(
                row_id=f"shabalina_row{i}",
                gene=row["Gene"],
                accession=row["Accession_number"],
                guide_seq=guide_seq,
                duplex_len=len(guide_seq),  # no overhang recorded in this source (19-mer only)
                mrna_window=mrna_window,
                site_start=window_site_start,
                site_len=len(sense),
                has_flanking_context=has_flanking_context,
                label=100.0 - float(row["Activity_Remaining_Pct"]),
                technology="Heterogeneous compilation (Shabalina et al. 2006)",
                source="Shabalina2006",
            )
        )
    return records


def _load_martinelli_records(flank_nt: int, data_dir: Path | None = None) -> list[SiRNARecord]:
    """Martinelli et al. / sirna-repro (bioRxiv, siRNAmod-derived): a 907-row
    chemically-modified siRNA corpus that ships with NO gene-identity column
    at all (only sequence, per-molecule modification descriptor, PCT, and a
    source PMID/patent ID) -- unusable as-is for leave-one-gene-out
    evaluation. 221 of the 907 rows were traced back to a real target by
    reading each source document's methods section for its stated target,
    then confirming that assignment by exact 19nt substring match of the
    sense strand's non-overhang core against the real target transcript (6
    distinct targets recovered: EGFP, ACP5, APOB, Luciferase_firefly,
    Luciferase_renilla, NPY). The remaining ~686 rows -- most notably a
    single patent that is 60% of the corpus -- could not be resolved this
    way and are not included here; see data/DATA_SOURCES.md for per-target
    provenance and data/DATA_SOURCE_LEDGER.md for what's still unresolved
    and why.

    Unlike every other loader in this file, both strands are taken directly
    from the source (already paired, including each row's own 3' overhang)
    rather than reverse-complementing one from the other -- the source gives
    real per-strand modification detail that a derived antisense would
    lose. Modification descriptors here are per-molecule free text (e.g.
    "locked nucleic acid", "2-fluoro* 2-o-methyl"), not per-position, so
    `sense_modifications`/`antisense_modifications` stay None; only the
    coarser `modification_chemistry` summary is populated.
    """
    data_dir = data_dir or DATA_DIR
    csv_path = data_dir / "martinelli_extra.csv"
    fasta_path = data_dir / "martinelli_transcripts.fasta"
    if not csv_path.exists() or not fasta_path.exists():
        return []

    df = pd.read_csv(csv_path)
    transcripts = {acc: _dna_to_rna(seq) for acc, seq in read_fasta(fasta_path).items()}

    records: list[SiRNARecord] = []
    for _, row in df.iterrows():
        sense = str(row["Sequence"]).strip().upper().replace("T", "U")
        # Only the first 19nt ("core") is expected to match the transcript --
        # the trailing 1-2nt are each row's own synthetic 3' overhang, not
        # genomic sequence.
        sense_core = sense[:19]
        guide_seq = str(row["Sequence_antisense"]).strip().upper().replace("T", "U")
        transcript = transcripts.get(row["Accession_number"])
        mrna_window, window_site_start, has_flanking_context = _locate_window(
            sense_core, transcript, flank_nt
        )

        sense_mod = str(row["Modification_sense"]).strip()
        antisense_mod = str(row["Modification_antisense"]).strip()
        mods = sorted({m for m in (sense_mod, antisense_mod) if m and m != "0"})
        is_modified = len(mods) > 0
        chemistry = (
            " / ".join(mods) + " (per-molecule, Martinelli/sirna-repro)"
            if is_modified
            else None
        )

        records.append(
            SiRNARecord(
                row_id=f"martinelli_{row['Experiment_ID']}",
                gene=row["Gene"],
                accession=row["Accession_number"],
                guide_seq=guide_seq,
                duplex_len=min(19, len(guide_seq)),
                mrna_window=mrna_window,
                site_start=window_site_start,
                site_len=len(sense_core),
                has_flanking_context=has_flanking_context,
                label=float(row["PCT"]),
                technology=(
                    "Reporter/qPCR knockdown assay (chemically modified siRNA; "
                    "Martinelli et al./sirna-repro)"
                ),
                source="Martinelli_sirna_repro",
                is_modified=is_modified,
                modification_chemistry=chemistry,
            )
        )
    return records


# --- CMsiRNAdb (He et al. 2026, BMC Bioinformatics) ------------------------
#
# CMsiRNAdb is CC BY-NC-ND 4.0. The "ND" (No Derivatives) term means we can
# redistribute the original, unmodified download, but NOT a filtered/
# curated/collapsed adaptation of it -- so unlike every other source in this
# file, there is no pre-made "cmsirnadb_*_extra.csv" shipped in data/raw/.
# Instead, `cmsirnadb_full_raw.tsv` (the untouched original patent-literature
# master, 43,153 rows across 13 genes) is the only CMsiRNAdb artifact this
# repo ships, and the two loaders below do all filtering/collapsing here, in
# code, at load time -- every caller reproduces their own local copy of the
# derived data instead of downloading an adaptation from us. See
# data/DATA_SOURCES.md and data/CMSIRNADB_FULL_RETRIEVAL.md.
#
# (The two transcript FASTAs these loaders also read, cmsirnadb_transcripts
# .fasta / cmsirnadb_full_transcripts.fasta, are independently fetched from
# NCBI RefSeq -- public domain, not CMsiRNAdb material -- so they're fine to
# ship as-is.)

CMSIRNADB_CORE_LEN = 19  # duplex core length; extra 3' bases are overhang
CMSIRNADB_MOUSE_ACCESSION = "NM_153565.2"  # Mus musculus Pcsk9
CMSIRNADB_NONHUMAN_CELL_TYPES = frozenset({"Mus musculus", "Non-human hepatocytes"})
CMSIRNADB_PCSK9_ACCESSION = "NM_174936.4"  # canonical human PCSK9 coding transcript
CMSIRNADB_INCLISIRAN_CORE = "AAGCAAAACAGGUCUAGAA"  # LEQVIO's real target -- keep unseen
CMSIRNADB_INHIBITION_RANGE = (-50.0, 100.0)  # outside this: raw data-entry outliers


def _read_cmsirnadb_master(data_dir: Path | None = None) -> pd.DataFrame | None:
    """Read the untouched original CMsiRNAdb patent-literature TSV, if
    present, with %inhibition pre-parsed to numeric (a handful of rows carry
    corrupt values, e.g. -7,103,597 -- handled by CMSIRNADB_INHIBITION_RANGE
    downstream, not here)."""
    data_dir = data_dir or DATA_DIR
    tsv_path = data_dir / "cmsirnadb_full_raw.tsv"
    if not tsv_path.exists():
        return None
    df = pd.read_csv(tsv_path, sep="\t", dtype=str)
    df["_inhibition"] = pd.to_numeric(df["Inhibition"], errors="coerce")
    return df


def _cmsirnadb_locate_core(
    sense_full: str, transcript: str | None, core_len: int = CMSIRNADB_CORE_LEN) -> str | None:
    """CMsiRNAdb's raw Sense_seqence column runs longer than the real
    duplex core for some rows (extra flanking bases from the reported
    chemistry, not always trimmed from the 3' end only). Tries every
    `core_len`-nt window of the sense sequence against `transcript` and
    returns the first one that's a real substring match; falls back to the
    5' `core_len` nt (this project's usual "3' overhang stripped"
    convention) if no window matches or no transcript is available.

    A small fraction of raw rows (~3% of the non-PCSK9 genes) have
    modification-notation characters (parentheses, ambiguity codes like
    'N'/'B'/'V', stray digits) bleeding into this column instead of clean
    sequence -- a data-entry issue in the source, not something introduced
    here. Returns None for any row whose sequence isn't a clean A/C/G/U
    string of at least `core_len` nt, so the caller can skip it the same
    way outlier %inhibition rows are skipped.
    """
    sense_full = _dna_to_rna(sense_full)
    if len(sense_full) < core_len or set(sense_full) - set("ACGU"):
        return None
    if transcript:
        for i in range(len(sense_full) - core_len + 1):
            window = sense_full[i : i + core_len]
            if window in transcript:
                return window
    return sense_full[:core_len]


_CMSIRNADB_BARE_BASES = frozenset({"A", "C", "G", "U", "T"})


def _cmsirnadb_parse_modification_types(
    types_field: str, seq_len: int) -> tuple[str | None, ...] | None:
    """Parse a raw `Modification_Types_{Sense,Antisense}_strand` cell, e.g.
    `"1*2'-O-Methylcytidine || 2*2'-O-Methyladenosine || 3*G || ..."`, into a
    per-position tuple aligned to the FULL raw strand (not any windowed
    core). A bare base letter at a position (the source's own sentinel for
    "no annotation there") maps to None (unmodified); anything else is kept
    as the modified-nucleoside name verbatim, exactly as CMsiRNAdb gives it.

    Returns None (not a tuple of Nones) if the field is missing/empty or
    doesn't cleanly parse into exactly `seq_len` positions numbered 1..N --
    treated as "no usable modification annotation for this row" rather than
    guessed at, matching this module's usual when-the-raw-data-has-a-
    data-entry-issue-on-this-row: skip/fall back, don't guess.
    """
    if not isinstance(types_field, str) or not types_field.strip():
        return None
    parsed: dict[int, str] = {}
    for part in types_field.split("||"):
        part = part.strip()
        pos_str, sep, value = part.partition("*")
        if not sep or not pos_str.strip().isdigit():
            return None
        parsed[int(pos_str.strip())] = value.strip()
    if set(parsed) != set(range(1, seq_len + 1)):
        return None
    return tuple(
        None if parsed[i].upper() in _CMSIRNADB_BARE_BASES else parsed[i]
        for i in range(1, seq_len + 1)
    )


def _cmsirnadb_align_modifications(
    seq_full_rna: str, target: str, types_field: str) -> tuple[str | None, ...] | None:
    """Locate `target` (the core/guide sequence this project actually
    stores) as a substring of `seq_full_rna` (the raw strand the source's
    position-indexed modification annotation is keyed to), and slice the
    parsed per-position modification list down to that window so it lines
    up 1:1 with `target`. Returns None if the annotation doesn't parse or
    `target` isn't found in `seq_full_rna` (e.g. a data-entry mismatch
    between the sequence and modification columns) -- same graceful,
    don't-guess fallback as everywhere else in this loader.
    """
    parsed = _cmsirnadb_parse_modification_types(types_field, len(seq_full_rna))
    if parsed is None:
        return None
    offset = seq_full_rna.find(target)
    if offset == -1:
        return None
    return parsed[offset : offset + len(target)]


def _cmsirnadb_chemistry_summary(
    sense_mods: tuple[str | None, ...] | None,
    antisense_mods: tuple[str | None, ...] | None,
) -> tuple[bool, str | None]:
    """Roll a record's per-position modification tuples up into
    (is_modified, short human-readable chemistry summary). Classifies each
    distinct modified-nucleoside name into a coarse family by substring
    match (2'-OMe, 2'-F, 2'-deoxy, phosphorothioate backbone, 5' vinyl-
    phosphonate cap, inverted-abasic cap, lipid/GalNAc conjugate) so the
    summary stays short and readable even though CMsiRNAdb's raw chemistry
    names are quite granular (e.g. distinguishing all 4 bases per family)."""
    names = [n for n in (*(sense_mods or ()), *(antisense_mods or ())) if n]
    if not names:
        return False, None
    text = " ".join(names).lower()
    families = []
    if "2'-o-methyl" in text:
        families.append("2'-OMe")
    if "2'-fluoro" in text:
        families.append("2'-F")
    if "2'-deoxy" in text:
        families.append("2'-deoxy")
    if "phosphorothioate" in text:
        families.append("PS-backbone")
    if "vinyl phosphonate" in text:
        families.append("5'-VP cap")
    if "inverted abasic" in text:
        families.append("abasic cap")
    if "hexadecyl" in text or "cholesterol" in text or "galnac" in text:
        families.append("lipid/GalNAc-conjugate")
    summary = "/".join(dict.fromkeys(families)) if families else "modified (uncategorized)"
    return True, f"{summary} (per-position, CMsiRNAdb)"


def _load_cmsirnadb_records(flank_nt: int, data_dir: Path | None = None) -> list[SiRNARecord]:
    """CMsiRNAdb, human PCSK9 subset -- derived at load time from the raw
    master TSV (see the module-level CMsiRNAdb note above for why). Only
    PCSK9 rows are used here; the other 12 genes are handled by
    _load_cmsirnadb_full_records.

    - **Species filter**: the raw PCSK9 rows mix human and non-human
      (mostly mouse) data under the same gene label -- excludes rows on
      accession NM_153565.2 (Mus musculus Pcsk9) and rows with
      Cell_Type "Mus musculus" or "Non-human hepatocytes".
    - **Inclisiran exclusion**: PCSK9 is also the target of LEQVIO
      (inclisiran), one of the FDA-approved drugs some downstream
      consumers use for external validation. Any row whose antisense
      sequence contains inclisiran's real 19nt target core is dropped, so
      that drug stays genuinely unseen training data.
    - Every surviving row is located against the ONE canonical human
      coding transcript (NM_174936.4) regardless of its own stated
      accession -- several accessions in the raw data are non-coding
      RefSeq variants that don't contain the real target site (see
      data/DATA_SOURCES.md).
    - Not deduplicated/collapsed (unlike _load_cmsirnadb_full_records):
      repeated measurements of the same duplex are kept as separate rows.
    - **Chemical modification**: the raw table's per-position modification
      columns (real 2'-O-Me/2'-F/2'-deoxy/phosphorothioate/vinyl-phosphonate
      -cap/lipid-conjugate chemistry, not just a "this source is modified"
      flag) are parsed and attached via `is_modified`/`modification_chemistry`
      /`sense_modifications`/`antisense_modifications` when a row's
      annotation cleanly aligns with the sequence we located -- see
      `_cmsirnadb_align_modifications` and data/DATA_SOURCES.md.
    """
    data_dir = data_dir or DATA_DIR
    fasta_path = data_dir / "cmsirnadb_transcripts.fasta"
    df = _read_cmsirnadb_master(data_dir)
    if df is None or not fasta_path.exists():
        return []
    transcripts = {acc: _dna_to_rna(seq) for acc, seq in read_fasta(fasta_path).items()}
    canonical_transcript = transcripts.get(CMSIRNADB_PCSK9_ACCESSION)

    pcsk9 = df[df["Target_Gene"] == "PCSK9"]
    is_nonhuman = (pcsk9["Accession_number"] == CMSIRNADB_MOUSE_ACCESSION) | (
        pcsk9["Cell_Type"].isin(CMSIRNADB_NONHUMAN_CELL_TYPES)
    )
    pcsk9 = pcsk9[~is_nonhuman]

    records: list[SiRNARecord] = []
    for i, row in pcsk9.iterrows():
        antisense_rna = _dna_to_rna(str(row["Antisense_seqence"]))
        if CMSIRNADB_INCLISIRAN_CORE in antisense_rna:
            continue

        sense = _cmsirnadb_locate_core(row["Sense_seqence"], canonical_transcript)
        if sense is None:
            continue
        guide_seq = _revcomp(sense)
        mrna_window, window_site_start, has_flanking_context = _locate_window(
            sense, canonical_transcript, flank_nt
        )
        sense_mods = _cmsirnadb_align_modifications(
            _dna_to_rna(str(row["Sense_seqence"])), sense,
            str(row.get("Modification_Types_Sense_strand", "")),
        )
        antisense_mods = _cmsirnadb_align_modifications(
            _dna_to_rna(str(row["Antisense_seqence"])), guide_seq,
            str(row.get("Modification_Types_Antisense_strand", "")),
        )
        is_modified, chemistry = _cmsirnadb_chemistry_summary(sense_mods, antisense_mods)
        records.append(
            SiRNARecord(
                row_id=f"cmsirnadb_row{i}",
                gene="PCSK9",
                accession=CMSIRNADB_PCSK9_ACCESSION,
                guide_seq=guide_seq,
                duplex_len=len(guide_seq),  # overhang length not separately modeled
                mrna_window=mrna_window,
                site_start=window_site_start,
                site_len=len(sense),
                has_flanking_context=has_flanking_context,
                label=float(row["_inhibition"]),
                technology=f"CMsiRNAdb patent-derived, chemically modified ({row['Cell_Type']})",
                source="CMsiRNAdb",
                is_modified=is_modified,
                modification_chemistry=chemistry,
                sense_modifications=sense_mods,
                antisense_modifications=antisense_mods,
            )
        )
    return records


def _load_cmsirnadb_full_records(
    flank_nt: int,
    existing_sequences: frozenset[str] = frozenset(),
    data_dir: Path | None = None,
) -> list[SiRNARecord]:
    """CMsiRNAdb, the other 12 genes (AGT, ANGPTL3, APP, CTNNB1, HSD17B13,
    INHBE, LPA, MAPT, MARC1, MSTN, PLN, PNPLA3) -- derived at load time from
    the same raw master TSV as _load_cmsirnadb_records (see the module-level
    CMsiRNAdb note above). PCSK9 is deliberately excluded here so
    _load_cmsirnadb_records stays the single source of PCSK9 rows and its
    inclisiran held-out logic isn't bypassed.

    - **Outlier removal**: a handful of raw rows have corrupt %inhibition
      values (e.g. -7,103,597); anything outside CMSIRNADB_INHIBITION_RANGE
      is dropped as a data-entry error, not a real measurement.
    - **Dedup against sequences already loaded from other sources**
      (`existing_sequences`, checked strand-agnostically): skips rows that
      would just duplicate a guide/sense sequence already present elsewhere
      in the dataset, so this addition is "net-new" data rather than
      padding out an existing example.
    - **Collapsed to one row per unique (gene, accession, guide/sense
      duplex)**: repeated measurements of the same duplex (different dose/
      time/cell conditions) are combined; the label is the MEDIAN
      %inhibition across a duplex's replicates (robust to the
      noise-around-zero that produces small negative readings).
      `technology` is taken from one representative row's Cell_Type.
    - Sequence = the 19nt sense core (3' overhangs stripped), located
      against the row's own stated-accession transcript where we have it
      fetched, else falls back to duplex-only context like every other
      source's unmapped rows.
    """
    data_dir = data_dir or DATA_DIR
    fasta_path = data_dir / "cmsirnadb_full_transcripts.fasta"
    df = _read_cmsirnadb_master(data_dir)
    if df is None or not fasta_path.exists():
        return []
    transcripts = {acc: _dna_to_rna(seq) for acc, seq in read_fasta(fasta_path).items()}

    others = df[df["Target_Gene"] != "PCSK9"].copy()
    lo, hi = CMSIRNADB_INHIBITION_RANGE
    others = others[others["_inhibition"].between(lo, hi)]

    groups: dict[tuple[str, str, str], list[float]] = {}
    representative: dict[tuple[str, str, str], pd.Series] = {}
    for _, row in others.iterrows():
        transcript = transcripts.get(row["Accession_number"])
        sense = _cmsirnadb_locate_core(row["Sense_seqence"], transcript)
        if sense is None:
            continue
        guide_seq = _revcomp(sense)
        if guide_seq in existing_sequences or sense in existing_sequences:
            continue
        key = (row["Target_Gene"], row["Accession_number"], sense)
        groups.setdefault(key, []).append(row["_inhibition"])
        representative.setdefault(key, row)

    records: list[SiRNARecord] = []
    for i, (key, inhibitions) in enumerate(groups.items()):
        gene, accession, sense = key
        row = representative[key]
        guide_seq = _revcomp(sense)
        transcript = transcripts.get(accession)
        mrna_window, window_site_start, has_flanking_context = _locate_window(
            sense, transcript, flank_nt
        )
        # Modification annotation comes from the representative row only
        # (same one `technology`'s Cell_Type is taken from) -- replicate
        # measurements of one duplex share the same chemical entity, so
        # this doesn't lose per-record information the way collapsing the
        # label itself would.
        sense_mods = _cmsirnadb_align_modifications(
            _dna_to_rna(str(row["Sense_seqence"])), sense,
            str(row.get("Modification_Types_Sense_strand", "")),
        )
        antisense_mods = _cmsirnadb_align_modifications(
            _dna_to_rna(str(row["Antisense_seqence"])), guide_seq,
            str(row.get("Modification_Types_Antisense_strand", "")),
        )
        is_modified, chemistry = _cmsirnadb_chemistry_summary(sense_mods, antisense_mods)
        records.append(
            SiRNARecord(
                row_id=f"cmsirnadb_full_row{i}",
                gene=gene,
                accession=accession,
                guide_seq=guide_seq,
                duplex_len=len(guide_seq),  # overhang length not separately modeled
                mrna_window=mrna_window,
                site_start=window_site_start,
                site_len=len(sense),
                has_flanking_context=has_flanking_context,
                label=float(median(inhibitions)),
                technology=f"CMsiRNAdb patent-derived, chemically modified ({row['Cell_Type']})",
                source="CMsiRNAdb_full",
                is_modified=is_modified,
                modification_chemistry=chemistry,
                sense_modifications=sense_mods,
                antisense_modifications=antisense_mods,
            )
        )
    return records


def load_records(
    csv_path: Path | None = None,
    fasta_path: Path | None = None,
    flank_nt: int = FLANK_NT,
    data_dir: Path | str | None = None,
    include_sirna_efficacy: bool = True,
    include_monopoli: bool = True,
    include_pdcd1: bool = True,
    include_shabalina: bool = True,
    include_martinelli: bool = True,
    include_cmsirnadb: bool = True,
    include_cmsirnadb_full: bool = True ) -> list[SiRNARecord]:
    """Load the full merged siRNA-efficacy dataset as a list of SiRNARecord.

    Each record pairs an siRNA guide sequence with the local mRNA window
    around its real target site (located by exact substring search in the
    full-length transcript, with graceful fallback to duplex-only context
    when a site can't be located) and its experimentally measured knockdown
    label. See data/DATA_SOURCES.md for full provenance of every source.

    Every source -- including the primary siRNAEfficacyDB set -- is gated
    behind its own `include_*` flag, so any combination of sources
    (including none) can be loaded; none is forced on unconditionally.

    `data_dir` points every source (except any of `csv_path`/`fasta_path`
    given explicitly, which still win for the primary source) at a specific
    directory of fetched files, e.g. `load_records(data_dir="./my_data")`.
    This is a plain function argument -- no `SIRNA_DATA_DIR` env var or
    export required. If omitted, falls back to `SIRNA_DATA_DIR` if set, else
    the package's default relative `data/raw/` location (see `DATA_DIR`).

    Licensing: this function's code is MIT licensed, but the DATA it returns
    is not -- most sources are Creative Commons Non-Commercial and restrict
    commercial use (see NOTICE.md). Calling this prints a one-time reminder
    to stderr (silence it with SIRNA_DATA_QUIET=1).
    """
    _maybe_show_license_notice()
    resolved_dir = Path(data_dir) if data_dir is not None else DATA_DIR
    csv_path = csv_path or resolved_dir / "sirna_efficacy.csv"
    fasta_path = fasta_path or resolved_dir / "mrna_transcripts.fasta"

    records: list[SiRNARecord] = []
    if include_sirna_efficacy:
        records += _load_sirnaefficacydb_records(csv_path, fasta_path, flank_nt)
    if include_monopoli:
        records += _load_monopoli_records(flank_nt, resolved_dir)
    if include_pdcd1:
        records += _load_pdcd1_records(flank_nt, resolved_dir)
    if include_shabalina:
        records += _load_shabalina_records(flank_nt, resolved_dir)
    if include_martinelli:
        records += _load_martinelli_records(flank_nt, resolved_dir)
    if include_cmsirnadb:
        records += _load_cmsirnadb_records(flank_nt, resolved_dir)
    if include_cmsirnadb_full:
        # Strand-agnostic sequence index of everything loaded so far, so the
        # 12-gene CMsiRNAdb addition only contributes genuinely new
        # sequences (see _load_cmsirnadb_full_records's docstring).
        existing_sequences: set[str] = set()
        for r in records:
            core = r.guide_seq[: r.duplex_len]
            existing_sequences.add(core)
            existing_sequences.add(_revcomp(core))
        records += _load_cmsirnadb_full_records(
            flank_nt, frozenset(existing_sequences), resolved_dir
        )
    return records


if __name__ == "__main__":
    recs = load_records()
    n_with_context = sum(r.has_flanking_context for r in recs)
    n_genes = len(set(r.gene for r in recs))
    print(
        f"Loaded {len(recs)} records across {n_genes} genes; "
        f"{n_with_context} with real flanking mRNA context"
    )
