# CMsiRNAdb full retrieval — the other 12 genes

Task: get MORE data, EXCLUSIVELY siRNA (synthetic siRNA duplexes with numeric knockdown).
No shRNA, no dsRNA, no miRNA, no ASO.

## Source
`https://www.cellknowledge.com.cn/CMsiRNAdb/` (Zhang lab, Chengdu Univ. of TCM / UESTC),
file `download/CMsiRNA_data_update.tsv`. CMsiRNAdb = Chemically-Modified siRNA database,
curated from patent literature. Contacts: yangzhang@cdutcm.edu.cn, zhy1001@alu.uestc.edu.cn.
The 13 per-gene `patent_dataset_*.tsv` files were also downloaded and confirmed to be an
exact subset of the master (0 additional rows), so only the master is kept.

**License: CC BY-NC-ND 4.0 — "No Derivatives".** This repo only ships the untouched
original master TSV, `data/raw/cmsirnadb_full_raw.tsv` (43,153 rows, all with numeric
% inhibition, 13 target genes: AGT, ANGPTL3, APP, CTNNB1, HSD17B13, INHBE, LPA, MAPT,
MARC1, MSTN, PCSK9, PLN, PNPLA3). Everything described below — outlier removal,
transcript location, deduplication, collapsing repeat measurements — happens at
*load time* in `_load_cmsirnadb_full_records()` (`src/sirna_data/raw_loader.py`), not
as a precomputed file, so no filtered/adapted CMsiRNAdb derivative is redistributed.
Anyone using `sirna_data` reproduces this derivation locally from the original download.

## What's covered vs. what's in the PCSK9 subset
This dataset separately integrates the **PCSK9 slice** of CMsiRNAdb (3,107 rows) via
`_load_cmsirnadb_records()` — see `DATA_SOURCES.md`. This document covers the other
**12 genes, 40,046 rows**, via `_load_cmsirnadb_full_records()`.

## Derivation pipeline (all in `_load_cmsirnadb_full_records()`)
1. **Start**: 40,046 raw non-PCSK9 rows.
2. **Outlier removal**: drop rows with `%inhibition` outside `[-50, 100]` (the source
   has a few corrupt values, e.g. `-7,103,597`) → **39,686 rows** (−360).
3. **Contamination filter**: drop rows whose sense sequence contains non-ACGU
   characters after RNA conversion (modification notation, ambiguity codes, or stray
   characters bleeding into the sequence column — parentheses, `I`/`N`/`B`/`V`, digits,
   newlines) → **38,345 rows** (−1,341).
4. **Cross-source dedup**: any row whose guide or sense sequence (strand-agnostic)
   already exists in the baseline dataset (siRNAEfficacyDB + Monopoli2023 +
   Shabalina2006 + CMsiRNAdb PCSK9) is skipped, so this addition only contributes
   genuinely new sequences. Verified: **0 rows excluded** here against the real
   baseline set — the 12-gene slice has no sequence overlap with the rest of the data.
5. **Transcript location**: each row's sense sequence is searched (sliding 19nt
   window) against its gene's NCBI RefSeq transcript
   (`data/raw/cmsirnadb_full_transcripts.fasta`, 34 records fetched via NCBI efetch);
   the located 19nt core is the target site. Rows with no window match fall back to
   the first 19nt of the sense sequence.
6. **Collapse to unique duplexes**: downstream feature engineering (e.g. graph
   construction) typically carries no dose/time/cell features, so raw patent
   measurements would otherwise become many near-identical training examples.
   Collapsed to one record per unique `(gene, accession, target site)`; label =
   **median % inhibition** across that duplex's replicate measurements (robust to
   the noise-around-zero that produces small negative readings) → **9,601 unique
   duplex records**.

**6,117 of 9,601 records (63.7%)** locate with full 30nt flanking mRNA context; the
rest fall back to duplex-only context, same mechanism as every other source.

Records by gene: PNPLA3 2,066 / HSD17B13 1,985 / APP 952 / AGT 872 / MARC1 823 /
INHBE 670 / MAPT 630 / LPA 556 / ANGPTL3 551 / CTNNB1 352 / PLN 135 / MSTN 9.

## siRNAEfficacyDB (also checked, NOT re-added)
Pulled `siRNAEfficacyDB/download/siRNA_all.txt` (the canonical source of the classic sets:
Huesken 2005 x2,431, Katoh 2007 x702, Reynolds x244, Vickers x76, Harborth x44, Ui-Tei x37,
Khvorova x10 = 3,544 numeric rows). Dedup showed **only 12 new rows** —
`sirna_efficacy.csv` here already IS this dataset. Nothing added.

## Trainability note
These are synthetic siRNA duplexes with directly-measured % inhibition — the SAME target
quantity as this dataset's existing numeric-%KD data. Unlike the DEMETER2/shRNA screens,
this IS a drop-in extension of the supervised set. Caveat: heavily chemically modified and
patent-derived (dose/assay conditions vary), and highly redundant at the sequence level
(40,046 raw rows -> 9,601 unique duplexes) — group-aware CV by sequence is advisable to
avoid leakage. `technology` is tagged the same way as the PCSK9 subset:
`"CMsiRNAdb patent-derived, chemically modified (<cell type>)"`.

## Files
- `data/raw/cmsirnadb_full_raw.tsv` — full 43,153-row CMsiRNAdb master (the only
  CMsiRNAdb artifact shipped; PCSK9 + these 12 genes are both derived from it at load
  time — see the ND-compliance note above)
- `data/raw/cmsirnadb_full_transcripts.fasta` — 34 mRNA transcripts (NCBI efetch)
- `data/cmsirnadb_new_sirna.png` — %inhibition + per-gene figure

## Integration into the loader
Wired into `load_records()` via `include_cmsirnadb_full=True`
(`src/sirna_data/raw_loader.py`, `_load_cmsirnadb_full_records()`), which builds a
strand-agnostic index of every sequence already loaded from the other sources before
running, so this addition only ever contributes new sequences (step 4 above).

Effect on the dataset (verified by running `load_records()`):
- Without: **6,577 records / 87 genes**
- With: **16,178 records / 97 genes** (+9,601 records, +10 new genes:
  AGT ANGPTL3 CTNNB1 HSD17B13 INHBE LPA MARC1 MSTN PLN PNPLA3; APP & MAPT deepened)

Anyone consuming this data with a cached/precomputed downstream representation
(e.g. a graph cache) should rebuild it after pulling this update, since the
underlying record set changed.
