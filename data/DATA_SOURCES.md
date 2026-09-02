# Data sources

## Summary

Audit of every siRNA-efficacy source integrated into this dataset: whether it provided
**trainable data** (siRNA sequence paired with a numeric knockdown/efficacy value
suitable for a regression/classification target), and how much. For the broader
landscape of sources considered — including ones never obtained or not pursued
further — see `POTENTIAL_DATA_SOURCES.md`.

- **Trainable records currently integrated: 7,505** across **103 genes**, all with a
  numeric %-knockdown label. (18,072 records / 113 genes if the CMsiRNAdb full-database
  retrieval and Davis2025 below are also included — both on by default via
  `include_cmsirnadb_full=True` / `include_davis2025=True`; see their own sections.)
- **7 sources** supply that headline 7,505; **siRNAEfficacyDB (3,532)** and the
  **CMsiRNAdb PCSK9 subset (2,756)** are 84% of it.
- **siRecords** was recovered (3,117 rated records, ~1,400 new gene accessions) but is
  **NOT trainable as-is** — its label is a 4-level ordinal rating, not a numeric
  percentage. See "Considered, not integrated: siRecords" below.
- **Chemical modification data**: `SiRNARecord` carries `is_modified`/
  `modification_chemistry`/`sense_modifications`/`antisense_modifications` — see
  "Chemical modification data" below.

### Integrated — trainable data (7,505 records)

| Source | Records | New genes | Metric | Notes |
|---|---|---|---|---|
| **siRNAEfficacyDB** (Zhang 2024) | **3,532** | 41 (baseline) | numeric %Inhibition | Primary source; `sirna_efficacy.csv`. |
| **CMsiRNAdb — PCSK9 subset** (He 2026) | **2,756** | +1 | numeric inhibition | Patent-derived, chemically modified; derived at load time from `cmsirnadb_full_raw.tsv` (CC BY-NC-ND, no derivative file shipped). |
| **Shabalina 2006** | **269** | +41 | numeric (100−Activity) | 269 new after dedup vs 653 in paper; `shabalina_extra.csv`. |
| **Monopoli 2023** | **20** | +4 (APP/MAPT/BACE1/SNCA) | numeric (100−reporter) | Modified sdRNA; `monopoli_extra.csv`. |
| **PDCD1 panel (Xu/Zhao 2024, siRNABERT)** | **8** | +1 (PDCD1) | numeric (qPCR knockdown %) | Recovered from deleted file in repo's git history, verified against NM_005018.3; license unresolved but `pdcd1_extra.csv` is committed at the user's explicit request — see `NOTICE.md`. |
| **Martinelli 2023 / sirna-repro** | **577** | +12 (EGFP, ACP5, APOB, Luciferase_firefly, Luciferase_renilla, NPY, VEGFA, KAZRIN, MIR155HG, CDKN1B, SOD2_chimp, CASR_rhesus) | numeric (PCT) | 577 of 907 rows resolved — 253 by source-document tracing, 324 more by brute-force sequence matching. See the section below. |
| **OligoGraph repo (Sciabola 2013 + Harborth 2001)** | **343** | +3 (HIF1A, HK2, HPSE; the 43 Harborth-traced rows merge into the existing Lamin A group) | `label * 100`, **interpreted as %KD** (not independently verified — see the section below) | Gene identity independently verified by exact 19nt match; the metric is a stated interpretation of an undocumented normalized value, not a confirmed true value. |
| **TOTAL** | **7,505** | **103 genes** | — | — |

### Obtained, not integrated

| Source | Records | Why not trainable | Path to make it trainable |
|---|---|---|---|
| **siRecords** (05 release) | 3,295 new by seq / 3,117 rated | Efficacy is **ordinal** (Very high/High/Medium/Low), not numeric %. Also: redistribution rights unresolved — see `NOTICE.md` and the dedicated section below | Trace the per-record PubMed IDs to primary papers and extract reported %KD |

### Definition of "trainable" used here

A record is trainable here if it pairs an siRNA **sequence** with a **numeric**
knockdown/inhibition value (continuous %), suitable as a regression or classification
target. Ordinal ratings (siRecords) and other proxy metrics are excluded on that
definition even though they encode efficacy information — using them would require
either a different model formulation or converting the metric first.

## siRNA sequences + knockdown efficacy: siRNAEfficacyDB

- Zhang et al. 2024, *IET Systems Biology*, "siRNAEfficacyDB: An experimentally
  supported small interfering RNA efficacy database".
  https://cellknowledge.com.cn/siRNAEfficacy
- 3,544 siRNA records across 42 human genes, compiled from published assays
  (Huesken et al. 2005 Nat. Biotechnol. and others).
- License: Creative Commons Attribution Non-Commercial (CC BY-NC) — free to use,
  distribute, and reproduce for non-commercial purposes with attribution to the
  original authors and underlying studies.
- Fetched by `sirna_data.fetch.sirna_efficacy` (`sirna-data-fetch`) into
  `data/raw/sirna_efficacy.csv`.

We deliberately did NOT use the merged CSVs from the OligoFormer GitHub repo
(github.com/lulab/OligoFormer), even though it bundles the same classic
Huesken/Reynolds/Vickers/etc. benchmarks. That repo ships under a proprietary
Tsinghua University license that prohibits redistribution and use in competing
products, so its prepared data files are not safe to vendor into this project.

## Full-length mRNA transcripts: NCBI Nucleotide (RefSeq/GenBank)

- Fetched by accession number (the `Accession_number` column in
  siRNAEfficacyDB) via NCBI E-utilities (`efetch`), a free public API.
  https://www.ncbi.nlm.nih.gov/books/NBK25501/
- Public-domain/freely reusable sequence records, no license restriction.
- Saved to `data/raw/mrna_transcripts.fasta`.

## Supplementary siRNA data: Monopoli et al. 2023

- Monopoli, Korkin & Khvorova 2023, *Molecular Therapy Nucleic Acids*,
  "Asymmetric trichotomous partitioning overcomes dataset limitations in
  building machine learning models for predicting siRNA efficacy"
  (doi:10.1016/j.omtn.2023.06.010). CC BY 4.0.
- Table S3: 20 siRNAs against 4 genes (`APP`, `MAPT`, `BACE1`, `SNCA`) not
  present in siRNAEfficacyDB — a real, if small, extension of gene coverage
  (41 → 45 genes) for leave-one-gene-out CV.
- Retrieved via Europe PMC's public `supplementaryFiles` REST API
  (`https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10338369/supplementaryFiles`),
  a legitimate, documented bulk-access endpoint — not scraping, no bot-detection
  involved. (Contrast with the ThermoFisher catalog we investigated and declined
  to use: that required replaying stolen Akamai bot-mitigation tokens against
  an internal API, which their Terms of Use and `robots.txt` both prohibit.)
  The supplementary PDF (`mmc1.pdf`) was parsed with `pdftotext -layout`; the
  20 rows were verified by exact substring match against the corresponding
  NCBI RefSeq transcripts before being transcribed into
  `sirna_data.fetch.monopoli`.
- **Important caveat**: these are not standard unmodified siRNA duplexes.
  They use a cholesterol-conjugated, asymmetric (15-nt sense / 20-nt
  antisense) "sdRNA" architecture with heavy 2'-fluoro/2'-O-methyl/
  phosphorothioate modification (from Shmushkovich et al. 2018, Nucleic
  Acids Research, the dataset Monopoli's model was trained on). This is
  now flagged on every record (`is_modified=True`, `modification_chemistry`
  set -- see "Chemical modification data" below), but RNAfold still has no
  model of these modifications, so we fold them as if they were plain
  unmodified RNA — a real approximation whose accuracy on this chemistry is
  unverified. `technology="Dual-luciferase reporter assay (modified sdRNA)"`
  falls into `graph_build.py`'s existing "other" one-hot bucket, which at
  least lets the model separate this subset's systematic effects from the
  primary dataset's.
- Label conversion: the source reports "% reporter expression remaining"
  (lower = more potent); we store `label = 100 - reporter_remaining_pct` to
  match `%Inhibition`'s convention (higher = more knockdown).
- We also identified Table S1 of the *primary* source (Shmushkovich et al.
  2018, 356 siRNAs, same chemistry, also CC BY / open access, retrieved the
  same way) but did not integrate it: that table has no gene-identity column
  at all (only an internal compound ID), so it cannot be gene-grouped for
  LOGO-CV and would only be usable as always-in-training augmentation data,
  requiring pipeline support we haven't built. Left for a future iteration
  if more training signal (as opposed to more evaluable genes) becomes the
  priority.

## Supplementary siRNA data: Davis et al. 2025 (966 net-new rows after dedup, 0 new genes)

- Davis, Hildebrand, MacMillan, Monopoli, ... Pai & Khvorova 2025, *Nucleic Acids
  Research* 53(12):gkaf479, "Systematic analysis of siRNA and mRNA features
  impacting fully chemically modified siRNA efficacy"
  (doi:10.1093/nar/gkaf479, PMC12205987). **CC BY 4.0** — the same lab as
  Monopoli 2023 (Kathryn R Monopoli is a co-author here too), and the same
  four genes (`APP`, `MAPT`, `BACE1`, `SNCA`).
- Supplemental Table S1: 1,248 raw rows, each a fully chemically modified
  siRNA (2'-OMe/2'-F, one of 3 scaffold architectures: "Blunt", "Asymmetric
  2'-OMe/-F", "Asymmetric 2'-OMe Rich") measured by both a native QuantiGene
  2.0 assay (always present) and a dual-glo luciferase reporter assay
  (present for 536 rows, not currently used here). Filtered to the paper's
  own "Included in Filtered Dataset" == "Yes" subset (1,011 rows) before
  loading — the excluded 237 target sites weren't confidently expressed in
  the tested SH-SY5Y cells per the paper's own RNA-seq/3P-seq analysis.
  `data/raw/davis2025_extra.csv` ships this 1,011-row filtered/derived CSV
  directly (CC BY 4.0 permits redistributing a derivative, unlike
  CMsiRNAdb's CC BY-NC-ND).
- **Obtaining this file was not straightforward to automate**: Oxford
  Academic's supplementary-files CDN link for this article sits behind a
  Cloudflare Turnstile ("Verify you are human") challenge — a real,
  interactive bot-check, not a `robots.txt`/ToS restriction like the
  ThermoFisher case under Monopoli 2023 above, so no automated fetch was
  attempted (this tooling won't solve CAPTCHAs). The user downloaded
  `gkaf479_supplemental_files.zip` manually from their own browser
  session and supplied `Supplemental Tables.xlsx` from it; `davis2025_extra.csv`
  was derived from that file's "Supplemental Table 1" sheet. Because
  acquisition required a human in the loop, there is no
  `sirna_data/fetch/davis2025.py` — same pattern as `pdcd1_extra.csv`/
  `martinelli_extra.csv`/`oligograph_extra.csv`, which also have no
  corresponding fetcher (see `src/sirna_data/fetch/__init__.py`).
- **No transcript FASTA for this source**: unlike every other loader, this
  source ships its own local 50nt mRNA window per row directly in the raw
  table (the target's "consensus sequence across mRNA variants expressed in
  SH-SY5Y cells", computed by the paper's own authors from RNA-seq + 3P-seq)
  — so `_load_davis2025_records` uses that window as-is instead of calling
  `_locate_window` against a fetched full-length transcript. Verified: the
  stored 20mer target site is an exact substring of the stored 50mer window,
  at a constant 15nt offset, for all 1,248 raw rows. 22 of the 1,011 kept
  rows have `?` placeholder characters in the window's flanking region (the
  source's own consensus-calling marking an ambiguous position across mRNA
  variants — never inside the 20mer target itself); those rows fall back to
  duplex-only context (`has_flanking_context=False`) rather than shipping a
  window with `?` in it.
- Label conversion: same convention as Monopoli/Shabalina —
  `label = 100 - Native Assay Average (% Untreated Control)`.
- **Sequence-identity caveat (read before trusting per-position chemistry
  from this source)**: the raw table's "Antisense/Sense Strand Sequence and
  Chemical Modification Scaffold" columns embed real per-position 2'-OMe/
  2'-F chemistry annotation, in principle parseable the same way
  CMsiRNAdb's `Modification_Types_*_strand` columns are. But cross-checking
  every one of the 1,248 rows found these two columns' embedded base calls
  are **not** simple reverse complements of the verified 20mer target site,
  or of each other, under any of the 16 reversal/complement orientations
  tried (average ~10-15 mismatches out of 20nt in every case — close to
  what unrelated random sequences would give). Whatever indexing or
  strand-orientation convention actually produced those two columns could
  not be confidently reconstructed here. Rather than guess, `guide_seq` is
  derived the same way every other source lacking an independently
  trustworthy antisense column handles it: `_revcomp()` of the verified
  target site. `is_modified=True` and `modification_chemistry` carry only
  the coarse scaffold-name tag parsed from `Compound Name` (dataset-level,
  not per-position — `sense_modifications`/`antisense_modifications` stay
  unset). A future revision that resolves the raw notation's real
  orientation (e.g. by reading the paper's Methods/figure legends more
  closely) could upgrade this to real per-position data.
- **Dedup**: Monopoli 2023 (same 4 genes, same lab) and CMsiRNAdb's full
  retrieval (which also covers APP/MAPT) plausibly overlap this set by
  sequence. `load_records()` loads this source last, deduping against every
  other source's sequences (strand-agnostic) — of 1,011 filtered rows, 966
  survive as genuinely new duplexes (45 dropped as exact-sequence repeats of
  something already loaded). All 4 genes were already covered by Monopoli/
  CMsiRNAdb_full, so this source adds 0 new genes — its value is depth
  (966 more measurements) on 4 genes this dataset already had only sparse
  coverage of, not breadth.

## Supplementary siRNA data: PDCD1 panel (Xu, Zhao et al. 2024 / siRNABERT)

- Xu, Xu, Xie, Zhao, Yu & Feng 2024, *GENE*, "BERT-siRNA: siRNA target
  prediction based on BERT pre-trained interpretable model" (DOI
  10.1016/j.gene.2024.148330), and the associated code repo
  github.com/ChengkuiZhao/siRNABERT.
- 8 siRNAs against `PDCD1` (PD-1, the immune-checkpoint gene) with a
  dual-readout luciferase-reporter + qPCR knockdown efficiency for each —
  the only source in this dataset with real, gene-identified, numeric
  efficacy data for PDCD1. +1 gene (97 → **98 genes**) for
  leave-one-gene-out CV.
- **Recovery**: the repo's own `data/ExperimentalData/PDCD1_8.csv` was
  deleted from `main` (commit `d2ad931`, "Delete data/ExperimentalData
  directory") but is still present in the parent commit `f3254ee` — pulled
  from there via `raw.githubusercontent.com/ChengkuiZhao/siRNABERT/f3254ee/
  data/ExperimentalData/PDCD1_8.csv`, not from the current `main` tree. A
  sibling file in the same deleted directory, `Covid19.csv` (11 candidate
  SARS-CoV-2-targeting siRNA sequences), was also recovered but has **no
  efficacy/label column at all** in the committed version — despite the
  repo's own `TestExperimentalData.py` code expecting one (`values[:,4]`) —
  so it is not usable as trainable data and was not integrated.
- **Verification**: all 8 sense-strand sequences were confirmed by exact
  substring match against the real NCBI RefSeq `NM_005018.3` (human PDCD1)
  mRNA transcript before being trusted — not assumed from the file's own
  column labels.
- Label: the source gives both a luciferase-reporter efficiency and a qPCR
  efficiency per row (both already expressed as 0–1 fractional knockdown,
  e.g. `0.972747`); we store `label = 100 * Efficiency_QPCR` (mRNA-level,
  matching this dataset's dominant `%Inhibition` convention) and keep the
  luciferase value alongside as `Efficiency_LUC_Pct` for reference.
- **License caveat, unresolved — committed anyway at the user's explicit
  request**: unlike every other source in this file, the siRNABERT repo
  ships no `LICENSE` file (all-rights-reserved by default under GitHub's
  terms), and the associated *GENE* (Elsevier) paper is not confirmed
  open-access. This is the same situation as siRecords below: obtained
  from a real but non-canonical channel (deleted-file git history rather
  than a publisher-sanctioned download). The gap was raised and the user
  chose to commit `data/raw/pdcd1_extra.csv` to this repo regardless — see
  `NOTICE.md` for the commercial-use caveat that travels with it. The
  transcript FASTA (`pdcd1_transcripts.fasta`) is pure NCBI RefSeq (public
  domain) and is committed as normal.
- Fetched by hand (no `sirna_data.fetch.*` module yet, given the small
  size) into `data/raw/pdcd1_extra.csv` and `data/raw/pdcd1_transcripts.fasta`;
  loaded by `_load_pdcd1_records` in `src/sirna_data/raw_loader.py`.

## Supplementary siRNA data: Shabalina, Spiridonov & Ogurtsov 2006

- Shabalina, Spiridonov & Ogurtsov 2006, *BMC Bioinformatics* 7:65,
  "Computational models with thermodynamic and composition features improve
  siRNA design" (doi:10.1186/1471-2105-7-65). CC BY 2.0.
- Additional File 4 ("TableS1A"): a 653-siRNA / 52-gene heterogeneous
  training set the authors compiled from multiple published assays to train
  their "ThermoComposition" method. Retrieved via Europe PMC's public
  `supplementaryFiles` REST API
  (`https://www.ebi.ac.uk/europepmc/webservices/rest/PMC1431570/supplementaryFiles`),
  same legitimate route as the Monopoli data above — not scraping.
- We investigated this initially hoping to add the *Matveeva et al. 2007*
  (PMID 17426130) 3336-siRNA compilation instead, but its only listed
  download link (a University of Utah personal page) has been dead for
  years, has no Wayback Machine snapshot, and isn't mirrored as journal
  supplementary data; 3 of its 4 constituent source papers are fully
  paywalled with no legitimate open copy anywhere we could find, and the 4th
  (Jagla et al. 2005, Sloan-Kettering, obtained directly from the user) turned
  out to never have published its per-siRNA sequence data at all — only
  aggregate rule-level statistics. This Shabalina et al. 2006 dataset was
  found while reading Matveeva's reference list; unlike Matveeva's, it is
  open access with the actual data still retrievable.
- **Deduplication**: this table is itself a compilation, and roughly half of
  its 653 rows turned out to be exact antisense-sequence duplicates of genes
  already in siRNAEfficacyDB (traced to the same underlying Khvorova et al.
  2003 and other classic assays siRNAEfficacyDB itself draws from). We kept
  only the rows targeting genes absent from our existing 45-gene set: **269
  rows across 41 new genes**, mostly Hsieh et al. 2004's PI3K-pathway siRNA
  library (`PTEN`/`TSC1`/`TSC2`/`AKT1`/`AKT2`/`IGF1R`/`MAPK14`/`GSK3A`/
  `GSK3B`/`MYC`/`RAB13`/`EIF4EBP1`/`CBL`/`CBLB`/`CSK`/`ILK`/`PIK3R1`/
  `PIK3R2`/`PIK3CA`/`IRS1`/`FOXO1`/`FOXO4`/`SKP1`/etc.), plus a handful of
  one-off genes from other sources (`HRAS`, `PSKH1`, `RB1`, `CDKN1A`, mouse
  `MyoD`, human/mouse tissue factor `F3_human`/`F3_mouse`). This takes
  leave-one-gene-out CV from 45 to **86 genes**.
  - The accession→gene mapping and every exclusion/correction were derived
    by hand from NCBI esummary/efetch lookups on all 52 accessions in the
    source table (the table itself only has accessions, not gene symbols) —
    see `sirna_data.fetch.shabalina` for the full mapping and the
    reasoning per exclusion. Notably: `NM_000314` (PTEN) was excluded as a
    duplicate of the already-present `MMAC1` gene under its old name;
    `NM_004351` is `CBLB`, a distinct paralog from `NM_005188`/`CBL`, not a
    duplicate despite the superficially similar name; `U47298` is the pGL3
    luciferase reporter vector backbone, i.e. the same "Firefly luciferase"
    gene already present; `M25346` (a puromycin-resistance marker, kept as
    gene `PAC`) and the tissue-factor orthologs are kept as legitimate
    distinct non-endogenous targets, the same way "Firefly luciferase",
    "SEAP", and "EGFP" already are in the primary dataset.
- Label conversion: the source's `Activ` column is % activity/expression
  remaining (lower = more potent, same convention as Monopoli's data above);
  we store `label = 100 - Activ` to match `%Inhibition`.
- No per-row assay/technology detail is given in the source table (unlike
  siRNAEfficacyDB), so all 269 rows are tagged
  `technology="Heterogeneous compilation (Shabalina et al. 2006)"`, which
  falls into `graph_build.py`'s existing "other" one-hot bucket.
- Target-site location: exact substring search of the derived sense sequence
  (reverse complement of the source's antisense 19-mer) against the fetched
  RefSeq/GenBank transcript, same as the primary dataset. 266/269 (98.9%)
  located successfully; the remaining 3 fall back to duplex-only context via
  `has_flanking_context`, same graceful degradation as the primary dataset.
- Fetched by `sirna_data.fetch.shabalina` (`sirna-data-fetch`) into
  `data/raw/shabalina_extra.csv` and `data/raw/shabalina_transcripts.fasta`.

## Supplementary siRNA data: Martinelli 2023 / sirna-repro (577 of 907 rows, 12 new genes/reporters)

- Martinelli 2023, bioRxiv preprint / `sirna-repro` (siRNAmod-derived
  reproduction dataset). 907 chemically-modified siRNAs pooled from 30
  source PMIDs/patents, kept locally (gitignored, not redistributed) at
  `data/raw/sirna_repro_martinelli_907.csv` — see `NOTICE.md`.
- **The blocker and how it was resolved**: as originally obtained, this
  table has no gene-identity column at all — only sequence, a per-molecule
  modification descriptor, `PCT` (%inhibition, already in this dataset's
  standard 0–100 convention), and a source PMID/patent ID. That alone made
  it unusable for leave-one-gene-out CV (see `POTENTIAL_DATA_SOURCES.md`'s
  original "not integrated" verdict). Rather than trust each source
  document's stated target from its abstract/title, every candidate
  gene/reporter assignment was **confirmed computationally**: the sense
  strand's first 19nt ("core", excluding each row's own 2nt 3' synthetic
  overhang) was checked for exact substring containment in the real target
  transcript, fetched fresh from NCBI/GenBank for that purpose.
- **Pass 1 — source-document tracing: 253 of 907 rows (27.9%) resolved this
  way, across 7 distinct real targets** — +7 genes/reporters (98 → **105
  genes**) for leave-one-gene-out CV:

  | Gene/reporter | Rows | Reference accession | Source PMIDs |
  |---|---|---|---|
  | `EGFP` | 74 | U55763.1 (pEGFP-C1 vector) | 19282453, 18575806, 12923253, 17363479, 22287630, 12408823 |
  | `Luciferase_firefly` | 58 | U47296.1 (pGL3-Control, *luc+*) | 15653644, 22411910, 25699137, 22260772, US 20080249039 A1, US 8653252 B2 |
  | `Luciferase_renilla` | 43 | AF025846.1 (pRL-TK) | 15653644, 17150641, 17924376, US 8653252 B2 |
  | `APOB` | 34 | NM_000384.3 | 21047800, 19917641 |
  | `ACP5` (TRACP) | 32 | NM_001111035.2 | 17511001 |
  | `NPY` | 8 | NM_012614.2 | 15653644, US 8653252 B2 |
  | `VEGFA` | 4 | NM_001025366.3 | 21985606 |

- **A second resolution pass found 32 more rows whose "sense strand" column
  was itself corrupted but whose antisense (guide) strand was intact and
  correct.** For 28 of these rows (all of PMID 19917641's 23 rows, plus 3
  from 17924376 and 2 from 25699137), the row's stated "sense" sequence did
  not match any known reference (the reason these were originally left
  unresolved) but was also not the reverse complement of the row's own
  antisense field — i.e., internally inconsistent, not just externally
  unverifiable. Computing the reverse complement of the antisense field
  and checking *that* against the reference transcripts found a clean
  19nt match every time: for the 3 rows from 17924376 and 2 from 25699137,
  the resulting corrected sequence differs from the stored "sense" field by
  only 1-2 characters (consistent with a data-entry/OCR typo, e.g. stored
  `CUUACGCUCUGUACUUCGA` vs. correct `CUUACGCUGAGUACUUCGA`); for all 23 rows
  of PMID 19917641 the stored "sense" field is the exact same,
  completely-unrelated 19-mer repeated across every row regardless of PMID
  19917641's real (varying) PCT values -- consistent with a single
  copy/fill-down error in the original spreadsheet compilation rather than
  23 independent scrambled controls (the original, more cautious read of
  this PMID's rows -- documented in an earlier revision of this file --
  is superseded by this finding: the antisense strand's guide-to-target
  relationship is confirmed correct against real human APOB,
  NM_000384.3). For these 28 rows, the CSV's `Sequence` column stores the
  *corrected* 19nt sense (derived from the verified-correct antisense), not
  the source table's original (unreliable) sense value; `Sequence_antisense`
  is unchanged, still taken directly from the source.
- A separate check confirmed **4 more rows (PMID 21985606) target VEGFA**
  (NM_001025366.3) directly and cleanly on both strands -- no correction
  needed, just a reference sequence (VEGFA) not previously in this
  project's Martinelli reference set.
- **Pass 2 — brute-force sequence matching (a follow-up pass over the
  ~654 rows pass 1 left unresolved): 324 more rows resolved, +5 more
  genes** — +12 genes/reporters total (98 → **110 genes**) for
  leave-one-gene-out CV. Rather than trying to trace each source
  document's stated target first, this pass fetched candidate reference
  transcripts and checked every still-unresolved sense strand against them
  by exact 19nt substring match directly — the same verification standard
  as pass 1, just reached from the reference sequence inward instead of
  from the citation outward:

  | Gene/reporter | Rows | Reference accession | Source |
  |---|---|---|---|
  | `Luciferase_firefly` (GL2 site) | 64 | X65324 (pGL2-Control) | PMID 17924376 (37), 15919084 (20), 17150641 (7) — all three papers used the classic Elbashir/Tuschl "GL2" anti-firefly-luciferase siRNA (`CGUACGCGGAAUACUUCGA`), a different target site from the pGL3-Control (`luc+`) site already covered by pass 1's `Luciferase_firefly` rows. Kept under the same `Luciferase_firefly` gene label since it's the same protein target — `Accession_number` distinguishes the two sites (U47296 vs X65324) per row. |
  | `KAZRIN` | 39 | NM_201628 | US20120088815A1/EP2415869A1 |
  | `MIR155HG` (BIC) | 42 | NR_001458 | US20120088815A1/EP2415869A1 |
  | `CDKN1B` | 60 | NM_004064 | US20120088815A1/EP2415869A1 |
  | `SOD2_chimp` | 51 | NM_001009022 (*Pan troglodytes*) | US20120088815A1/EP2415869A1 |
  | `CASR_rhesus` | 68 | XM_001111972 (*Macaca mulatta*, PREDICTED) | US20120088815A1/EP2415869A1 |

  The five patent-derived genes came from that patent's own "Example"
  sections, which give a fully-worked case per gene (e.g. "Human KAZRIN
  (Genbank Accession No. NM_201628) ... nucleotides 718-738 ... SEQ ID
  NO:2") — but rather than trust the patent's SEQ ID NO pairing (a small
  summarization pass over a document this large is not reliable enough to
  copy numbers from directly), each of the five genes' real RefSeq/GenBank
  transcript was fetched independently and matched against all 544
  then-unresolved patent rows by exact substring search, letting the
  sequence data itself confirm or refute every assignment. `SOD2_chimp`
  and `CASR_rhesus` are kept as species-tagged gene labels distinct from
  any future human `SOD2`/`CASR` addition, the same convention already
  used for Shabalina's `F3_human`/`F3_mouse` pair. 309 of these 324 rows
  (95.4%) locate with full 30nt flanking mRNA context; the rest fall back
  to duplex-only context. All 324 matches are forward-orientation (sense
  strand identical to the mRNA window, not its reverse complement) — no
  further sense/antisense correction needed, unlike pass 1's 28-row fix.
- **The remaining ~330 rows are NOT included** and remain
  unresolved:
  - The single largest chunk (**284 rows, ~86% of what's left**) is the
    remainder of patent US20120088815A1/EP2415869A1 (pass 2 above resolved
    260 of the patent's 544 rows to 5 genes found in the patent's own
    worked examples). The patent's own text says it tested "target genes
    in table 1-182" via a luciferase-reporter assay; the fetched portion of
    the document only reaches tables 1-51, so the genes covered by tables
    52-182 remain unidentified — reaching them would need either full
    access to the rest of this very large document or BLAST access,
    neither available in this environment.
  - PMID 16598842 (15 rows): not yet investigated (only found as a citation
    in other papers' reference lists, not fetched directly — the paper's
    full text sits behind an Elsevier paywall not reachable from this
    environment).
  - PMID 23820891 (6 rows), "5' Unlocked Nucleic Acid Modification Improves
    siRNA Targeting" (Snead et al. 2013, PMC3732871): confirmed via full
    text to target the HIV-1 transcript (an siRNA called "siH5", assayed
    via strand-specific dual-luciferase reporter, not a genomic human
    target) -- out of scope for this project's gene-symbol convention
    without adding an HIV reference sequence, not pursued.
  - PMID 22889374 (4 rows): confirmed via abstract to target Enterovirus 71
    (EV71)'s 5'UTR -- a viral, not human, target; not pursued (would need a
    viral genome reference).
  - PMID 21141919 (4 rows) and PMID 22982308 (2 rows): both titled/abstracted
    as anti-MDR1 (ABCB1) work and share a near-identical sequence across
    both papers, but neither the sense field nor the antisense-reverse-
    complement of any of these 6 rows matches human ABCB1 (NM_000927.5),
    nor mouse Abcb1a (NM_011076) or Abcb1b (NM_011075) -- left unresolved;
    the papers' own full text (not reachable from this environment) would
    be needed to confirm the actual species/isoform/paralog.
  - PMID 15653644 (2 rows), 25699137 (2 rows), 22260772 (2 rows): small
    remainders in three papers otherwise resolved above (EGFP/Luciferase/
    NPY/APOB). 15653644's own text confirms it also tested an siRNA
    against the SARS-CoV genome (a 5th target the paper describes but
    doesn't give in its own tables) -- almost certainly these 2 rows,
    viral and out of scope, not pursued further. 25699137's and 22260772's
    remaining sequences don't match any reference already fetched for this
    project and their full text wasn't reachable this round.
  - PMID 20005874 (3 rows): confirmed via abstract to target coxsackievirus
    B3 -- viral, not pursued.
  - A handful of smaller/single-row documents (~10 rows total) not yet
    investigated: PMID 22895883, 17616127, 17539595, 23682837.
- **Both strands are taken directly from the source** (except the 28
  sense-corrected rows above), not derived by reverse-complementing one
  from the other — unlike every other loader in this file, Martinelli's
  raw table gives independently-recorded sense and antisense sequences
  (each with its own overhang and modification descriptor), and using the
  source's real antisense preserves that detail instead of reconstructing
  an idealized one.
- **Chemical modification**: every one of the 577 rows carries a real,
  non-placeholder per-molecule modification descriptor — locked nucleic
  acid, hexitol nucleic acid, 2'-fluoro, 2'-O-methyl, 2'-deoxy, unlocked
  nucleic acid, and 4-thioribose all appear. `is_modified`/
  `modification_chemistry` are populated from this; like Monopoli2023 (and
  unlike CMsiRNAdb), the annotation is per-molecule, not per-position, so
  `sense_modifications`/`antisense_modifications` stay `None`.
- Loaded by `_load_martinelli_records` in `src/sirna_data/raw_loader.py`
  from `data/raw/martinelli_extra.csv` (the 577-row derivative with added
  gene/accession columns, and corrected sense sequences for the 28 pass-1
  rows noted above) and `data/raw/martinelli_transcripts.fasta` (the 13
  reference sequences above, 7 from pass 1 plus 6 from pass 2). **License**:
  the source `sirna-repro` dataset is CC BY-NC 4.0, which — unlike
  CMsiRNAdb's CC BY-NC-ND — permits derivatives, so this filtered/
  gene-annotated subset is committed to the repo like Monopoli/Shabalina,
  not gitignored like PDCD1/siRecords.

## Supplementary siRNA data: OligoGraph repo (343 rows, 3 new genes; label interpreted as %KD)

- github.com/drugparadigm/OligoGraph — a public training-data repo for an
  siRNA-efficacy GNN, cloned directly (its web UI is unreliable/blocked for
  this purpose — repeatedly misreported the repo's actual file structure).
  Ships 4 CSVs: `Hu.csv` (2,472 rows), `Mix.csv` (472 rows), `Taka.csv`
  (~3,412 rows), `Simone.csv` (322 rows), each with `siRNA`/`mRNA`/`label`
  columns (`mRNA` = sense/target-context strand, `siRNA` = its exact reverse
  complement — verified programmatically for 20 sample rows across files,
  matches this project's own loader convention of storing the sense strand
  and deriving `guide_seq` by reverse-complementing it).
- **Overlap check**: every OligoGraph sequence was checked against this
  project's full existing corpus by exact-core-sequence match. `Hu.csv`'s
  2,361 overlapping rows are all already present here (from siRNAEfficacyDB
  and others) and add nothing new. `Taka.csv` has no discoverable literature
  citation at all (checked several candidate papers, none matched) and was
  not pursued. `Mix.csv` and `Simone.csv` are where the genuinely new
  content is.
- **Gene identity — independently verified, not trusted from any file/table
  label**: 343 rows (300 from `Simone.csv`, 43 from `Mix.csv`) were matched
  by exact 19nt substring search of the sense strand against real NCBI
  RefSeq transcripts, fetched fresh for this purpose. All 343 located with
  full 30nt flanking mRNA context.

  | Gene | Rows | Reference accession | Traced to |
  |---|---|---|---|
  | `HIF1A` | 100 | NM_001530.4 | Sciabola et al. 2013 |
  | `HK2` | 95 | NM_000189.5 | Sciabola et al. 2013 |
  | `HPSE` | 105 | NM_006665.5 | Sciabola et al. 2013 |
  | `Lamin A` (merged into the existing gene group, not a new `LMNA` label) | 43 | NM_170707.4 | Harborth et al. 2001 (via Sciabola 2013's compilation) |

  `Simone.csv`'s 300 rows trace to Sciabola et al. 2013's "HUVK" siRNA
  training compilation (via its title/scope match plus a literature search;
  the paper's own Supplementary Table S4 could not be fetched to confirm
  row-for-row — see the metric caveat below). `Mix.csv`'s 43 lamin-A/C rows
  are labeled by Sciabola 2013 as citing Harborth et al. 2001 — but
  Harborth's own paper (read in full via an open academia.edu mirror) does
  not actually contain a 44-siRNA panel with individual values; it
  describes a single validated lamin A/C siRNA used as one tool within a
  21-gene essential-gene screen. So these 43 rows' real origin is most
  likely Huesken et al. 2005's large-scale tiling-array supplementary data
  (the systematic source Sciabola's compilation most plausibly drew the
  numeric values from), not Harborth's paper itself — not yet independently
  confirmed.
- **THE METRIC — read this before using this subset for anything where the
  label's exact scale matters.** OligoGraph does not ship the original
  papers' reported %-knockdown values, only its own `label` column (~0–1
  range), which is undocumented anywhere in the OligoGraph repo (its own
  `Preprocess2.py` preprocessing script and README were both checked —
  neither defines it). An initial attempt fit `label * 134.1 = true
  %Inhibition` from `Hu.csv`'s 2,361 overlapping rows (R²=1.0, seemingly
  exact) — but the same formula, tested against `Mix.csv`'s own 420
  overlapping rows, gave R²=0.82 with concrete counterexamples (the same
  `label=0.25` value mapping to four different true percentages: 62, 62, 0,
  -27). `Hu.csv`, `Mix.csv`, and `Taka.csv` each independently max out at
  exactly `label=1.0` — consistent with each source file being
  max-normalized to its own scale by whatever pipeline built these CSVs,
  not one shared constant across the repo. `Simone.csv` has zero rows
  overlapping this project's existing corpus, so the formula was never even
  checkable against it. Pulling the true values directly from the primary
  sources was attempted and blocked: Sciabola et al. 2013's Supplementary
  Table S4 returns HTTP 403 from Oxford's CDN; PubMed, PMC, and the journal
  page itself are reCAPTCHA/paywall-blocked; Europe PMC's API is
  rate-limited. **Given that, `label * 100` is used directly as %KD here per
  explicit instruction — a stated interpretation of an undocumented
  normalized value, not an independently verified true reported knockdown
  percentage.** Anyone relying on this subset for anything where the exact
  scale matters should treat it accordingly, and revisit if Huesken et al.
  2005's own supplementary data (the more plausible true source for the
  `Mix.csv` rows) or Sciabola 2013's Table S4 become reachable.
- Loaded by `_load_oligograph_records` in `src/sirna_data/raw_loader.py`
  from `data/raw/oligograph_extra.csv` (343-row derivative:
  `Sequence,Gene,Accession_number,Pct_Inhibition,Source_Paper`) and
  `data/raw/oligograph_transcripts.fasta` (4 reference sequences). **License**:
  OligoGraph's own repo carries no explicit license file — treated as
  all-rights-reserved by default, same convention as the PDCD1 panel; see
  `NOTICE.md`.

## Supplementary siRNA data: CMsiRNAdb, human PCSK9 subset (2,756 rows, 1 new gene)

- He et al. 2026, *BMC Bioinformatics* 27:33, "CMsiRNAdb: a database of
  chemically modified siRNA silencing efficiency for nucleic acid drug
  design" (DOI 10.1186/s12859-025-06359-y). **CC BY-NC-ND 4.0** -- unlike
  every other source in this file, the "ND" (No Derivatives) term means we
  can redistribute the *original, unmodified* download but not a
  filtered/curated/collapsed adaptation of it. Built by the same research
  group as siRNAEfficacyDB (our primary source, same
  `cellknowledge.com.cn` platform) as an explicit companion covering
  *chemically modified* siRNAs, which siRNAEfficacyDB doesn't. Live
  database at `cellknowledge.com.cn/CMsiRNAdb/`, with a no-login bulk TSV
  download -- not a paywalled supplement like the two papers investigated
  and declined above.
- **ND-compliant design**: this repo ships only the untouched original
  download, `data/raw/cmsirnadb_full_raw.tsv` (43,153 rows across 13
  genes: AGT, ANGPTL3, APP, CTNNB1, HSD17B13, INHBE, LPA, MAPT, MARC1,
  MSTN, PCSK9, PLN, PNPLA3). All filtering, species exclusion, and
  collapsing described below happens at *load time* in
  `src/sirna_data/raw_loader.py` (`_load_cmsirnadb_records`), not as a
  pre-computed file -- every caller reproduces their own local copy of the
  derived data instead of downloading an adaptation from us. (The
  transcript FASTA these loaders also read is independently fetched from
  NCBI RefSeq -- public domain, not CMsiRNAdb material -- so that ships
  as-is.) Fetched by `sirna_data.fetch.cmsirnadb` (`sirna-data-fetch`) into
  `data/raw/cmsirnadb_full_raw.tsv`, `data/raw/cmsirnadb_transcripts.fasta`,
  and `data/raw/cmsirnadb_full_transcripts.fasta`.
- Only the human PCSK9 subset is used by `_load_cmsirnadb_records` --
  **2,756 of 3,107 raw PCSK9 rows kept** -- takes leave-one-gene-out CV
  from 86 to **87 genes**. The other 12 genes are handled separately (see
  "CMsiRNAdb — full retrieval, the other 12 genes" below).
- **Species filtering**: the raw PCSK9 rows mix human and non-human data
  under the same gene label. Excluded outright: rows on accession
  `NM_153565.2` (*Mus musculus* Pcsk9) and rows with `Cell_Type` of
  `Mus musculus` or `Non-human hepatocytes` -- 321 rows total (3,107 ->
  2,786).
- **Accession quirk, resolved**: many raw rows cite `NR_110451.3`, a human
  but *non-coding* RefSeq transcript variant of PCSK9 (not the
  protein-coding mRNA). Rather than fetching and folding a second,
  non-coding reference sequence, every surviving row (regardless of its
  own stated accession) is located against the ONE canonical human coding
  transcript, **NM_174936.4**; rows whose site isn't found there fall back
  to duplex-only context via the existing `has_flanking_context`
  mechanism, exactly like every other source's unmapped rows.
- **Critical exclusion -- protecting external validation against LEQVIO**:
  PCSK9 is also the target of LEQVIO/inclisiran, one of the 5 FDA-approved
  drugs in the external-validation set below. Any row whose antisense
  sequence contains inclisiran's real 19nt target core
  (`AAGCAAAACAGGUCUAGAA`) is dropped -- 28 rows excluded (2,786 -> 2,758)
  -- so that drug stays genuinely unseen data for anyone using it for
  external validation. None of the other 4 external-validation drugs'
  target genes (TTR, ALAS1, HAO1) are among CMsiRNAdb's 13 genes, so no
  further exclusion is needed there.
- **Raw data-entry contamination**: a small number of rows (2 for PCSK9)
  have modification-notation characters (parentheses, ambiguity codes)
  bleeding into the sequence column instead of clean bases -- dropped as
  unusable (2,758 -> **2,756** final).
- 2,707/2,756 kept rows (98.2%) locate with full 30nt flanking context
  against NM_174936.4; the rest fall back to duplex-only context.
- **Known caveats, not fixed**: (1) All sequences here carry real chemical
  modifications (2'-O-methyl, phosphorothioate, 2'-fluoro, etc. -- 36
  types across the full database). Their identity and position ARE now
  captured -- see "Chemical modification data" below -- but folding/
  structure prediction (RNAfold) still has no model of modified bases and
  treats every sequence as plain unmodified RNA, same approximation
  already accepted for Monopoli2023, arguably on a wider and more
  heterogeneous set of chemistries here since this spans 11 different
  patents/filers rather than one lab's single design. (2) Patent-derived, not assay-paper-derived:
  quality/protocol consistency across 11 different patent filers is
  inherently more heterogeneous than the academic sources. (3) Real
  sequence duplication: some sequences are repeated across many
  concentration/cell-type combinations or independently claimed in
  multiple patents. Not deduplicated (unlike the other-12-genes addition
  below, which does collapse repeated measurements) -- a real imbalance
  worth knowing about. `technology` is tagged
  `"CMsiRNAdb patent-derived, chemically modified (<cell type>)"` per row,
  which falls into downstream feature-engineering's "other" bucket like
  the rest of this file's non-primary sources.

## CMsiRNAdb — full retrieval, the other 12 genes (9,601 rows, 10 new genes)

- Task: get more data, exclusively siRNA (synthetic siRNA duplexes with numeric
  knockdown) — no shRNA, no dsRNA, no miRNA, no ASO. Source:
  `https://www.cellknowledge.com.cn/CMsiRNAdb/` (Zhang lab, Chengdu Univ. of TCM /
  UESTC), file `download/CMsiRNA_data_update.tsv`. Contacts: yangzhang@cdutcm.edu.cn,
  zhy1001@alu.uestc.edu.cn. The 13 per-gene `patent_dataset_*.tsv` files were also
  downloaded and confirmed to be an exact subset of the master (0 additional rows), so
  only the master is kept.
- **License: CC BY-NC-ND 4.0 — "No Derivatives"**, same as the PCSK9 subset above.
  This repo only ships the untouched original master TSV,
  `data/raw/cmsirnadb_full_raw.tsv` (43,153 rows, all with numeric % inhibition, 13
  target genes: AGT, ANGPTL3, APP, CTNNB1, HSD17B13, INHBE, LPA, MAPT, MARC1, MSTN,
  PCSK9, PLN, PNPLA3). Everything below — outlier removal, transcript location,
  deduplication, collapsing repeat measurements — happens at *load time* in
  `_load_cmsirnadb_full_records()` (`src/sirna_data/raw_loader.py`), not as a
  precomputed file, so no filtered/adapted CMsiRNAdb derivative is redistributed. This
  section covers the other **12 genes, 40,046 rows**, separate from the PCSK9 slice
  (3,107 rows, via `_load_cmsirnadb_records()`) covered above.
- **Derivation pipeline** (all in `_load_cmsirnadb_full_records()`):
  1. **Start**: 40,046 raw non-PCSK9 rows.
  2. **Outlier removal**: drop rows with `%inhibition` outside `[-50, 100]` (the
     source has a few corrupt values, e.g. `-7,103,597`) → **39,686 rows** (−360).
  3. **Contamination filter**: drop rows whose sense sequence contains non-ACGU
     characters after RNA conversion (modification notation, ambiguity codes, or
     stray characters bleeding into the sequence column) → **38,345 rows** (−1,341).
  4. **Cross-source dedup**: any row whose guide or sense sequence (strand-agnostic)
     already exists in the baseline dataset (siRNAEfficacyDB + Monopoli2023 +
     Shabalina2006 + CMsiRNAdb PCSK9) is skipped. Verified: **0 rows excluded** here —
     the 12-gene slice has no sequence overlap with the rest of the data.
  5. **Transcript location**: each row's sense sequence is searched (sliding 19nt
     window) against its gene's NCBI RefSeq transcript
     (`data/raw/cmsirnadb_full_transcripts.fasta`, 34 records fetched via NCBI
     efetch); rows with no window match fall back to the first 19nt of the sense
     sequence.
  6. **Collapse to unique duplexes**: raw patent measurements would otherwise become
     many near-identical training examples, so rows are collapsed to one record per
     unique `(gene, accession, target site)`; label = **median % inhibition** across
     that duplex's replicate measurements → **9,601 unique duplex records**.
- **6,117 of 9,601 records (63.7%)** locate with full 30nt flanking mRNA context; the
  rest fall back to duplex-only context, same mechanism as every other source.
- Records by gene: PNPLA3 2,066 / HSD17B13 1,985 / APP 952 / AGT 872 / MARC1 823 /
  INHBE 670 / MAPT 630 / LPA 556 / ANGPTL3 551 / CTNNB1 352 / PLN 135 / MSTN 9.
- Also re-checked `siRNAEfficacyDB/download/siRNA_all.txt` (the canonical source of
  the classic sets: Huesken 2005 x2,431, Katoh 2007 x702, Reynolds x244, Vickers x76,
  Harborth x44, Ui-Tei x37, Khvorova x10 = 3,544 numeric rows) directly from
  CMsiRNAdb's sister site — dedup showed only 12 new rows against `sirna_efficacy.csv`,
  confirming nothing further to add from that source.
- **Trainability note**: these are synthetic siRNA duplexes with directly-measured %
  inhibition — the same target quantity as this dataset's existing numeric-%KD data.
  This is a drop-in extension of the supervised set. Caveat: heavily chemically
  modified and patent-derived (dose/assay conditions vary), and highly redundant at
  the sequence level (40,046 raw rows -> 9,601 unique duplexes) — group-aware CV by
  sequence is advisable to avoid leakage. `technology` is tagged the same way as the
  PCSK9 subset: `"CMsiRNAdb patent-derived, chemically modified (<cell type>)"`.
- **Files**: `data/raw/cmsirnadb_full_raw.tsv` (the full 43,153-row master; PCSK9 and
  these 12 genes are both derived from it at load time), `data/raw/cmsirnadb_full_transcripts.fasta`
  (34 mRNA transcripts, NCBI efetch), and `data/cmsirnadb_new_sirna.png` (%inhibition +
  per-gene figure).
- **Integration**: wired into `load_records()` via `include_cmsirnadb_full=True`
  (`src/sirna_data/raw_loader.py`, `_load_cmsirnadb_full_records()`), which builds a
  strand-agnostic index of every sequence already loaded from the other sources before
  running (step 4 above). Effect on the dataset (verified by running `load_records()`):
  without, **6,577 records / 87 genes**; with, **16,178 records / 97 genes** (+9,601
  records, +10 new genes: AGT ANGPTL3 CTNNB1 HSD17B13 INHBE LPA MARC1 MSTN PLN PNPLA3;
  APP & MAPT deepened). Anyone consuming this data with a cached/precomputed downstream
  representation (e.g. a graph cache) should rebuild it after pulling this update.

## Considered, not integrated: siRecords (Ren/Gong 2006–09)

- siRecords 04/28/05 release, retrieved from the Internet Archive (the live siRecords
  servers at `sirecords.umn.edu`, `c1.accurascience.com`, and `sirecords.biolead.org`
  are all defunct — HTTP 502).
- **License status — do not treat as "open"**: siRecords' own definitive writeup, Ren
  et al. 2009, *Nucleic Acids Research* 37 (Database issue) D146-D149, "siRecords: a
  database of mammalian RNAi experiments and efficacies" (doi:10.1093/nar/gkn817), is
  itself published under CC BY-NC 2.0 UK — but that license covers the *article*
  (text/figures), not a grant to redistribute the underlying bulk dataset. The paper's
  own DATA ACCESS section states the actual terms for the data itself: "The siRecords
  web site is publicly accessible through the URL http://siRecords.umn.edu/siRecords.
  Academic users can obtain a copy of the current release of the dataset by sending an
  email" to the corresponding author — a controlled, individual-request distribution
  model restricted to academic users, not a blanket open-data license. Separately, the
  database's later mirror host, AccuraScience (`c1.accurascience.com/siRecords/`),
  publishes a general Terms of Use for its site stating downloaded content is for
  "personal non-commercial use" only and may not be redistributed "for any other
  purpose whatsoever without the prior written permission" of AccuraScience. Neither of
  those covers how this repo obtained the data: `sirecords_efficacy.csv` was recovered
  from an Internet Archive snapshot of the live site, not through the authors'
  sanctioned academic-request channel — so even the narrow "academic users, on request"
  permission the original paper describes doesn't technically apply to this copy. Net
  finding: siRecords' data was never established as freely redistributable, and there
  is no license (CC or otherwise) that clearly covers bulk redistribution of it as done
  here. Treat `sirecords_efficacy.csv` and `sirecords_new_only.csv` as **unresolved
  license risk**, not merely "unverified" — see `NOTICE.md`. Both files are kept
  locally but excluded from git (`.gitignore`), same as the PDCD1 extra data above.
- **Overlap methodology**: matching replicates this dataset's own dedup rule — exact
  nucleotide-sequence identity, strand-agnostic. Every existing guide/target sequence
  (from `sirna_efficacy.csv` antisense+sense columns, plus the Shabalina/Monopoli/
  CMsiRNAdb extras) was normalised (U→T, uppercased) and indexed both as-is and as
  reverse-complement, down to overlapping 19-mers. A siRecords sequence counts as
  "already in this dataset" if it (or its reverse complement) shares a 19-mer with any
  existing sequence; sequences shorter than 19 nt were tested as exact substrings.
- **Headline** — of the 4,162 siRecords rows, 275 have no usable sequence (too short /
  blank) and are excluded from the record-level comparison:

  | Level | New | Already in this dataset | siRecords total (usable) |
  |---|---|---|---|
  | **siRNA records** (by exact sequence) | **3,295 (84.8%)** | 592 (15.2%) | 3,887 |
  | **Unique sequences** | 2,901 (84.1%) | 550 (15.9%) | 3,451 |
  | **Target-gene accessions** | 1,534 | 49 | 1,583 |

  The 592 overlapping records trace to the existing sources exactly as expected —
  siRNAEfficacyDB (primary) 446, Shabalina 2006 145, CMsiRNAdb 1 — because
  siRNAEfficacyDB and Shabalina 2006 are themselves compilations of the same classic
  mid-2000s assays (Huesken 2005, Reynolds, Khvorova, Vickers, Hsieh) that siRecords
  aggregates, so the shared core overlaps.
- **What is genuinely new — and the critical caveat**: 3,295 records (3,117 carrying
  an efficacy rating), spanning ~1,400 new target-gene accessions, are not currently
  in this dataset by sequence — on its face a large potential expansion of gene
  coverage. However, this "new" data is **not directly train-ready**, for one decisive
  reason: siRecords efficacy is a **4-level ORDINAL rating** (Very high / High /
  Medium / Low), not a numeric `%Inhibition`. Every integrated source in this dataset
  stores a continuous knockdown/inhibition percentage (a regression/
  classification-ready target); siRecords does not provide per-record percentages,
  only the coarse bin. The new-only rated breakdown is: Very high 1,145 · High 1,003 ·
  Medium 485 · Low 484. To use the new siRecords records you would have to either (a)
  train/evaluate on the ordinal label directly (a different target from the numeric
  percentage this dataset otherwise provides), or (b) go back to the **PubMed IDs**
  siRecords provides (present for essentially every row) and extract the reported
  numeric knockdown from the primary papers — the same provenance-tracing approach
  used elsewhere in this file (e.g. Martinelli 2023 above).
- Secondary caveats: siRecords predates and overlaps the existing academic sources, so
  the 15% that overlaps is redundant and should be dropped on integration; cell line /
  assay / concentration metadata is present but formatted differently from this
  dataset's schema and would need mapping to the `Technology` one-hot buckets used by
  downstream feature engineering; sequence lengths are heterogeneous (12–64 nt; median
  19) where the rest of this dataset assumes ~19–21-mers.
- **Files**: `data/raw/sirecords_efficacy.csv` (full 4,162-record siRecords release,
  all fields) and `data/raw/sirecords_new_only.csv` (the 3,295 records whose sequence
  is not already in this dataset — 3,117 with an efficacy rating — the candidate
  extension set, overlap already removed).

## External validation: five FDA-approved siRNA drug sequences (not training data)

- Real-world antisense-strand sequences (with 2'-F/2'-O-Me/phosphorothioate
  chemical modification notation) and reported clinical inhibition for 5
  FDA-approved siRNA drugs — AMVUTTRA (vutrisiran/TTR), GIVLAARI
  (givosiran/ALAS1), LEQVIO (inclisiran/PCSK9), OXLUMO (lumasiran/HAO1), and
  ONPATTRO (patisiran/TTR) — came from AttSiOff (Liu, Yuan, Pan, Shen & Jin
  2024, *Med-X* 2:5, "AttSiOff: a self-attention-based approach on siRNA
  design with inhibition and off-target effect prediction", DOI:
  10.1007/s44258-024-00019-1), Supplementary Table S1, fetched via Europe
  PMC/Springer's public supplementary-materials link
  (static-content.springer.com/esm/art%3A10.1007%2Fs44258-024-00019-1/...).
- Deliberately not part of any train/val/test split — included here purely
  as a held-out, real-world validation set for anyone benchmarking a model
  trained on this dataset against approved drugs it never saw in training.
- Chemical-modification notation was parsed down to plain bases, then each
  drug's real antisense/sense assignment and the exact 3' overhang boundary
  were **verified computationally** (not assumed from table column order,
  which turned out unreliable — ONPATTRO's antisense/sense were transposed
  in a naive read of the source table) by exact substring search of every
  strand/orientation/trim combination against the real NCBI RefSeq
  transcript, keeping only the confirmed longest match. TTR: NM_000371.4.
  ALAS1: NM_000688.6 (not the first hit a naive gene-symbol lookup returns,
  an unrelated PREDICTED XM_ isoform that doesn't contain the target site —
  fetched the canonical NM_ record directly instead). PCSK9: NM_001407247.1.
  HAO1: NM_017545.3.
- AttSiOff's own predicted inhibition + ranking-percentile for these same
  five drugs (their Supplementary Table S3) is available for comparison if
  you want to benchmark against it — note AttSiOff itself has no runnable
  code/weights published anywhere, and its DH/DR/DT training data is not
  publicly downloadable (only "available from the corresponding author upon
  reasonable request" per the paper), so only its self-reported numbers can
  be used, not a re-run.

## Chemical modification data

Most of this dataset is standard/unmodified synthetic siRNA, but a growing
minority is chemically modified (2'-O-methyl, 2'-fluoro, 2'-deoxy,
phosphorothioate backbone, lipid/GalNAc conjugates, etc.) -- the kind of
stabilizing chemistry used in real therapeutic siRNAs. Until now that
distinction existed only in prose (this file's per-source notes, e.g.
"chemically modified" in CMsiRNAdb's `technology` string); `SiRNARecord`
itself had no field for it, so every downstream feature-engineering step
was silently treating every record as plain unmodified RNA regardless of
what was actually assayed.

`SiRNARecord` now carries four modification fields:

- `is_modified: bool` -- whether this specific measured molecule carries
  any known chemical modification. `False` (the default) for every source
  unless a loader explicitly sets otherwise.
- `modification_chemistry: str | None` -- a short human-readable summary of
  the chemistry class, e.g. `"2'-OMe/2'-F/PS-backbone (per-position,
  CMsiRNAdb)"` or a dataset-level note when no per-position detail exists.
- `sense_modifications` / `antisense_modifications: tuple[str | None, ...] | None`
  -- per-position modified-nucleoside name (e.g. `"2'-O-Methylcytidine"`)
  or `None` (confirmed unmodified/natural ribonucleotide) at each index,
  aligned 1:1 with the corresponding strand's stored sequence. The whole
  field is `None` (not a tuple of `None`s) when no per-position annotation
  is available at all for that record -- distinct from "we checked, it's
  unmodified everywhere here".

**Coverage, by source:**

- **CMsiRNAdb** (PCSK9 subset + the other-12-genes addition, ~12,357
  records) -- the raw `cmsirnadb_full_raw.tsv` already carries real
  per-position modification columns (`Modification_Types_{Sense,
  Antisense}_strand`, one `position*chemistry-name` entry per nucleotide of
  the full raw strand) that were previously read for sequence/label only
  and otherwise discarded. `_cmsirnadb_align_modifications` in
  `src/sirna_data/raw_loader.py` parses these and slices them to line up
  with whichever window of the raw strand this project actually located
  and stored, so the modification tuple's indices match `sense`/`guide_seq`
  exactly. On the live dataset, **~64% of CMsiRNAdb records resolve real
  per-position chemistry** (the rest either have no annotation for that
  row, or the annotation didn't align cleanly with the located
  sequence -- left as `None` rather than guessed at); common chemistry
  classes found include 2'-OMe/2'-F/phosphorothioate-backbone (the most
  common single combination), 2'-OMe-only, and a smaller GalNAc/lipid-
  conjugate and 5'-vinyl-phosphonate-cap tail matching modern
  ESC-platform-style therapeutic chemistry. For the collapsed
  `_load_cmsirnadb_full_records` groups (multiple raw rows sharing one
  duplex, different assay conditions), modification data comes from the
  same representative row `technology`'s `Cell_Type` is already taken
  from -- replicate measurements of one duplex share one chemical entity.
- **Monopoli2023** (20 records) -- `is_modified=True` with a fixed,
  dataset-level `modification_chemistry` summary ("sdRNA: heavy
  2'-F/2'-OMe/phosphorothioate, cholesterol-conjugated"), since the paper
  describes one architecture applied uniformly to all 20 rows but doesn't
  give a per-position map the way CMsiRNAdb's raw table does --
  `sense_modifications`/`antisense_modifications` stay `None`.
- **Shmushkovich et al. 2018** (356, not integrated -- see
  `POTENTIAL_DATA_SOURCES.md`) -- same chemistry class as Monopoli
  (Monopoli's model was trained on it), but blocked from integration by
  its own no-gene-identity problem, unrelated to modification data.
- **Martinelli 2023 / sirna-repro** (577 of 907 rows integrated -- see the
  dedicated section above) -- a genuine per-siRNA (not per-position)
  modification-type column (e.g. "hexitol nucleic acid", "2-fluoro
  2-O-methyl", or "0" for unmodified); every one of the 577 gene-resolved
  rows carries real modification data, wired up via `is_modified`/
  `modification_chemistry` like Monopoli2023. The other ~330 rows remain
  blocked, mostly by the same still-unreached portion (tables 52-182) of
  patent US20120088815A1, plus a handful of smaller papers whose full text
  wasn't reachable this round or that turned out to be viral rather than
  human targets.
- **5 FDA-approved drugs external-validation set** -- the AttSiOff source
  table had 2'-F/2'-O-Me/phosphorothioate notation, but it was parsed down
  to plain bases during the original verification work (see that section
  above) and, separately, no code/data file for this set exists in the
  repo yet at all -- reconstructing its modification data would mean
  starting over from the source supplement, not just wiring up something
  already committed.
- **Every other source** (siRNAEfficacyDB, Shabalina 2006, siRecords, PDCD1
  panel) -- standard/unmodified synthetic siRNA; `is_modified=False` by
  the schema's default, nothing to wire up.

## Known data-quality caveats (do not "fix" silently — filtered/flagged instead)

- 12 rows have `Accession_number == "-"` (a Renilla luciferase assay control,
  not an endogenous gene) — dropped during acquisition.
- The "Takayuki" EGFP-reporter subset (702 rows, a 1-nt-resolution tiling
  screen across the EGFP CDS) is mapped in siRNAEfficacyDB to accession
  `NZ_CP024869`, whose current RefSeq record is a ~3.7 Mb *Dietzia* sp.
  bacterial genome, not a human reporter plasmid. This looked like a bad
  accession at first, but all 702 target sites (verified) match a single
  tightly-clustered ~700 bp window at one locus of that assembly containing
  the exact canonical EGFP coding sequence — almost certainly lab-plasmid
  contamination baked into that particular genome assembly (a documented,
  recurring issue in NCBI genome submissions), not a citation/mapping error.
  It happens to give correct, self-consistent local sequence context for
  this subset, so `raw_loader.py` uses it as-is.
- ~49 rows (out of 3,532) fail the exact-substring target-site search against
  their nominal transcript (concentrated in accessions `M15077`, `XM_214061`,
  `NM_012864`, `NM_014501`, and a handful of others) — most likely isoform,
  UTR, or cloning-construct differences between the assay's original mRNA
  and the current RefSeq record. These fall back to duplex-only context (see
  `has_flanking_context` in the processed dataset) rather than being dropped.
- A few accessions (e.g. `XM_214061`, `XM_371822`) are old, sometimes
  superseded RefSeq predictions; NCBI still resolves them today but that is
  not guaranteed to remain true indefinitely.
- `%Inhibition` values are real experimental measurements and are noisy: the
  range in the raw data is roughly -27.8 to 134.1 (i.e., below 0% or above
  100%), which is expected assay noise, not a bug. We do not clip these by
  default; see `src/sirna_data/raw_loader.py` for an optional clip.
