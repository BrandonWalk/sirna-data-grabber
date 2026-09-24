# Data sources

## Summary

Audit of every siRNA-efficacy source integrated into this dataset: whether it provided
**trainable data** (siRNA sequence paired with a numeric knockdown/efficacy value
suitable for a regression/classification target), and how much. For the broader
landscape of sources considered — including ones never obtained or not pursued
further — see `POTENTIAL_DATA_SOURCES.md`.

- **Trainable records currently integrated: 7,518** across **108 genes**, all with a
  numeric %-knockdown label. (18,085 records / 117 genes if the CMsiRNAdb full-database
  retrieval and Davis2025 below are also included — both on by default via
  `include_cmsirnadb_full=True` / `include_davis2025=True`; see their own sections.)
- **7 sources** supply that headline 7,510; **siRNAEfficacyDB (3,532)** and the
  **CMsiRNAdb PCSK9 subset (2,756)** are 84% of it.
- **Chemical modification data**: `SiRNARecord` carries `is_modified`/
  `modification_chemistry`/`sense_modifications`/`antisense_modifications` — see
  "Chemical modification data" below.

### Integrated — trainable data (7,518 records)

| Source | Records | New genes | Metric | Notes |
|---|---|---|---|---|
| **siRNAEfficacyDB** (Zhang 2024) | **3,488** | 40 (baseline) | numeric %Inhibition | Primary source; `sirna_efficacy.csv`. Its 44-row `Lamin A` block is corrupt (44 copies of one duplex) and is superseded by the Harborth 2003 section below — 3,532 rows load when `include_harborth2003=False`. |
| **CMsiRNAdb — PCSK9 subset** (He 2026) | **2,756** | +1 | numeric inhibition | Patent-derived, chemically modified; derived at load time from `cmsirnadb_full_raw.tsv` (CC BY-NC-ND, no derivative file shipped). |
| **Shabalina 2006** | **269** | +41 | numeric (100−Activity) | 269 new after dedup vs 653 in paper; `shabalina_extra.csv`. |
| **Monopoli 2023** | **20** | +4 (APP/MAPT/BACE1/SNCA) | numeric (100−reporter) | Modified sdRNA; `monopoli_extra.csv`. |
| **Martinelli 2023 / sirna-reproduction** | **577** | +12 (EGFP, ACP5, APOB, Luciferase_firefly, Luciferase_renilla, NPY, VEGFA, KAZRIN, MIR155HG, CDKN1B, SOD2_chimp, CASR_rhesus) | numeric (PCT) | 577 of 907 rows resolved — 253 by source-document tracing, 324 more by brute-force sequence matching. See the section below. |
| **Sciabola et al. 2013 in-house panel** (NAR Supp. Tables S3/S4) | **356** | +4 (BIRC5, EZH2, MTOR, BRAF; also re-sources HIF1A/HK2/HPSE and adds rows to MYC/CTNNB1/PIK3CA) | numeric %Inhibition (mean of the doses each row was screened at) | The paper's own Hep3B/QuantiGene 2.0 measurements, CC BY-NC 3.0. `sciabola2013_extra.csv`. |
| **Harborth et al. 2003 lamin A/C panel** (via Ichihara et al. 2007) | **44** | +0 (replaces the corrupt `Lamin A` block in the primary source) | numeric % Inhibition, **protein-level** (immunoblot, HeLa) | The real 44-duplex panel, from Ichihara 2007's CC BY-NC 2.0 UK republication. Supersedes both second-hand copies the corpus used to carry. `harborth2003_extra.csv`. |
| **TOTAL** | **7,518** | **108 genes** | — | — |

### Obtained, not integrated

| Source | Records | Why not trainable | Path to make it trainable |
|---|---|---|---|

### Definition of "trainable" used here

A record is trainable here if it pairs an siRNA **sequence** with a **numeric**
knockdown/inhibition value (continuous %), suitable as a regression or classification
target. Ordinal ratings and other proxy metrics are excluded on that
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

We deliberately did NOT use the merged CSVs from the Model B GitHub repo
(github.com/lulab/Model B), even though it bundles the same classic
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
  was derived from that file's "Supplemental Table 1" sheet. **That manual step
  is no longer needed**: Europe PMC's public supplementaryFiles API serves the
  same `gkaf479_supplemental_files.zip` with no challenge, so
  `sirna_data.fetch.davis2025` now fetches and filters it in-process (the
  workbook is read by `sirna_data/fetch/_xlsx.py`, no openpyxl required).
  Verified: the fetcher reproduces the committed `davis2025_extra.csv` row for
  row, every cell.
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

## Supplementary siRNA data: Martinelli 2023 / sirna-reproduction (577 of 907 rows, 12 new genes/reporters)

- Martinelli 2023, bioRxiv preprint / `sirna-reproduction` (siRNAmod-derived
  reproduction dataset). 907 chemically-modified siRNAs pooled from 30
  source PMIDs/patents; only the 577 rows whose targets could be confirmed
  are integrated here.
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
- Fetched by `sirna_data.fetch.martinelli` (`sirna-data-fetch`) from the
  upstream 907-row table on `raw.githubusercontent.com`, with gene identity
  recovered mechanically: the longest window of either strand that occurs in
  exactly ONE of the 13 traced transcripts decides the row's gene, and a
  window matching two (the firefly luciferase reporters overlap) is dropped
  rather than guessed. Verified: **a re-fetch reproduces the committed
  `martinelli_extra.csv` byte for byte.**
- **Canonical form**: the hand-built file originally mixed conventions --
  DNA overhangs written `T` on some rows and `U` on others, `PCT` as both
  `0` and `0.0`, rows in no particular order -- so a re-fetch would have
  rewritten 473 of 577 rows. `canonical_row` in the fetcher defines the form
  (sequences in the RNA alphabet, since the chemistry lives in the
  modification columns rather than the letters; percentages as numbers; rows
  ordered by experiment number), `fetch()` writes it, and
  `python -m sirna_data.fetch.martinelli --normalize data/raw/martinelli_extra.csv
  [--transcripts data/raw/martinelli_transcripts.fasta]` rewrites an existing
  file into it. That command has already been run over the committed file,
  which is why the two now agree exactly; it is idempotent, so running it
  again is a no-op.
- **What the stored target site is**: for a row placed by its sense strand,
  the verified window spliced back into that strand (so the source's 2nt
  overhang survives). For the 28 rows placed by their antisense strand --
  where the source's sense column does not match any transcript -- the
  reverse complement's last 19 bases, since the 2nt reverse-complementing
  the guide's 3' overhang are not part of the target.
- Loaded by `_load_martinelli_records` in `src/sirna_data/raw_loader.py`
  from `data/raw/martinelli_extra.csv` (the 577-row derivative with added
  gene/accession columns, and corrected sense sequences for the 28 pass-1
  rows noted above) and `data/raw/martinelli_transcripts.fasta` (the 13
  reference sequences above, 7 from pass 1 plus 6 from pass 2). **License**:
  the source `sirna-reproduction` dataset is CC BY-NC 4.0, which — unlike
  CMsiRNAdb's CC BY-NC-ND — permits derivatives, so this filtered/
  gene-annotated subset is committed to the repo like Monopoli/Shabalina.

## Supplementary siRNA data: Harborth et al. 2003 lamin A/C panel (44 rows, via Ichihara et al. 2007)

- **The experiments**: Harborth, Elbashir, Vandenburgh, Manninga, Scaringe,
  Weber & Tuschl 2003, *Antisense Nucleic Acid Drug Development* 13:83-105,
  "Sequence, chemical, and structural variation of small interfering RNAs and
  short hairpin RNAs and the effect on mammalian gene silencing"
  (doi:10.1089/108729003321629638, PMID 12804036). 44 standard 21-nt duplexes
  tiling lamin A/C, read out as **lamin A/C protein** by immunoblot /
  immunofluorescence in **human HeLa** cells (mouse SW3T3 served as a second,
  ~2x less responsive cell type).
- **The values**: Harborth 2003 is **not** open access (Mary Ann Liebert; no
  PMC copy, subscription archive) and carries no reuse license. The numbers
  loaded here come instead from its republication in **Ichihara, Naito, Ui-Tei,
  Sakuma, Juni, Ueda & Saigo 2007**, *Nucleic Acids Research* 35(18):e123,
  "Thermodynamic instability of siRNA duplex is a prerequisite for dependable
  prediction of siRNA activities" (doi:10.1093/nar/gkm699, PMC2094068),
  **CC BY-NC 2.0 UK**: "© 2007 The Author(s) This is an Open Access article
  distributed under the terms of the Creative Commons Attribution
  Non-Commercial License (http://creativecommons.org/licenses/by-nc/2.0/uk/)
  which permits unrestricted non-commercial use, distribution, and
  reproduction in any medium, provided the original work is properly cited."
  Its supplementary workbook (`nar_gkm699_...File007.xls`, fetched from PMC)
  tabulates every dataset the i-Score model was built on with a per-row
  `Authors` column; the 44 rows whose author is `Harborth` carry gene,
  accession, antisense 21-mer, sense 19-mer and % inhibition. **Cite both
  papers** when using this subset: Harborth for the experiments, Ichihara for
  the compilation these numbers were taken from.
- **Why this replaced two other copies of the same panel.** Before this
  source existed, the corpus carried the panel twice, badly, and neither copy
  was licensed:
  1. **siRNAEfficacyDB's `Lamin A` block is corrupt.** All 44 of its rows are
     byte-identical: same sequence pair (`GAGCUCCUGCAGGUCCUCCuu` /
     `GGAGGACCUGCAGGAGCUC` -- this panel's B1), same 83.0% label, same cell
     line, dose and timepoint. Checked against the whole file: every other
     duplex in `sirna_efficacy.csv` appears exactly **once** (3,488 duplexes
     at x1, this one at x44), and no duplex anywhere in that file carries two
     different labels, so this is not a replicate-measurement pattern -- rows
     that should hold B2-B44 hold B1's data. siRNAEfficacyDB's own accounting
     lists "Harborth x44", so the intent was clearly this panel. Those 44 rows
     are dropped in favour of this source (`_load_sirnaefficacydb_records`'s
     `superseded_genes`).
- **Label**: `% Inhibition` as reported -- **protein-level** knockdown in
  HeLa, not the mRNA-level %inhibition most of this corpus carries. Treat it
  as a different measurement family (same caveat class as Martinelli's
  PCT/reporter rows). Sanity check against the paper's own abstract ("26 of
  44 tested standard 21-nt siRNA duplexes reduced the protein expression by at
  least 90%"): these rows give 25 at >=90 and 1 below 50, the one-row
  difference being a boundary/rounding call in the abstract's count.
- **Accession**: the paper and Ichihara's table cite `AH001498`, a segmented
  GenBank gene record; only 42 of the 44 sense strands locate in it, while all
  44 locate in the RefSeq mRNA `NM_170707`, so records carry `NM_170707` and
  the raw CSV preserves the table's own accession as `Source_Accession`. This
  also means the `Lamin A` gene group no longer spans two accessions.
- **Verification**: both strands come from the raw table rather than being
  derived by revcomp, and were cross-checked against each other --
  `Antisense_21mer[:19]` is the exact reverse complement of `Sense_19mer` for
  all 44 rows -- and every sense strand was located by exact substring match
  in `NM_170707`.
- Fetched by `sirna_data.fetch.harborth2003` (`sirna-data-fetch`) from
  Europe PMC's supplementaryFiles API for PMC2094068: the workbook is read
  in-process (`sirna_data/fetch/_ole.py`, no xlrd/LibreOffice needed), the
  rows whose `Authors` column is `Harborth` are kept, and every sense strand
  is checked against `NM_170707` before anything is written. Verified: the
  fetcher reproduces the committed CSV row for row.
- Loaded by `_load_harborth2003_records` in `src/sirna_data/raw_loader.py`
  from `data/raw/harborth2003_extra.csv` (44 rows) and
  `data/raw/harborth2003_transcripts.fasta` (`NM_170707`, pure NCBI RefSeq,
  public domain). Tests: `tests/test_harborth2003.py`,
  `tests/test_fetch_harborth2003.py`, `tests/test_fetch_ole.py`.

## Supplementary siRNA data: Sciabola et al. 2013 in-house panel (356 rows, 4 new genes)

- Sciabola, Cao, Orozco, Faustino & Stanton 2013, *Nucleic Acids Research*
  41(3):1383-1394, "Improved nucleic acid descriptors for siRNA efficacy
  prediction" (doi:10.1093/nar/gks1191, PMID 23241392, PMC3561943).
  **CC BY-NC 3.0** -- the article's own stated terms, printed on PMC's
  copyright block: "This is an Open Access article distributed under the
  terms of the Creative Commons Attribution License
  (http://creativecommons.org/licenses/by-nc/3.0/), which permits
  non-commercial reuse, distribution, and reproduction in any medium,
  provided the original work is properly cited." Supplementary data rides on
  the article's license, so the derived CSV is redistributable here with
  attribution, non-commercially -- same footing as siRNAEfficacyDB.
- **What it is**: the paper's own in-house panel, run at the Pfizer
  Oligonucleotide Therapeutic Unit -- not a literature compilation. It
  designed 21-nt siRNAs ("19p2": 19nt duplex core + 2nt target-matching 3'
  overhangs, per the paper's Figure 4) against ten hepatocellular-carcinoma
  relevant genes and screened them at up to four concentrations.
  Supplementary **Table S3** carries the sequences, Supplementary **Table
  S4** the measured % inhibition per dose.
- **Assay**: "Hep3B cells (American Type Culture Collection) were grown in
  EMEM (ATCC) supplemented with 10% fetal calf serum", transfected with
  Lipofectamine RNAiMAX, and "QuantiGene 2.0 assay (Affymetrix Inc. Santa
  Clara, CA) was used to measure the expression level of target genes before
  and after knockdown in Hep3B cell lines" 48 h post-transfection. So every
  row is mRNA-level knockdown in Hep3B -- the same readout family as
  siRNAEfficacyDB's qPCR rows and Davis 2025's QuantiGene rows.
- **Obtaining the file**: Oxford's own supplementary CDN link for this
  article returns HTTP 403, but PMC hosts the same file --
  `supp_gks1191_nar-01814-n-2012-File002.doc`, linked from PMC3561943, a
  21-page Word document containing Tables S1-S4. It was fetched through a
  browser (PMC gates supplementary downloads behind its own client-side
  check; no CAPTCHA was solved and nothing was scraped around it), converted
  browser and, at first, converted with LibreOffice. That manual step is no
  longer needed: `sirna_data.fetch.sciabola2013` (`sirna-data-fetch --only
  sciabola2013`) now pulls the same file through Europe PMC's public
  supplementaryFiles API and parses the Word tables directly, with no
  dependency beyond pandas and the standard library (see
  `sirna_data/fetch/_ole.py`). Verified: the fetcher reproduces the committed
  `sciabola2013_extra.csv` row for row, so no local copy of the source
  document is kept in this repo.
- **Rows**: Table S3/S4 list 361 21-mers -- HIF1A 100, HK2 99, HPSE 110, and
  seven follow-up genes at 5-10 each (MTOR/`FRAP1` 10, BIRC5/`SURVIVIN` 9,
  MYC/`c-Myc` 9, EZH2 9, BRAF/`bRAF` 5, CTNNB1 5, PIK3CA 5). 356 are
  integrated. The table's 591 R-/L-Dicer substrates (25/27-mers, a different
  duplex architecture) are **not** loaded -- a candidate for a future
  `include_*` subset, not part of this addition.
- **Gene identity, verified per row, and a strand surprise**: Table S3's
  sequence column turned out to MIX orientations -- 333 of the 361 rows are
  the antisense (guide) strand and 23 are the sense strand (all 10 MTOR
  rows, all 5 PIK3CA, 5 of 9 BIRC5, 2 of 9 MYC, 1 of 5 BRAF). Every row was
  therefore tested in both orientations against a freshly fetched NCBI
  RefSeq transcript for its gene, and `sciabola2013_extra.csv` stores the
  normalized **sense** strand (`Sense_21mer`) so the loader follows this
  repo's usual convention (store sense, derive the guide by revcomp), with
  the table's own string and the orientation it proved to be preserved
  alongside as `Table_S3_Sequence` / `Table_S3_Strand`. 355 of the 356 kept
  rows locate as a full 21nt exact match (the design's overhangs are
  target-matching); 1 locates by its 19nt core.
- **Dropped**: 5 of 361 rows (HK2 `19p2_103`/`19p2_133`, HPSE
  `19p2_240`/`19p2_241`/`19p2_242`) match no checked transcript in either
  orientation -- 6-8 mismatches at their best alignment, so not a
  transcription typo -- and were dropped under the same verify-or-drop rule
  every other source here follows. HPSE variants `NM_001098540.3` /
  `NM_001166498.3` and the older `NM_006665.2` / `NM_000189.3` records were
  checked too. One HPSE row (`19p2_244`) locates only in `NM_001098540.3`
  and carries that accession, which is why HPSE shows two accessions in the
  README gene table.
- **Gene labels are normalized** to current official symbols so they group
  with the rest of this corpus instead of fragmenting it (`HIF1a` -> `HIF1A`,
  `SURVIVIN` -> `BIRC5`, `c-Myc` -> `MYC`, `FRAP1` -> `MTOR`, `bRAF` ->
  `BRAF`). This is a deliberate exception to the "raw gene strings are not
  normalized" rule that keeps `EGFP`/`EGFP ` apart: there the two spellings
  are genuinely different source entries, whereas here they are the same
  genes other sources already contribute under their official symbols. Each
  row's original spelling is kept as `Source_Gene_Label`.
- **LABEL -- the one thing to read before using these labels.** Which doses
  exist differs by gene: HIF1A/HK2 at 0.08/0.4/2/10 nM, HPSE at 0.08/0.4/2 nM
  (the paper stopped at 2 nM for HPSE -- "the hit rate was high enough at
  2 nM to not require a higher concentration"), and the seven follow-up genes
  at 10 nM only. `label` is the **mean of the doses present for that row**,
  so it is a dose-averaged potency for the three tiled genes but a single
  10 nM reading for the follow-up genes -- not an apples-to-apples quantity
  across genes, and systematically lower than a top-dose label for
  HIF1A/HK2/HPSE (e.g. `19p2_1`: 12.0 as a 3-dose mean vs 33.0 at 10 nM).
  All four dose columns are shipped verbatim in the raw CSV
  (`Pct_Inhibition_{008,04,2,10}nM`), and `_load_sciabola2013_records`
  computes the mean in code, so switching to a single-dose or
  top-dose-only label is a one-line change rather than a re-extraction.
  Negative values are the source's own (measured expression above untreated
  control) and are kept as reported.
- Loaded by `_load_sciabola2013_records` in `src/sirna_data/raw_loader.py`
  from `data/raw/sciabola2013_extra.csv` (356-row derivative) and
  `data/raw/sciabola2013_transcripts.fasta` (11 RefSeq transcripts, pure
  NCBI RefSeq and therefore public domain). Tests: `tests/test_sciabola2013.py`.

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
  these 12 genes are both derived from it at load time) and `data/raw/cmsirnadb_full_transcripts.fasta`
  (34 mRNA transcripts, NCBI efetch).
- **Integration**: wired into `load_records()` via `include_cmsirnadb_full=True`
  (`src/sirna_data/raw_loader.py`, `_load_cmsirnadb_full_records()`), which builds a
  strand-agnostic index of every sequence already loaded from the other sources before
  running (step 4 above). Effect on the dataset (verified by running `load_records()`):
  without, **6,577 records / 87 genes**; with, **16,178 records / 97 genes** (+9,601
  records, +10 new genes: AGT ANGPTL3 CTNNB1 HSD17B13 INHBE LPA MARC1 MSTN PLN PNPLA3;
  APP & MAPT deepened). Anyone consuming this data with a cached/precomputed downstream
  representation (e.g. a graph cache) should rebuild it after pulling this update.

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
- **Martinelli 2023 / sirna-reproduction** (577 of 907 rows integrated -- see the
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
- **Every other source** (siRNAEfficacyDB, Shabalina 2006) -- standard/
  unmodified synthetic siRNA; `is_modified=False` by
  the schema's default, nothing to wire up.

## Known data-quality caveats (do not "fix" silently — filtered/flagged instead)

- **siRNAEfficacyDB's `Lamin A` block is corrupt (44 identical rows).** All
  44 rows of that gene are byte-identical — same sequence pair, same 83.0%
  label, same cell line, dose and timepoint — i.e. 44 copies of one duplex
  (Harborth et al. 2003's B1) where rows B2–B44 of that panel should be. This
  is the only such block in the file: a duplicate-key scan of all 3,532 rows
  finds 3,488 duplexes appearing exactly once and this one appearing 44
  times, with no duplex anywhere carrying two different labels (so it is not
  a replicate-measurement pattern). Not silently deduped: the whole block is
  superseded by the real panel, loaded from Ichihara et al. 2007 — see the
  Harborth 2003 section above. `include_harborth2003=False` loads the corrupt
  rows again, unchanged.
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
