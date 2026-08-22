# Data-source ledger — trainable data audit

Audit of every siRNA-efficacy source considered for this dataset, recording whether it
provided **trainable data** (siRNA sequence paired with a numeric knockdown/efficacy
value suitable for a regression/classification target), and how much.

Companion machine-readable table: `data_source_ledger.csv`.
Source landscape: `POTENTIAL_DATA_SOURCES.md` · Integrated-source detail: `DATA_SOURCES.md`.

## Bottom line

- **Trainable records currently integrated: 6,838** across **95 genes**, all with a numeric %-knockdown label. (16,439 records / 105 genes if the CMsiRNAdb full-database retrieval is also included — see `CMSIRNADB_FULL_RETRIEVAL.md`.)
- **6 sources** supply that data; **siRNAEfficacyDB (3,532)** and the **CMsiRNAdb PCSK9 subset (2,756)** are 92% of it.
- **siRecords** was recovered (3,117 rated records, ~1,400 new gene accessions) but is **NOT trainable as-is** — its label is a 4-level ordinal rating, not a numeric percentage. It becomes trainable only by tracing its PubMed IDs to the primary papers.
- Everything else is either already embedded (the classic benchmarks), a proxy metric (functional-genomics screens), unavailable/defunct, or reagent-only.
- **Chemical modification data**: `SiRNARecord` now carries `is_modified`/`modification_chemistry`/`sense_modifications`/`antisense_modifications`. CMsiRNAdb's raw table already had real per-position modification chemistry that was previously parsed for sequence/label only and discarded -- now wired up (~64% of CMsiRNAdb's ~12,357 records resolve real per-position data). Monopoli2023 gets a coarse dataset-level tag. See "Chemical modification data" in `DATA_SOURCES.md`.

## Integrated — trainable data in this dataset (6,838 records)

| Source | Records | New genes | Metric | Notes |
|---|---|---|---|---|
| **siRNAEfficacyDB** (Zhang 2024) | **3,532** | 41 (baseline) | numeric %Inhibition | Primary source; `sirna_efficacy.csv`. |
| **CMsiRNAdb — PCSK9 subset** (He 2026) | **2,756** | +1 | numeric inhibition | Patent-derived, chemically modified; derived at load time from `cmsirnadb_full_raw.tsv` (CC BY-NC-ND, no derivative file shipped — see `DATA_SOURCES.md`). |
| **Shabalina 2006** | **269** | +41 | numeric (100−Activity) | 269 new after dedup vs 653 in paper; `shabalina_extra.csv`. |
| **Monopoli 2023** | **20** | +4 (APP/MAPT/BACE1/SNCA) | numeric (100−reporter) | Modified sdRNA; `monopoli_extra.csv`. |
| **PDCD1 panel (Xu/Zhao 2024, siRNABERT)** | **8** | +1 (PDCD1) | numeric (qPCR knockdown %) | Recovered from deleted file in repo's git history, verified against NM_005018.3; license unresolved but `pdcd1_extra.csv` is committed at the user's explicit request — see `DATA_SOURCES.md`/`NOTICE.md`. |
| **Martinelli 2023 / sirna-repro** | **253** | +7 (EGFP, ACP5, APOB, Luciferase_firefly, Luciferase_renilla, NPY, VEGFA) | numeric (PCT) | 253 of 907 rows traced from no-gene-identity to a real target: 221 by tracing each source PMID/patent and verifying by exact 19nt substring match, 28 more by discovering their stored "sense" field was corrupted (typo, or a single copy/fill-down error across 23 rows) and deriving the correct sense from their verified-correct antisense strand instead, 4 more (VEGFA) by direct match. ~654 rows remain unresolved (60% is one still-unidentified patent; several newly-identified documents are viral, not human, targets). Chemically modified. CC BY-NC 4.0 permits derivatives — `martinelli_extra.csv` committed. See `DATA_SOURCES.md`. |
| **TOTAL** | **6,838** | **95 genes** | — | — |

## Obtained but NOT integrated

| Source | Records | Why not trainable | Path to make it trainable |
|---|---|---|---|
| **siRecords** (05 release) | 3,295 new by seq / 3,117 rated | Efficacy is **ordinal** (Very high/High/Medium/Low), not numeric %. Also: redistribution rights unresolved — see `NOTICE.md` / `sirecords_overlap_analysis.md` | Trace the per-record PubMed IDs to primary papers and extract reported %KD |

## Candidate — trainable in principle, not pursued

| Source | Est. records | Metric | Blocker |
|---|---|---|---|
| CMsiRNAdb — full (other 12 genes) | ~40,000 | numeric | Highly redundant, patent-derived, all modified chemistry; PCSK9 already taken |
| Katoh / Ichihara i-Score | small | numeric (renormalized) | Recompilation of benchmarks already embedded |
| GEO / ArrayExpress | unbounded | supp. tables | Needs per-experiment manual extraction |
| PubMed / PMC supplements | unbounded | per-paper | Manual extraction; also the route for siRecords |

## Rejected / already embedded / unavailable

| Source | Trainable? | Reason |
|---|---|---|
| Huesken 2005 (~2,431) | redundant | Core of siRNAEfficacyDB — already in |
| siRNABERT test data (katoh.csv 662 + reynold.csv 253 + Covid19.csv 11) | **no (99% redundant / no label)** | katoh.csv + reynold.csv: sequence-dedup vs. baseline shows 902/915 already present; remaining 12 have no gene-identity column. Covid19.csv (recovered from deleted git history, commit `f3254ee`): 11 candidate sequences with no efficacy/label column in the committed version. The genuinely new, usable part of this repo — the `PDCD1_8.csv` panel, also recovered from the same deleted history — **is integrated**; see the Integrated table above. |
| Reynolds / Khvorova / Ui-Tei / Vickers | redundant | Cross-embedded in curated DBs |
| Shmushkovich 2018 (356) | **no** | No gene-identity column → cannot group for leave-one-gene-out CV |
| Martinelli 2023 / sirna-repro, remaining ~654 of 907 rows | **not yet** | Same disqualifier as Shmushkovich 2018 for the rows that couldn't be traced to a real gene — no gene-identity column in the source, and per-document tracing didn't resolve them (60% is one still-unidentified 182-gene patent; several other documents identified as viral, not human, targets -- HIV, EV71, coxsackievirus B3; see `DATA_SOURCES.md`). The 253 rows that WERE traced and verified are integrated — see the Integrated table above. |
| Matveeva 2007 (~3,336) | unknown | Download dead, no Wayback capture, sub-sources paywalled |
| siRNA_Design (Neurotechnology-at-ETH-Zurich) | **no (nothing new)** | Tools/scraper pipeline repo, not a data source. Its 4 bundled datasets are Huesken 2005 (already embedded), Ichihara 2007 = Reynolds/Vickers/Haborth/Ui-Tei/Khvorova compilation (already embedded), a 38-record Mysara subset of that (redundant), and a 476-record Fellmann shRNA subset (excluded — shRNA ≠ siRNA, same as DepMap/TRC). No new primary data. |
| siRNAdb / HuSiDa (legacy) | unknown | Defunct/superseded; not verifiable this session |
| Vendor design datasets | **no** | Proprietary; raw data not public (ThermoFisher declined per ToS) |
| DepMap / DRIVE / Achilles, TRC/pLKO | **proxy only** | shRNA depletion, not reporter %KD; shRNA ≠ siRNA potency; Geff is modeled, not measured. Retrieved and briefly loadable, then removed entirely — see `FUNCTIONAL_GENOMICS_SCREENS.md`. |
| GenomeRNAi | **proxy only** | phenotype calls, not reporter %KD |
| Addgene | **no** | Construct sequences only, no efficacy |
| 5 FDA-approved siRNA drugs | **no (held out)** | External validation set, deliberately not training data |

## Definition of "trainable" used here

A record is trainable here if it pairs an siRNA **sequence** with a **numeric**
knockdown/inhibition value (continuous %), suitable as a regression or classification
target. Ordinal ratings (siRecords), depletion scores (DepMap), and phenotype calls
(GenomeRNAi) are excluded on that definition even though they encode efficacy
information — using them would require either a different model formulation or
converting the metric first.

_Note: legacy-database live-status checks (HuSiDa, siRNAdb, GenomeRNAi, full CMsiRNAdb) were
blocked by the sandbox network allowlist this session, so their "unavailable" status reflects
prior findings and documented history, not a fresh fetch._
