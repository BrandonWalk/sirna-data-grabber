# Data-source ledger — trainable data audit

Audit of every siRNA-efficacy source considered for this dataset, recording whether it
provided **trainable data** (siRNA sequence paired with a numeric knockdown/efficacy
value suitable for a regression/classification target), and how much.

Companion machine-readable table: `data_source_ledger.csv`.
Source landscape: `POTENTIAL_DATA_SOURCES.md` · Integrated-source detail: `DATA_SOURCES.md`.

## Bottom line

- **Trainable records currently integrated: 6,577** across **87 genes**, all with a numeric %-knockdown label. (16,178 records / 97 genes if the CMsiRNAdb full-database retrieval is also included — see `CMSIRNADB_FULL_RETRIEVAL.md`.)
- **4 sources** supply that data; **siRNAEfficacyDB (3,532)** and the **CMsiRNAdb PCSK9 subset (2,756)** are 96% of it.
- **siRecords** was recovered (3,117 rated records, ~1,400 new gene accessions) but is **NOT trainable as-is** — its label is a 4-level ordinal rating, not a numeric percentage. It becomes trainable only by tracing its PubMed IDs to the primary papers.
- Everything else is either already embedded (the classic benchmarks), a proxy metric (functional-genomics screens), unavailable/defunct, or reagent-only.

## Integrated — trainable data in this dataset (6,577 records)

| Source | Records | New genes | Metric | Notes |
|---|---|---|---|---|
| **siRNAEfficacyDB** (Zhang 2024) | **3,532** | 41 (baseline) | numeric %Inhibition | Primary source; `sirna_efficacy.csv`. |
| **CMsiRNAdb — PCSK9 subset** (He 2026) | **2,756** | +1 | numeric inhibition | Patent-derived, chemically modified; derived at load time from `cmsirnadb_full_raw.tsv` (CC BY-NC-ND, no derivative file shipped — see `DATA_SOURCES.md`). |
| **Shabalina 2006** | **269** | +41 | numeric (100−Activity) | 269 new after dedup vs 653 in paper; `shabalina_extra.csv`. |
| **Monopoli 2023** | **20** | +4 (APP/MAPT/BACE1/SNCA) | numeric (100−reporter) | Modified sdRNA; `monopoli_extra.csv`. |
| **TOTAL** | **6,577** | **87 genes** | — | — |

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
| Reynolds / Khvorova / Ui-Tei / Vickers | redundant | Cross-embedded in curated DBs |
| Shmushkovich 2018 (356) | **no** | No gene-identity column → cannot group for leave-one-gene-out CV |
| Matveeva 2007 (~3,336) | unknown | Download dead, no Wayback capture, sub-sources paywalled |
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
