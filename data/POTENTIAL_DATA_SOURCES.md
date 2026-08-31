# Potential data sources for siRNA knockdown efficacy

Catalogue of candidate sources for siRNA sequence → measured knockdown/efficacy data,
compiled during data-source review for this dataset. Once a source is fully obtained
and integrated it is removed from this file — it's documented in `DATA_SOURCES.md`
instead. This file stays the broader landscape of what's still only partly obtained,
never pursued, or unavailable.

Legend for **Status**: 🟡 obtained, not integrated · 🔎 candidate · ⛔ unavailable/defunct

## 1. Curated siRNA-efficacy databases

| Source | Content | Efficacy metric | Status | Notes |
|---|---|---|---|---|
| **siRecords** (Ren/Gong 2006–09) | ~4,200 records / ~1,580 accessions | **4-level ordinal** (Very high/High/Medium/Low) | 🟡 this review | live site defunct; 05 release recovered from Internet Archive. Ordinal label, not %; overlaps existing sources ~15% |
| **siRNAdb / HuSiDa** | older curated validated-siRNA sets | mixed | 🔎 | dated; largely superseded |
| Vendor design datasets (Dharmacon/Thermo, Qiagen, Sigma) | proprietary training sets behind design tools | n/a | ⛔ | tools usable, raw data not public; ThermoFisher catalog declined (bot-mitigation/ToS) |

## 2. Classic benchmark assay datasets (the ML staples)

These are the foundational sequence-vs-efficacy sets; note they are largely **already embedded**
in siRNAEfficacyDB / Shabalina2006, so integrating them again mostly adds redundancy.

| Source | Content | Status | Notes |
|---|---|---|---|
| **Katoh & Suzuki / Ichihara (i-Score)** | re-normalized recompilations | 🔎 | documented normalization; heavy overlap with above |
| **siRNABERT** (Xu/Zhao 2024, GENE; github.com/ChengkuiZhao/siRNABERT) | code repo; `data/testdata/{katoh,reynold}.csv` (662+253 rows, current tree) + `data/ExperimentalData/Covid19.csv` (deleted from `main` in commit `d2ad931`, still recoverable from parent commit `f3254ee`) | ⛔ nothing further | katoh.csv: 662/662 already present. reynold.csv: 240/253 already present, remaining 12 have no gene-identity column. Covid19.csv (recovered from git history): 11 candidate sequences but no efficacy/label column in the committed version. (This repo's PDCD1_8.csv panel was integrated — see `DATA_SOURCES.md`.) |
| **Matveeva et al. 2007** (3,336-siRNA compilation) | large compilation | ⛔ | only download link dead, no Wayback, 3/4 sources paywalled |
| **OligoGraph repo** (github.com/drugparadigm/OligoGraph) | third-party GNN repo; ships `Hu/Mix/Taka/Simone.csv` (3,857 rows total, no gene column) | 🟢 343/343 integrated (metric is an interpretation, not a verified true value) | Gene identity for 343 rows (300 `Simone.csv` + 43 `Mix.csv`) was resolved and independently verified by exact 19nt substring match against real NCBI reference transcripts, then integrated — see `DATA_SOURCES.md`'s dedicated OligoGraph section for the full writeup. **The metric caveat still stands**: OligoGraph ships only its own undocumented `label` column (~0–1), not the source papers' true reported %-knockdown. A `label*134.1=true-%` fit was exact on `Hu.csv` (R²=1.0) but failed on `Mix.csv` (R²=0.82, same label value mapping to four different true percentages) — each source file appears to be independently max-normalized to its own scale, not one shared constant. Sciabola et al. 2013's Supplementary Table S4 (blocked: Oxford CDN 403, PubMed/PMC reCAPTCHA, Europe PMC 429) and Harborth et al. 2001's paper (read in full — doesn't actually contain the claimed 44-siRNA panel at all) were both pursued directly and neither yielded true values; Huesken et al. 2005's own supplementary data is the more plausible true source for the Mix.csv/lamin-A-C rows specifically but hasn't been attempted. Given that, `label*100` is used directly as %KD **per explicit instruction, as a stated interpretation** — not an independently verified true value. `Hu.csv`'s 2,361 rows all overlap this dataset already; `Taka.csv` (702 rows) has no discoverable citation and wasn't pursued. |
| **siRNA_Design** (Neurotechnology-at-ETH-Zurich, GitHub) | Design-tool scraping + scoring pipeline repo; ships 4 "validation" datasets in `sirna_datasets_analysis/datasets_excel/` | ⛔ nothing new | Dataset A = Huesken 2005 (2,431, already embedded via siRNAEfficacyDB). Dataset B = Ichihara 2007, itself a compilation of Reynolds/Vickers/Haborth/Ui-Tei/Khvorova (419, already embedded). Dataset C = Mysara 2012, a 38-record subset of B (redundant subset). Dataset D = Fellmann 2011 via Mysara (476) — **shRNA, not siRNA**, same disqualifier as shRNA screens below (shRNA ≠ siRNA potency). The repo's own in-vitro validation (3 candidates, >90% knockdown) is narrative only, not structured data. |

## 3. Chemically-modified / therapeutic siRNA data (drug-relevant)

| Source | Content | Status | Notes |
|---|---|---|---|
| **Shmushkovich et al. 2018** (NAR) | 356 siRNAs, same chemistry | 🔎 | no gene-identity column → augmentation-only, not gene-groupable |
| **Martinelli 2023 / sirna-repro** (bioRxiv, siRNAmod-derived) | 907 chemically-modified siRNAs, 30 source documents | 🟡 577/907 integrated, ~330 remain blocked | 577 of 907 rows have now been resolved to a real gene/reporter target and integrated — see `DATA_SOURCES.md`. Pass 1 traced 253 rows to their source documents; Pass 2 resolved another 324 rows by brute-force exact 19nt substring matching against real fetched reference transcripts (no trust placed in any source table's own labels): 64 rows reused the classic GL2 anti-luciferase site (X65324) already present under a different reporter accession, and 260 rows matched five genes/orthologs (KAZRIN, MIR155HG/BIC, CDKN1B, chimp SOD2, rhesus CASR) named in the worked examples of patent US20120088815A1. The remaining ~330 rows are still blocked; the single largest chunk (284 rows, 86% of what's left) is the rest of that same patent — its own text confirms it tested up to 182 real genes across tables 1–182, but WebFetch can only surface the reachable "Example" tables (52 of them) and cannot paginate the full document, so the remaining genes' names are unreachable in this environment. A handful of smaller documents (PMID 16598842, 25699137, 22260772, and the MDR1/ABCB1 pair 21141919+22982308) remain unresolved behind paywalls or reCAPTCHA-blocked full text. Several other newly-identified documents (HIV, Enterovirus 71, coxsackievirus B3) are viral rather than human targets and weren't pursued further. Full raw 907-row CSV kept locally, gitignored (not redistributed) — see `data/raw/sirna_repro_martinelli_907.csv`. **Checked the live upstream directly (crdd.osdd.net/servers/sirnamod/)**: its own bulk CSV export currently returns only 414 rows (all one PMID) against the site's claimed 4,894 — looks broken/stale (old 2016 OSDD/IMTECH server). Not that it matters: the search page's full field list (ID, PMID, per-strand modification name/position/component, activity, assay/cell/transfection) confirms siRNAmod has no gene/target-identity field at the schema level, so a working export would hit the exact same disqualifier. |
| Alnylam / Arrowhead patents & pubs | modification-vs-activity | 🔎 | scattered in med-chem literature & patent filings |

## 4. Large-scale functional-genomics screens (indirect, proxy efficacy)

Genome-scale shRNA/siRNA reagent screens — relative depletion/activity, not clean reporter %KD.
Useful only if proxy on-target efficacy is acceptable.

| Source | Content | Status | Notes |
|---|---|---|---|
| **DepMap / Project DRIVE / Achilles (RNAi)** | genome-scale shRNA dropout, 100s of cell lines | 🔎 | DEMETER/DEMETER2 estimate reagent efficacy; shRNA ≠ siRNA potency |
| **The RNAi Consortium (TRC) / pLKO** | shRNA library sequences + screen performance | 🔎 | reagent definitions + phenotype |
| **GenomeRNAi** | aggregated RNAi screen phenotypes (human, Drosophila) | 🔎 | phenotype-level, not %KD |

## 5. General repositories to mine for raw measurements

| Source | Content | Status | Notes |
|---|---|---|---|
| **GEO / ArrayExpress (BioStudies)** | knockdown-validation experiments, supp. tables | 🔎 | programmatically queryable |
| **PubMed / PMC supplementary tables** | per-paper %KD in supplements | 🔎 | siRecords' PubMed IDs are an index into this |
| **Addgene** | shRNA/siRNA construct sequences | 🔎 | reagent/sequence side |

## Recurring caveats when merging any of these

1. **Metric mismatch** — %mRNA KD (qPCR) ≠ %protein KD (Western) ≠ reporter-normalized inhibition ≠ ordinal rating (siRecords). Harmonize before pooling. This is the single biggest pitfall.
2. **Assay conditions** — concentration, timepoint, cell line, transfection/delivery all shift measured KD. Curated benchmark sets fix these; aggregated DBs do not.
3. **Unmodified vs chemically modified** — a hard domain gap; models trained on unmodified academic data generalize poorly to modified therapeutic siRNAs (folded as plain RNA here — see DATA_SOURCES.md).
4. **shRNA ≠ siRNA potency** — screen scores fold in Dicer processing and expression level.
5. **Redundancy** — the classic benchmarks are heavily cross-embedded in the curated DBs; dedup on exact (strand-agnostic) sequence before adding, as this project already does.
6. **Licensing** — verify redistribution terms (OligoFormer's prepared CSVs, e.g., are proprietary-licensed and were deliberately not vendored).
