# Potential data sources for siRNA knockdown efficacy

Catalogue of candidate sources for siRNA sequence → measured knockdown/efficacy data,
compiled during data-source review for this dataset. Sources already integrated are
documented in `DATA_SOURCES.md`; this file is the broader landscape
of what exists and what is realistically obtainable.

Legend for **Status**: ✅ integrated · 🟡 obtained, not integrated · 🔎 candidate · ⛔ unavailable/defunct

## 1. Curated siRNA-efficacy databases

| Source | Content | Efficacy metric | Status | Notes |
|---|---|---|---|---|
| **siRNAEfficacyDB** (Zhang 2024) | 3,544 records / 42 genes | numeric %Inhibition | ✅ primary | `cellknowledge.com.cn/siRNAEfficacy`, CC BY-NC |
| **CMsiRNAdb** (He 2026) | 43,153 entries / 13 genes, chemically modified | numeric inhibition | ✅ PCSK9 subset only | patent-derived; same group as siRNAEfficacyDB |
| **siRecords** (Ren/Gong 2006–09) | ~4,200 records / ~1,580 accessions | **4-level ordinal** (Very high/High/Medium/Low) | 🟡 this review | live site defunct; 05 release recovered from Internet Archive. Ordinal label, not %; overlaps existing sources ~15% |
| **siRNAdb / HuSiDa** | older curated validated-siRNA sets | mixed | 🔎 | dated; largely superseded |
| Vendor design datasets (Dharmacon/Thermo, Qiagen, Sigma) | proprietary training sets behind design tools | n/a | ⛔ | tools usable, raw data not public; ThermoFisher catalog declined (bot-mitigation/ToS) |

## 2. Classic benchmark assay datasets (the ML staples)

These are the foundational sequence-vs-efficacy sets; note they are largely **already embedded**
in siRNAEfficacyDB / Shabalina2006, so integrating them again mostly adds redundancy.

| Source | Content | Status | Notes |
|---|---|---|---|
| **Huesken et al. 2005** (Nat. Biotechnol.) | ~2,431 siRNAs / 34 genes, normalized inhibition | ✅ via siRNAEfficacyDB | the BIOPREDsi training set; single most-reused benchmark |
| **Reynolds et al. 2004** | 180 siRNAs, 2 genes | ✅ via compilations | origin of "Reynolds rules" |
| **Khvorova / Ui-Tei / Amarzguioui / Vickers** | small foundational sets | ✅ via compilations | mid-2000s design-rule papers |
| **Shabalina, Spiridonov & Ogurtsov 2006** (BMC Bioinf.) | 653-siRNA / 52-gene compilation | ✅ 269 new rows / 41 genes | CC BY 2.0; Europe PMC supplementary API |
| **Katoh & Suzuki / Ichihara (i-Score)** | re-normalized recompilations | 🔎 | documented normalization; heavy overlap with above |
| **siRNABERT** (Xu/Zhao 2024, GENE; github.com/ChengkuiZhao/siRNABERT) | code repo; `data/testdata/{katoh,reynold}.csv` (662+253 rows, current tree) + `data/ExperimentalData/{Covid19,PDCD1_8}.csv` (deleted from `main` in commit `d2ad931`, still recoverable from parent commit `f3254ee`) | ✅ PDCD1_8 integrated, rest redundant/unusable | katoh.csv: 662/662 already present. reynold.csv: 240/253 already present, remaining 12 have no gene-identity column. Covid19.csv (recovered from git history): 11 candidate sequences but no efficacy/label column in the committed version. **PDCD1_8.csv (recovered from git history): 8 real siRNAs against `PDCD1`, verified against NM_005018.3, integrated** — +1 gene (97→98). License unresolved (no LICENSE file in source repo); committed as `pdcd1_extra.csv` at the user's explicit request despite this gap — see `data/DATA_SOURCES.md` and `NOTICE.md`. |
| **Matveeva et al. 2007** (3,336-siRNA compilation) | large compilation | ⛔ | only download link dead, no Wayback, 3/4 sources paywalled |
| **siRNA_Design** (Neurotechnology-at-ETH-Zurich, GitHub) | Design-tool scraping + scoring pipeline repo; ships 4 "validation" datasets in `sirna_datasets_analysis/datasets_excel/` | ⛔ nothing new | Dataset A = Huesken 2005 (2,431, already embedded via siRNAEfficacyDB). Dataset B = Ichihara 2007, itself a compilation of Reynolds/Vickers/Haborth/Ui-Tei/Khvorova (419, already embedded). Dataset C = Mysara 2012, a 38-record subset of B (redundant subset). Dataset D = Fellmann 2011 via Mysara (476) — **shRNA, not siRNA**, same disqualifier already applied to DepMap/TRC (see `FUNCTIONAL_GENOMICS_SCREENS.md`). The repo's own in-vitro validation (3 candidates, >90% knockdown) is narrative only, not structured data. |

## 3. Chemically-modified / therapeutic siRNA data (drug-relevant)

| Source | Content | Status | Notes |
|---|---|---|---|
| **Monopoli, Korkin & Khvorova 2023** (MTNA) | 20 siRNAs / 4 CNS genes, sdRNA chemistry | ✅ | CC BY 4.0; cholesterol-conjugated modified sdRNA |
| **Shmushkovich et al. 2018** (NAR) | 356 siRNAs, same chemistry | 🔎 | no gene-identity column → augmentation-only, not gene-groupable |
| **Martinelli 2023 / sirna-repro** (bioRxiv, siRNAmod-derived) | 907 chemically-modified siRNAs, 30 source documents | ✅ 253 of 907 rows | originally blocked by the same disqualifier as Shmushkovich 2018 (no gene-identity column, only sequence + modification + PCT + source PMID/patent-ID). Resolved for 253 rows (7 genes/reporters: EGFP, ACP5, APOB, Luciferase_firefly, Luciferase_renilla, NPY, VEGFA) by tracing each source PMID/patent's stated target and confirming it with exact 19nt substring match against the real transcript, including 28 rows whose stored "sense" field was itself corrupted (fixed by deriving the correct sense from their verified-correct antisense) — see `DATA_SOURCES.md`. The remaining ~654 rows are still blocked; the single largest chunk (544 rows, 60%) is one patent, US20120088815A1, whose own text confirms it tested up to 182 real genes but doesn't name them in the fetched text, and resolving bare 19-mers to genes would need BLAST access unavailable in this environment. Several other newly-identified documents (HIV, Enterovirus 71, coxsackievirus B3) are viral rather than human targets and weren't pursued further. Full raw 907-row CSV kept locally, gitignored (not redistributed) — see `data/raw/sirna_repro_martinelli_907.csv`; the 253-row gene-annotated subset is committed as `data/raw/martinelli_extra.csv` (CC BY-NC 4.0 permits derivatives). **Checked the live upstream directly (crdd.osdd.net/servers/sirnamod/)**: its own bulk CSV export currently returns only 414 rows (all one PMID) against the site's claimed 4,894 — looks broken/stale (old 2016 OSDD/IMTECH server). Not that it matters: the search page's full field list (ID, PMID, per-strand modification name/position/component, activity, assay/cell/transfection) confirms siRNAmod has no gene/target-identity field at the schema level, so a working export would hit the exact same disqualifier. |
| **CMsiRNAdb** (full, 13 genes) | chemically modified, patent-derived | ✅ all 13 genes integrated | see `CMSIRNADB_FULL_RETRIEVAL.md`; on by default (`include_cmsirnadb_full=True`) |
| Alnylam / Arrowhead patents & pubs | modification-vs-activity | 🔎 | scattered in med-chem literature & patent filings |
| **5 FDA-approved siRNA drugs** (via AttSiOff Suppl.) | vutrisiran, givosiran, inclisiran, lumasiran, patisiran | ✅ external validation only | held-out; NOT training data |

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
