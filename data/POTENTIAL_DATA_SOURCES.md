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
| **Matveeva et al. 2007** (3,336-siRNA compilation) | large compilation | ⛔ | only download link dead, no Wayback, 3/4 sources paywalled |

## 3. Chemically-modified / therapeutic siRNA data (drug-relevant)

| Source | Content | Status | Notes |
|---|---|---|---|
| **Monopoli, Korkin & Khvorova 2023** (MTNA) | 20 siRNAs / 4 CNS genes, sdRNA chemistry | ✅ | CC BY 4.0; cholesterol-conjugated modified sdRNA |
| **Shmushkovich et al. 2018** (NAR) | 356 siRNAs, same chemistry | 🔎 | no gene-identity column → augmentation-only, not gene-groupable |
| **CMsiRNAdb** (full, 13 genes) | chemically modified, patent-derived | 🟡 12 genes not yet pursued | see row in §1 |
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
