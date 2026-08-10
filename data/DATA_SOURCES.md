# Data sources

## siRNA sequences + knockdown efficacy: siRNAEfficacyDB

- Zhang et al. 2024, *IET Systems Biology*, "siRNAEfficacyDB: An experimentally
  supported small interfering RNA efficacy database".
  https://cellknowledge.com.cn/siRNAEfficacy
- 3,544 siRNA records across 42 human genes, compiled from published assays
  (Huesken et al. 2005 Nat. Biotechnol. and others).
- License: Creative Commons Attribution Non-Commercial (CC BY-NC) — free to use,
  distribute, and reproduce for non-commercial purposes with attribution to the
  original authors and underlying studies.
- Fetched by `scripts/download_data.py` into `data/raw/sirna_efficacy.csv`.

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
  `scripts/download_monopoli_data.py`.
- **Important caveat**: these are not standard unmodified siRNA duplexes.
  They use a cholesterol-conjugated, asymmetric (15-nt sense / 20-nt
  antisense) "sdRNA" architecture with heavy 2'-fluoro/2'-O-methyl/
  phosphorothioate modification (from Shmushkovich et al. 2018, Nucleic
  Acids Research, the dataset Monopoli's model was trained on). RNAfold has
  no model of these modifications, so we fold them as if they were plain
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
    see `scripts/download_shabalina_data.py` for the full mapping and the
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
- Fetched by `scripts/download_shabalina_data.py` into
  `data/raw/shabalina_extra.csv` and `data/raw/shabalina_transcripts.fasta`.

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
  as-is.)
- Only the human PCSK9 subset is used by `_load_cmsirnadb_records` --
  **2,756 of 3,107 raw PCSK9 rows kept** -- takes leave-one-gene-out CV
  from 86 to **87 genes**. The other 12 genes are handled separately (see
  `data/CMSIRNADB_FULL_RETRIEVAL.md`).
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
  types across the full database) that RNAfold has no model of, folded as
  plain unmodified RNA -- same approximation already accepted for
  Monopoli2023, arguably on a wider and more heterogeneous set of
  chemistries here since this spans 11 different patents/filers rather
  than one lab's single design. (2) Patent-derived, not assay-paper-derived:
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
