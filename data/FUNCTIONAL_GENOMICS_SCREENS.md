# Functional-genomics RNAi screens — investigated, not included

Follow-up to the data-source ledger: actually pulling the large-scale screen data
(DepMap/DRIVE/Achilles, TRC/pLKO, GenomeRNAi) rather than just noting them as proxy sources.

## DepMap DEMETER2 (Achilles + DRIVE + Marcotte) — retrieved, evaluated, removed

Source: DepMap `DEMETER2 data (v6)`, figshare article 6025238 (CC BY 4.0).
Combines the Broad **Achilles** and Novartis **Project DRIVE** genome-scale shRNA
screens plus Marcotte 2016, all reprocessed through the DEMETER2 model
(McFarland et al. 2018, *Nat. Commun.*).

**224,722 shRNA hairpins**, mapped to **17,527 genes**, was pulled into
`data/raw/depmap_demeter2_shRNA_efficacy.csv` and briefly had a loader
(`load_demeter2_records()`) that returned a distinct `ShRNARecord` type,
explicitly kept separate from `load_records()`/`SiRNARecord`.

**Both the data and the loader were removed.** This is shRNA, not siRNA, and
its `Geff` column is a *modeled* latent efficacy score (inferred from
multi-cell-line dropout, scaled 0-1) -- not a directly-measured wet-lab
knockdown percentage. Every other record in this repo carries an actual
experimental %knockdown/%inhibition measurement; this source has no analog to
that at all, only a model's estimate of one. On review, keeping data with no
real measured-knockdown value in a project whose whole schema
(`SiRNARecord.label`) promises exactly that isn't worth the confusion risk,
even behind a separate type and a well-documented boundary -- see the
trainability reasons below for the full case. The 26MB of source CSVs
were deleted from `data/raw/`; `data/DEMETER2_README.txt` (the upstream
release notes) is kept for reference since it's tiny and carries no data.

### TRC / pLKO shRNA library — same disposition

The Broad GPP portal (`portals.broadinstitute.org/gpp`) serves the TRC library only
through an interactive JavaScript clone-search with no bulk download. This was moot
either way: the Achilles screens **are** TRC-library screens, so the TRC/pLKO hairpin
sequences were already present as the 228,003 unique barcodes in DEMETER2's
`shRNA-mapping.csv` (`depmap_shRNA_gene_map.csv`) -- also removed along with the
efficacy file, for the same reason.

### GenomeRNAi — never retrieved (broken backend) + not trainable anyway

Two access paths tried:
- `www.genomernai.org` — HTTP blocked at the sandbox allowlist; HTTPS returns a TLS
  `TLSV1_UNRECOGNIZED_NAME` error (the server rejects its own hostname in the SNI
  handshake — a server-side misconfiguration no client can work around).
- **`genomernai.dkfz.de`** — the DKFZ mirror IS alive and serves pages, but its bulk
  "Download all experiments" endpoint returns **HTTP 500**, and the two static download
  links (`frequentHitters`, `ArchivedScreens`) return **404**. The download layer is
  non-functional.

**More importantly, GenomeRNAi would not be trainable for this dataset's purposes even if downloaded.**
It is a *phenotype* database: its records are screen-level phenotype calls (hit / no-hit
and free-text phenotype descriptions) indexed by gene, aggregated across heterogeneous
cell-based and in-vivo RNAi screens. It does not store per-reagent numeric knockdown
percentages, and many of its reagents are Drosophila/C. elegans dsRNAs or undisclosed
siRNA pools without a clean sequence->%KD pairing. It is a phenotype-discovery resource,
not a knockdown-efficacy training source. See contacts below if the phenotype layer is
wanted for a different purpose.

## Why DEMETER2 was ultimately excluded

This is **proxy efficacy, not the same target as this dataset's numeric %KD** --
and unlike every other source in this repo, there is no real measured knockdown
value anywhere in it:

1. **shRNA ≠ siRNA.** Hairpins are transcribed and Dicer-processed; potency folds in
   expression level and processing efficiency, so Geff is not interchangeable with a
   transfected-siRNA reporter knockdown.
2. **Geff is a modeled latent variable**, inferred from multi-cell-line dropout, not a
   wet-lab percentage. Scale (0–1) is not a knockdown fraction, and no wet-lab
   percentage exists anywhere in the source data to fall back on.
3. **Guide vs full hairpin.** The 21-nt barcode is the hairpin stem; deriving the mature
   guide requires the pLKO loop/processing model.

If a future need arises for a large auxiliary/pretraining corpus or an
independent gene-coverage benchmark, this is re-fetchable from the figshare
link above (CC BY 4.0, no access restriction) -- nothing here prevents
revisiting it, but it isn't part of this repo's data or package by default.

## Files
- `data/DEMETER2_README.txt` — upstream release notes (kept for reference; the CSVs themselves were removed)

## Contacts (for GenomeRNAi and for anything behind a portal)

- **GenomeRNAi** — Boutros lab, German Cancer Research Center (DKFZ) / Heidelberg.
  Michael Boutros, `m.boutros@dkfz.de` (lab: `http://www.dkfz.de/signaling`). GenomeRNAi
  was developed by Thomas Horn / Esther Schmidt / Florian Wilsch in that group; the DKFZ
  signaling-and-functional-genomics group is the maintainer of record.
- **Broad DepMap / GPP (Achilles, TRC/pLKO)** — Cancer Data Science, Broad Institute.
  General: `depmap@broadinstitute.org`. GPP/TRC reagent questions: `genetic-perturbation@broadinstitute.org`.
- **Project DRIVE** — originated at Novartis Institutes for BioMedical Research (McDonald
  et al. 2017); the reprocessed data is distributed through DepMap, so route requests via
  `depmap@broadinstitute.org`.

_Contact addresses are institutional/lab addresses from the projects' public pages and
publications; verify current validity before relying on them._
