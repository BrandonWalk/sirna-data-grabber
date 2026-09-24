# sirna-data-grabber

[![Tests](https://github.com/BrandonWalk/sirna-data-grabber/actions/workflows/tests.yml/badge.svg)](https://github.com/BrandonWalk/sirna-data-grabber/actions/workflows/tests.yml)
[![PyPI version](https://img.shields.io/pypi/v/sirna-data-grabber.svg)](https://pypi.org/project/sirna-data-grabber/)
[![Python versions](https://img.shields.io/pypi/pyversions/sirna-data-grabber.svg)](https://pypi.org/project/sirna-data-grabber/)
[![License](https://img.shields.io/pypi/l/sirna-data-grabber.svg)](LICENSE)

A standalone siRNA knockdown-efficacy dataset: the raw data files, full
provenance/license documentation, and a small reusable Python package
(`sirna_data`) for loading it -- and, since `pip install sirna-data-grabber`
alone can't ship most of this non-commercial data, a bundled `sirna-data-fetch`
command that re-fetches it from its original sources. Any project that wants
this dataset can depend on this repo (or just the PyPI package) rather than
vendoring a copy of the data or the loading code.

**Currently: 18,077 siRNA records across 116 genes** (`load_records()`
default). Every source is individually toggleable via its own `include_<data-source>`
flag -- see [`data/DATA_SOURCES.md`](data/DATA_SOURCES.md) for the
per-source breakdown and audit.

## Data sources at a glance

| Source | Published | siRNAs | Genes | License |
|---|---|---|---|---|
| [siRNAEfficacyDB](https://cellknowledge.com.cn/siRNAEfficacy) (Zhang et al.) | 2024 | 3,488 | 40 | CC BY-NC |
| [CMsiRNAdb](https://cellknowledge.com.cn/CMsiRNAdb/) (He et al.) | 2026 | 12,357 | 13 | CC BY-NC-ND 4.0[^nd] |
| Shabalina, Spiridonov & Ogurtsov | 2006 | 269 | 41 | CC BY 2.0 |
| Martinelli / sirna-reproduction | 2023 | 577 | 12 | CC BY-NC 4.0 |
| Davis, Monopoli et al. (NAR gkaf479) | 2025 | 966 | 4 | CC BY 4.0 |
| Monopoli, Korkin & Khvorova | 2023 | 20 | 4 | CC BY 4.0 |
| [Sciabola et al. in-house panel](https://doi.org/10.1093/nar/gks1191) (NAR Supp. Tables S3/S4) | 2013 | 356 | 10 | CC BY-NC 3.0 |
| [Harborth et al. lamin A/C panel](https://doi.org/10.1093/nar/gkm699) (via Ichihara et al. 2007) | 2003 | 44 | 1 | CC BY-NC 2.0 UK |
| **Total** | | **18,077** | **116** | |

"Published" is the year of the paper/database each source comes from, not
when it was added here -- see [`data/DATA_SOURCES.md`](data/DATA_SOURCES.md)
for full citations, license terms, and how each source's data was verified.
"Genes" is how many distinct genes/reporters that source contributes to this
dataset; some genes (e.g. `APP`, `MAPT`) are covered by more than one
source, so the per-source counts don't sum to the 116 total. CMsiRNAdb's
count combines its PCSK9 subset and the other-12-genes addition (same
underlying paper) -- see the gene-level table below for the split. Davis
2025 and Monopoli 2023 are from the same lab and cover the same 4 genes
(`APP`/`MAPT`/`BACE1`/`SNCA`) -- Davis 2025's 966 is already net of
deduping against Monopoli 2023 and CMsiRNAdb (full) by exact sequence.

[^nd]: "CC BY-NC-ND" means no derivatives may be redistributed, which is
    why only CMsiRNAdb's untouched original TSV is committed and its
    per-gene subsets are derived at load time -- see
    [`NOTICE.md`](NOTICE.md).

## Genes in this dataset

All 116 genes currently in `load_records()`'s default output, the source
dataset(s) each came from, how many siRNA records target that gene, and the
length of the real mRNA/GenBank transcript its target sites were located in,
and the cell line(s) it was tested in.
Computed directly from the fetched `data/raw/` files, not hand-maintained --
for genes with more than one distinct transcript accession across records
(marked [^multi]), the length shown is for the one used by the most records.

<details>
<summary>Show all 116 genes</summary>

| Gene | Source dataset | siRNAs | Transcript length (nt) | Cell line / test system [^cells] |
|---|---|---|---|---|
| ACP5 | Martinelli 2023 / sirna-reproduction | 32 | 1,683 | not recorded |
| AGT | CMsiRNAdb (full) | 872 | 2,148 [^multi] | Hep3B (638), HepG2 (98), transgenic mice (in vivo) (60), +6 more |
| AKT1 | Shabalina 2006 | 5 | 3,008 | not recorded |
| AKT2 | Shabalina 2006 | 4 | 5,250 | not recorded |
| ALPG | Shabalina 2006 | 11 | 2,492 | not recorded |
| ANGPTL3 | CMsiRNAdb (full) | 551 | 2,926 [^multi] | Hep3B (222), primary cyno hepatocytes (184), transgenic mice (in vivo) (96), +6 more |
| APOB | Martinelli 2023 / sirna-reproduction | 34 | 14,121 | not recorded |
| APP | CMsiRNAdb (full) + Davis 2025 + Monopoli 2023 | 1,244 | 3,358 [^multi] | BE(2)-C (353), SH-SY5Y (284), primary cyno hepatocytes (269), +5 more |
| BACE1 | Davis 2025 + Monopoli 2023 | 183 | 5,835 | SH-SY5Y (180), not recorded (3) |
| BIRC5 | Sciabola 2013 in-house panel (NAR Supp. S3/S4) | 9 | 2,574 | Hep3B |
| BRAF | Sciabola 2013 in-house panel (NAR Supp. S3/S4) | 5 | 6,459 | Hep3B |
| C6orf110 | siRNAEfficacyDB | 145 | 3,465 | HeLa |
| Cacnb1 | siRNAEfficacyDB | 46 | 3,393 | HeLa |
| CASR_rhesus | Martinelli 2023 / sirna-reproduction | 68 | 3,144 | not recorded |
| CBL | Shabalina 2006 | 5 | 11,168 | not recorded |
| CBLB | Shabalina 2006 | 5 | 3,354 | not recorded |
| CDC34 | siRNAEfficacyDB | 57 | 1,418 | HeLa |
| CDKN1A | Shabalina 2006 | 5 | 2,117 | not recorded |
| CDKN1B | Martinelli 2023 / sirna-reproduction | 60 | 2,410 | not recorded |
| CSK | Shabalina 2006 | 5 | 2,743 | not recorded |
| CTNNB1 | CMsiRNAdb (full) + Sciabola 2013 in-house panel (NAR Supp. S3/S4) | 357 | 3,488 [^multi] | Hep3B (350), HeLa (7) |
| Cyclophilin B | siRNAEfficacyDB | 90 | 851 | HEK293 |
| DAD1 | Shabalina 2006 | 5 | 684 | not recorded |
| DBI | siRNAEfficacyDB | 9 | 675 | HEK293 |
| EGFP | Martinelli 2023 / sirna-reproduction | 74 | 1,470 | not recorded |
| EGFP[^egfp2] | siRNAEfficacyDB | 702 | N/A [^egfp] | HeLa |
| EIF4EBP1 | Shabalina 2006 | 4 | 827 | not recorded |
| EZH2 | Sciabola 2013 in-house panel (NAR Supp. S3/S4) | 9 | 2,654 | Hep3B |
| F3_human | Shabalina 2006 | 14 | 2,104 | not recorded |
| F3_mouse | Shabalina 2006 | 10 | 1,821 | not recorded |
| Firefly luciferase | siRNAEfficacyDB | 87 | 2,387 | HEK293 |
| FireflyLuc | siRNAEfficacyDB | 46 | 2,387 [^multi] | avg. of CHO-K1/HeLa/E14TG2a (37), HEK293 (9) |
| FLJ11011 | siRNAEfficacyDB | 78 | 8,412 | HeLa |
| FLJ16071 | Shabalina 2006 | 14 | 2,773 | not recorded |
| FOXO1 | Shabalina 2006 | 5 | 5,779 | not recorded |
| FOXO4 | Shabalina 2006 | 5 | 3,644 | not recorded |
| Fxyd6 | siRNAEfficacyDB | 72 | 1,766 | HeLa |
| FYN | Shabalina 2006 | 5 | 3,628 | not recorded |
| GAPDH | siRNAEfficacyDB | 20 | 1,285 | HEK293 |
| GSK3A | Shabalina 2006 | 5 | 2,193 | not recorded |
| GSK3B | Shabalina 2006 | 5 | 7,782 | not recorded |
| HIF1A | Sciabola 2013 in-house panel (NAR Supp. S3/S4) | 100 | 3,946 | Hep3B |
| HIP2 | siRNAEfficacyDB | 79 | 5,153 | HeLa |
| HK2 | Sciabola 2013 in-house panel (NAR Supp. S3/S4) | 97 | 5,626 | Hep3B |
| HPSE | Sciabola 2013 in-house panel (NAR Supp. S3/S4) | 107 | 4,652 [^multi] | Hep3B |
| HRAS | Shabalina 2006 | 10 | 570 | not recorded |
| HSD17B13 | CMsiRNAdb (full) | 1,985 | 2,260 [^multi] | COS7 (639), primary human hepatocytes (549), primary cyno hepatocytes (476), +7 more |
| HSPC150 | siRNAEfficacyDB | 77 | 878 | HeLa |
| ICAM-1 | siRNAEfficacyDB | 40 | 2,986 | T24 |
| IGF1R | Shabalina 2006 | 21 | 12,235 | not recorded |
| ILK | Shabalina 2006 | 5 | 1,759 | not recorded |
| INHBE | CMsiRNAdb (full) | 670 | 2,460 [^multi] | Hep3B (595), primary cyno hepatocytes (75) |
| IRS1 | Shabalina 2006 | 5 | 9,771 | not recorded |
| ITGB1 | Shabalina 2006 | 5 | 3,735 | not recorded |
| KAZRIN | Martinelli 2023 / sirna-reproduction | 39 | 2,641 | not recorded |
| Lamin A | Harborth 2003 (via Ichihara 2007) | 44 | 3,178 | HeLa |
| LPA | CMsiRNAdb (full) | 556 | 6,431 [^multi] | Huh7 (376), Hep3B (111), transgenic mice (in vivo) (56), +1 more |
| Luciferase_firefly | Martinelli 2023 / sirna-reproduction | 122 | 6,047 [^multi] | not recorded |
| Luciferase_renilla | Martinelli 2023 / sirna-reproduction | 43 | 1,969 | not recorded |
| LYPD1 | Shabalina 2006 | 14 | 3,458 | not recorded |
| MAPK14 | Shabalina 2006 | 8 | 4,222 | not recorded |
| MAPT | CMsiRNAdb (full) + Davis 2025 + Monopoli 2023 | 917 | 6,816 [^multi] | T98G (379), SH-SY5Y (292), BE(2)-C (213), +4 more |
| MARC1 | CMsiRNAdb (full) | 823 | 1,020 [^multi] | Hep3B (341), Huh7 (259), transgenic mice (in vivo) (174), +3 more |
| MIR155HG | Martinelli 2023 / sirna-reproduction | 42 | 1,500 | not recorded |
| MMAC1 | siRNAEfficacyDB | 36 | 3,160 | T24 |
| Mmp7 | siRNAEfficacyDB | 150 | 1,043 | HeLa |
| MSTN | CMsiRNAdb (full) | 9 | 2,705 [^multi] | CD-1 mouse muscle (in vivo) (2), Hepa1-6 (2), chicken embryo myoblasts (2), +2 more |
| MTOR | Sciabola 2013 in-house panel (NAR Supp. S3/S4) | 10 | 8,721 | Hep3B |
| MYC | Sciabola 2013 in-house panel (NAR Supp. S3/S4) + Shabalina 2006 | 14 | 3,721 | Hep3B (9), not recorded (5) |
| MyoD | Shabalina 2006 | 5 | 1,833 | not recorded |
| NOG | siRNAEfficacyDB | 71 | 1,913 | HeLa |
| NPY | Martinelli 2023 / sirna-reproduction | 8 | 567 | not recorded |
| P2rx2 | siRNAEfficacyDB | 77 | 1,833 | HeLa |
| P2RX3 | siRNAEfficacyDB | 90 | 3,792 | HeLa |
| PAC | Shabalina 2006 | 10 | 906 | not recorded |
| PCSK9 | CMsiRNAdb (PCSK9) | 2,756 | 3,637 | HeLa (1,825), HepG2 (379), HEK293A (258), +4 more |
| PDPK1 | Shabalina 2006 | 5 | 7,184 | not recorded |
| PIK3CA | Shabalina 2006 + Sciabola 2013 in-house panel (NAR Supp. S3/S4) | 10 | 9,259 | not recorded (5), Hep3B (5) |
| PIK3R1 | Shabalina 2006 | 5 | 3,371 | not recorded |
| PIK3R2 | Shabalina 2006 | 5 | 3,980 | not recorded |
| PLK | siRNAEfficacyDB | 10 | 2,123 | HEK293 |
| PLN | CMsiRNAdb (full) | 135 | 2,480 | Hepa1-6 |
| PNPLA3 | CMsiRNAdb (full) | 2,066 | 2,753 [^multi] | Hep3B (680), COS7 (624), primary cyno hepatocytes (298), +6 more |
| PSKH1 | Shabalina 2006 | 4 | 3,460 | not recorded |
| RAB13 | Shabalina 2006 | 5 | 1,164 | not recorded |
| RAB6IP1 | siRNAEfficacyDB | 126 | 4,991 | HeLa |
| RB1 | Shabalina 2006 | 5 | 4,768 | not recorded |
| RPS6 | Shabalina 2006 | 5 | 1,369 | not recorded |
| RPS6KA1 | Shabalina 2006 | 5 | 3,192 | not recorded |
| RPS6KA3 | Shabalina 2006 | 5 | 7,987 | not recorded |
| SEAP | siRNAEfficacyDB | 17 | 2,754 [^multi] | HEK293 |
| SEPTIN2 | Shabalina 2006 | 5 | 3,251 | not recorded |
| SKP1 | Shabalina 2006 | 5 | 2,616 | not recorded |
| SNCA | Davis 2025 + Monopoli 2023 | 224 | 3,177 | SH-SY5Y (220), not recorded (4) |
| SOD2_chimp | Martinelli 2023 / sirna-reproduction | 51 | 600 | not recorded |
| SOST | siRNAEfficacyDB | 75 | 2,296 | HeLa |
| TC10 | siRNAEfficacyDB | 67 | 4,780 | HeLa |
| TCAP | siRNAEfficacyDB | 144 | 1,532 | HeLa |
| TSC1 | Shabalina 2006 | 5 | 8,598 | not recorded |
| TSC2 | Shabalina 2006 | 5 | 6,415 | not recorded |
| UBE2B | siRNAEfficacyDB | 79 | 2,241 | HeLa |
| UBE2C | siRNAEfficacyDB | 76 | 777 | HeLa |
| UBE2D3 | siRNAEfficacyDB | 78 | 3,976 | HeLa |
| UBE2E3 | siRNAEfficacyDB | 79 | 1,555 | HeLa |
| UBE2G1 | siRNAEfficacyDB | 79 | 4,167 | HeLa |
| UBE2H | siRNAEfficacyDB | 70 | 5,162 | HeLa |
| UBE2I | siRNAEfficacyDB | 64 | 2,850 | HeLa |
| UBE2J1 | siRNAEfficacyDB | 49 | 4,164 | HeLa |
| UBE2L3 | siRNAEfficacyDB | 53 | 2,861 | HeLa |
| UBE2L6 | siRNAEfficacyDB | 72 | 1,219 | HeLa |
| UBE2M | siRNAEfficacyDB | 76 | 1,159 | HeLa |
| UBE2N | siRNAEfficacyDB | 79 | 4,877 | HeLa |
| UBE2S | siRNAEfficacyDB | 79 | 2,559 | HeLa |
| UBE2V1 | siRNAEfficacyDB | 74 | 2,539 | HeLa |
| Ufc1 | siRNAEfficacyDB | 70 | 888 | HeLa |
| VEGFA | Martinelli 2023 / sirna-reproduction | 4 | 3,660 | not recorded |

</details>

[^multi]: This gene has more than one distinct transcript accession across
its records in the raw data (different isoforms/predicted RefSeq entries
used for different rows) -- the length shown is for the accession used by
the largest number of records, not necessarily all of them.
[^egfp]: This `EGFP` row's 702 rows are mapped in siRNAEfficacyDB to
accession `NZ_CP024869`, which currently resolves to a ~3.7 Mb bacterial
genome assembly, not the actual EGFP transcript -- almost certainly
lab-plasmid contamination baked into that assembly (see "Known data-quality
caveats" in [`data/DATA_SOURCES.md`](data/DATA_SOURCES.md)). All 702 target
sites still verify correctly against a small window of that assembly, so
it's usable for target-site context, but its full length is not a
meaningful "EGFP transcript length" and is omitted here rather than shown as
3,720,309 nt. `Firefly luciferase` and `FireflyLuc` are also two separate
string labels in the source data for what is conceptually the same
reporter, kept distinct here since that's how `load_records()` actually
groups them.
[^egfp2]: This dataset has two textually-distinct `EGFP` gene entries: the
Martinelli row uses a clean `"EGFP"` gene string, while siRNAEfficacyDB's
own `Gene` column has a trailing space (`"EGFP "`) -- a pre-existing
data-entry quirk in that source, not introduced by adding Martinelli. They
group separately here and in `load_records()` because that's how the raw
gene strings actually compare, not silently merged.

[^cells]: The cell line (or in vivo model) each gene's knockdown was measured
in, taken from the raw data wherever the source records it:
siRNAEfficacyDB's `Cell` column and CMsiRNAdb's `Cell_Type` column (names
lightly normalized, e.g. `Hela` -> `HeLa`, `Be(2)C cell line` -> `BE(2)-C`;
"in vivo" marks animal-model rows). Davis 2025 rows are listed as SH-SY5Y,
the cell line that paper's native assay was run in (see
[`data/DATA_SOURCES.md`](data/DATA_SOURCES.md)). Sciabola 2013's rows are
Hep3B (QuantiGene 2.0 mRNA, 48 h) and the Harborth 2003 lamin-A/C panel is
HeLa (lamin A/C protein by immunoblot) -- both from those papers' own
Methods rather than a column in the data. Shabalina 2006, Martinelli 2023 and Monopoli 2023
don't include a cell-line field in the data shipped here and their papers
weren't checked for this column, so their rows show "not recorded" rather
than a value guessed at. Where a gene was tested in more
than one system, the top three are shown with their siRNA counts.

## License

**The code in this repo (`sirna_data`, `tests/`) is MIT licensed** — see
[`LICENSE`](LICENSE). Use it, modify it, ship it commercially, whatever you
want.

**The data in `data/raw/` is NOT covered by that license.** It's redistributed
under each original source's own terms, and most of those sources are
**non-commercial only** (CC BY-NC / CC BY-NC-ND). Loading the data with this
permissively-licensed code does not lift those restrictions — you still have
to comply with them separately. See [`NOTICE.md`](NOTICE.md) for the
per-source summary and [`data/DATA_SOURCES.md`](data/DATA_SOURCES.md) for
full terms before using the data itself, especially commercially.

The same per-source table is available in code (`sirna_data.list_sources()`),
and `load_records(licenses=[...])` loads only the sources carrying the
licenses you name -- see [Loading only the licenses you can
use](#loading-only-the-licenses-you-can-use).

## What's here

```
LICENSE                    MIT license -- covers the code only, not data/raw/
NOTICE.md                  per-source data license summary (see License section above)
data/
  raw/                       fetched CSVs + FASTA transcripts (the actual dataset)
  DATA_SOURCES.md            full provenance, license terms, and trainable-data
                              audit for every source (including CMsiRNAdb full-
                              database retrieval)
  POTENTIAL_DATA_SOURCES.md  landscape of sources investigated
  data_source_ledger.csv     machine-readable companion to DATA_SOURCES.md
src/sirna_data/
  raw_loader.py               load + merge every source into SiRNARecord rows
  genes.py                     list_genes / describe_genes -- what genes are available to load
  licenses.py                  machine-readable per-source data-license table (load_records(licenses=...))
  ncbi_fetch.py                fetch a gene's RefSeq mRNA transcript by symbol
  sequence_utils.py            DNA/RNA sequence helpers (to_rna, to_dna, transcribe_template_to_mrna)
  splitting.py                 train_test_split / leave_n_genes_out dataset splitters
  evaluation.py                evaluate_predictions + PredictionMetrics/GeneCorrelation
  rank_confidence.py           probability/confidence model for "how many top-K predictions to check"
  rank_confidence_cli.py       `sirna-rank-confidence` entry point ([project.scripts])
  rank_confidence_plot.py      optional matplotlib plotting for rank_confidence (requires [plot] extra)
  rank_confidence_plot_cli.py  `sirna-rank-confidence-plot` entry point ([project.scripts], requires [plot] extra)
  __init__.py                  public API
  fetch/                       sirna-data-fetch CLI + per-source fetchers (see Install below)
    cli.py                       `sirna-data-fetch` entry point ([project.scripts])
    sirna_efficacy.py            siRNAEfficacyDB + NCBI -> sirna_efficacy.csv, mrna_transcripts.fasta
    monopoli.py                  Monopoli et al. 2023 supplementary data -> monopoli_*
    shabalina.py                 Shabalina et al. 2006 supplementary data -> shabalina_*
    cmsirnadb.py                 CMsiRNAdb + NCBI -> cmsirnadb_full_raw.tsv, cmsirnadb*_transcripts.fasta
tests/
  test_raw_loader.py          unit tests for raw_loader.py (fixtures, no real data needed)
  test_genes.py               unit tests for genes.py
  test_licenses.py            unit tests for licenses.py (incl. registry/loader sync check)
  test_sciabola2013.py        unit tests for the Sciabola 2013 loader
  test_harborth2003.py        unit tests for the Harborth 2003 loader + its supersede rule
  test_ncbi_fetch.py          unit tests for ncbi_fetch.py (mocked HTTP calls)
  test_fetch_cli.py           unit tests for fetch/cli.py
  test_sequence_utils.py      unit tests for sequence_utils.py
  test_splitting.py           unit tests for splitting.py
  test_evaluation.py          unit tests for evaluation.py
  test_rank_confidence.py     unit tests for rank_confidence.py
  test_rank_confidence_cli.py unit tests for rank_confidence_cli.py
  test_rank_confidence_plot.py unit tests for rank_confidence_plot.py (skipped without the [plot] extra)
  test_rank_confidence_plot_cli.py unit tests for rank_confidence_plot_cli.py (skipped without the [plot] extra)
  conftest.py                 shared pytest fixtures
```

Start with [`data/DATA_SOURCES.md`](data/DATA_SOURCES.md) for what's in the
dataset, where it came from, and the bottom-line audit (7,510 trainable
records across 107 genes, 7 sources — 18,077 records / 116 genes if the
optional CMsiRNAdb full-database retrieval and Davis2025 are also
included). Primary source is **siRNAEfficacyDB** (Zhang et al. 2024, CC
BY-NC); see the docs for the rest and their individual license terms
before reusing this data outside this project.

## Install

`sirna-data-grabber` is [on PyPI](https://pypi.org/project/sirna-data-grabber/),
so most users just need:

```
pip install sirna-data-grabber
```

That installs the `sirna_data` package plus the `sirna-data-fetch` command
(no extras needed). Since the PyPI package can't ship most of this
non-commercial data, use `sirna-data-fetch` to reconstruct it from its
original sources into a local directory:

```
sirna-data-fetch --dest ./my_data
```

Then point `sirna_data` at that directory. Two equivalent ways to do this --
pass it directly, no env var needed:

```python
from sirna_data import load_records
records = load_records(data_dir="./my_data")
```

or export it once as `SIRNA_DATA_DIR` and call `load_records()` with no
arguments:

```
export SIRNA_DATA_DIR=./my_data
```

`sirna-data-fetch --only sirna_efficacy monopoli` fetches a subset instead of
all nine sources; see `sirna-data-fetch --help`. **Every source
`load_records()` reads has a fetcher**, so a bare `pip install` plus one
fetch reconstructs the whole dataset quoted at the top of this file --
nothing here depends on having the git checkout.
`tests/test_fetch_coverage.py` enforces that: adding a loadable source
without a fetcher fails the suite.

### From a git checkout

If you're working from this repo instead (e.g. to browse `data/raw/` and the
provenance docs alongside the code, or to contribute):

```
python3 -m venv .venv && source .venv/bin/activate
pip install -e .
```

This installs `sirna_data` in editable mode, so it resolves `data/raw/`
relative to the checkout automatically -- no `sirna-data-fetch`,
`SIRNA_DATA_DIR`, or `data_dir` needed if `data/raw/` already has the files.
If you copy the `data/` folder somewhere else, point at it with either
`data_dir` or `SIRNA_DATA_DIR` as shown above.

## Usage

```python
from sirna_data import load_records, fetch_mrna_by_gene

records = load_records()  # reads from data_dir / SIRNA_DATA_DIR / default data/raw/, in that order
print(len(records), "records across", len({r.gene for r in records}), "genes")

r = records[0]
r.guide_seq       # siRNA antisense strand
r.mrna_window      # local mRNA context around the real target site
r.label            # experimental %knockdown / %inhibition
r.source           # provenance, e.g. "siRNAEfficacyDB"

# Chemical modification (most records are standard/unmodified; a minority
# -- currently CMsiRNAdb, Monopoli2023, Martinelli_sirna_reproduction, and
# Davis2025 -- are chemically modified):
r.is_modified              # bool
r.modification_chemistry   # short summary, e.g. "2'-OMe/2'-F/PS-backbone (per-position, CMsiRNAdb)"
r.sense_modifications       # per-position modified-nucleoside name or None; CMsiRNAdb only
r.antisense_modifications   # same, for the guide strand

# Look up any gene's RefSeq transcript live from NCBI:
transcript = fetch_mrna_by_gene("TP53")
transcript.accession, transcript.sequence
```

`load_records()` takes `include_sirna_efficacy` / `include_monopoli` /
`include_shabalina` / `include_martinelli` /
`include_harborth2003` / `include_sciabola2013` / `include_cmsirnadb` /
`include_cmsirnadb_full` / `include_davis2025` flags
(all default `True`) to include or exclude any individual source.

`data_dir` (a `Path` or `str`) points every source at a specific directory of
fetched files, as a plain function argument -- no `SIRNA_DATA_DIR` export
required. It falls back to `SIRNA_DATA_DIR` if set, then the package's
default relative `data/raw/` location, in that order.

### Listing the genes available to load

```python
from sirna_data import list_genes, describe_genes

list_genes()          # ['ACP5', 'AGT', 'AKT1', ...] -- all 116, sorted
len(list_genes())     # 116

# same thing with per-gene detail:
for info in describe_genes():
    info.gene         # 'APP'
    info.n_records    # 1244
    info.sources      # ('CMsiRNAdb_full', 'Davis2025', 'Monopoli2023')
    info.licenses     # ('CC BY 4.0', 'CC BY-NC-ND 4.0')  -- mixed: see below
    info.accessions   # ('NM_000484', 'NM_001198823', ...)
```

Both are computed from the actual fetched data files rather than a
hardcoded table, so they can't go stale, and both accept every
`load_records()` argument -- `list_genes(licenses=["CC BY 4.0"])`,
`list_genes(include_cmsirnadb_full=False)`, `list_genes(data_dir="./my_data")`
-- so you can ask what a given subset contains before loading it. Since
that means a full load under the hood, pass records you already have to
skip it: `list_genes(records)`.

Gene strings are the raw sources' own, not normalized -- `"EGFP"` and
`"EGFP "` stay two entries, for the reason in the gene-table footnote above.

### Loading only the licenses you can use

The data's licenses vary per source and most are non-commercial (see
[License](#license) below). That table is also available in code, so
"load only what my project is allowed to use" doesn't mean maintaining your
own list of `include_*` flags:

```python
from sirna_data import load_records, list_licenses, list_sources, license_for_source

list_licenses()
# ['CC BY 2.0', 'CC BY 4.0', 'CC BY-NC', 'CC BY-NC 2.0 UK', 'CC BY-NC 3.0', 'CC BY-NC 4.0', 'CC BY-NC-ND 4.0', 'unresolved']

# only the permissively-licensed sources (Shabalina 2006, Monopoli 2023, Davis 2025):
records = load_records(licenses=["CC BY 4.0", "CC BY 2.0"])
len(records), len({r.gene for r in records})   # 1255, 45

for source in list_sources():
    source.license_id                   # 'CC BY-NC-ND 4.0'
    source.commercial_use               # False  (None == unresolved, NOT "probably fine")
    source.derivatives_redistributable  # False  (the "ND" term)
    source.notes                        # the caveat worth reading

license_for_source(records[0].source).license_id   # per-record provenance -> license
```

License ids are matched case- and punctuation-insensitively (`"cc-by-4.0"`
works) but are version-specific, since two sources here state two different
things (`"CC BY-NC"` vs `"CC BY-NC 4.0"`). An id no source carries raises
`ValueError` listing the valid ones, rather than quietly loading nothing.
`licenses=` intersects with the `include_*` flags, and the two
unresolved-license sources are never selected by a CC license -- only by
asking for `"unresolved"` explicitly.

This registry is a machine-readable copy of [`NOTICE.md`](NOTICE.md)'s
table, not legal advice; read the real terms before relying on it.

### Splitting into train/test

```python
from sirna_data import load_records, train_test_split, leave_n_genes_out

records = load_records()

# sklearn-style train_test_split, but grouped by gene by default so no gene
# straddles both splits (see by_gene below for why this matters):
train, test = train_test_split(records, test_size=0.2, random_state=0)

# leave-N-genes-out cross-validation: a generator yielding one (train, test)
# fold per group of N genes, until every gene has been held out exactly once
for train, test in leave_n_genes_out(records, n=5, random_state=0):
    ...  # train + evaluate a model on this fold
```

`train_test_split` mirrors `sklearn.model_selection.train_test_split`'s name
and parameters (`test_size`, `random_state`) -- the one addition is
`by_gene` (default `True`), which sklearn has no equivalent for. With
`by_gene=True`, every record for a given gene goes entirely into train or
entirely into test, so a model can't partly "solve" a test siRNA just by
having seen another siRNA against the same gene during training. Pass
`by_gene=False` for a plain per-record random split with no regard for gene.

`leave_n_genes_out(records, n, random_state=None)` generalizes leave-one-
gene-out cross-validation: it shuffles the distinct genes once, partitions
them into consecutive groups of `n`, and yields one `(train, test)` fold per
group -- so every gene appears in exactly one test fold across the full
iteration (`n=1` reproduces classic leave-one-gene-out CV). If the gene
count isn't evenly divisible by `n`, the last fold holds out fewer than `n`
genes.

### Rank confidence: how many top predictions do you need to check?

```python
from sirna_data import min_top_k_for_confidence, probability_true_top_in_predicted_top_k

# Given only a correlation between a model's predicted and true rankings of
# 4561 candidate items, how many of the top-predicted items do you need to
# check to be 95% confident the true best one is among them?
min_top_k_for_confidence(n_items=4561, confidence=0.95, pcc=0.3686)

# Or ask it the other way: given you check the top 50, how confident can
# you be that the true best item is in there?
probability_true_top_in_predicted_top_k(50, 4561, pcc=0.3686)
```

Both take the correlation as either `pcc` (Pearson's r, used directly) or
`spcc` (Spearman's rho, converted internally) -- exactly one of the two.
`top_n` (default 1) generalizes the question from "is the single true best
item captured" to "is at least one of the true top `top_n` items captured"
-- pass e.g. `top_n=10` to ask about catching any of the top 10, which
needs a smaller K for the same confidence. See `sirna_data.rank_confidence`'s
module docstring for the full model and its caveats -- this is a planning
heuristic (generally conservative), not a certified statistical bound.

Also installed: the `sirna-rank-confidence` CLI --
`sirna-rank-confidence --pcc 0.3686 --n-items 4561 --confidence 0.99 0.95 0.9`.

#### Comparing multiple models at once

```python
from sirna_data import min_top_k_for_confidence_multi, probability_curves_for_pccs

pccs = [0.2, 0.4, 0.6]  # one Pearson correlation per model to compare

# {pcc: min top-K needed for 95% confidence}, one entry per model
min_top_k_for_confidence_multi(pccs, n_items=4561, confidence=0.95)

# {pcc: [probability at each K in a default spread of K's]}, one entry per model
probability_curves_for_pccs(pccs, n_items=4561)
```

Both run the single-model function above once per PCC in the list --
`min_top_k_for_confidence_multi` for a straight side-by-side "tests needed"
comparison, `probability_curves_for_pccs` for the full probability-vs-K
curve each model traces out (this is what the plotting function below
draws). Pass `k_values` to either the fixed set of K's you want the
comparison at instead of the default spread.

#### Plotting probability vs. number of tests

```python
from sirna_data.rank_confidence_plot import plot_probability_vs_num_tests

plot_probability_vs_num_tests(pccs, n_items=4561, save_path="curves.png")
```

One curve per PCC, x-axis is K (number of top-predicted items checked),
y-axis is the probability of capturing at least one true top-`top_n` item
at that K -- lets you see at a glance how the number of tests needed
relates to each model's correlation. Requires the optional `plot` extra
(`pip install sirna-data-grabber[plot]`) for matplotlib -- not installed by
the core package, and this function lives in its own
`sirna_data.rank_confidence_plot` module (not `sirna_data`'s top-level
import) specifically so nothing else in this package needs matplotlib.
Returns the `matplotlib.axes.Axes` for further customization; pass an
existing `ax=` to draw on it instead of creating a new figure.

Per-point markers, marker size, and line style are all configurable via
`marker` / `markersize` / `linestyle` (each forwarded straight to
`Axes.plot`) -- e.g. `marker=None` for plain lines with no dots, useful
once `k_values` gets dense enough that individual markers just clutter the
curve:

```python
# Dots (default):
plot_probability_vs_num_tests(pccs, n_items=4561, save_path="curves.png")

# Plain lines, no per-point markers:
plot_probability_vs_num_tests(pccs, n_items=4561, marker=None, save_path="curves.png")

# Dashed lines with square markers:
plot_probability_vs_num_tests(
    pccs, n_items=4561, marker="s", linestyle="--", save_path="curves.png"
)
```

Also installed: the `sirna-rank-confidence-plot` CLI, a thin wrapper
around the same function that writes straight to a file --

```
sirna-rank-confidence-plot --pcc 0.2 0.4 0.6 --n-items 4561 \
    --marker none --save-path curves.png
```

`--marker`/`--marker-size`/`--linestyle` mirror the Python function's
`marker`/`markersize`/`linestyle` (pass `--marker none` or `--linestyle
none` for no markers / no connecting line, respectively); `--labels` sets
the legend text per `--pcc` value; `--k-max`/`--num-points`/`--k-values`
control which K's get plotted. Run `sirna-rank-confidence-plot --help` for
the full option list.

## Using this from another project

Install as a sibling checkout in editable mode:

```
pip install -e ../sirna-data-grabber
```

That gives you `import sirna_data` with no other coupling — this repo only
depends on pandas and requests, and knows nothing about any particular
downstream model or feature-engineering pipeline.

## Tests

```
pip install -e ".[test]"
pytest
```

Tests run entirely against small in-memory/tmp-dir fixtures (see
`tests/conftest.py`) and mocked HTTP calls, so they don't touch the real
dataset or the network.

## Linting and type checking

```
pip install -e ".[lint]"
ruff check .
mypy
```

Both run in CI on every pull request (`.github/workflows/tests.yml`), alongside
the test matrix.
