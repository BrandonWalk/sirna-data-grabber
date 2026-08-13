# sirna-data-grabber

A standalone siRNA knockdown-efficacy dataset: the raw data files, full
provenance/license documentation, and a small reusable Python package
(`sirna_data`) for loading it -- and, since `pip install sirna-data-grabber`
alone can't ship most of this non-commercial data, a bundled `sirna-data-fetch`
command that re-fetches it from its original sources. Any project that wants
this dataset can depend on this repo (or just the PyPI package) rather than
vendoring a copy of the data or the loading code.

**Currently: 16,178 siRNA records across 97 genes** (`load_records()`
default). Every source is individually toggleable via its own `include_*`
flag -- see [`data/DATA_SOURCE_LEDGER.md`](data/DATA_SOURCE_LEDGER.md) for
the per-source breakdown.

## Genes in this dataset

All 97 genes currently in `load_records()`'s default output, the source
dataset(s) each came from, how many siRNA records target that gene, and the
length of the real mRNA/GenBank transcript its target sites were located in.
Computed directly from the fetched `data/raw/` files, not hand-maintained --
for genes with more than one distinct transcript accession across records
(marked [^multi]), the length shown is for the one used by the most records.

<details>
<summary>Show all 97 genes</summary>

| Gene | Source dataset | siRNAs | Transcript length (nt) |
|---|---|---|---|
| AGT | CMsiRNAdb (full) | 872 | 2,148 [^multi] |
| AKT1 | Shabalina 2006 | 5 | 3,008 |
| AKT2 | Shabalina 2006 | 4 | 5,250 |
| ALPG | Shabalina 2006 | 11 | 2,492 |
| ANGPTL3 | CMsiRNAdb (full) | 551 | 2,926 [^multi] |
| APP | CMsiRNAdb (full) + Monopoli 2023 | 960 | 3,358 [^multi] |
| BACE1 | Monopoli 2023 | 3 | 5,835 |
| C6orf110 | siRNAEfficacyDB | 145 | 3,465 |
| Cacnb1 | siRNAEfficacyDB | 46 | 3,393 |
| CBL | Shabalina 2006 | 5 | 11,168 |
| CBLB | Shabalina 2006 | 5 | 3,354 |
| CDC34 | siRNAEfficacyDB | 57 | 1,418 |
| CDKN1A | Shabalina 2006 | 5 | 2,117 |
| CSK | Shabalina 2006 | 5 | 2,743 |
| CTNNB1 | CMsiRNAdb (full) | 352 | 3,488 |
| Cyclophilin B | siRNAEfficacyDB | 90 | 851 |
| DAD1 | Shabalina 2006 | 5 | 684 |
| DBI | siRNAEfficacyDB | 9 | 675 |
| EGFP | siRNAEfficacyDB | 702 | N/A [^egfp] |
| EIF4EBP1 | Shabalina 2006 | 4 | 827 |
| F3_human | Shabalina 2006 | 14 | 2,104 |
| F3_mouse | Shabalina 2006 | 10 | 1,821 |
| Firefly luciferase | siRNAEfficacyDB | 87 | 2,387 |
| FireflyLuc | siRNAEfficacyDB | 46 | 2,387 [^multi] |
| FLJ11011 | siRNAEfficacyDB | 78 | 8,412 |
| FLJ16071 | Shabalina 2006 | 14 | 2,773 |
| FOXO1 | Shabalina 2006 | 5 | 5,779 |
| FOXO4 | Shabalina 2006 | 5 | 3,644 |
| Fxyd6 | siRNAEfficacyDB | 72 | 1,766 |
| FYN | Shabalina 2006 | 5 | 3,628 |
| GAPDH | siRNAEfficacyDB | 20 | 1,285 |
| GSK3A | Shabalina 2006 | 5 | 2,193 |
| GSK3B | Shabalina 2006 | 5 | 7,782 |
| HIP2 | siRNAEfficacyDB | 79 | 5,153 |
| HRAS | Shabalina 2006 | 10 | 570 |
| HSD17B13 | CMsiRNAdb (full) | 1,985 | 2,260 [^multi] |
| HSPC150 | siRNAEfficacyDB | 77 | 878 |
| ICAM-1 | siRNAEfficacyDB | 40 | 2,986 |
| IGF1R | Shabalina 2006 | 21 | 12,235 |
| ILK | Shabalina 2006 | 5 | 1,759 |
| INHBE | CMsiRNAdb (full) | 670 | 2,460 |
| IRS1 | Shabalina 2006 | 5 | 9,771 |
| ITGB1 | Shabalina 2006 | 5 | 3,735 |
| Lamin A | siRNAEfficacyDB | 44 | 9,756 |
| LPA | CMsiRNAdb (full) | 556 | 6,431 |
| LYPD1 | Shabalina 2006 | 14 | 3,458 |
| MAPK14 | Shabalina 2006 | 8 | 4,222 |
| MAPT | CMsiRNAdb (full) + Monopoli 2023 | 635 | 6,816 [^multi] |
| MARC1 | CMsiRNAdb (full) | 823 | 1,020 [^multi] |
| MMAC1 | siRNAEfficacyDB | 36 | 3,160 |
| Mmp7 | siRNAEfficacyDB | 150 | 1,043 |
| MSTN | CMsiRNAdb (full) | 9 | 2,705 [^multi] |
| MYC | Shabalina 2006 | 5 | 3,721 |
| MyoD | Shabalina 2006 | 5 | 1,833 |
| NOG | siRNAEfficacyDB | 71 | 1,913 |
| P2rx2 | siRNAEfficacyDB | 77 | 1,833 |
| P2RX3 | siRNAEfficacyDB | 90 | 3,792 |
| PAC | Shabalina 2006 | 10 | 906 |
| PCSK9 | CMsiRNAdb (PCSK9) | 2,756 | 3,637 |
| PDPK1 | Shabalina 2006 | 5 | 7,184 |
| PIK3CA | Shabalina 2006 | 5 | 9,259 |
| PIK3R1 | Shabalina 2006 | 5 | 3,371 |
| PIK3R2 | Shabalina 2006 | 5 | 3,980 |
| PLK | siRNAEfficacyDB | 10 | 2,123 |
| PLN | CMsiRNAdb (full) | 135 | 2,480 |
| PNPLA3 | CMsiRNAdb (full) | 2,066 | 2,753 [^multi] |
| PSKH1 | Shabalina 2006 | 4 | 3,460 |
| RAB13 | Shabalina 2006 | 5 | 1,164 |
| RAB6IP1 | siRNAEfficacyDB | 126 | 4,991 |
| RB1 | Shabalina 2006 | 5 | 4,768 |
| RPS6 | Shabalina 2006 | 5 | 1,369 |
| RPS6KA1 | Shabalina 2006 | 5 | 3,192 |
| RPS6KA3 | Shabalina 2006 | 5 | 7,987 |
| SEAP | siRNAEfficacyDB | 17 | 2,754 [^multi] |
| SEPTIN2 | Shabalina 2006 | 5 | 3,251 |
| SKP1 | Shabalina 2006 | 5 | 2,616 |
| SNCA | Monopoli 2023 | 4 | 3,177 |
| SOST | siRNAEfficacyDB | 75 | 2,296 |
| TC10 | siRNAEfficacyDB | 67 | 4,780 |
| TCAP | siRNAEfficacyDB | 144 | 1,532 |
| TSC1 | Shabalina 2006 | 5 | 8,598 |
| TSC2 | Shabalina 2006 | 5 | 6,415 |
| UBE2B | siRNAEfficacyDB | 79 | 2,241 |
| UBE2C | siRNAEfficacyDB | 76 | 777 |
| UBE2D3 | siRNAEfficacyDB | 78 | 3,976 |
| UBE2E3 | siRNAEfficacyDB | 79 | 1,555 |
| UBE2G1 | siRNAEfficacyDB | 79 | 4,167 |
| UBE2H | siRNAEfficacyDB | 70 | 5,162 |
| UBE2I | siRNAEfficacyDB | 64 | 2,850 |
| UBE2J1 | siRNAEfficacyDB | 49 | 4,164 |
| UBE2L3 | siRNAEfficacyDB | 53 | 2,861 |
| UBE2L6 | siRNAEfficacyDB | 72 | 1,219 |
| UBE2M | siRNAEfficacyDB | 76 | 1,159 |
| UBE2N | siRNAEfficacyDB | 79 | 4,877 |
| UBE2S | siRNAEfficacyDB | 79 | 2,559 |
| UBE2V1 | siRNAEfficacyDB | 74 | 2,539 |
| Ufc1 | siRNAEfficacyDB | 70 | 888 |

</details>

[^multi]: This gene has more than one distinct transcript accession across
its records in the raw data (different isoforms/predicted RefSeq entries
used for different rows) -- the length shown is for the accession used by
the largest number of records, not necessarily all of them.
[^egfp]: `EGFP`'s 702 rows are mapped in siRNAEfficacyDB to accession
`NZ_CP024869`, which currently resolves to a ~3.7 Mb bacterial genome
assembly, not the actual EGFP transcript -- almost certainly lab-plasmid
contamination baked into that assembly (see "Known data-quality caveats" in
[`data/DATA_SOURCES.md`](data/DATA_SOURCES.md)). All 702 target sites still
verify correctly against a small window of that assembly, so it's usable for
target-site context, but its full length is not a meaningful "EGFP
transcript length" and is omitted here rather than shown as 3,720,309 nt.
`Firefly luciferase` and `FireflyLuc` are also two separate string labels in
the source data for what is conceptually the same reporter, kept distinct
here since that's how `load_records()` actually groups them.

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

## What's here

```
LICENSE                    MIT license -- covers the code only, not data/raw/
NOTICE.md                  per-source data license summary (see License section above)
data/
  raw/                       fetched CSVs + FASTA transcripts (the actual dataset)
  DATA_SOURCES.md            full provenance + license terms for every source
  DATA_SOURCE_LEDGER.md      audit: what's trainable, what's not, and why
  CMSIRNADB_FULL_RETRIEVAL.md   detail on the CMsiRNAdb full-database retrieval
  DEMETER2_README.txt        upstream release notes for DepMap DEMETER2 (investigated, not included -- see FUNCTIONAL_GENOMICS_SCREENS.md)
  FUNCTIONAL_GENOMICS_SCREENS.md   notes on functional-genomics screen sources considered
  POTENTIAL_DATA_SOURCES.md  landscape of sources investigated
  sirecords_overlap_analysis.md    siRecords overlap/dedup analysis
  data_source_ledger.csv     machine-readable companion to DATA_SOURCE_LEDGER.md
  *.png                      figures referenced by the docs above
src/sirna_data/
  raw_loader.py               load + merge every source into SiRNARecord rows
  ncbi_fetch.py                fetch a gene's RefSeq mRNA transcript by symbol
  __init__.py                  public API
  fetch/                       sirna-data-fetch CLI + per-source fetchers (see Install below)
    cli.py                       `sirna-data-fetch` entry point ([project.scripts])
    sirna_efficacy.py            siRNAEfficacyDB + NCBI -> sirna_efficacy.csv, mrna_transcripts.fasta
    monopoli.py                  Monopoli et al. 2023 supplementary data -> monopoli_*
    shabalina.py                 Shabalina et al. 2006 supplementary data -> shabalina_*
    cmsirnadb.py                 CMsiRNAdb + NCBI -> cmsirnadb_full_raw.tsv, cmsirnadb*_transcripts.fasta
tests/
  test_raw_loader.py          unit tests for raw_loader.py (fixtures, no real data needed)
  test_ncbi_fetch.py          unit tests for ncbi_fetch.py (mocked HTTP calls)
  conftest.py                 shared pytest fixtures
```

Start with [`data/DATA_SOURCES.md`](data/DATA_SOURCES.md) for what's in the
dataset and where it came from; [`data/DATA_SOURCE_LEDGER.md`](data/DATA_SOURCE_LEDGER.md)
for the bottom-line audit (6,577 trainable records across 87 genes, 4
sources — 16,178 records / 97 genes if the optional CMsiRNAdb full-database
retrieval is also included). Primary source is **siRNAEfficacyDB** (Zhang
et al. 2024, CC BY-NC); see the docs for the rest and their individual
license terms before reusing this data outside this project.

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
all four sources; see `sirna-data-fetch --help`.

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

# Look up any gene's RefSeq transcript live from NCBI:
transcript = fetch_mrna_by_gene("TP53")
transcript.accession, transcript.sequence
```

`load_records()` takes `include_sirna_efficacy` / `include_monopoli` /
`include_shabalina` / `include_cmsirnadb` / `include_cmsirnadb_full` flags
(all default `True`) to include or exclude any individual source, including
the primary siRNAEfficacyDB set -- no source is loaded unconditionally.

`data_dir` (a `Path` or `str`) points every source at a specific directory of
fetched files, as a plain function argument -- no `SIRNA_DATA_DIR` export
required. It falls back to `SIRNA_DATA_DIR` if set, then the package's
default relative `data/raw/` location, in that order.

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
