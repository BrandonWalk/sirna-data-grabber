# sirna-data-grabber

A standalone siRNA knockdown-efficacy dataset: the raw data files, full
provenance/license documentation, and a small reusable Python package
(`sirna_data`) for loading it -- and, since `pip install sirna-data-grabber`
alone can't ship most of this non-commercial data, a bundled `sirna-data-fetch`
command that re-fetches it from its original sources. Any project that wants
this dataset can depend on this repo (or just the PyPI package) rather than
vendoring a copy of the data or the loading code.

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

```
python3 -m venv .venv && source .venv/bin/activate
pip install -e .
```

This installs the `sirna_data` package in editable mode, so it resolves
`data/raw/` relative to this checkout automatically. If you copy the `data/`
folder somewhere else, point at it explicitly instead:

```
export SIRNA_DATA_DIR=/path/to/data/raw
```

Re-fetching the raw data from scratch (not required if `data/raw/` already
has the files) uses the `sirna-data-fetch` command, installed automatically
with the package -- no extras needed:

```
sirna-data-fetch --dest data/raw
```

This also means a plain `pip install sirna-data-grabber` from PyPI (with no
git checkout at all) can reconstruct the full dataset itself:

```
pip install sirna-data-grabber
sirna-data-fetch --dest ./my_data
export SIRNA_DATA_DIR=./my_data
```

`sirna-data-fetch --only sirna_efficacy monopoli` fetches a subset instead of
all four sources; see `sirna-data-fetch --help`.

## Usage

```python
from sirna_data import load_records, fetch_mrna_by_gene

records = load_records()  # list[SiRNARecord]
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
