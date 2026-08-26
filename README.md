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

**Currently: 16,439 siRNA records across 105 genes** (`load_records()`
default). Every source is individually toggleable via its own `include_*`
flag -- see [`data/DATA_SOURCE_LEDGER.md`](data/DATA_SOURCE_LEDGER.md) for
the per-source breakdown.

## Data sources at a glance

| Source | Published | siRNAs | Genes |
|---|---|---|---|
| [siRNAEfficacyDB](https://cellknowledge.com.cn/siRNAEfficacy) (Zhang et al.) | 2024 | 3,532 | 41 |
| [CMsiRNAdb](https://cellknowledge.com.cn/CMsiRNAdb/) (He et al.) | 2026 | 12,357 | 13 |
| Shabalina, Spiridonov & Ogurtsov | 2006 | 269 | 41 |
| Martinelli / sirna-repro | 2023 | 253 | 7 |
| Monopoli, Korkin & Khvorova | 2023 | 20 | 4 |
| PDCD1 panel (Xu, Zhao et al. / siRNABERT) | 2024 | 8 | 1 |
| **Total** | | **16,439** | **105** |

"Published" is the year of the paper/database each source comes from, not
when it was added here -- see [`data/DATA_SOURCES.md`](data/DATA_SOURCES.md)
for full citations, license terms, and how each source's data was verified.
"Genes" is how many distinct genes/reporters that source contributes to this
dataset; some genes (e.g. `APP`, `MAPT`) are covered by more than one
source, so the per-source counts don't sum to the 105 total. CMsiRNAdb's
count combines its PCSK9 subset and the other-12-genes addition (same
underlying paper) -- see the gene-level table below for the split.

## Genes in this dataset

All 105 genes currently in `load_records()`'s default output, the source
dataset(s) each came from, how many siRNA records target that gene, and the
length of the real mRNA/GenBank transcript its target sites were located in.
Computed directly from the fetched `data/raw/` files, not hand-maintained --
for genes with more than one distinct transcript accession across records
(marked [^multi]), the length shown is for the one used by the most records.

<details>
<summary>Show all 105 genes</summary>

| Gene | Source dataset | siRNAs | Transcript length (nt) |
|---|---|---|---|
| ACP5 | Martinelli 2023 / sirna-repro | 32 | 1,683 |
| AGT | CMsiRNAdb (full) | 872 | 2,148 [^multi] |
| AKT1 | Shabalina 2006 | 5 | 3,008 |
| AKT2 | Shabalina 2006 | 4 | 5,250 |
| ALPG | Shabalina 2006 | 11 | 2,492 |
| ANGPTL3 | CMsiRNAdb (full) | 551 | 2,926 [^multi] |
| APOB | Martinelli 2023 / sirna-repro | 34 | 14,121 |
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
| EGFP | Martinelli 2023 / sirna-repro | 74 | 1,470 |
| EGFP[^egfp2] | siRNAEfficacyDB | 702 | N/A [^egfp] |
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
| INHBE | CMsiRNAdb (full) | 670 | 2,460 [^multi] |
| IRS1 | Shabalina 2006 | 5 | 9,771 |
| ITGB1 | Shabalina 2006 | 5 | 3,735 |
| Lamin A | siRNAEfficacyDB | 44 | 9,756 |
| LPA | CMsiRNAdb (full) | 556 | 6,431 [^multi] |
| Luciferase_firefly | Martinelli 2023 / sirna-repro | 58 | 1,932 |
| Luciferase_renilla | Martinelli 2023 / sirna-repro | 43 | 1,969 |
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
| NPY | Martinelli 2023 / sirna-repro | 8 | 567 |
| P2rx2 | siRNAEfficacyDB | 77 | 1,833 |
| P2RX3 | siRNAEfficacyDB | 90 | 3,792 |
| PAC | Shabalina 2006 | 10 | 906 |
| PCSK9 | CMsiRNAdb (PCSK9) | 2,756 | 3,637 |
| PDCD1 | siRNABERT PDCD1 panel (Xu/Zhao 2024) | 8 | 2,097 |
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
| VEGFA | Martinelli 2023 / sirna-repro | 4 | 3,660 |

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
dataset and where it came from; [`data/DATA_SOURCE_LEDGER.md`](data/DATA_SOURCE_LEDGER.md)
for the bottom-line audit (6,838 trainable records across 95 genes, 6
sources — 16,439 records / 105 genes if the optional CMsiRNAdb full-database
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

# Chemical modification (most records are standard/unmodified; a minority
# -- currently CMsiRNAdb, Monopoli2023, and Martinelli_sirna_repro -- are
# chemically modified):
r.is_modified              # bool
r.modification_chemistry   # short summary, e.g. "2'-OMe/2'-F/PS-backbone (per-position, CMsiRNAdb)"
r.sense_modifications       # per-position modified-nucleoside name or None; CMsiRNAdb only
r.antisense_modifications   # same, for the guide strand

# Look up any gene's RefSeq transcript live from NCBI:
transcript = fetch_mrna_by_gene("TP53")
transcript.accession, transcript.sequence
```

`load_records()` takes `include_sirna_efficacy` / `include_monopoli` /
`include_pdcd1` / `include_shabalina` / `include_martinelli` /
`include_cmsirnadb` / `include_cmsirnadb_full` flags (all default `True`) to
include or exclude any individual source, including the primary
siRNAEfficacyDB set -- no source is loaded unconditionally.

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
