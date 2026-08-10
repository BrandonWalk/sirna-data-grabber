# siRecords vs. this dataset — new-data overlap analysis

_Generated during data-source review. siRecords 04/28/05 release retrieved from the
Internet Archive (the live siRecords servers at `sirecords.umn.edu`,
`c1.accurascience.com`, and `sirecords.biolead.org` are all defunct — HTTP 502)._

_Snapshot note: this analysis was run against the dataset as it existed at the time,
which still included a small user-supplied data source (91 records) since removed.
The overlap results and per-source breakdown below are as computed then; the
"currently spans 127 genes / 128 accessions" figure below is stale relative to the
current dataset (see `DATA_SOURCE_LEDGER.md` for current totals) but the dedup
methodology and findings are otherwise unaffected._

## License status (researched after the fact — do not treat as "open")

siRecords' own definitive writeup, Ren et al. 2009, *Nucleic Acids Research* 37
(Database issue) D146-D149, "siRecords: a database of mammalian RNAi experiments
and efficacies" (doi:10.1093/nar/gkn817), is itself published under CC BY-NC 2.0 UK
— but that license covers the *article* (text/figures), not a grant to redistribute
the underlying bulk dataset. The paper's own **DATA ACCESS** section states the
actual terms for the data itself:

> "The siRecords web site is publicly accessible through the URL
> http://siRecords.umn.edu/siRecords. **Academic users can obtain a copy of the
> current release of the dataset by sending an email** to [the corresponding
> author]."

That is a controlled, individual-request distribution model restricted to
academic users — not a blanket open-data license. Separately, the database's
later mirror host, AccuraScience (`c1.accurascience.com/siRecords/`), publishes a
general Terms of Use for its site stating downloaded content is for "personal
non-commercial use" only and may not be "copie[d], broadcast, download[ed],
store[d] ... transmit[ted] ... for any other purpose whatsoever without the prior
written permission" of AccuraScience.

**Neither of those covers how this repo obtained the data.** `sirecords_efficacy.csv`
here was recovered from an Internet Archive snapshot of the live site, not through
the authors' sanctioned academic-request channel — so even the narrow "academic
users, on request" permission the original paper describes doesn't technically
apply to this copy. Net finding: siRecords' data was never established as freely
redistributable, and there is no license (CC or otherwise) that clearly covers
bulk redistribution of it as done here. Treat `sirecords_efficacy.csv` and
`sirecords_new_only.csv` as **unresolved license risk**, not merely "unverified" —
see `NOTICE.md`.

## Headline

| Level | New | Already in this dataset | siRecords total (usable) |
|---|---|---|---|
| **siRNA records** (by exact sequence) | **3295 (84.8%)** | 592 (15.2%) | 3,887 |
| **Unique sequences** | 2,901 (84.1%) | 550 (15.9%) | 3,451 |
| **Target-gene accessions** | 1,534 | 49 | 1,583 |

Of the 4,162 siRecords rows, 275 have no usable sequence (too short / blank) and are excluded from the record-level comparison.

## How overlap was determined

Matching replicates this dataset's own dedup rule: **exact nucleotide-sequence identity**, strand-agnostic.
- Every existing guide/target sequence (from `sirna_efficacy.csv` antisense+sense columns, plus `user_provided`, `shabalina`, `monopoli`, `cmsirnadb` extras) was normalised (U→T, uppercased) and indexed both as-is and as reverse-complement, down to overlapping 19-mers.
- A siRecords sequence counts as "already in this dataset" if it (or its reverse complement) shares a 19-mer with any existing sequence. Sequences shorter than 19 nt were tested as exact substrings.

The 592 overlapping records trace to the existing sources exactly as expected from `DATA_SOURCES.md`:

| Traces to existing source | overlapping records |
|---|---|
| siRNAEfficacyDB (primary) | 446 |
| Shabalina 2006 | 145 |
| CMsiRNAdb | 1 |

This is the anticipated result: siRNAEfficacyDB and Shabalina2006 are themselves compilations of the same classic mid-2000s assays (Huesken 2005, Reynolds, Khvorova, Vickers, Hsieh) that siRecords aggregates, so the shared core overlaps.

## What is genuinely new — and the critical caveat

**3295 records (3117 carrying an efficacy rating), spanning ~1,400 new target-gene accessions, are not currently in this dataset by sequence.** On its face this is a large potential expansion of gene coverage (this dataset currently spans 127 genes / 128 accessions; siRecords references 1,583 accessions).

**However, this "new" data is not directly train-ready, for one decisive reason:**

- **siRecords efficacy is a 4-level ORDINAL rating (Very high / High / Medium / Low), not a numeric `%Inhibition`.** Every existing source in this dataset stores a continuous knockdown/inhibition percentage (a regression/classification-ready target). siRecords does **not** provide per-record percentages — only the coarse bin. The new-only rated breakdown is: Very high 1,145 · High 1,003 · Medium 485 · Low 484.

To use the new siRecords records you would have to either (a) train/evaluate on the ordinal label directly (a different target from the numeric percentage this dataset otherwise provides), or (b) go back to the **PubMed IDs** siRecords provides (present for essentially every row) and extract the reported numeric knockdown from the primary papers — the same provenance-tracing approach `DATA_SOURCES.md` already documents for other sources.

Secondary caveats:
- siRecords predates and overlaps the existing academic sources, so the 15% that overlaps is redundant and should be dropped on integration (the new-only file below already excludes it).
- Cell line / assay / concentration metadata is present but formatted differently from this dataset's schema and would need mapping to the `Technology` one-hot buckets used by downstream feature engineering.
- Sequence lengths are heterogeneous (12–64 nt; median 19); the rest of this dataset assumes ~19–21-mers.

## Files

- `data/raw/sirecords_efficacy.csv` — full 4,162-record siRecords release (all fields).
- `data/raw/sirecords_new_only.csv` — the 3295 records whose sequence is **not** already in this dataset (3117 with an efficacy rating) — the candidate extension set, overlap already removed.
- `data/POTENTIAL_DATA_SOURCES.md` — catalogue of candidate siRNA-efficacy data sources (siRecords and others).
