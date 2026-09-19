# Contributing to sirna-data-grabber

Thanks for your interest. This repo is two things at once -- a small MIT
Python package (`sirna_data`) and a **provenance-audited dataset** assembled
from eight separately-licensed sources -- and the second half is what makes
contributing here a little different from a normal Python project. Most of
the rules below exist because a wrong step on the data side is a licensing
problem, not just a bug.

Most useful contributions, roughly in order:

1. **A new data source.** The biggest lever. See [Adding a data
   source](#adding-a-data-source) -- it is a checklist, not a one-file
   change.
2. **Data-quality findings.** A mis-mapped accession, a mislabeled gene, a
   duplicate that slipped past dedup. See [Reporting a data
   problem](#reporting-a-data-problem).
3. **Resolving an unresolved license.** One source (the REMOVED panel)
   ships no `LICENSE` file, so it is all-rights-reserved by default and its
   derivative file is kept out of git. Another route works too: when the
   same measurements exist in an open-access paper, sourcing them from there
   beats chasing permission. If you can get an explicit license statement from those
   authors, that is a genuinely valuable contribution -- open an issue
   with the correspondence.
4. **Loader, splitting, and docs improvements.**

## Development setup

```bash
git clone https://github.com/BrandonWalk/sirna-data-grabber.git
cd sirna-data-grabber
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[test,lint,plot]"
```

Editable installs resolve `data/raw/` relative to the checkout, so if the
raw files are present, `load_records()` works with no `data_dir` or
`SIRNA_DATA_DIR`. `data/raw/` is not fully committed (see
[Licensing](#licensing-the-part-you-cannot-skip)) -- reconstruct it from the
original sources with:

```bash
sirna-data-fetch --dest data/raw
```

Requires network access to `cellknowledge.com.cn`, Europe PMC, GitHub raw
and NCBI E-utilities. `sirna-data-fetch` covers every loadable source, so a
bare `pip install` reconstructs the full dataset;
`tests/test_fetch_coverage.py` fails if a source is added without a
fetcher, which is the check that keeps the README's totals honest.

## Running the checks

CI runs exactly three things; run them locally before opening a PR:

```bash
ruff check .        # lint (E, F, I, UP, B; line-length 100)
mypy                # type check -- config pins python_version = 3.12
pytest              # unit tests
```

The test matrix covers Python 3.10, 3.11, 3.12 and 3.13, so keep syntax and
typing compatible with 3.10 (the `requires-python` floor) even though mypy
itself runs under 3.12 -- see the comment in `pyproject.toml` for why those
differ.

**Tests must not touch the network and must not require `data/raw/`.**
Every existing test builds its input from fixtures in `tests/conftest.py` or
mocks the HTTP call (`tests/test_ncbi_fetch.py` is the model). A test that
silently passes because a maintainer happens to have the raw files locally
is worse than no test. Tests requiring the `[plot]` extra should skip
cleanly without it, as `test_rank_confidence_plot.py` does.

## Licensing: the part you cannot skip

The code is MIT. **The data is not**, and the sources do not share terms:
some are CC BY, several are non-commercial (CC BY-NC), one is
no-derivatives (CC BY-NC-ND 4.0), and two are unresolved. Three rules
follow from that, and PRs that break them cannot be merged regardless of
how good the data is:

1. **Never commit data whose license does not permit redistribution.**
   Unresolved-license sources are kept locally, excluded via `.gitignore`,
   and documented in [`NOTICE.md`](NOTICE.md). If you are not sure a source
   may be redistributed, the answer is no until it is established in
   writing.
2. **No-derivatives sources are filtered at load time, never shipped
   derived.** CMsiRNAdb's untouched original TSV is what is committed; its
   per-gene subsets are computed inside `raw_loader.py` when the user calls
   `load_records()`. Do not add a pre-derived CSV for an ND source.
3. **Every new source needs its license researched and recorded before its
   code is written.** "The paper is open access" is not a data license, and
   a repo with no `LICENSE` file is all-rights-reserved by default under
   GitHub's terms -- not public domain.

If a license question is genuinely ambiguous, say so explicitly in the docs
rather than picking the convenient reading. `SourceLicense` has a `None`
state for `commercial_use` / `derivatives_redistributable` precisely so
"unknown" can be represented honestly, and the docstring's warning applies
to contributors too: do not read `None` as "probably fine."

## Adding a data source

Work through this in order. Steps 1-2 come before any code.

**1. Research the license.** Find the actual terms covering the *data*
(supplementary files, database download, repo contents) -- not just the
article. Record the URL and the exact stated license id.

**2. Open an issue** with the source, its size, its genes, the efficacy
metric it reports, and what you found in step 1. This is where "is this
worth integrating and may we redistribute it" gets settled.

**3. Write a fetcher** at `src/sirna_data/fetch/<source>.py` exposing
`fetch(dest: Path) -> None`, which downloads from the original location and
writes cleaned file(s) into `dest`. Follow `fetch/sirna_efficacy.py`: stdlib
plus pandas only (no new dependencies -- the fetch CLI is installed with the
base package), print progress, rate-limit polite (NCBI allows 3 req/s
without an API key), and write only files that are meant to persist.

**4. Register it** in `fetch/cli.py`'s `SOURCES` dict and mention it in
`fetch/__init__.py`'s docstring.

**5. Write the loader** as `_load_<source>_records(flank_nt, data_dir)` in
`raw_loader.py`, returning `list[SiRNARecord]`. Things to get right:

- `guide_seq` is the **antisense** strand, RNA alphabet. If the source
  publishes sense strands, reverse-complement them (`_revcomp`) and say so
  in the docs.
- `mrna_window` must be real transcript context located in the actual
  full-length sequence by exact substring search (`_locate_window`), with
  the documented fallback to duplex-only context when a site cannot be
  located. Set `has_flanking_context` accordingly. Do not fabricate
  context.
- `label` is normalized to **percent knockdown / inhibition**. If the
  source reports remaining expression, convert it and document the
  conversion.
- `source` is a stable provenance string (e.g. `"Davis2025"`). Use a
  `_`-suffixed discriminator when one source file contains several
  underlying studies (`"<registered source>_<study>"`); `license_for_source`
  resolves those by prefix.
- Chemical modification fields stay at their defaults unless the raw data
  carries real annotation. `sense_modifications` /
  `antisense_modifications` are aligned 1:1 by position or left `None`
  entirely -- not padded or guessed.
- If the source overlaps existing data, dedup by exact sequence against
  what is already loaded, the way `_load_davis2025_records` and
  `_load_cmsirnadb_full_records` do, and state in the docs what it was
  deduped against.

**6. Add an `include_<source>: bool = True` flag** to `load_records()` and
gate the loader behind it. Every source is individually toggleable; none is
forced on.

**7. Add a `SourceLicense` entry** to `SOURCE_LICENSES` in `licenses.py`,
keyed by the flag suffix, with `record_sources` listing every `source`
string the loader emits.

**8. Update the documentation.** All four, kept consistent:

- [`NOTICE.md`](NOTICE.md) -- the per-source license table (mirrors
  `SOURCE_LICENSES`; keep the two in sync).
- [`data/DATA_SOURCES.md`](data/DATA_SOURCES.md) -- a full section: citation,
  license terms, what was fetched, how records were verified, what was
  dropped and why, and any caveats.
- [`data/data_source_ledger.csv`](data/data_source_ledger.csv) -- the
  machine-readable row.
- [`README.md`](README.md) -- the "Data sources at a glance" and gene tables,
  plus the record/gene totals in the intro. The gene table is computed from
  the fetched files rather than hand-maintained; regenerate it rather than
  editing rows by hand.

**9. Add tests** using fixtures, covering at minimum: the loader's happy
path, a record that cannot be located in its transcript, and the
`include_*` flag turning the source off.

A source that arrives without steps 1, 7 and 8 will get sent back for them,
so it is easier to do them as you go.

## Reporting a data problem

Open an issue with the `row_id`, the source, what the raw file says, what
the loader produced, and what you believe is correct.

The governing policy is in
[`data/DATA_SOURCES.md`](data/DATA_SOURCES.md) -- "Known data-quality
caveats (do not 'fix' silently)". **Upstream errors are documented and
flagged, not quietly patched.** The dataset's value is that it reports what
the sources actually say; a silent correction makes it irreproducible
against the original. The `EGFP` / `NZ_CP024869` accession is the worked
example: it is almost certainly lab-plasmid contamination in a public
assembly, and it is documented as such and excluded from the transcript-
length column -- not rewritten.

Two textually-distinct gene strings that mean the same thing (`"EGFP"` vs
`"EGFP "` with a trailing space) are likewise kept distinct, because that is
how the raw data actually compares. Normalizing them is a behavior change
that needs discussing in an issue first.

## Code style

- Ruff with `E, F, I, UP, B`, line length 100. Run `ruff check .` (and
  `ruff format` if you use it) before committing.
- `from __future__ import annotations` at the top of every module; modern
  built-in generics (`list[str]`, `X | None`) throughout.
- Public functions carry type annotations and docstrings. `mypy` runs with
  `check_untyped_defs`, `warn_unused_ignores` and `warn_redundant_casts` --
  a `# type: ignore` that is not needed is an error.
- This codebase comments *why*, not *what*, and existing comments explaining
  a non-obvious decision (a dedup rule, a license constraint, a version pin)
  are load-bearing. Please write in the same register, and do not strip
  them.

## Pull requests

- Branch off `main`; keep each PR to one concern.
- Say what you changed and why in the description; for data changes,
  include the before/after record and gene counts.
- All three CI checks must pass.
- New behavior needs tests. Behavior changes need the docs updated in the
  same PR.

## Releases

Maintainer-only:

1. Bump `version` in `pyproject.toml`.
2. Update the record/gene counts in `README.md` if the dataset changed.
3. Publish a GitHub Release. `.github/workflows/publish.yml` builds and
   uploads to PyPI via Trusted Publishing (OIDC) -- no token is stored in
   the repo.

## Questions

Open an issue at
https://github.com/BrandonWalk/sirna-data-grabber/issues. For anything about
data licensing or reuse, please read [`NOTICE.md`](NOTICE.md) first -- it
answers most of it, and it is the authoritative version.

By contributing, you agree that your code contributions are licensed under
the repo's [MIT license](LICENSE). Data contributions remain under their
original sources' terms, which the MIT license does not and cannot override.
