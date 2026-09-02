# NOTICE — data licenses

This repo's code (`sirna_data`, `tests/`) is MIT licensed (see
[`LICENSE`](LICENSE)) and may be used, modified, and redistributed freely,
including commercially.

**The datasets in `data/raw/` are separate from the code and are NOT MIT
licensed.** Each was fetched from its original publisher and is redistributed
here under that publisher's own terms. The MIT license on the code does not
extend to the data, and using this permissively-licensed loader to read the
data does not lift the data's own restrictions. Most sources below are
**non-commercial only** — read this table before using the data for anything
beyond non-commercial research, and see
[`data/DATA_SOURCES.md`](data/DATA_SOURCES.md) for full terms and provenance
per source.

## Sources loaded by `load_records()`

| Source | File(s) | License | Commercial use? |
|---|---|---|---|
| siRNAEfficacyDB (Zhang et al. 2024) | `sirna_efficacy.csv` | CC BY-NC | **No** — non-commercial only |
| Monopoli et al. 2023 | `monopoli_extra.csv` | CC BY 4.0 | Yes, with attribution |
| PDCD1 panel (Xu/Zhao et al. 2024, siRNABERT repo) | `pdcd1_extra.csv` | **Unresolved — no LICENSE file in the source repo, paper not confirmed open-access.** Committed at the user's explicit request despite this; see below. | Unclear — do not assume |
| Shabalina et al. 2006 | `shabalina_extra.csv` | CC BY 2.0 | Yes, with attribution |
| OligoGraph repo (Sciabola et al. 2013 + Harborth et al. 2001) | `oligograph_extra.csv` | **Unresolved — no LICENSE file in the source repo.** | Unclear — do not assume |
| CMsiRNAdb (He et al. 2026) | `cmsirnadb_full_raw.tsv` | CC BY-NC-ND 4.0 | **No** — non-commercial only, and the "ND" term means only the original unmodified file may be redistributed (see below) |
| Davis et al. 2025 (NAR gkaf479) | `davis2025_extra.csv` | CC BY 4.0 | Yes, with attribution |
| NCBI RefSeq/GenBank transcripts | `*_transcripts.fasta` | Public domain | Yes, unrestricted |

**PDCD1 panel's unresolved license**: recovered from a deleted file in
github.com/ChengkuiZhao/siRNABERT's git history (commit `f3254ee`, before
its removal in `d2ad931`) — a real but non-canonical channel, and that repo
carries no `LICENSE` file (all-rights-reserved by default). Loaded by
`load_records()` when present locally (`include_pdcd1=True` by default,
gracefully returns nothing if the file is absent). Unlike siRecords below,
`pdcd1_extra.csv` **is committed to this repo** — at the user's explicit
request, made and confirmed with the license gap called out beforehand. If
you're relying on this data outside non-commercial research, verify
redistribution rights directly with the siRNABERT authors first; this
repo's own permissive (MIT) license does not extend to it. The companion
transcript FASTA is pure NCBI RefSeq (public domain) and is committed
normally. See `data/DATA_SOURCES.md`.

**CMsiRNAdb's "No Derivatives" term**: this repo ships only the untouched
original `cmsirnadb_full_raw.tsv` download. All filtering, collapsing, and
transformation happens in code at load time
(`_load_cmsirnadb_records`/`_load_cmsirnadb_full_records` in
`src/sirna_data/raw_loader.py`), not as a precomputed derivative file — so no
adaptation of CMsiRNAdb's data is redistributed, only the original plus code
that anyone can run themselves. See `data/DATA_SOURCES.md`.

## Other files present in `data/raw/` but not used by anything in `sirna_data`

These were investigated as candidate sources (see
[`data/DATA_SOURCES.md`](data/DATA_SOURCES.md)) and kept for
reference/provenance, but nothing in `sirna_data` reads them:

| File(s) | Source | License |
|---|---|---|
| `sirecords_efficacy.csv`, `sirecords_new_only.csv` | siRecords (04/28/05 release, via Internet Archive) | **Restricted, not established for redistribution.** Researched: the database's own paper (Ren et al. 2009, NAR 37:D146-D149) is CC BY-NC licensed as an *article*, but its DATA ACCESS section says bulk copies were only ever given to "academic users" who emailed the authors directly — not published as an open download. The data here was recovered from an Internet Archive snapshot, not that channel, so no license actually covers this copy. See `data/DATA_SOURCES.md` for the full writeup and sources. Kept locally but excluded from git (`.gitignore`). |

## If you're not sure whether your use is covered

None of the above is legal advice. If your use case isn't clearly
non-commercial research, check the original source's license directly (links
in `data/DATA_SOURCES.md`) or contact the original authors before relying on
this data.
