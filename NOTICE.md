# NOTICE — data licenses

This repo's code (`sirna_data`, `tests/`) is MIT licensed (see
[`LICENSE`](LICENSE)) and may be used, modified, and redistributed freely,
including commercially.

The datasets in `data/raw/` are separate from the code and are NOT MIT
licensed. Each was fetched from its original publisher and is redistributed
here under that publisher's own terms. The MIT license on the code does not
extend to the data, and using this permissively-licensed loader to read the
data does not lift the data's own restrictions. Most sources below are
non-commercial only — read this table before using the data for anything
beyond non-commercial research, and see
[`data/DATA_SOURCES.md`](data/DATA_SOURCES.md) for full terms and provenance
per source.

This table is also available in code as
`sirna_data.list_sources()` / `sirna_data.list_licenses()`, and
`load_records(licenses=[...])` loads only the sources carrying the licenses
you name (see the README's "Loading only the licenses you can use"). The
code copy lives in `src/sirna_data/licenses.py` and is kept in sync with
this file — it is a convenience, not legal advice.

## Sources loaded by `load_records()`

| Source | File(s) | License | Commercial use? |
|---|---|---|---|
| siRNAEfficacyDB (Zhang et al. 2024) | `sirna_efficacy.csv` | CC BY-NC | No — non-commercial only |
| Monopoli et al. 2023 | `monopoli_extra.csv` | CC BY 4.0 | Yes, with attribution |
| Shabalina et al. 2006 | `shabalina_extra.csv` | CC BY 2.0 | Yes, with attribution |
| Martinelli et al. 2023 / `sirna-reproduction` | `martinelli_extra.csv` | CC BY-NC 4.0 — the bioRxiv preprint's stated terms; attribution required, derivatives allowed (this file is a derivative: 577 of the source's 907 rows, with gene identity recovered computationally and corrupted sense strands corrected here, not supplied by the source — see `data/DATA_SOURCES.md`). | No — non-commercial only |
| Harborth et al. 2003 lamin A/C panel, via Ichihara et al. 2007 (NAR supplement) | `harborth2003_extra.csv` | CC BY-NC 2.0 UK — Ichihara 2007's terms (Harborth 2003 itself is not open access and carries no reuse license; these values come from Ichihara's republication). Cite both papers. | No — non-commercial only |
| Sciabola et al. 2013 in-house panel (NAR Supp. Tables S3/S4) | `sciabola2013_extra.csv` | CC BY-NC 3.0 — the article's own stated terms; attribution required, derivatives allowed. | No — non-commercial only |
| CMsiRNAdb (He et al. 2026) | `cmsirnadb_full_raw.tsv` | CC BY-NC-ND 4.0 | No — non-commercial only, and the "ND" term means only the original unmodified file may be redistributed (see below) |
| Davis et al. 2025 (NAR gkaf479) | `davis2025_extra.csv` | CC BY 4.0 | Yes, with attribution |
| NCBI RefSeq/GenBank transcripts | `*_transcripts.fasta` | Public domain | Yes, unrestricted |

CMsiRNAdb's "No Derivatives" term: this repo ships only the untouched
original `cmsirnadb_full_raw.tsv` download. All filtering, collapsing, and
transformation happens in code at load time
(`_load_cmsirnadb_records`/`_load_cmsirnadb_full_records` in
`src/sirna_data/raw_loader.py`), not as a precomputed derivative file — so no
adaptation of CMsiRNAdb's data is redistributed, only the original plus code
that anyone can run themselves. See `data/DATA_SOURCES.md`.

## If you're not sure whether your use is covered

None of the above is legal advice. If your use case isn't clearly
non-commercial research, check the original source's license directly (links
in `data/DATA_SOURCES.md`) or contact the original authors before relying on
this data.
