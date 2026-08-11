"""Command-line entry point: `sirna-data-fetch` (registered via
[project.scripts] in pyproject.toml). Fetches every source
`sirna_data.load_records()` reads, from its original location, into a local
directory.

    sirna-data-fetch --dest ./my_data
    export SIRNA_DATA_DIR=./my_data

This is installed automatically with `pip install sirna-data-grabber` --
no extras needed, since every fetcher here only needs pandas (already a
core dependency) plus the standard library.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

from . import cmsirnadb, monopoli, shabalina, sirna_efficacy

# Order matters only for display; each fetcher is independent.
SOURCES = {
    "sirna_efficacy": sirna_efficacy.fetch,
    "monopoli": monopoli.fetch,
    "shabalina": shabalina.fetch,
    "cmsirnadb": cmsirnadb.fetch,
}


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="sirna-data-fetch",
        description=(
            "Fetch the siRNA knockdown-efficacy dataset (siRNAEfficacyDB + "
            "supplementary sources) from its original sources into a local "
            "directory, for use with sirna_data.load_records()."
        ),
    )
    parser.add_argument(
        "--dest",
        type=Path,
        default=Path(os.environ.get("SIRNA_DATA_DIR", "data/raw")),
        help=(
            "Directory to write the fetched files into "
            "(default: $SIRNA_DATA_DIR if set, else ./data/raw)."
        ),
    )
    parser.add_argument(
        "--only",
        nargs="+",
        choices=sorted(SOURCES),
        metavar="SOURCE",
        help=f"Fetch only these sources instead of all four ({', '.join(sorted(SOURCES))}).",
    )
    args = parser.parse_args()

    dest: Path = args.dest.resolve()
    dest.mkdir(parents=True, exist_ok=True)
    names = args.only or sorted(SOURCES)

    for name in names:
        print(f"=== {name} ===")
        SOURCES[name](dest)
        print()

    print(f"Done. To use this data:\n  export SIRNA_DATA_DIR={dest}")


if __name__ == "__main__":
    main()
