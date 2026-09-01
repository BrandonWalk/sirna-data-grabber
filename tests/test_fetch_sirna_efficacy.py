from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

from sirna_data.fetch import sirna_efficacy


def _fasta_bytes(entries: dict[str, str]) -> bytes:
    """entries: accession (with optional version) -> sequence."""
    lines = []
    for acc, seq in entries.items():
        lines.append(f">{acc} some description")
        lines.append(seq)
    return ("\n".join(lines) + "\n").encode()


def _mock_urlopen(payload: bytes):
    cm = MagicMock()
    cm.__enter__.return_value.read.return_value = payload
    return cm


def test_download_efficacy_table_reads_tsv_and_reports_counts(tmp_path: Path, capsys):
    tsv_path = tmp_path / "source.tsv"
    tsv_path.write_text("Gene\tAccession_number\nTP53\tNM_000546\nTP53\tNM_000546\nBRCA1\t-\n")

    def fake_urlretrieve(url, filename):
        Path(filename).write_bytes(tsv_path.read_bytes())

    with patch.object(
        sirna_efficacy.urllib.request, "urlretrieve", side_effect=fake_urlretrieve
    ) as mock_retrieve:
        df = sirna_efficacy.download_efficacy_table()

    mock_retrieve.assert_called_once()
    assert mock_retrieve.call_args[0][0] == sirna_efficacy.EFFICACY_DB_URL
    assert len(df) == 3
    assert set(df["Gene"]) == {"TP53", "BRCA1"}
    out = capsys.readouterr().out
    assert "Downloaded 3 siRNA records covering 2 genes" in out


def test_fetch_mrna_fasta_single_batch_parses_headers_and_dedupes():
    with patch.object(
        sirna_efficacy.urllib.request,
        "urlopen",
        return_value=_mock_urlopen(_fasta_bytes({"NM_000001.2": "ACGT", "NM_000002.1": "TTGG"})),
    ) as mock_urlopen, patch.object(sirna_efficacy.time, "sleep") as mock_sleep:
        result = sirna_efficacy.fetch_mrna_fasta(["NM_000002", "NM_000001", "NM_000001"])

    assert mock_urlopen.call_count == 1
    # version suffix stripped from FASTA headers when building keys
    assert result == {"NM_000001": "ACGT", "NM_000002": "TTGG"}
    mock_sleep.assert_called_once_with(sirna_efficacy.REQUEST_DELAY_S)


def test_fetch_mrna_fasta_batches_large_accession_lists():
    accessions = [f"NM_{i:06d}" for i in range(sirna_efficacy.BATCH_SIZE + 5)]
    batch1 = {acc: "ACGT" for acc in accessions[: sirna_efficacy.BATCH_SIZE]}
    batch2 = {acc: "TTGG" for acc in accessions[sirna_efficacy.BATCH_SIZE :]}
    with patch.object(
        sirna_efficacy.urllib.request,
        "urlopen",
        side_effect=[_mock_urlopen(_fasta_bytes(batch1)), _mock_urlopen(_fasta_bytes(batch2))],
    ) as mock_urlopen, patch.object(sirna_efficacy.time, "sleep"):
        result = sirna_efficacy.fetch_mrna_fasta(accessions)

    assert mock_urlopen.call_count == 2
    assert len(result) == len(accessions)


def test_fetch_mrna_fasta_empty_input_makes_no_request():
    with patch.object(sirna_efficacy.urllib.request, "urlopen") as mock_urlopen:
        result = sirna_efficacy.fetch_mrna_fasta([])
    mock_urlopen.assert_not_called()
    assert result == {}


def test_fetch_writes_csv_and_fasta_dropping_control_rows(tmp_path: Path, capsys):
    fake_df_rows = {
        "Gene": ["TP53", "TP53", "CTRL"],
        "Accession_number": ["NM_000546", "NM_000546", "-"],
    }
    import pandas as pd

    fake_df = pd.DataFrame(fake_df_rows)

    with patch.object(
        sirna_efficacy, "download_efficacy_table", return_value=fake_df
    ) as mock_download, patch.object(
        sirna_efficacy,
        "fetch_mrna_fasta",
        return_value={"NM_000546": "ACGUACGU"},
    ) as mock_fetch_fasta:
        sirna_efficacy.fetch(tmp_path)

    mock_download.assert_called_once()
    mock_fetch_fasta.assert_called_once_with(["NM_000546"])

    csv_path = tmp_path / "sirna_efficacy.csv"
    fasta_path = tmp_path / "mrna_transcripts.fasta"
    assert csv_path.exists() and fasta_path.exists()

    written = pd.read_csv(csv_path)
    assert len(written) == 2  # control row (Accession_number == "-") dropped
    assert "CTRL" not in set(written["Gene"])

    fasta_text = fasta_path.read_text()
    assert ">NM_000546" in fasta_text
    assert "ACGUACGU" in fasta_text

    out = capsys.readouterr().out
    assert "Dropped 1 control rows without a gene accession" in out


def test_fetch_warns_about_unresolved_accessions(tmp_path: Path, capsys):
    import pandas as pd

    fake_df = pd.DataFrame(
        {"Gene": ["TP53", "BRCA1"], "Accession_number": ["NM_000546", "NM_007294"]}
    )
    with (
        patch.object(sirna_efficacy, "download_efficacy_table", return_value=fake_df),
        patch.object(sirna_efficacy, "fetch_mrna_fasta", return_value={"NM_000546": "ACGT"}),
    ):
        sirna_efficacy.fetch(tmp_path)

    out = capsys.readouterr().out
    assert "WARNING: 1 accessions did not resolve" in out
    assert "NM_007294" in out
