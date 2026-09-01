from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd

from sirna_data.fetch import cmsirnadb


def _fasta_bytes(entries: dict[str, str]) -> bytes:
    lines = []
    for acc, seq in entries.items():
        lines.append(f">{acc} description")
        lines.append(seq)
    return ("\n".join(lines) + "\n").encode()


def _mock_urlopen(payload: bytes):
    cm = MagicMock()
    cm.__enter__.return_value.read.return_value = payload
    return cm


def test_download_cmsirnadb_table_writes_raw_tsv_and_reads_as_str(tmp_path: Path, capsys):
    def fake_urlretrieve(url, filename):
        Path(filename).write_text(
            "Target_Gene\tAccession_number\nPCSK9\tNM_174936.4\nAGT\tNM_000029.4\n"
        )

    with patch.object(
        cmsirnadb.urllib.request, "urlretrieve", side_effect=fake_urlretrieve
    ) as mock_retrieve:
        df = cmsirnadb.download_cmsirnadb_table(tmp_path)

    mock_retrieve.assert_called_once()
    assert mock_retrieve.call_args[0][0] == cmsirnadb.CMSIRNADB_URL
    assert mock_retrieve.call_args[0][1] == tmp_path / "cmsirnadb_full_raw.tsv"
    assert (tmp_path / "cmsirnadb_full_raw.tsv").exists()
    assert len(df) == 2
    # dtype=str: accession strings must round-trip exactly, not get coerced
    assert df["Accession_number"].tolist() == ["NM_174936.4", "NM_000029.4"]

    out = capsys.readouterr().out
    assert "Downloaded 2 CMsiRNAdb rows covering 2 genes" in out


def test_fetch_transcript_fasta_keeps_full_versioned_accession():
    payload = _fasta_bytes({"NM_174936.4": "ACGT"})
    with patch.object(
        cmsirnadb.urllib.request, "urlopen", return_value=_mock_urlopen(payload)
    ) as mock_urlopen, patch.object(cmsirnadb.time, "sleep") as mock_sleep:
        result = cmsirnadb.fetch_transcript_fasta(["NM_174936.4"])

    mock_urlopen.assert_called_once()
    mock_sleep.assert_called_once_with(cmsirnadb.REQUEST_DELAY_S)
    # unlike sirna_efficacy.fetch_mrna_fasta, the version suffix is KEPT here
    assert result == {"NM_174936.4": "ACGT"}


def test_fetch_transcript_fasta_batches_large_lists():
    accessions = [f"NM_{i:06d}.1" for i in range(cmsirnadb.BATCH_SIZE + 3)]
    batch1 = {acc: "AAAA" for acc in accessions[: cmsirnadb.BATCH_SIZE]}
    batch2 = {acc: "TTTT" for acc in accessions[cmsirnadb.BATCH_SIZE :]}
    with patch.object(
        cmsirnadb.urllib.request,
        "urlopen",
        side_effect=[_mock_urlopen(_fasta_bytes(batch1)), _mock_urlopen(_fasta_bytes(batch2))],
    ) as mock_urlopen, patch.object(cmsirnadb.time, "sleep"):
        result = cmsirnadb.fetch_transcript_fasta(accessions)

    assert mock_urlopen.call_count == 2
    assert len(result) == len(accessions)


def test_write_fasta_writes_one_record_per_entry(tmp_path: Path):
    path = tmp_path / "out.fasta"
    cmsirnadb.write_fasta(path, {"ACC1": "ACGT", "ACC2": "TTTT"})
    text = path.read_text()
    assert text == ">ACC1\nACGT\n>ACC2\nTTTT\n"


def test_fetch_splits_pcsk9_from_other_genes_into_separate_fasta_files(tmp_path: Path, capsys):
    fake_df = pd.DataFrame(
        {
            "Target_Gene": ["PCSK9", "PCSK9", "AGT", "MSTN"],
            "Accession_number": [
                "NM_000000.1",  # ignored for PCSK9 rows -- canonical accession used instead
                None,
                "NM_000029.4",
                "NM_005259.3",
            ],
        }
    )
    pcsk9_result = {cmsirnadb.PCSK9_ACCESSION: "PCSK9SEQ"}
    others_result = {"NM_000029.4": "AGTSEQ", "NM_005259.3": "MSTNSEQ"}

    with patch.object(
        cmsirnadb, "download_cmsirnadb_table", return_value=fake_df
    ), patch.object(
        cmsirnadb, "fetch_transcript_fasta", side_effect=[pcsk9_result, others_result]
    ) as mock_fetch:
        cmsirnadb.fetch(tmp_path)

    first_call, second_call = mock_fetch.call_args_list
    assert first_call.args[0] == [cmsirnadb.PCSK9_ACCESSION]
    assert sorted(second_call.args[0]) == ["NM_000029.4", "NM_005259.3"]

    pcsk9_fasta = (tmp_path / "cmsirnadb_transcripts.fasta").read_text()
    assert ">NM_174936.4" in pcsk9_fasta
    full_fasta = (tmp_path / "cmsirnadb_full_transcripts.fasta").read_text()
    assert ">NM_000029.4" in full_fasta and ">NM_005259.3" in full_fasta

    out = capsys.readouterr().out
    assert "Fetching PCSK9 canonical transcript" in out
    assert "Fetching 2 unique transcripts for the other 12 genes" in out


def test_fetch_warns_when_pcsk9_or_other_accessions_unresolved(tmp_path: Path, capsys):
    fake_df = pd.DataFrame(
        {"Target_Gene": ["PCSK9", "AGT"], "Accession_number": [None, "NM_000029.4"]}
    )
    with patch.object(
        cmsirnadb, "download_cmsirnadb_table", return_value=fake_df
    ), patch.object(cmsirnadb, "fetch_transcript_fasta", side_effect=[{}, {}]):
        cmsirnadb.fetch(tmp_path)

    out = capsys.readouterr().out
    assert f"WARNING: {cmsirnadb.PCSK9_ACCESSION} did not resolve" in out
    assert "WARNING: 1 accessions did not resolve" in out
