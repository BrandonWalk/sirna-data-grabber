from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

from sirna_data.fetch import monopoli


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


def test_fetch_transcripts_requests_all_four_genes_and_strips_version():
    payload = _fasta_bytes(
        {
            "NM_000484.4": "AAAA",
            "NM_005910.6": "CCCC",
            "NM_012104.5": "GGGG",
            "NM_000345.4": "UUUU",
        }
    )
    with patch.object(
        monopoli.urllib.request, "urlopen", return_value=_mock_urlopen(payload)
    ) as mock_urlopen:
        result = monopoli.fetch_transcripts()

    (call_args,), _ = mock_urlopen.call_args
    for acc in monopoli.GENE_ACCESSION.values():
        assert acc in call_args
    assert result == {
        "NM_000484": "AAAA",
        "NM_005910": "CCCC",
        "NM_012104": "GGGG",
        "NM_000345": "UUUU",
    }


def test_fetch_writes_table_s3_csv_and_transcripts_fasta(tmp_path: Path, capsys):
    fake_sequences = {acc: "ACGU" for acc in monopoli.GENE_ACCESSION.values()}
    with patch.object(monopoli, "fetch_transcripts", return_value=fake_sequences):
        monopoli.fetch(tmp_path)

    csv_path = tmp_path / "monopoli_extra.csv"
    fasta_path = tmp_path / "monopoli_transcripts.fasta"
    assert csv_path.exists() and fasta_path.exists()

    lines = csv_path.read_text().strip().splitlines()
    assert lines[0] == "Sequence,Gene,Accession_number,Reporter_Remaining_Pct,SD_Pct"
    assert len(lines) == 1 + len(monopoli.TABLE_S3)
    # spot-check one known row round-trips with its gene's accession
    first_seq, first_gene, first_remaining, first_sd = monopoli.TABLE_S3[0]
    assert lines[1] == (
        f"{first_seq},{first_gene},{monopoli.GENE_ACCESSION[first_gene]},"
        f"{first_remaining},{first_sd}"
    )

    fasta_text = fasta_path.read_text()
    for acc in monopoli.GENE_ACCESSION.values():
        assert f">{acc}" in fasta_text

    out = capsys.readouterr().out
    assert f"Wrote {csv_path} ({len(monopoli.TABLE_S3)} rows)" in out
    assert f"Wrote {fasta_path} ({len(fake_sequences)} transcripts)" in out
