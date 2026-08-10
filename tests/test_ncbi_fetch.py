from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest
import requests

from sirna_data.ncbi_fetch import GeneNotFoundError, fetch_mrna_by_gene


def _mock_response(*, json_data=None, text="", status_ok=True):
    resp = MagicMock()
    resp.json.return_value = json_data or {}
    resp.text = text
    if status_ok:
        resp.raise_for_status.return_value = None
    else:
        resp.raise_for_status.side_effect = requests.HTTPError("boom")
    return resp


def test_empty_gene_raises_without_any_http_call():
    with patch("sirna_data.ncbi_fetch.requests.get") as mock_get:
        with pytest.raises(GeneNotFoundError):
            fetch_mrna_by_gene("   ")
        mock_get.assert_not_called()


def test_esearch_no_hits_raises_and_skips_efetch():
    esearch_resp = _mock_response(json_data={"esearchresult": {"idlist": []}})
    with patch("sirna_data.ncbi_fetch.requests.get", return_value=esearch_resp) as mock_get:
        with pytest.raises(GeneNotFoundError, match="no RefSeq mRNA found"):
            fetch_mrna_by_gene("NOTAREALGENE")
        assert mock_get.call_count == 1  # never reached efetch


def test_successful_fetch_converts_dna_to_rna_and_parses_header():
    esearch_resp = _mock_response(json_data={"esearchresult": {"idlist": ["999"]}})
    efetch_resp = _mock_response(
        text=">NM_000123.4 Homo sapiens TEST gene, mRNA\nACGTACGT\nACGT\n"
    )
    with patch(
        "sirna_data.ncbi_fetch.requests.get", side_effect=[esearch_resp, efetch_resp]
    ) as mock_get:
        result = fetch_mrna_by_gene("TEST")

    assert mock_get.call_count == 2
    assert result.accession == "NM_000123.4"
    assert result.description == "Homo sapiens TEST gene, mRNA"
    assert result.sequence == "ACGUACGUACGU"  # T -> U, multi-line concatenated


def test_fetch_uses_esearch_id_and_passes_organism():
    esearch_resp = _mock_response(json_data={"esearchresult": {"idlist": ["42"]}})
    efetch_resp = _mock_response(text=">ACC1 desc\nACGT\n")
    with patch(
        "sirna_data.ncbi_fetch.requests.get", side_effect=[esearch_resp, efetch_resp]
    ) as mock_get:
        fetch_mrna_by_gene("MYGENE", organism="Mus musculus")

    esearch_call, efetch_call = mock_get.call_args_list
    assert "Mus musculus" in esearch_call.kwargs["params"]["term"]
    assert esearch_call.kwargs["params"]["db"] == "nuccore"
    assert efetch_call.kwargs["params"]["id"] == "42"
    assert efetch_call.kwargs["params"]["rettype"] == "fasta"


def test_malformed_fasta_response_raises():
    esearch_resp = _mock_response(json_data={"esearchresult": {"idlist": ["1"]}})
    efetch_resp = _mock_response(text="not a fasta file")
    with patch(
        "sirna_data.ncbi_fetch.requests.get", side_effect=[esearch_resp, efetch_resp]
    ):
        with pytest.raises(GeneNotFoundError, match="unexpected response"):
            fetch_mrna_by_gene("TEST")


def test_empty_sequence_response_raises():
    esearch_resp = _mock_response(json_data={"esearchresult": {"idlist": ["1"]}})
    efetch_resp = _mock_response(text=">ACC1 header only, no sequence lines\n")
    with patch(
        "sirna_data.ncbi_fetch.requests.get", side_effect=[esearch_resp, efetch_resp]
    ):
        with pytest.raises(GeneNotFoundError, match="empty sequence"):
            fetch_mrna_by_gene("TEST")


def test_http_error_on_esearch_propagates():
    esearch_resp = _mock_response(status_ok=False)
    with patch("sirna_data.ncbi_fetch.requests.get", return_value=esearch_resp):
        with pytest.raises(requests.HTTPError):
            fetch_mrna_by_gene("TEST")


def test_http_error_on_efetch_propagates():
    esearch_resp = _mock_response(json_data={"esearchresult": {"idlist": ["1"]}})
    efetch_resp = _mock_response(status_ok=False)
    with patch(
        "sirna_data.ncbi_fetch.requests.get", side_effect=[esearch_resp, efetch_resp]
    ):
        with pytest.raises(requests.HTTPError):
            fetch_mrna_by_gene("TEST")
