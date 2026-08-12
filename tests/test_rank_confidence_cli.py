from __future__ import annotations

import pytest

from sirna_data.rank_confidence_cli import DEFAULT_CONFIDENCE_LEVELS, main


def test_pcc_prints_one_row_per_confidence_level(monkeypatch: pytest.MonkeyPatch, capsys):
    monkeypatch.setattr(
        "sys.argv",
        ["sirna-rank-confidence", "--pcc", "0.3686", "--n-items", "4561", "--confidence", "0.9"],
    )
    main()
    out = capsys.readouterr().out
    assert "pcc = 0.3686, n_items = 4561" in out
    assert "0.9" in out


def test_spcc_flag_labelled_correctly(monkeypatch: pytest.MonkeyPatch, capsys):
    monkeypatch.setattr(
        "sys.argv",
        ["sirna-rank-confidence", "--spcc", "0.5", "--n-items", "1000", "--confidence", "0.8"],
    )
    main()
    out = capsys.readouterr().out
    assert "spcc = 0.5, n_items = 1000" in out


def test_default_confidence_levels_used_when_omitted(monkeypatch: pytest.MonkeyPatch, capsys):
    monkeypatch.setattr(
        "sys.argv", ["sirna-rank-confidence", "--pcc", "0.4", "--n-items", "500"]
    )
    main()
    out = capsys.readouterr().out
    # one printed row per default confidence level
    for level in DEFAULT_CONFIDENCE_LEVELS:
        assert f"{level:>10.4g}" in out


def test_multiple_confidence_levels_all_printed(monkeypatch: pytest.MonkeyPatch, capsys):
    monkeypatch.setattr(
        "sys.argv",
        [
            "sirna-rank-confidence",
            "--pcc",
            "0.4",
            "--n-items",
            "500",
            "--confidence",
            "0.99",
            "0.5",
        ],
    )
    main()
    out = capsys.readouterr().out
    assert out.count("\n") >= 4  # header line + 2 rows + blank line, at least


def test_requires_exactly_one_of_spcc_or_pcc(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr("sys.argv", ["sirna-rank-confidence", "--n-items", "100"])
    with pytest.raises(SystemExit):
        main()

    monkeypatch.setattr(
        "sys.argv",
        ["sirna-rank-confidence", "--spcc", "0.5", "--pcc", "0.5", "--n-items", "100"],
    )
    with pytest.raises(SystemExit):
        main()


def test_requires_n_items(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr("sys.argv", ["sirna-rank-confidence", "--pcc", "0.5"])
    with pytest.raises(SystemExit):
        main()


def test_rejects_out_of_range_pcc_via_underlying_function(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(
        "sys.argv",
        ["sirna-rank-confidence", "--pcc", "1.5", "--n-items", "100", "--confidence", "0.9"],
    )
    with pytest.raises(ValueError):
        main()
