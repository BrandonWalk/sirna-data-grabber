from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from sirna_data.fetch import cli


@pytest.fixture
def mock_sources(monkeypatch: pytest.MonkeyPatch):
    """Replace every real fetcher with a no-op mock -- these tests only
    exercise the CLI's argument parsing / dest resolution / --only
    filtering, never real network calls."""
    mocks = {name: MagicMock() for name in cli.SOURCES}
    monkeypatch.setattr(cli, "SOURCES", mocks)
    return mocks


def test_default_dest_is_data_raw(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mock_sources
):
    resolved_tmp = tmp_path.resolve()
    monkeypatch.delenv("SIRNA_DATA_DIR", raising=False)
    monkeypatch.chdir(resolved_tmp)
    monkeypatch.setattr("sys.argv", ["sirna-data-fetch"])

    cli.main()

    expected = resolved_tmp / "data" / "raw"
    assert expected.is_dir()
    for mock in mock_sources.values():
        mock.assert_called_once_with(expected)


def test_sirna_data_dir_env_var_used_as_default(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mock_sources
):
    env_dest = (tmp_path / "from_env").resolve()
    monkeypatch.setenv("SIRNA_DATA_DIR", str(env_dest))
    monkeypatch.setattr("sys.argv", ["sirna-data-fetch"])

    cli.main()

    assert env_dest.is_dir()
    for mock in mock_sources.values():
        mock.assert_called_once_with(env_dest)


def test_explicit_dest_overrides_env_var(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mock_sources
):
    monkeypatch.setenv("SIRNA_DATA_DIR", str(tmp_path / "ignored"))
    explicit_dest = (tmp_path / "explicit").resolve()
    monkeypatch.setattr("sys.argv", ["sirna-data-fetch", "--dest", str(explicit_dest)])

    cli.main()

    assert explicit_dest.is_dir()
    for mock in mock_sources.values():
        mock.assert_called_once_with(explicit_dest)


def test_only_filters_to_requested_sources(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mock_sources
):
    resolved_tmp = tmp_path.resolve()
    monkeypatch.setattr(
        "sys.argv",
        ["sirna-data-fetch", "--dest", str(resolved_tmp), "--only", "monopoli", "shabalina"],
    )

    cli.main()

    mock_sources["monopoli"].assert_called_once_with(resolved_tmp)
    mock_sources["shabalina"].assert_called_once_with(resolved_tmp)
    mock_sources["sirna_efficacy"].assert_not_called()
    mock_sources["cmsirnadb"].assert_not_called()


def test_only_rejects_unknown_source(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mock_sources
):
    monkeypatch.setattr(
        "sys.argv", ["sirna-data-fetch", "--dest", str(tmp_path), "--only", "not_a_source"]
    )

    with pytest.raises(SystemExit):
        cli.main()

    for mock in mock_sources.values():
        mock.assert_not_called()
