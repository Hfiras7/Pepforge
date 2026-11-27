"""Tests de l'interface en ligne de commande."""

import pytest

from pepforge.cli.main import main


def test_parse_command(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["parse", "Ac-RGDS-NH2"]) == 0
    out = capsys.readouterr().out
    assert "RGDS" in out
    assert "amidée" in out


def test_parse_warns_on_default_termini(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["parse", "RGDS"]) == 0
    assert "valeur par défaut" in capsys.readouterr().out


def test_parse_error_returns_1(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["parse", "RGBX"]) == 1
    assert "Erreur" in capsys.readouterr().err


def test_config_template(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["config-template"]) == 0
    assert "champ de force" in capsys.readouterr().out


def test_parse_warns_on_ambiguous_reading(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["parse", "Ac-ALA-NHMe"]) == 0
    out = capsys.readouterr().out
    assert "ALA-LEU-ALA" in out
    assert "lecture à confirmer" in out
