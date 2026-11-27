"""Tests de l'interface en ligne de commande."""

import pytest

from pepforge.cli.main import main


def test_parse_command(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["parse", "--notation", "1", "Ac-RGDS-NH2"]) == 0
    out = capsys.readouterr().out
    assert "code à une lettre (monolitère)" in out
    assert "Ac-RGDS-NH2" in out
    assert "Ac-Arg-Gly-Asp-Ser-NH2" in out
    assert "amidée" in out


def test_parse_warns_on_default_termini(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["parse", "-n", "1", "RGDS"]) == 0
    assert "valeur par défaut" in capsys.readouterr().out


def test_parse_error_returns_1(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["parse", "-n", "1", "RGBX"]) == 1
    assert "Erreur" in capsys.readouterr().err


def test_config_template(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["config-template"]) == 0
    assert "champ de force" in capsys.readouterr().out


# --- Mode de lecture déclaré au départ (revue du module 1) ----------------------


@pytest.mark.parametrize("command", ["parse", "convert"])
def test_the_reading_mode_is_required(command: str, capsys: pytest.CaptureFixture[str]) -> None:
    assert main([command, "Ac-RGDS-NH2"]) == 1
    err = capsys.readouterr().err
    assert "Mode de lecture non précisé" in err
    assert "--notation 1" in err
    assert "--notation 3" in err
    assert "Traceback" not in err


@pytest.mark.parametrize("command", ["parse", "convert"])
def test_an_unknown_reading_mode_is_refused(
    command: str, capsys: pytest.CaptureFixture[str]
) -> None:
    assert main([command, "-n", "2", "Ac-RGDS-NH2"]) == 1
    assert "Mode de lecture inconnu : '2'" in capsys.readouterr().err


@pytest.mark.parametrize(
    ("notation", "expected"),
    [
        ("1", "Code à trois lettres : Ac-Ala-Leu-Ala-NHMe"),
        ("3", "Code à trois lettres : Ac-Ala-NHMe"),
    ],
)
def test_parse_reads_in_the_declared_mode(
    notation: str, expected: str, capsys: pytest.CaptureFixture[str]
) -> None:
    assert main(["parse", "-n", notation, "Ac-ALA-NHMe"]) == 0
    assert expected in capsys.readouterr().out


@pytest.mark.parametrize(
    ("notation", "text", "converted"),
    [
        ("1", "Ac-RGDS-NH2", "Ac-Arg-Gly-Asp-Ser-NH2"),
        ("3", "Ac-Arg-Gly-Asp-Ser-NH2", "Ac-RGDS-NH2"),
        ("one_letter", "NH2-RGDS-CO-NH2", "H-Arg-Gly-Asp-Ser-NH2"),
        ("three_letter", "Ac-Ala-NHMe", "Ac-A-NHMe"),
    ],
)
def test_convert_switches_to_the_other_code(
    notation: str, text: str, converted: str, capsys: pytest.CaptureFixture[str]
) -> None:
    assert main(["convert", "-n", notation, text]) == 0
    captured = capsys.readouterr()
    assert captured.out == converted + "\n"
    assert captured.err == ""


def test_convert_keeps_warnings_off_the_result(capsys: pytest.CaptureFixture[str]) -> None:
    """Seule la séquence convertie va sur la sortie standard, pour pouvoir la réutiliser."""
    assert main(["convert", "-n", "1", "RGDS"]) == 0
    captured = capsys.readouterr()
    assert captured.out == "Arg-Gly-Asp-Ser\n"
    assert "extrémité N non déclarée" in captured.err
    assert "extrémité C non déclarée" in captured.err


def test_convert_error_returns_1(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["convert", "-n", "1", "RGDfK"]) == 1
    err = capsys.readouterr().err
    assert "Erreur" in err
    assert "acide aminé D" in err
