"""Tests de l'analyse des séquences et des extrémités."""

import pytest

from pepforge.exceptions import SequenceParseError, UnknownResidueError, UnknownTerminusError
from pepforge.io.sequence import parse_sequence


def test_one_letter_without_termini_uses_defaults() -> None:
    pep = parse_sequence("RGDS")
    assert pep.one_letter_sequence == "RGDS"
    assert pep.n_terminus.key == "free_amine"
    assert pep.c_terminus.key == "free_acid"
    assert "n_terminus_source" in pep.metadata
    assert "c_terminus_source" in pep.metadata


@pytest.mark.parametrize(
    ("text", "n_key", "c_key"),
    [
        ("Ac-RGDS-NH2", "acetyl", "amide"),
        ("NH2-RGDS-CO-NH2", "free_amine", "amide"),
        ("H-RGDS-OH", "free_amine", "free_acid"),
        ("H2N-RGDS-COOH", "free_amine", "free_acid"),
        ("Ac-RGDS-CO-NH-CH3", "acetyl", "n_methylamide"),
    ],
)
def test_termini_notations(text: str, n_key: str, c_key: str) -> None:
    pep = parse_sequence(text)
    assert pep.one_letter_sequence == "RGDS"
    assert pep.n_terminus.key == n_key
    assert pep.c_terminus.key == c_key


def test_amide_is_not_n_methylamide() -> None:
    """-CO-NH2 et -CO-NH-CH3 sont deux extrémités distinctes (NHE ≠ NME en AMBER)."""
    assert parse_sequence("RGDS-CO-NH2").c_terminus != parse_sequence("RGDS-CO-NH-CH3").c_terminus


def test_three_letter_notation() -> None:
    pep = parse_sequence("Ac-Arg-Gly-Asp-Ser-NH2")
    assert pep.three_letter_sequence == "ARG-GLY-ASP-SER"
    assert pep.composition() == {"ARG": 1, "GLY": 1, "ASP": 1, "SER": 1}


def test_single_title_case_token_is_three_letter() -> None:
    assert parse_sequence("Ala").one_letter_sequence == "A"
    assert parse_sequence("ALA").one_letter_sequence == "ALA"


def test_whitespace_is_ignored() -> None:
    assert parse_sequence(" Ac - RGD S - NH2 ").one_letter_sequence == "RGDS"


def test_histidine_alone_is_a_residue() -> None:
    pep = parse_sequence("H")
    assert pep.three_letter_sequence == "HIS"


@pytest.mark.parametrize("text", ["", "   ", "RG--DS", "-RGDS"])
def test_malformed_sequences(text: str) -> None:
    with pytest.raises(SequenceParseError):
        parse_sequence(text)


def test_unknown_residue() -> None:
    with pytest.raises(UnknownResidueError):
        parse_sequence("RGBX")
    with pytest.raises(UnknownResidueError):
        parse_sequence("Arg-Xyz-Gly")


def test_unknown_terminus_lookup() -> None:
    from pepforge.residues import load_default_library

    with pytest.raises(UnknownTerminusError):
        load_default_library().c_terminus("CO-SH")
