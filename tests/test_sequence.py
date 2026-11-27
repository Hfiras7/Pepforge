"""Tests de l'analyse des séquences et des extrémités."""

import dataclasses

import pytest

from pepforge.core import AMBIGUOUS_READING, Residue
from pepforge.exceptions import SequenceParseError, UnknownResidueError, UnknownTerminusError
from pepforge.io.sequence import parse_sequence
from pepforge.residues import ResidueLibrary, load_default_library


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


# --- Notation des non standards entre crochets (étape 1) ------------------------


def _library_with_aib() -> ResidueLibrary:
    """Bibliothèque par défaut plus un non standard, l'acide α-aminoisobutyrique."""
    library = load_default_library()
    residues = dict(library.residues)
    residues["AIB"] = Residue(code="AIB", one_letter=None, standard=False)
    return dataclasses.replace(library, residues=residues)


@pytest.mark.parametrize(
    ("text", "codes"),
    [
        ("Ac-RG[Aib]S-NH2", ("ARG", "GLY", "AIB", "SER")),
        ("[Aib]RG", ("AIB", "ARG", "GLY")),
        ("RG[Aib]", ("ARG", "GLY", "AIB")),
        ("[Aib][aib]", ("AIB", "AIB")),
        ("[Aib]", ("AIB",)),
        ("Ac-Arg-Gly-[Aib]-Ser-NH2", ("ARG", "GLY", "AIB", "SER")),
        ("Ac-Arg-Gly-Aib-Ser-NH2", ("ARG", "GLY", "AIB", "SER")),
        ("Ac - RG [ Aib ] S - NH2", ("ARG", "GLY", "AIB", "SER")),
    ],
)
def test_bracket_notation(text: str, codes: tuple[str, ...]) -> None:
    peptide = parse_sequence(text, _library_with_aib())
    assert tuple(r.code for r in peptide.residues) == codes


def test_bracketed_residue_keeps_its_library_entry() -> None:
    peptide = parse_sequence("Ac-RG[Aib]S-NH2", _library_with_aib())
    assert peptide.residue_at(3).standard is False
    assert peptide.one_letter_sequence == "RGXS"
    assert peptide.n_terminus.key == "acetyl"
    assert peptide.c_terminus.key == "amide"


def test_bracketed_standard_code_is_the_same_residue() -> None:
    assert parse_sequence("Ac-RG[Asp]S-NH2") == parse_sequence("Ac-RGDS-NH2")


def test_both_notations_give_the_same_peptide() -> None:
    library = _library_with_aib()
    assert parse_sequence("Ac-RG[Aib]S-NH2", library) == parse_sequence(
        "Ac-Arg-Gly-[Aib]-Ser-NH2", library
    )


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ("RG[AibS", "jamais fermé"),
        ("RGAib]S", "sans crochet ouvrant"),
        ("RG[]S", "Crochets vides"),
        ("RG[[Aib]]S", "imbriqués"),
        ("RG[Ai-b]S", "non refermé avant le tiret"),
        ("Ac-RG[Aib-NH2", "non refermé avant le tiret"),
        ("Arg-Gly[Aib]-Ser", "occupe un bloc entier"),
        ("Arg-[Aib][Aib]-Ser", "occupe un bloc entier"),
    ],
)
def test_malformed_brackets(text: str, message: str) -> None:
    with pytest.raises(SequenceParseError, match=message):
        parse_sequence(text, _library_with_aib())


def test_bracketed_code_absent_from_the_library() -> None:
    """La bibliothèque par défaut ne contient que les vingt standards."""
    with pytest.raises(UnknownResidueError, match="Résidu inconnu : 'Aib'"):
        parse_sequence("Ac-RG[Aib]S-NH2")


def test_residues_are_a_tuple() -> None:
    assert isinstance(parse_sequence("RGDS").residues, tuple)


# --- Lecture ambiguë signalée (revue du module 1) --------------------------------


@pytest.mark.parametrize(
    ("text", "codes"),
    [
        ("ALA", ("ALA", "LEU", "ALA")),
        ("Ac-ALA-NHMe", ("ALA", "LEU", "ALA")),
        ("ala", ("ALA", "LEU", "ALA")),
        ("H-CYS-OH", ("CYS", "TYR", "SER")),
    ],
)
def test_ambiguous_block_is_read_in_one_letter_codes_and_flagged(
    text: str, codes: tuple[str, ...]
) -> None:
    peptide = parse_sequence(text)
    assert tuple(r.code for r in peptide.residues) == codes
    note = peptide.metadata[AMBIGUOUS_READING]
    assert "-".join(codes) in note
    assert f"écrire « {codes[0].title()} »" in note


@pytest.mark.parametrize(
    "text", ["Ala", "Ac-Ala-NHMe", "RGDS", "Ac-RGDS-NH2", "Ac-ALA-GLY-NH2", "AL[Asp]"]
)
def test_unambiguous_sequences_are_not_flagged(text: str) -> None:
    assert AMBIGUOUS_READING not in parse_sequence(text).metadata


def test_unreadable_ambiguous_block_says_how_to_write_the_residue() -> None:
    with pytest.raises(UnknownResidueError, match="écrire « Glu »"):
        parse_sequence("GLU")
