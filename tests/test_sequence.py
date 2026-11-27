"""Tests de la lecture et de l'écriture des séquences et de leurs extrémités."""

import dataclasses
import inspect

import pytest

from pepforge.core import Link, LinkKind, Modification, Peptide, Residue
from pepforge.exceptions import (
    SequenceParseError,
    UnknownResidueError,
    UnknownTerminusError,
    UnsupportedFeatureError,
)
from pepforge.io.sequence import (
    SEQUENCE_NOTATION,
    SequenceNotation,
    convert_sequence,
    format_sequence,
    parse_sequence,
)
from pepforge.residues import ResidueLibrary, load_default_library

ONE = SequenceNotation.ONE_LETTER
THREE = SequenceNotation.THREE_LETTER


def test_one_letter_without_termini_uses_defaults() -> None:
    pep = parse_sequence("RGDS", notation=ONE)
    assert pep.one_letter_sequence == "RGDS"
    assert pep.n_terminus.key == "free_amine"
    assert pep.c_terminus.key == "free_acid"
    assert "n_terminus_source" in pep.metadata
    assert "c_terminus_source" in pep.metadata


@pytest.mark.parametrize(
    ("text", "notation", "n_key", "c_key"),
    [
        ("Ac-RGDS-NH2", ONE, "acetyl", "amide"),
        ("NH2-RGDS-CO-NH2", ONE, "free_amine", "amide"),
        ("H-RGDS-OH", ONE, "free_amine", "free_acid"),
        ("H2N-RGDS-COOH", ONE, "free_amine", "free_acid"),
        ("Ac-RGDS-CO-NH-CH3", ONE, "acetyl", "n_methylamide"),
        ("Ac-Arg-Gly-Asp-Ser-NH2", THREE, "acetyl", "amide"),
        ("NH2-Arg-Gly-Asp-Ser-CO-NH2", THREE, "free_amine", "amide"),
        ("H-Arg-Gly-Asp-Ser-OH", THREE, "free_amine", "free_acid"),
        ("Ac-Arg-Gly-Asp-Ser-CO-NH-CH3", THREE, "acetyl", "n_methylamide"),
    ],
)
def test_termini_notations(text: str, notation: SequenceNotation, n_key: str, c_key: str) -> None:
    pep = parse_sequence(text, notation=notation)
    assert pep.one_letter_sequence == "RGDS"
    assert pep.n_terminus.key == n_key
    assert pep.c_terminus.key == c_key


def test_amide_is_not_n_methylamide() -> None:
    """-CO-NH2 et -CO-NH-CH3 sont deux extrémités distinctes (NHE ≠ NME en AMBER)."""
    amide = parse_sequence("RGDS-CO-NH2", notation=ONE)
    methylamide = parse_sequence("RGDS-CO-NH-CH3", notation=ONE)
    assert amide.c_terminus != methylamide.c_terminus


def test_three_letter_notation() -> None:
    pep = parse_sequence("Ac-Arg-Gly-Asp-Ser-NH2", notation=THREE)
    assert pep.three_letter_sequence == "ARG-GLY-ASP-SER"
    assert pep.composition() == {"ARG": 1, "GLY": 1, "ASP": 1, "SER": 1}


def test_whitespace_is_ignored() -> None:
    assert parse_sequence(" Ac - RGD S - NH2 ", notation=ONE).one_letter_sequence == "RGDS"


def test_histidine_alone_is_a_residue() -> None:
    pep = parse_sequence("H", notation=ONE)
    assert pep.three_letter_sequence == "HIS"


@pytest.mark.parametrize("notation", [ONE, THREE])
@pytest.mark.parametrize("text", ["", "   ", "RG--DS", "-RGDS"])
def test_malformed_sequences(text: str, notation: SequenceNotation) -> None:
    with pytest.raises(SequenceParseError):
        parse_sequence(text, notation=notation)


def test_unknown_residue() -> None:
    with pytest.raises(UnknownResidueError):
        parse_sequence("RGBX", notation=ONE)
    with pytest.raises(UnknownResidueError):
        parse_sequence("Arg-Xyz-Gly", notation=THREE)


def test_unknown_terminus_lookup() -> None:
    with pytest.raises(UnknownTerminusError):
        load_default_library().c_terminus("CO-SH")


# --- Mode de lecture déclaré au départ (revue du module 1) ----------------------


def test_the_reading_mode_has_no_default() -> None:
    """Rien n'est deviné : sans mode de lecture, l'appel lui-même est refusé."""
    parameter = inspect.signature(parse_sequence).parameters["notation"]
    assert parameter.default is inspect.Parameter.empty
    assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
    with pytest.raises(TypeError):
        parse_sequence("RGDS")  # type: ignore[call-arg]


@pytest.mark.parametrize(
    ("text", "notation", "codes"),
    [
        ("ALA", ONE, ("ALA", "LEU", "ALA")),
        ("Ac-ALA-NHMe", ONE, ("ALA", "LEU", "ALA")),
        ("H-CYS-OH", ONE, ("CYS", "TYR", "SER")),
        ("ALA", THREE, ("ALA",)),
        ("Ala", THREE, ("ALA",)),
        ("ala", THREE, ("ALA",)),
        ("Ac-ALA-NHMe", THREE, ("ALA",)),
        ("H-CYS-OH", THREE, ("CYS",)),
    ],
)
def test_the_declared_mode_decides_the_reading(
    text: str, notation: SequenceNotation, codes: tuple[str, ...]
) -> None:
    """`ALA` est Ala-Leu-Ala en code à une lettre, l'alanine en code à trois lettres."""
    peptide = parse_sequence(text, notation=notation)
    assert tuple(r.code for r in peptide.residues) == codes


@pytest.mark.parametrize(("notation", "value"), [(ONE, "one_letter"), (THREE, "three_letter")])
def test_the_declared_mode_is_recorded(notation: SequenceNotation, value: str) -> None:
    text = "RGDS" if notation is ONE else "Arg-Gly"
    assert parse_sequence(text, notation=notation).metadata[SEQUENCE_NOTATION] == value
    assert parse_sequence(text, notation=value).metadata[SEQUENCE_NOTATION] == value


@pytest.mark.parametrize("value", ["1", "auto", "ONE_LETTER", ""])
def test_an_unknown_mode_is_refused(value: str) -> None:
    with pytest.raises(SequenceParseError, match="Mode de lecture inconnu"):
        parse_sequence("RGDS", notation=value)


def test_three_letter_codes_ignore_case() -> None:
    peptide = parse_sequence("ac-arg-GLY-Asp-nh2", notation=THREE)
    assert peptide.three_letter_sequence == "ARG-GLY-ASP"
    assert (peptide.n_terminus.key, peptide.c_terminus.key) == ("acetyl", "amide")


@pytest.mark.parametrize(("text", "position"), [("RGDfK", 4), ("rgds", 1), ("Ac-RG[Asp]fK-NH2", 4)])
def test_one_letter_refuses_lower_case(text: str, position: int) -> None:
    """Une minuscule désigne d'ordinaire un acide aminé D : jamais lue comme un L."""
    with pytest.raises(SequenceParseError, match=f"en position {position} .*acide aminé D"):
        parse_sequence(text, notation=ONE)


def test_one_letter_refuses_several_blocks() -> None:
    with pytest.raises(SequenceParseError, match=r"un seul bloc.*Arg, Gly, Asp"):
        parse_sequence("Arg-Gly-Asp", notation=ONE)


@pytest.mark.parametrize(
    ("text", "error"),
    [("Ala", SequenceParseError), ("GLU", UnknownResidueError), ("Pro", SequenceParseError)],
)
def test_one_letter_error_points_to_the_three_letter_mode(
    text: str, error: type[Exception]
) -> None:
    with pytest.raises(error, match="choisir le code à trois lettres"):
        parse_sequence(text, notation=ONE)


@pytest.mark.parametrize("text", ["RGDS", "Ac-RGDS-NH2", "rgd"])
def test_three_letter_error_points_to_the_one_letter_mode(text: str) -> None:
    with pytest.raises(UnknownResidueError, match="code à une lettre, choisir ce mode"):
        parse_sequence(text, notation=THREE)


def test_a_one_letter_body_spelled_like_a_terminus_fails_loudly() -> None:
    """`AC` (Ala-Cys) écrit seul devant `-NH2` se lit comme l'extrémité acétylée.

    Limitation connue : la lecture échoue alors bruyamment, jamais en silence, et
    déclarer l'extrémité N suffit à lever l'ambiguïté.
    """
    with pytest.raises(UnknownResidueError):
        parse_sequence("AC-NH2", notation=ONE)
    peptide = parse_sequence("H-AC-NH2", notation=ONE)
    assert peptide.three_letter_sequence == "ALA-CYS"
    assert (peptide.n_terminus.key, peptide.c_terminus.key) == ("free_amine", "amide")


# --- Notation des non standards entre crochets (étape 1) ------------------------


def _library_with_aib() -> ResidueLibrary:
    """Bibliothèque par défaut plus un non standard, l'acide α-aminoisobutyrique."""
    library = load_default_library()
    residues = dict(library.residues)
    residues["AIB"] = Residue(code="AIB", one_letter=None, standard=False)
    return dataclasses.replace(library, residues=residues)


@pytest.mark.parametrize(
    ("text", "notation", "codes"),
    [
        ("Ac-RG[Aib]S-NH2", ONE, ("ARG", "GLY", "AIB", "SER")),
        ("[Aib]RG", ONE, ("AIB", "ARG", "GLY")),
        ("RG[Aib]", ONE, ("ARG", "GLY", "AIB")),
        ("[Aib][aib]", ONE, ("AIB", "AIB")),
        ("[Aib]", ONE, ("AIB",)),
        ("[Aib]", THREE, ("AIB",)),
        ("Aib", THREE, ("AIB",)),
        ("Ac-Arg-Gly-[Aib]-Ser-NH2", THREE, ("ARG", "GLY", "AIB", "SER")),
        ("Ac-Arg-Gly-Aib-Ser-NH2", THREE, ("ARG", "GLY", "AIB", "SER")),
        ("Ac - RG [ Aib ] S - NH2", ONE, ("ARG", "GLY", "AIB", "SER")),
    ],
)
def test_bracket_notation(text: str, notation: SequenceNotation, codes: tuple[str, ...]) -> None:
    peptide = parse_sequence(text, _library_with_aib(), notation=notation)
    assert tuple(r.code for r in peptide.residues) == codes


def test_bracketed_residue_keeps_its_library_entry() -> None:
    peptide = parse_sequence("Ac-RG[Aib]S-NH2", _library_with_aib(), notation=ONE)
    assert peptide.residue_at(3).standard is False
    assert peptide.one_letter_sequence == "RGXS"
    assert peptide.n_terminus.key == "acetyl"
    assert peptide.c_terminus.key == "amide"


def test_bracketed_standard_code_is_the_same_residue() -> None:
    assert parse_sequence("Ac-RG[Asp]S-NH2", notation=ONE) == parse_sequence(
        "Ac-RGDS-NH2", notation=ONE
    )


def test_both_notations_give_the_same_peptide() -> None:
    library = _library_with_aib()
    assert parse_sequence("Ac-RG[Aib]S-NH2", library, notation=ONE) == parse_sequence(
        "Ac-Arg-Gly-[Aib]-Ser-NH2", library, notation=THREE
    )


@pytest.mark.parametrize(
    ("text", "notation", "message"),
    [
        ("RG[AibS", ONE, "jamais fermé"),
        ("RGAib]S", ONE, "sans crochet ouvrant"),
        ("RG[]S", ONE, "Crochets vides"),
        ("RG[[Aib]]S", ONE, "imbriqués"),
        ("RG[Ai-b]S", ONE, "non refermé avant le tiret"),
        ("Ac-RG[Aib-NH2", ONE, "non refermé avant le tiret"),
        ("Arg-Gly[Aib]-Ser", THREE, "occupe un bloc entier"),
        ("Arg-[Aib][Aib]-Ser", THREE, "occupe un bloc entier"),
        ("Arg[Aib]", THREE, "occupe un bloc entier"),
    ],
)
def test_malformed_brackets(text: str, notation: SequenceNotation, message: str) -> None:
    with pytest.raises(SequenceParseError, match=message):
        parse_sequence(text, _library_with_aib(), notation=notation)


def test_bracketed_code_absent_from_the_library() -> None:
    """La bibliothèque par défaut ne contient que les vingt standards."""
    with pytest.raises(UnknownResidueError, match="Résidu inconnu : 'Aib'"):
        parse_sequence("Ac-RG[Aib]S-NH2", notation=ONE)


def test_residues_are_a_tuple() -> None:
    assert isinstance(parse_sequence("RGDS", notation=ONE).residues, tuple)


# --- Conversion entre les deux codes, à la manière de PyMOL (revue du module 1) --

_STANDARD_20_THREE = (
    "Ala-Cys-Asp-Glu-Phe-Gly-His-Ile-Lys-Leu-Met-Asn-Pro-Gln-Arg-Ser-Thr-Val-Trp-Tyr"
)

_PAIRS = [
    ("Ac-RGDS-NH2", "Ac-Arg-Gly-Asp-Ser-NH2"),
    ("H-RGDS-OH", "H-Arg-Gly-Asp-Ser-OH"),
    ("RGDS", "Arg-Gly-Asp-Ser"),
    ("Ac-RGDS", "Ac-Arg-Gly-Asp-Ser"),
    ("RGDS-NHMe", "Arg-Gly-Asp-Ser-NHMe"),
    ("Ac-ALA-NHMe", "Ac-Ala-Leu-Ala-NHMe"),
    ("Ac-A-NHMe", "Ac-Ala-NHMe"),
    ("ACDEFGHIKLMNPQRSTVWY", _STANDARD_20_THREE),
]


@pytest.mark.parametrize(("one_letter", "three_letter"), _PAIRS)
def test_conversion_in_both_directions(one_letter: str, three_letter: str) -> None:
    assert convert_sequence(one_letter, ONE, THREE) == three_letter
    assert convert_sequence(three_letter, THREE, ONE) == one_letter


@pytest.mark.parametrize(
    ("text", "source", "target", "expected"),
    [
        ("NH2-RGDS-CO-NH2", ONE, THREE, "H-Arg-Gly-Asp-Ser-NH2"),
        ("H2N-Arg-Gly-CO-OH", THREE, ONE, "H-RG-OH"),
        ("Ac-RGDS-CO-NH-CH3", ONE, ONE, "Ac-RGDS-NHMe"),
        ("ac-arg-GLY-nh2", THREE, THREE, "Ac-Arg-Gly-NH2"),
        (" Ac - RGD S - NH2 ", ONE, ONE, "Ac-RGDS-NH2"),
        ("Ac-RG[Asp]S-NH2", ONE, ONE, "Ac-RGDS-NH2"),
    ],
)
def test_conversion_writes_the_canonical_form(
    text: str, source: SequenceNotation, target: SequenceNotation, expected: str
) -> None:
    assert convert_sequence(text, source, target) == expected


def test_conversion_of_a_non_standard_residue() -> None:
    library = _library_with_aib()
    three = convert_sequence("Ac-RG[Aib]S-NH2", ONE, THREE, library)
    assert three == "Ac-Arg-Gly-[AIB]-Ser-NH2"
    assert convert_sequence(three, THREE, ONE, library) == "Ac-RG[AIB]S-NH2"


@pytest.mark.parametrize(
    ("text", "flags"),
    [
        ("RGDS", {"n_terminus_source", "c_terminus_source"}),
        ("Ac-RGDS", {"c_terminus_source"}),
        ("RGDS-NH2", {"n_terminus_source"}),
        ("Ac-RGDS-NH2", set()),
    ],
)
def test_an_undeclared_terminus_stays_undeclared(text: str, flags: set[str]) -> None:
    """Convertir ne doit pas transformer une extrémité supposée en extrémité choisie."""
    three = convert_sequence(text, ONE, THREE)
    reread = parse_sequence(three, notation=THREE)
    assert {k for k in ("n_terminus_source", "c_terminus_source") if k in reread.metadata} == flags


@pytest.mark.parametrize(
    ("text", "notation"),
    [(one, ONE) for one, _ in _PAIRS] + [(three, THREE) for _, three in _PAIRS],
)
@pytest.mark.parametrize("target", [ONE, THREE])
def test_round_trip_gives_back_the_same_peptide(
    text: str, notation: SequenceNotation, target: SequenceNotation
) -> None:
    peptide = parse_sequence(text, notation=notation)
    assert parse_sequence(format_sequence(peptide, target), notation=target) == peptide


def test_notation_other_and_labels() -> None:
    assert ONE.other is THREE
    assert THREE.other is ONE
    assert "monolitère" in ONE.label
    assert "trilitère" in THREE.label


@pytest.mark.parametrize(
    ("peptide", "lost"),
    [
        (
            parse_sequence("Ac-CRGDC-NH2", notation=ONE).with_link(Link(LinkKind.DISULFIDE, 1, 5)),
            "pont disulfure 1–5",
        ),
        (
            parse_sequence("H-RGDS-OH", notation=ONE).with_link(Link(LinkKind.HEAD_TO_TAIL, 1, 4)),
            "cyclisation tête-queue 1–4",
        ),
        (
            parse_sequence("Ac-RGDS-NH2", notation=ONE).with_modification(
                Modification("phosphorylation", 4)
            ),
            "phosphorylation en position 4",
        ),
        (
            parse_sequence("Ac-AHA-NH2", notation=ONE).with_protonation(2, "his_delta"),
            "his_delta en position 2",
        ),
    ],
)
@pytest.mark.parametrize("target", [ONE, THREE])
def test_writing_refuses_what_a_sequence_cannot_carry(
    peptide: Peptide, lost: str, target: SequenceNotation
) -> None:
    with pytest.raises(UnsupportedFeatureError, match=lost):
        format_sequence(peptide, target)


def test_writing_refuses_a_text_that_does_not_read_back() -> None:
    """His amidée sans extrémité N déclarée s'écrirait `H-NH2`, illisible : refus."""
    peptide = parse_sequence("His-NH2", notation=THREE)
    assert format_sequence(peptide, THREE) == "His-NH2"
    with pytest.raises(UnsupportedFeatureError, match="ne se relit pas"):
        format_sequence(peptide, ONE)


def test_writing_refuses_a_residue_the_library_would_read_differently() -> None:
    """Le texte produit est relu et comparé : un écart est refusé, jamais rendu."""
    odd_arginine = Residue(code="ARG", one_letter=None, standard=True)
    peptide = parse_sequence("Ac-RGDS-NH2", notation=ONE).mutate(1, odd_arginine)
    with pytest.raises(UnsupportedFeatureError, match="se relit comme un autre peptide"):
        format_sequence(peptide, THREE)


def test_writing_an_empty_peptide_is_refused() -> None:
    library = load_default_library()
    empty = Peptide(residues=(), n_terminus=library.n_terminus(), c_terminus=library.c_terminus())
    with pytest.raises(UnsupportedFeatureError, match="sans résidu"):
        format_sequence(empty, ONE)
