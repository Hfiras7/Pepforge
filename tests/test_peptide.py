"""Tests du modèle de données du peptide (étape 1 du plan de construction)."""

from __future__ import annotations

import copy
import dataclasses
import pickle
import re
from pathlib import Path

import pytest

from pepforge.core import (
    C_TERMINUS_DEFAULTED,
    N_TERMINUS_DEFAULTED,
    Link,
    LinkKind,
    Modification,
    Peptide,
    Residue,
    Terminus,
)
from pepforge.exceptions import PeptideModelError
from pepforge.io import parse_sequence
from pepforge.residues import load_default_library

CORE = Path(__file__).resolve().parents[1] / "src" / "pepforge" / "core"

# Noms de résidus, d'extrémités et de patches propres à un champ de force. Le
# plan exige les six premiers ; les suivants couvrent les variantes AMBER de
# protonation et les noms CHARMM, champ de force retenu pour la thèse.
FORCE_FIELD_NAMES = (
    "NHE",
    "NME",
    "HIE",
    "HID",
    "HIP",
    "CYX",
    "ACE",
    "CYM",
    "ASH",
    "GLH",
    "LYN",
    "HSD",
    "HSE",
    "HSP",
    "NTER",
    "CTER",
    "CT2",
    "CT3",
    "DISU",
    "ASPP",
    "GLUP",
    "LSN",
)
_FORBIDDEN = re.compile(r"\b(" + "|".join(FORCE_FIELD_NAMES) + r")\b", re.IGNORECASE)


def _forbidden_names(text: str) -> list[str]:
    return [match.group(0) for match in _FORBIDDEN.finditer(text)]


def _rgds() -> Peptide:
    return parse_sequence("Ac-RGDS-NH2")


def _cys_peptide() -> Peptide:
    """Ac-CRGDC-NH2 avec un pont disulfure entre les positions 1 et 5."""
    peptide = parse_sequence("Ac-CRGDC-NH2")
    return peptide.with_link(Link(LinkKind.DISULFIDE, 1, 5))


# --- Aucun nom de champ de force dans core -----------------------------------


def test_core_contains_no_force_field_name() -> None:
    offenders = {
        str(path.relative_to(CORE)): names
        for path in sorted(CORE.rglob("*.py"))
        if (names := _forbidden_names(path.read_text(encoding="utf-8")))
    }
    assert not offenders, f"Noms de champ de force dans core : {offenders}"


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("le résidu `NHE` d'AMBER", ["NHE"]),
        ("tautomère hie", ["hie"]),
        ("patch CT2", ["CT2"]),
        ("clé neutre his_delta", []),
        ("extrémité n_methylamide", []),
        ("une surface active", []),
    ],
)
def test_force_field_name_scan_has_teeth(text: str, expected: list[str]) -> None:
    """Le test précédent n'est pas vide de sens : il attrape bien ce qu'il cherche."""
    assert _forbidden_names(text) == expected


# --- Immuabilité ---------------------------------------------------------------


def test_peptide_fields_cannot_be_reassigned() -> None:
    peptide = _rgds()
    with pytest.raises(dataclasses.FrozenInstanceError):
        peptide.residues = ()  # type: ignore[misc]
    with pytest.raises(dataclasses.FrozenInstanceError):
        peptide.n_terminus = peptide.c_terminus  # type: ignore[misc]


def test_containers_are_immutable() -> None:
    peptide = _rgds()
    assert isinstance(peptide.residues, tuple)
    assert isinstance(peptide.links, tuple)
    assert isinstance(peptide.modifications, tuple)
    with pytest.raises(TypeError):
        peptide.metadata["input"] = "autre chose"  # type: ignore[index]


def test_lists_given_to_the_constructor_are_frozen() -> None:
    library = load_default_library()
    residues = [library.by_code("ARG"), library.by_code("GLY")]
    peptide = Peptide(
        residues=residues,  # type: ignore[arg-type]
        n_terminus=library.n_terminus(),
        c_terminus=library.c_terminus(),
    )
    residues.append(library.by_code("ASP"))
    assert isinstance(peptide.residues, tuple)
    assert len(peptide) == 2


def test_metadata_is_copied_at_construction() -> None:
    library = load_default_library()
    metadata = {"input": "RG"}
    peptide = Peptide(
        residues=(library.by_code("ARG"), library.by_code("GLY")),
        n_terminus=library.n_terminus(),
        c_terminus=library.c_terminus(),
        metadata=metadata,
    )
    metadata["input"] = "modifié après coup"
    assert peptide.metadata["input"] == "RG"


def test_peptide_survives_copy_pickle_and_asdict() -> None:
    """Une immuabilité qui empêcherait la copie ou la sérialisation serait un piège."""
    peptide = _cys_peptide().with_protonation(2, "arg_protonated")
    assert copy.deepcopy(peptide) == peptide
    restored = pickle.loads(pickle.dumps(peptide))
    assert restored == peptide
    assert restored.metadata == peptide.metadata
    assert dataclasses.asdict(peptide)["links"][0]["first"] == 1


def test_peptide_is_hashable() -> None:
    assert len({_rgds(), _rgds()}) == 1


# --- Égalité ---------------------------------------------------------------------


def test_notations_of_the_same_molecule_are_equal() -> None:
    one_letter = parse_sequence("Ac-RGDS-NH2")
    three_letters = parse_sequence("Ac-Arg-Gly-Asp-Ser-NH2")
    assert one_letter.metadata["input"] != three_letters.metadata["input"]
    assert one_letter == three_letters
    assert hash(one_letter) == hash(three_letters)


def test_declaration_order_does_not_change_equality() -> None:
    peptide = parse_sequence("Ac-CCGCC-NH2")
    first = peptide.with_link(Link(LinkKind.DISULFIDE, 1, 4)).with_link(
        Link(LinkKind.DISULFIDE, 2, 5)
    )
    second = peptide.with_link(Link(LinkKind.DISULFIDE, 5, 2)).with_link(
        Link(LinkKind.DISULFIDE, 4, 1)
    )
    assert first == second
    assert first.disulfides == ((1, 4), (2, 5))


def test_protonation_is_part_of_identity() -> None:
    peptide = parse_sequence("Ac-AHA-NH2")
    assert peptide.with_protonation(2, "his_delta") != peptide.with_protonation(2, "his_epsilon")


# --- Positions et mutation -----------------------------------------------------


def test_positions_start_at_one() -> None:
    peptide = _rgds()
    assert peptide.residue_at(1).code == "ARG"
    assert peptide.residue_at(4).code == "SER"


@pytest.mark.parametrize("position", [0, -1, 5, True])
def test_out_of_range_position_is_refused(position: int) -> None:
    with pytest.raises(PeptideModelError, match="positions valides de 1 à 4"):
        _rgds().residue_at(position)


def test_mutation_returns_a_new_peptide_and_keeps_the_original() -> None:
    original = _rgds()
    alanine = load_default_library().by_code("ALA")
    mutated = original.mutate(3, alanine)
    assert mutated is not original
    assert mutated.one_letter_sequence == "RGAS"
    assert original.one_letter_sequence == "RGDS"
    assert mutated.n_terminus == original.n_terminus
    assert mutated.c_terminus == original.c_terminus
    assert mutated.metadata == original.metadata


def test_mutation_out_of_range_is_refused() -> None:
    with pytest.raises(PeptideModelError):
        _rgds().mutate(5, load_default_library().by_code("ALA"))


def test_mutation_does_not_carry_the_old_protonation() -> None:
    library = load_default_library()
    peptide = parse_sequence("Ac-AHA-NH2").with_protonation(2, "his_delta")
    mutated = peptide.mutate(2, library.by_code("HIS"))
    assert mutated.residue_at(2).protonation is None


def test_mutation_is_refused_on_a_disulfide_cysteine() -> None:
    serine = load_default_library().by_code("SER")
    with pytest.raises(PeptideModelError, match=r"pont disulfure 1–5"):
        _cys_peptide().mutate(5, serine)


def test_mutation_is_refused_on_a_modified_residue() -> None:
    peptide = _rgds().with_modification(Modification("phosphorylation", 4))
    with pytest.raises(PeptideModelError, match="phosphorylation"):
        peptide.mutate(4, load_default_library().by_code("ALA"))


def test_mutation_is_allowed_after_removing_the_declaration() -> None:
    link = Link(LinkKind.DISULFIDE, 1, 5)
    serine = load_default_library().by_code("SER")
    peptide = _cys_peptide().without_link(link).mutate(5, serine)
    assert peptide.links == ()
    assert peptide.one_letter_sequence == "CRGDS"


def test_head_to_tail_does_not_block_mutating_its_ends() -> None:
    library = load_default_library()
    cyclic = parse_sequence("H-RGDS-OH").with_link(Link(LinkKind.HEAD_TO_TAIL, 1, 4))
    mutated = cyclic.mutate(1, library.by_code("LYS"))
    assert mutated.is_cyclic
    assert mutated.one_letter_sequence == "KGDS"


# --- Extrémités ------------------------------------------------------------------


def test_changing_a_terminus_makes_it_declared() -> None:
    library = load_default_library()
    defaulted = parse_sequence("RGDS")
    assert N_TERMINUS_DEFAULTED in defaulted.metadata
    assert C_TERMINUS_DEFAULTED in defaulted.metadata
    capped = defaulted.with_n_terminus(library.n_terminus("Ac")).with_c_terminus(
        library.c_terminus("NH2")
    )
    assert capped.n_terminus.key == "acetyl"
    assert capped.c_terminus.key == "amide"
    assert N_TERMINUS_DEFAULTED not in capped.metadata
    assert C_TERMINUS_DEFAULTED not in capped.metadata
    # L'original n'a pas bougé.
    assert defaulted.n_terminus.key == "free_amine"
    assert N_TERMINUS_DEFAULTED in defaulted.metadata


def test_setting_one_terminus_leaves_the_other_flagged() -> None:
    peptide = parse_sequence("RGDS").with_n_terminus(load_default_library().n_terminus("Ac"))
    assert N_TERMINUS_DEFAULTED not in peptide.metadata
    assert C_TERMINUS_DEFAULTED in peptide.metadata


def test_terminus_key_must_be_internal() -> None:
    with pytest.raises(PeptideModelError, match="clé interne mal formée"):
        Terminus(key="Free Amine", label="amine libre")


# --- Protonation -----------------------------------------------------------------


def test_protonation_is_set_and_cleared() -> None:
    peptide = parse_sequence("Ac-AHA-NH2")
    protonated = peptide.with_protonation(2, "his_protonated")
    assert protonated.residue_at(2).protonation == "his_protonated"
    assert peptide.residue_at(2).protonation is None
    assert protonated.with_protonation(2, None) == peptide


@pytest.mark.parametrize("key", ["HisDelta", "his-delta", "1his", "", "his delta"])
def test_protonation_key_must_be_internal(key: str) -> None:
    with pytest.raises(PeptideModelError, match="clé interne mal formée"):
        parse_sequence("Ac-AHA-NH2").with_protonation(2, key)


# --- Résidus ---------------------------------------------------------------------


@pytest.mark.parametrize("code", ["", "Arg", "A-B", "[AIB]", "AI B"])
def test_residue_code_must_be_well_formed(code: str) -> None:
    with pytest.raises(PeptideModelError, match="Code de résidu mal formé"):
        Residue(code=code)


@pytest.mark.parametrize("letter", ["", "r", "RG", "1"])
def test_one_letter_code_must_be_one_capital(letter: str) -> None:
    with pytest.raises(PeptideModelError, match="code à une lettre"):
        Residue(code="ARG", one_letter=letter)


# --- Liaisons ----------------------------------------------------------------------


def test_link_positions_are_put_in_order() -> None:
    link = Link(LinkKind.DISULFIDE, 7, 2)
    assert link.positions == (2, 7)
    assert link == Link(LinkKind.DISULFIDE, 2, 7)


def test_link_kind_given_as_text_is_converted() -> None:
    link = Link("disulfide", 1, 2)  # type: ignore[arg-type]
    assert link.kind is LinkKind.DISULFIDE


def test_unknown_link_kind_is_refused() -> None:
    with pytest.raises(PeptideModelError, match="Type de liaison inconnu"):
        Link("thioether", 1, 2)  # type: ignore[arg-type]


@pytest.mark.parametrize(("first", "second"), [(0, 2), (2, -1), (3, 3)])
def test_invalid_link_positions_are_refused(first: int, second: int) -> None:
    with pytest.raises(PeptideModelError):
        Link(LinkKind.DISULFIDE, first, second)


def test_link_beyond_the_chain_is_refused() -> None:
    with pytest.raises(PeptideModelError, match="hors du peptide"):
        _rgds().with_link(Link(LinkKind.DISULFIDE, 2, 5))


def test_head_to_tail_must_join_the_two_ends() -> None:
    with pytest.raises(PeptideModelError, match="relie la position 1"):
        _rgds().with_link(Link(LinkKind.HEAD_TO_TAIL, 1, 3))
    assert _rgds().with_link(Link(LinkKind.HEAD_TO_TAIL, 1, 4)).is_cyclic


def test_a_position_belongs_to_one_disulfide_only() -> None:
    peptide = parse_sequence("Ac-CCC-NH2").with_link(Link(LinkKind.DISULFIDE, 1, 2))
    with pytest.raises(PeptideModelError, match="deux ponts disulfure"):
        peptide.with_link(Link(LinkKind.DISULFIDE, 2, 3))


def test_duplicate_link_is_refused() -> None:
    with pytest.raises(PeptideModelError, match="déclarée deux fois"):
        _cys_peptide().with_link(Link(LinkKind.DISULFIDE, 5, 1))


def test_removing_an_undeclared_link_is_refused() -> None:
    with pytest.raises(PeptideModelError, match="n'est pas déclarée"):
        _rgds().without_link(Link(LinkKind.DISULFIDE, 1, 2))


# --- Modifications -------------------------------------------------------------------


def test_modifications_are_declared_and_removed() -> None:
    modification = Modification("phosphorylation", 4)
    modified = _rgds().with_modification(modification)
    assert modified.modifications == (modification,)
    assert modified.without_modification(modification) == _rgds()


def test_modification_beyond_the_chain_is_refused() -> None:
    with pytest.raises(PeptideModelError, match="hors du peptide"):
        _rgds().with_modification(Modification("phosphorylation", 5))


def test_duplicate_modification_is_refused() -> None:
    peptide = _rgds().with_modification(Modification("phosphorylation", 4))
    with pytest.raises(PeptideModelError, match="déclarée deux fois"):
        peptide.with_modification(Modification("phosphorylation", 4))


def test_removing_an_undeclared_modification_is_refused() -> None:
    with pytest.raises(PeptideModelError, match="n'est pas déclarée"):
        _rgds().without_modification(Modification("phosphorylation", 4))


@pytest.mark.parametrize(("key", "position"), [("Phospho", 1), ("phosphorylation", 0)])
def test_invalid_modification_is_refused(key: str, position: int) -> None:
    with pytest.raises(PeptideModelError):
        Modification(key, position)


def test_modifications_are_sorted() -> None:
    peptide = (
        _rgds()
        .with_modification(Modification("phosphorylation", 4))
        .with_modification(Modification("methylation", 1))
    )
    assert [m.position for m in peptide.modifications] == [1, 4]
