"""Tests des groupes ionisables, du choix des états et de la charge nette (étape 2)."""

from __future__ import annotations

import copy
import math
import re
from typing import Any

import pytest

from pepforge.core import (
    IonizationForm,
    Link,
    LinkKind,
    Modification,
    Peptide,
    Residue,
    Terminus,
)
from pepforge.exceptions import (
    ChemistryError,
    ConfigurationError,
    PeptideModelError,
    ReferenceDataError,
    UnsupportedFeatureError,
)
from pepforge.io.sequence import parse_sequence
from pepforge.residues import (
    acid_fraction,
    assign_protonation,
    load_default_library,
    load_ionizable_groups,
    parse_ionizable_groups,
    theoretical_net_charge,
)
from pepforge.residues.library import read_data_file
from pepforge.residues.protonation import DATA_FILE, check_ph, major_form

INTERNAL_KEY = re.compile(r"[a-z][a-z0-9_]*")


def _raw() -> dict[str, Any]:
    return copy.deepcopy(read_data_file(DATA_FILE))


def _group(data: dict[str, Any], key: str) -> dict[str, Any]:
    return next(g for g in data["groups"] if g["key"] == key)


# --- Le fichier de données ----------------------------------------------------


def test_data_file_loads_and_declares_its_sources() -> None:
    table = load_ionizable_groups()
    assert table.version == 1
    assert table.near_pka_window == 1.0
    for group in [*table.by_residue.values(), *table.by_terminus.values()]:
        assert group.pka_source in table.sources
        assert table.sources[group.pka_source].strip()


def test_every_standard_residue_and_terminus_is_classified() -> None:
    """Aucun résidu ni extrémité n'est supposé neutre faute d'être cité."""
    table = load_ionizable_groups()
    library = load_default_library()
    residues = set(table.by_residue) | table.non_ionizable_residues
    assert residues == set(library.residues)
    termini = set(table.by_terminus) | table.non_ionizable_termini
    assert termini == set(library.n_termini) | set(library.c_termini)


def test_reference_pka_values_are_those_of_the_cited_table() -> None:
    """Valeurs du tableau I de Grimsley 2009, colonne des pentapeptides d'alanine."""
    table = load_ionizable_groups()
    expected = {"ASP": 3.9, "GLU": 4.3, "HIS": 6.5, "CYS": 8.6, "TYR": 9.8, "LYS": 10.4}
    assert {code: table.by_residue[code].pka for code in expected} == expected
    assert table.by_terminus["free_amine"].pka == 8.0
    assert table.by_terminus["free_acid"].pka == 3.7


def test_arginine_pka_is_missing_until_validated() -> None:
    assert load_ionizable_groups().by_residue["ARG"].pka is None


def test_state_keys_are_neutral_internal_keys() -> None:
    table = load_ionizable_groups()
    for group in [*table.by_residue.values(), *table.by_terminus.values()]:
        for state in group.states:
            assert INTERNAL_KEY.fullmatch(state.key), state.key


def test_histidine_keys_match_those_announced_in_module_1() -> None:
    his = load_ionizable_groups().by_residue["HIS"]
    assert {s.key for s in his.states} == {"his_delta", "his_epsilon", "his_protonated"}


def test_protonated_form_carries_one_more_charge() -> None:
    table = load_ionizable_groups()
    for group in [*table.by_residue.values(), *table.by_terminus.values()]:
        acid = group.charge_of(IonizationForm.ACID)
        assert acid - group.charge_of(IonizationForm.BASE) == 1


def _mutate_missing_source(d: dict[str, Any]) -> None:
    _group(d, "asp_side_chain")["pka_source"] = "inconnue"


def _mutate_both_sites(d: dict[str, Any]) -> None:
    _group(d, "asp_side_chain")["terminus"] = "free_acid"


def _mutate_no_pka_field(d: dict[str, Any]) -> None:
    del _group(d, "lys_side_chain")["pka"]


def _mutate_text_pka(d: dict[str, Any]) -> None:
    _group(d, "lys_side_chain")["pka"] = "10,4"


def _mutate_tautomer_charge(d: dict[str, Any]) -> None:
    _group(d, "his_side_chain")["states"][1]["charge"] = 1


def _mutate_two_protons(d: dict[str, Any]) -> None:
    _group(d, "lys_side_chain")["states"][0]["charge"] = 2


def _mutate_no_base(d: dict[str, Any]) -> None:
    _group(d, "glu_side_chain")["states"][1]["form"] = "acid"


def _mutate_unknown_form(d: dict[str, Any]) -> None:
    _group(d, "glu_side_chain")["states"][1]["form"] = "zwitterion"


def _mutate_duplicate_state(d: dict[str, Any]) -> None:
    _group(d, "glu_side_chain")["states"][0]["key"] = "asp_protonated"


def _mutate_duplicate_residue(d: dict[str, Any]) -> None:
    _group(d, "glu_side_chain")["residue"] = "ASP"


def _mutate_both_classes(d: dict[str, Any]) -> None:
    d["non_ionizable_residues"].append("LYS")


def _mutate_window(d: dict[str, Any]) -> None:
    d["near_pka_window"]["value"] = 0


def _mutate_window_unjustified(d: dict[str, Any]) -> None:
    del d["near_pka_window"]["justification"]


def _mutate_float_charge(d: dict[str, Any]) -> None:
    _group(d, "lys_side_chain")["states"][0]["charge"] = 1.0


@pytest.mark.parametrize(
    "mutation",
    [
        _mutate_missing_source,
        _mutate_both_sites,
        _mutate_no_pka_field,
        _mutate_text_pka,
        _mutate_tautomer_charge,
        _mutate_two_protons,
        _mutate_no_base,
        _mutate_unknown_form,
        _mutate_duplicate_state,
        _mutate_duplicate_residue,
        _mutate_both_classes,
        _mutate_window,
        _mutate_window_unjustified,
        _mutate_float_charge,
    ],
)
def test_inconsistent_data_file_is_refused(mutation: Any) -> None:
    data = _raw()
    mutation(data)
    with pytest.raises(ReferenceDataError):
        parse_ionizable_groups(data)


# --- Fonctions élémentaires ---------------------------------------------------


@pytest.mark.parametrize(
    ("pka", "ph", "expected"),
    [(7.0, 7.0, 0.5), (4.0, 5.0, 1 / 11), (4.0, 3.0, 10 / 11), (4.0, 6.0, 1 / 101)],
)
def test_acid_fraction_follows_henderson_hasselbalch(
    pka: float, ph: float, expected: float
) -> None:
    assert acid_fraction(pka, ph) == pytest.approx(expected)


def test_major_form() -> None:
    assert major_form(6.5, 5.0) is IonizationForm.ACID
    assert major_form(6.5, 7.4) is IonizationForm.BASE
    assert major_form(6.5, 6.5) is None


@pytest.mark.parametrize("ph", [-0.1, 14.1, 74.0, math.nan, math.inf, "7.4", True, None])
def test_invalid_ph_is_refused(ph: object) -> None:
    with pytest.raises(ConfigurationError):
        check_ph(ph)


@pytest.mark.parametrize("ph", [0, 0.0, 7.4, 14])
def test_valid_ph_is_accepted(ph: float) -> None:
    assert check_ph(ph) == float(ph)


# --- Extrémités ---------------------------------------------------------------


def test_blocked_termini_carry_no_charge() -> None:
    result = assign_protonation(parse_sequence("Ac-GGG-NH2"), 7.4)
    assert result.assignments == ()
    assert result.net_charge == 0
    assert result.theoretical_net_charge == 0.0
    assert result.peptide.n_terminus.protonation is None
    assert result.peptide.c_terminus.protonation is None


def test_n_methylamide_carries_no_charge() -> None:
    result = assign_protonation(parse_sequence("Ac-GGG-NHMe"), 7.4)
    assert result.net_charge == 0


@pytest.mark.parametrize(
    ("ph", "n_state", "c_state", "charge"),
    [
        (2.0, "nterm_protonated", "cterm_protonated", 1),
        (7.4, "nterm_protonated", "cterm_deprotonated", 0),
        (11.0, "nterm_deprotonated", "cterm_deprotonated", -1),
    ],
)
def test_free_termini_follow_the_ph(ph: float, n_state: str, c_state: str, charge: int) -> None:
    result = assign_protonation(parse_sequence("H-GGG-OH"), ph)
    assert result.peptide.n_terminus.protonation == n_state
    assert result.peptide.c_terminus.protonation == c_state
    assert result.net_charge == charge


def test_terminus_protonation_keeps_the_defaulted_mark() -> None:
    peptide = parse_sequence("GGG")
    assert "n_terminus_source" in peptide.metadata
    result = assign_protonation(peptide, 7.4)
    assert "n_terminus_source" in result.peptide.metadata
    assert "c_terminus_source" in result.peptide.metadata


# --- Chaînes latérales --------------------------------------------------------


@pytest.mark.parametrize(
    ("code", "ph", "state", "charge"),
    [
        ("D", 3.0, "asp_protonated", 0),
        ("D", 7.4, "asp_deprotonated", -1),
        ("E", 3.0, "glu_protonated", 0),
        ("E", 7.4, "glu_deprotonated", -1),
        ("H", 5.0, "his_protonated", 1),
        ("C", 7.4, "cys_protonated", 0),
        ("C", 10.0, "cys_deprotonated", -1),
        ("Y", 7.4, "tyr_protonated", 0),
        ("Y", 12.0, "tyr_deprotonated", -1),
        ("K", 7.4, "lys_protonated", 1),
        ("K", 12.0, "lys_deprotonated", 0),
    ],
)
def test_each_side_chain_group(code: str, ph: float, state: str, charge: int) -> None:
    result = assign_protonation(parse_sequence(f"Ac-G{code}G-NH2"), ph)
    assert result.peptide.residue_at(2).protonation == state
    (assignment,) = result.assignments
    assert (assignment.site, assignment.state, assignment.charge) == ("2", state, charge)
    assert not assignment.imposed
    assert result.net_charge == charge


def test_non_ionizable_residues_get_no_state() -> None:
    result = assign_protonation(parse_sequence("Ac-ANQGILMFPSTWV-NH2"), 7.4)
    assert result.assignments == ()
    assert all(r.protonation is None for r in result.peptide.residues)


def test_neutral_histidine_without_tautomer_is_refused() -> None:
    with pytest.raises(ChemistryError, match="his_delta ou his_epsilon"):
        assign_protonation(parse_sequence("Ac-GHG-NH2"), 7.4)


@pytest.mark.parametrize("tautomer", ["his_delta", "his_epsilon"])
def test_imposed_histidine_tautomer_is_kept(tautomer: str) -> None:
    peptide = parse_sequence("Ac-GHG-NH2").with_protonation(2, tautomer)
    result = assign_protonation(peptide, 7.4)
    (assignment,) = result.assignments
    assert assignment.state == tautomer
    assert assignment.imposed
    assert result.net_charge == 0


def test_ph_equal_to_pka_requires_an_imposed_state() -> None:
    with pytest.raises(ChemistryError, match="égal au pKa"):
        assign_protonation(parse_sequence("Ac-GKG-NH2"), 10.4)
    imposed = parse_sequence("Ac-GKG-NH2").with_protonation(2, "lys_protonated")
    assert assign_protonation(imposed, 10.4).assignments[0].major_form is None


def test_imposed_minority_state_is_kept() -> None:
    peptide = parse_sequence("Ac-GDG-NH2").with_protonation(2, "asp_protonated")
    (assignment,) = assign_protonation(peptide, 7.4).assignments
    assert assignment.state == "asp_protonated"
    assert assignment.imposed
    assert assignment.form is IonizationForm.ACID
    assert assignment.major_form is IonizationForm.BASE


def test_state_of_another_group_is_refused() -> None:
    peptide = parse_sequence("Ac-GDG-NH2").with_protonation(2, "his_delta")
    with pytest.raises(ChemistryError, match="his_delta"):
        assign_protonation(peptide, 7.4)


def test_state_on_a_non_ionizable_residue_is_refused() -> None:
    """Limite 1 du module 1 : `his_delta` sur une alanine est maintenant refusé."""
    peptide = parse_sequence("Ac-GAG-NH2").with_protonation(2, "his_delta")
    with pytest.raises(ChemistryError, match="sans groupe ionisable"):
        assign_protonation(peptide, 7.4)


def test_state_on_a_blocked_terminus_is_refused() -> None:
    peptide = parse_sequence("Ac-GGG-NH2").with_n_terminus_protonation("nterm_protonated")
    with pytest.raises(ChemistryError, match="extrémité"):
        assign_protonation(peptide, 7.4)


def test_disulfide_cysteines_are_not_ionizable() -> None:
    peptide = parse_sequence("Ac-CGC-NH2").with_link(Link(LinkKind.DISULFIDE, 1, 3))
    result = assign_protonation(peptide, 10.0)
    assert result.assignments == ()
    assert result.net_charge == 0


def test_state_on_a_bridged_cysteine_is_refused() -> None:
    peptide = (
        parse_sequence("Ac-CGC-NH2")
        .with_link(Link(LinkKind.DISULFIDE, 1, 3))
        .with_protonation(1, "cys_protonated")
    )
    with pytest.raises(ChemistryError, match="pont disulfure"):
        assign_protonation(peptide, 7.4)


def test_arginine_is_refused_while_its_pka_is_missing() -> None:
    with pytest.raises(ReferenceDataError, match="arg_side_chain"):
        assign_protonation(parse_sequence("Ac-RGDS-NH2"), 7.4)
    with pytest.raises(ReferenceDataError):
        theoretical_net_charge(parse_sequence("Ac-RGDS-NH2"), 7.4)


def test_unclassified_residue_is_refused() -> None:
    peptide = parse_sequence("Ac-GGG-NH2").mutate(2, Residue("AIB", standard=False))
    with pytest.raises(ReferenceDataError, match="AIB"):
        assign_protonation(peptide, 7.4)


def test_cyclic_peptide_is_refused() -> None:
    peptide = parse_sequence("H-GKG-OH").with_link(Link(LinkKind.HEAD_TO_TAIL, 1, 3))
    with pytest.raises(UnsupportedFeatureError, match="cyclisation"):
        assign_protonation(peptide, 7.4)


def test_modified_peptide_is_refused() -> None:
    peptide = parse_sequence("Ac-GSG-NH2").with_modification(Modification("phosphorylation", 2))
    with pytest.raises(UnsupportedFeatureError, match="phosphorylation"):
        assign_protonation(peptide, 7.4)


# --- Résultat -----------------------------------------------------------------


def test_input_peptide_is_unchanged() -> None:
    peptide = parse_sequence("H-KDG-OH")
    before = copy.deepcopy(peptide)
    result = assign_protonation(peptide, 7.4)
    assert peptide == before
    assert result.peptide != peptide


def test_assignments_follow_the_chain_order() -> None:
    result = assign_protonation(parse_sequence("H-KDG-OH"), 7.4)
    assert [a.site for a in result.assignments] == ["N", "1", "2", "C"]
    assert [a.position for a in result.assignments] == [None, 1, 2, None]


def test_assignment_records_pka_and_source() -> None:
    (assignment,) = assign_protonation(parse_sequence("Ac-GKG-NH2"), 7.4).assignments
    assert assignment.pka == 10.4
    assert assignment.pka_source == "grimsley2009_table1"
    assert assignment.group == "lys_side_chain"


def test_result_names_its_data_set() -> None:
    result = assign_protonation(parse_sequence("Ac-GKG-NH2"), 7.4)
    assert result.data_source == "ionizable_groups.yaml v1 (pace_alanine_pentapeptides_v1)"
    assert result.near_pka_window == 1.0


def test_theoretical_charge_ignores_imposed_states() -> None:
    plain = parse_sequence("Ac-GDG-NH2")
    imposed = plain.with_protonation(2, "asp_protonated")
    assert theoretical_net_charge(plain, 7.4) == theoretical_net_charge(imposed, 7.4)


def test_protonated_peptide_can_be_reassigned_identically() -> None:
    """Relancer l'étape sur sa propre sortie ne change rien (états vus comme imposés)."""
    first = assign_protonation(parse_sequence("H-KDG-OH"), 7.4)
    second = assign_protonation(first.peptide, 7.4)
    assert second.peptide == first.peptide
    assert second.net_charge == first.net_charge


def test_terminus_protonation_is_part_of_identity() -> None:
    peptide = parse_sequence("H-GGG-OH")
    assert peptide.with_n_terminus_protonation("nterm_protonated") != peptide
    assert isinstance(peptide.with_c_terminus_protonation(None), Peptide)


def test_terminus_protonation_key_must_be_neutral() -> None:
    with pytest.raises(PeptideModelError):
        Terminus("free_amine", "amine libre", protonation="NH3+")
