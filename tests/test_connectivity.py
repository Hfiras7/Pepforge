"""Tests de l'étape C : liaisons covalentes déclarées et possibles, composants rattachés.

Toutes les structures sont construites ici, à géométrie connue : chaque distance
est choisie, donc chaque classement est attendu exactement. Rien n'est tiré d'un
fichier réel (ceux-là sont dans `validation/etape_c_connectivity.py`).

Cas couverts, comme demandé par Firas le 28/11/2026 : liaisons intra-chaîne et
entre chaînes, chélateurs (N-terminal, C-terminal, chaîne latérale), linkers,
ligands non liés, ions, liaisons déclarées (LINK, CONECT) ou non, cas ambigus
(bornes, conformations, recouvrement), et prévention des pertes silencieuses.
"""

from __future__ import annotations

import ast
import math
from dataclasses import replace
from itertools import pairwise
from pathlib import Path

import pytest
import yaml

from pepforge.core.peptide import LinkKind
from pepforge.core.structure import Atom, Connection, ConnectionKind, ResidueId, Structure
from pepforge.exceptions import ConfigurationError, ExtractionRefusedError, ReferenceDataError
from pepforge.io import (
    AttachmentSite,
    BondEvidence,
    BondStatus,
    Certainty,
    Decision,
    ExtractionDecisions,
    analyze_chain,
    extract_peptides,
    read_structure,
    write_structure,
)
from pepforge.io.connectivity import assess_declared, find_contacts
from pepforge.io.peptide_extraction import STRUCTURE_SET_ASIDE, ComponentCategory
from pepforge.io.structure_formats import CONECT_UNRESOLVED
from pepforge.report.journal import Journal
from pepforge.residues.structure_components import CovalentCriteria, load_covalent_criteria

CRITERIA = load_covalent_criteria()
STEP = 3.33
"""Pas entre deux résidus : N(i+1) à 1,33 Å du C(i)."""


# --- Construction à géométrie connue --------------------------------------------------------


def _x(number: int) -> float:
    return (number - 1) * STEP


def _chain(
    names: str = "ALA LYS GLY CYS SER",
    chain: str = "A",
    side: dict[int, tuple[str, str]] | None = None,
) -> list[Atom]:
    """Résidus 1..n le long de x. Liaison C(i)–N(i+1) de 1,33 Å, O(i)–N(i+1) de 2,25 Å.

    `side` : numéro → (atome, élément) d'une chaîne latérale, placé en (x+1, −3, 0).
    """
    atoms: list[Atom] = []
    serial = 0
    side = side or {}
    for number, name in enumerate(names.split(), start=1):
        x = _x(number)
        positions = {
            "N": (x, 0.0),
            "CA": (x + 1.0, 0.5),
            "C": (x + 2.0, 0.0),
            "O": (x + 1.33, 1.03),
        }
        for atom, (px, py) in positions.items():
            serial += 1
            atoms.append(Atom(serial, atom, atom[0], name, number, chain, px, py, 0.0))
        if number in side:
            atom, element = side[number]
            serial += 1
            atoms.append(Atom(serial, atom, element, name, number, chain, x + 1.0, -3.0, 0.0))
    return atoms


def _het(
    name: str,
    number: int,
    start: tuple[float, float, float],
    direction: tuple[float, float, float],
    count: int = 3,
    chain: str = "A",
    element: str = "C",
    serial: int = 900,
) -> list[Atom]:
    """Composant HETATM : `count` atomes alignés, à 1,5 Å l'un de l'autre."""
    return [
        Atom(
            serial + k,
            f"{element}{k + 1}",
            element,
            name,
            number,
            chain,
            start[0] + 1.5 * k * direction[0],
            start[1] + 1.5 * k * direction[1],
            start[2] + 1.5 * k * direction[2],
            hetero=True,
        )
        for k in range(count)
    ]


def _rid(number: int, name: str, chain: str = "A") -> ResidueId:
    return ResidueId(chain, number, "", name)


def _component(report: object, label: str) -> object:
    (found,) = [c for c in report.components if c.residue == label]  # type: ignore[attr-defined]
    return found


LYS_NZ = {2: ("NZ", "N")}
NZ = (_x(2) + 1.0, -3.0, 0.0)
"""Position de l'azote NZ de la lysine 2."""

BELOW_NZ = (NZ[0], NZ[1] - 1.45, 0.0)
"""1,45 Å sous NZ : distance d'une liaison C–N (fenêtre N–C : 0,63 à 1,92 Å)."""


# --- Critères : données sourcées, jamais supposées -----------------------------------------


def test_covalent_criteria_come_from_the_data_file_with_their_sources() -> None:
    assert CRITERIA.radii["C"] == 0.76 and CRITERIA.radii["N"] == 0.71
    assert CRITERIA.radii["S"] == 1.05 and CRITERIA.radii["O"] == 0.66
    assert CRITERIA.tolerance == 0.45 and CRITERIA.min_distance == 0.63
    assert CRITERIA.max_distance("N", "C") == pytest.approx(1.92)
    assert CRITERIA.max_distance("S", "S") == pytest.approx(2.55)
    assert "10.1039/B801115J" in CRITERIA.radii_source
    assert "Open Babel" in CRITERIA.tolerance_source
    assert "0.63" in CRITERIA.min_distance_source
    assert not CRITERIA.metals & set(CRITERIA.radii)
    assert {"ZN", "GA", "LU", "RE", "TC", "CU", "NA"} <= CRITERIA.metals


def test_an_element_without_a_sourced_radius_is_refused_not_guessed() -> None:
    with pytest.raises(ReferenceDataError, match="aucun rayon covalent sourcé"):
        CRITERIA.max_distance("C", "X")


# --- Lecture : CONECT, provenance, fusion avec LINK ----------------------------------------


PDB = """\
{link}ATOM      1  N   LYS A   1       0.000   0.000   0.000  1.00 10.00           N
ATOM      2  CA  LYS A   1       1.000   0.500   0.000  1.00 10.00           C
ATOM      3  C   LYS A   1       2.000   0.000   0.000  1.00 10.00           C
ATOM      4  O   LYS A   1       1.330   1.030   0.000  1.00 10.00           O
ATOM      5  NZ  LYS A   1       1.000  -3.000   0.000  1.00 10.00           N
HETATM    6  C1  DOT A 101       1.000  -4.450   0.000  1.00 10.00           C
HETATM    7  C2  DOT A 101       1.000  -5.950   0.000  1.00 10.00           C
{conect}END
"""
LINK = "LINK         NZ  LYS A   1                 C1  DOT A 101     1555   1555  1.45  \n"
CONECT = "CONECT    5    6\nCONECT    6    5    7\nCONECT    7    6\n"


def _pdb(tmp_path: Path, link: str = "", conect: str = "", name: str = "x.pdb") -> Path:
    path = tmp_path / name
    path.write_text(PDB.format(link=link, conect=conect), encoding="utf-8")
    return path


def test_conect_alone_is_read_as_a_declared_bond_of_unknown_kind(tmp_path: Path) -> None:
    journal = Journal()
    structure = read_structure(_pdb(tmp_path, conect=CONECT), journal=journal)
    (bond,) = structure.connections
    assert (bond.first.label(), bond.first_atom, bond.second.label(), bond.second_atom) == (
        "A:LYS1",
        "NZ",
        "A:DOT101",
        "C1",
    )
    assert bond.kind == ConnectionKind.UNKNOWN and bond.records == ("CONECT",)
    read = journal.entries[0].parameters
    assert read["conect_pairs"] == 2 and read["conect_internal_bonds"] == 1
    assert read["conect_unresolved"] == []


def test_link_and_conect_of_the_same_bond_give_one_bond_with_both_sources(
    tmp_path: Path,
) -> None:
    structure = read_structure(_pdb(tmp_path, link=LINK, conect=CONECT), journal=Journal())
    (bond,) = structure.connections
    assert bond.kind == ConnectionKind.COVALENT
    assert bond.records == ("LINK", "CONECT")
    assert bond.reported_distance == pytest.approx(1.45)


def test_unresolvable_conect_is_reported_and_blocks_the_extraction(tmp_path: Path) -> None:
    structure = read_structure(
        _pdb(tmp_path, conect="CONECT    5   99\n" + CONECT), journal=Journal()
    )
    assert "99 absent du modèle" in structure.metadata[CONECT_UNRESOLVED]
    report = analyze_chain(structure, "A")
    assert any("connectivité du fichier incertaine" in r for r in report.refusals)


def test_conect_serials_are_resolved_in_the_model_read(tmp_path: Path) -> None:
    """Dans un ensemble RMN, un numéro de série vaut pour le modèle qui le porte."""
    body = PDB.format(link="", conect="").replace("END\n", "")
    model_2 = body.replace("C1  DOT A 101", "C1  DOT A 102").replace(
        "C2  DOT A 101", "C2  DOT A 102"
    )
    text = f"MODEL        1\n{body}ENDMDL\nMODEL        2\n{model_2}ENDMDL\n{CONECT}END\n"
    path = tmp_path / "ensemble.pdb"
    path.write_text(text, encoding="utf-8")
    first = read_structure(path, journal=Journal(), model=1)
    second = read_structure(path, journal=Journal(), model=2)
    assert [c.second.label() for c in first.connections] == ["A:DOT101"]
    assert [c.second.label() for c in second.connections] == ["A:DOT102"]


def test_conect_only_bond_is_written_back_as_conect(tmp_path: Path) -> None:
    structure = read_structure(_pdb(tmp_path, conect=CONECT), journal=Journal())
    target = write_structure(structure, tmp_path / "out.pdb", "pdb", journal=Journal())
    text = target.read_text(encoding="utf-8")
    assert "CONECT" in text and "LINK " not in text
    (bond,) = read_structure(target, journal=Journal()).connections
    assert bond.records == ("CONECT",) and bond.kind == ConnectionKind.UNKNOWN
    mmcif = write_structure(structure, tmp_path / "out.cif", "mmcif", journal=Journal())
    (again,) = read_structure(mmcif, journal=Journal()).connections
    assert again.records == ("struct_conn",)  # mmCIF n'a pas de CONECT : provenance changée


def test_link_to_a_symmetry_mate_is_not_measured(tmp_path: Path) -> None:
    link = LINK.replace("1555   1555", "1555   2555")
    structure = read_structure(_pdb(tmp_path, link=link), journal=Journal())
    (bond,) = structure.connections
    assert bond.symmetry_mate
    assessment = assess_declared(bond, _by_residue(structure.atoms), CRITERIA)
    assert assessment.status is BondStatus.INSUFFICIENT_DATA
    assert assessment.certainty is Certainty.DECLARED_UNVERIFIED


# --- Liaison déclarée confrontée aux coordonnées ---------------------------------------------


def _by_residue(atoms: tuple[Atom, ...] | list[Atom]) -> dict[ResidueId, list[Atom]]:
    table: dict[ResidueId, list[Atom]] = {}
    for atom in atoms:
        table.setdefault(atom.residue_id, []).append(atom)
    return table


def _nz_link(**kwargs: object) -> Connection:
    return Connection(
        ConnectionKind.COVALENT,
        _rid(2, "LYS"),
        "NZ",
        _rid(101, "DOT"),
        "C1",
        records=("LINK",),
        **kwargs,  # type: ignore[arg-type]
    )


def test_declared_bond_at_a_bonding_distance_is_confirmed() -> None:
    atoms = _chain(side=LYS_NZ) + _het("DOT", 101, BELOW_NZ, (0, -1, 0))
    assessment = assess_declared(_nz_link(), _by_residue(atoms), CRITERIA)
    assert assessment.status is BondStatus.ACCEPTED
    assert assessment.certainty is Certainty.CONFIRMED
    assert assessment.distance == pytest.approx(1.45)
    assert assessment.max_distance == pytest.approx(1.92)


def test_stale_link_is_contradicted_and_refused() -> None:
    far = (NZ[0], NZ[1] - 4.0, 0.0)
    atoms = _chain(side=LYS_NZ) + _het("DOT", 101, far, (0, -1, 0))
    structure = Structure(atoms=tuple(atoms), connections=(_nz_link(),))
    report = analyze_chain(structure, "A")
    (bond,) = [b for b in report.bonds if b.evidence is BondEvidence.DECLARED]
    assert bond.status is BondStatus.INCOMPATIBLE and bond.certainty is Certainty.CONTRADICTED
    assert any("contredite par les coordonnées" in r for r in report.refusals)


def test_declared_bond_to_an_absent_atom_is_insufficient_data() -> None:
    atoms = _chain(side=LYS_NZ) + _het("DOT", 101, BELOW_NZ, (0, -1, 0))
    bond = replace(_nz_link(), second_atom="C9")
    assessment = assess_declared(bond, _by_residue(atoms), CRITERIA)
    assert assessment.status is BondStatus.INSUFFICIENT_DATA
    assert "absent" in assessment.reason


def test_conformations_that_disagree_make_a_declared_bond_insufficient() -> None:
    atoms = (
        _chain(side=LYS_NZ)
        + [replace(a, alt_loc="A") for a in _het("DOT", 101, BELOW_NZ, (0, -1, 0))]
        + [
            replace(a, alt_loc="B", serial=a.serial + 50)
            for a in _het("DOT", 101, (NZ[0], -8.0, 0.0), (0, -1, 0))
        ]
    )
    assessment = assess_declared(_nz_link(), _by_residue(atoms), CRITERIA)
    assert assessment.status is BondStatus.INSUFFICIENT_DATA
    assert "conformations se contredisent" in assessment.reason


def test_declared_bond_to_a_metal_is_never_assessed_as_covalent() -> None:
    zinc = Atom(950, "ZN", "ZN", "ZN", 201, "A", NZ[0], NZ[1] - 2.1, 0.0, hetero=True)
    bond = Connection(
        ConnectionKind.UNKNOWN, _rid(2, "LYS"), "NZ", _rid(201, "ZN"), "ZN", records=("CONECT",)
    )
    structure = Structure(atoms=(*_chain(side=LYS_NZ), zinc), connections=(bond,))
    report = analyze_chain(structure, "A")
    assert _component(report, "A:ZN201").category is ComponentCategory.ION  # type: ignore[attr-defined]
    (declared,) = [b for b in report.bonds if b.evidence is BondEvidence.DECLARED]
    assert declared.status is BondStatus.NOT_ASSESSED
    assert report.accepted


# --- Détection géométrique : bornes, hydrogènes, conformations, recouvrement --------------


def _pair(
    distance: float,
    first: str = "N",
    second: str = "C",
    alt: tuple[str, str] = ("", ""),
    criteria: CovalentCriteria = CRITERIA,
) -> list:  # type: ignore[type-arg]
    a = Atom(1, first, first, "LYS", 1, "A", 0.0, 0.0, 0.0, alt_loc=alt[0] or None)
    b = Atom(
        2,
        second + "1",
        second,
        "LIG",
        2,
        "A",
        distance,
        0.0,
        0.0,
        alt_loc=alt[1] or None,
        hetero=True,
    )
    return find_contacts([a], [a, b], criteria)


@pytest.mark.parametrize(
    ("distance", "expected"),
    [
        (0.629, BondStatus.INCOMPATIBLE),
        (0.631, BondStatus.POSSIBLE),
        (1.919, BondStatus.POSSIBLE),
        (1.921, None),
    ],
)
def test_contact_window_bounds_are_exact(distance: float, expected: BondStatus | None) -> None:
    found = _pair(distance)
    if expected is None:
        assert found == []
    else:
        (contact,) = found
        assert contact.status is expected
        assert contact.evidence is BondEvidence.GEOMETRIC
        assert contact.certainty is (
            Certainty.OVERLAP if expected is BondStatus.INCOMPATIBLE else Certainty.POSSIBLE
        )


def test_hydrogens_are_never_used_for_detection() -> None:
    assert _pair(1.0, first="N", second="H") == []


def test_atoms_of_incompatible_conformations_are_not_paired() -> None:
    assert _pair(1.4, alt=("A", "B")) == []
    assert len(_pair(1.4, alt=("A", ""))) == 1


def test_contact_with_a_metal_is_undetermined_never_covalent() -> None:
    (contact,) = _pair(2.1, second="ZN")
    assert contact.evidence is BondEvidence.UNDETERMINED
    assert contact.status is BondStatus.NOT_ASSESSED
    assert contact.certainty is Certainty.UNDETERMINED


# --- Chélateurs : sans liaison déclarée, jamais perdus -------------------------------------


@pytest.mark.parametrize(
    ("start", "direction", "site", "target"),
    [
        ((-1.40, 0.0, 0.0), (-1, 0, 0), AttachmentSite.N_TERMINAL, "A:ALA1:N"),
        ((_x(5) + 2.0 + 1.40, 0.0, 0.0), (1, 0, 0), AttachmentSite.C_TERMINAL, "A:SER5:C"),
        (BELOW_NZ, (0, -1, 0), AttachmentSite.SIDE_CHAIN, "A:LYS2:NZ"),
    ],
)
def test_chelator_without_declared_bond_is_possibly_bound_and_refused(
    start: tuple[float, float, float],
    direction: tuple[float, float, float],
    site: AttachmentSite,
    target: str,
) -> None:
    atoms = _chain(side=LYS_NZ) + _het("DOT", 101, start, direction, count=5)
    report = analyze_chain(Structure(atoms=tuple(atoms)), "A")
    chelator = _component(report, "A:DOT101")
    assert chelator.category is ComponentCategory.POSSIBLY_BOUND  # type: ignore[attr-defined]
    assert chelator.decision is Decision.REJECT  # type: ignore[attr-defined]
    assert chelator.attachment is site  # type: ignore[attr-defined]
    (contact,) = [b for b in report.bonds if b.evidence is BondEvidence.GEOMETRIC]
    assert {contact.first.label(), contact.second.label()} == {target, "A:DOT101:C1"}
    assert contact.certainty is Certainty.POSSIBLE
    assert not report.accepted
    assert any("A:DOT101" in r and "jamais exclu en silence" in r for r in report.refusals)


def test_chelator_with_declared_link_is_a_modification_with_its_site() -> None:
    atoms = _chain(side=LYS_NZ) + _het("DOT", 101, BELOW_NZ, (0, -1, 0), count=5)
    structure = Structure(atoms=tuple(atoms), connections=(_nz_link(),))
    report = analyze_chain(structure, "A")
    chelator = _component(report, "A:DOT101")
    assert chelator.category is ComponentCategory.SIDE_CHAIN_MODIFICATION  # type: ignore[attr-defined]
    assert chelator.attachment is AttachmentSite.SIDE_CHAIN  # type: ignore[attr-defined]
    assert not [b for b in report.bonds if b.evidence is BondEvidence.GEOMETRIC]
    assert not report.accepted


def test_possibly_bound_component_can_be_excluded_only_by_name(tmp_path: Path) -> None:
    atoms = _chain(side=LYS_NZ) + _het("DOT", 101, BELOW_NZ, (0, -1, 0))
    structure = Structure(atoms=tuple(atoms))
    decisions = ExtractionDecisions(components={"A:DOT101": "exclude_as_unbound"})
    journal = Journal()
    result = extract_peptides(structure, "A", journal=journal, decisions=decisions)
    chelator = _component(result.report, "A:DOT101")
    assert chelator.decision is Decision.EXCLUDE_AS_UNBOUND  # type: ignore[attr-defined]
    (peptide,) = result.peptides
    assert peptide.metadata[STRUCTURE_SET_ASIDE] == "A:DOT101=exclude_as_unbound"
    logged = journal.entries[-1].parameters
    assert any(c["decision"] == "exclude_as_unbound" for c in logged["components"])
    assert any(b["certainty"] == "possible" for b in logged["bonds"])


def test_declared_component_can_be_removed_only_by_name() -> None:
    atoms = _chain(side=LYS_NZ) + _het("NAG", 101, BELOW_NZ, (0, -1, 0))
    bond = replace(_nz_link(), second=_rid(101, "NAG"))
    structure = Structure(atoms=tuple(atoms), connections=(bond,))
    decisions = ExtractionDecisions(components={"A:NAG101": "remove_covalent"})
    result = extract_peptides(structure, "A", journal=Journal(), decisions=decisions)
    assert result.peptides[0].metadata[STRUCTURE_SET_ASIDE] == "A:NAG101=remove_covalent"
    assert [b.certainty for b in result.report.bonds] == [Certainty.CONFIRMED]


def test_a_decision_must_match_the_category_of_its_component() -> None:
    atoms = _chain(side=LYS_NZ) + _het("DOT", 101, BELOW_NZ, (0, -1, 0))
    declared = Structure(atoms=tuple(atoms), connections=(_nz_link(),))
    with pytest.raises(ConfigurationError, match="side_chain_modification"):
        analyze_chain(
            declared,
            "A",
            decisions=ExtractionDecisions(components={"A:DOT101": "exclude_as_unbound"}),
        )
    geometric = Structure(atoms=tuple(atoms))
    with pytest.raises(ConfigurationError, match="possibly_bound"):
        analyze_chain(
            geometric,
            "A",
            decisions=ExtractionDecisions(components={"A:DOT101": "remove_covalent"}),
        )
    with pytest.raises(ConfigurationError, match="n'est pas un composant"):
        analyze_chain(
            geometric,
            "A",
            decisions=ExtractionDecisions(components={"A:DOT999": "exclude_as_unbound"}),
        )


def test_component_of_another_chain_is_examined() -> None:
    atoms = _chain(side=LYS_NZ) + _het("DOT", 101, BELOW_NZ, (0, -1, 0), chain="Z")
    report = analyze_chain(Structure(atoms=tuple(atoms)), "A")
    chelator = _component(report, "Z:DOT101")
    assert chelator.category is ComponentCategory.POSSIBLY_BOUND  # type: ignore[attr-defined]
    assert not report.accepted


def test_declared_bond_to_a_component_of_another_chain_attaches_it() -> None:
    atoms = _chain(side=LYS_NZ) + _het("DOT", 101, BELOW_NZ, (0, -1, 0), chain="Z")
    bond = replace(_nz_link(), second=_rid(101, "DOT", "Z"))
    report = analyze_chain(Structure(atoms=tuple(atoms), connections=(bond,)), "A")
    assert _component(report, "Z:DOT101").category is ComponentCategory.SIDE_CHAIN_MODIFICATION  # type: ignore[attr-defined]


# --- Linkers ---------------------------------------------------------------------------------


def test_linker_of_several_residues_is_followed_to_its_end() -> None:
    atoms = (
        _chain(side=LYS_NZ)
        + _het("LK1", 101, BELOW_NZ, (0, -1, 0), count=2, serial=900)
        + _het("LK2", 102, (NZ[0], BELOW_NZ[1] - 3.0, 0.0), (0, -1, 0), count=2, serial=910)
        + _het("LK3", 103, (NZ[0], BELOW_NZ[1] - 6.0, 0.0), (0, -1, 0), count=2, serial=920)
    )
    report = analyze_chain(Structure(atoms=tuple(atoms)), "A")
    sites = {c.residue: (c.category, c.attachment) for c in report.components}
    assert sites == {
        "A:LK1101": (ComponentCategory.POSSIBLY_BOUND, AttachmentSite.SIDE_CHAIN),
        "A:LK2102": (ComponentCategory.POSSIBLY_BOUND, AttachmentSite.COMPONENT_ONLY),
        "A:LK3103": (ComponentCategory.POSSIBLY_BOUND, AttachmentSite.COMPONENT_ONLY),
    }


def test_component_touching_two_residues_is_bridging() -> None:
    side = {2: ("NZ", "N"), 3: ("NZ", "N")}
    a = (_x(2) + 1.0, -4.45, 0.0)
    atoms = [
        *_chain("ALA LYS LYS GLY", side=side),
        Atom(900, "C1", "C", "BRG", 101, "A", a[0], a[1], 0.0, hetero=True),
        Atom(901, "C2", "C", "BRG", 101, "A", a[0] + STEP / 2, a[1], 0.0, hetero=True),
        Atom(902, "C3", "C", "BRG", 101, "A", a[0] + STEP, a[1], 0.0, hetero=True),
    ]
    report = analyze_chain(Structure(atoms=tuple(atoms)), "A")
    assert _component(report, "A:BRG101").attachment is AttachmentSite.BRIDGING  # type: ignore[attr-defined]


# --- Ligands non liés et ions ----------------------------------------------------------------


def test_ligand_out_of_bonding_distance_stays_unbound_without_question() -> None:
    near = (NZ[0], NZ[1] - 2.5, 0.0)  # 2,5 Å : contact, pas liaison (fenêtre N–C 1,92 Å)
    atoms = _chain(side=LYS_NZ) + _het("LIG", 101, near, (0, -1, 0))
    result = extract_peptides(Structure(atoms=tuple(atoms)), "A", journal=Journal())
    ligand = _component(result.report, "A:LIG101")
    assert ligand.category is ComponentCategory.UNBOUND_LIGAND  # type: ignore[attr-defined]
    assert ligand.min_distance_to_chain == pytest.approx(2.5)  # type: ignore[attr-defined]
    assert result.report.bonds == ()


def test_ion_near_the_chain_stays_an_ion_and_its_contact_is_reported() -> None:
    zinc = Atom(950, "ZN", "ZN", "ZN", 201, "A", NZ[0], NZ[1] - 2.1, 0.0, hetero=True)
    result = extract_peptides(Structure(atoms=(*_chain(side=LYS_NZ), zinc)), "A", journal=Journal())
    assert _component(result.report, "A:ZN201").category is ComponentCategory.ION  # type: ignore[attr-defined]
    (contact,) = result.report.bonds
    assert contact.certainty is Certainty.UNDETERMINED


# --- Liaisons intra-chaîne et entre chaînes --------------------------------------------------


def _undeclared_disulfide() -> Structure:
    side = {2: ("SG", "S"), 5: ("SG", "S")}
    atoms = _chain("ALA CYS GLY SER CYS", side=side)
    moved = [
        replace(a, x=_x(2) + 1.0, y=-3.0, z=2.04) if (a.residue_number, a.name) == (5, "SG") else a
        for a in atoms
    ]
    return Structure(atoms=tuple(moved))


def test_undeclared_intra_chain_contact_is_refused() -> None:
    report = analyze_chain(_undeclared_disulfide(), "A")
    (contact,) = [b for b in report.bonds if b.evidence is BondEvidence.GEOMETRIC]
    assert contact.label() == "A:CYS2:SG–A:CYS5:SG"
    assert contact.distance == pytest.approx(2.04)
    assert any("dans la chaîne" in r and "A:CYS2–A:CYS5" in r for r in report.refusals)


def test_undeclared_disulfide_enters_the_peptide_by_an_explicit_decision() -> None:
    decisions = ExtractionDecisions(contacts={"A:CYS2–A:CYS5": "disulfide"})
    result = extract_peptides(_undeclared_disulfide(), "A", journal=Journal(), decisions=decisions)
    (peptide,) = result.peptides
    assert [(link.kind, link.positions) for link in peptide.links] == [(LinkKind.DISULFIDE, (2, 5))]


def test_not_bonded_decision_accepts_the_contact_without_a_link() -> None:
    decisions = ExtractionDecisions(contacts={"A:CYS2–A:CYS5": "not_bonded"})
    (peptide,) = extract_peptides(
        _undeclared_disulfide(), "A", journal=Journal(), decisions=decisions
    ).peptides
    assert peptide.links == ()


def test_disulfide_decision_on_a_contact_that_is_not_sg_sg_is_refused() -> None:
    side = {2: ("NZ", "N"), 5: ("OG", "O")}
    atoms = _chain("ALA LYS GLY SER SER", side=side)
    moved = tuple(
        replace(a, x=_x(2) + 1.0, y=-3.0, z=1.4) if (a.residue_number, a.name) == (5, "OG") else a
        for a in atoms
    )
    decisions = ExtractionDecisions(contacts={"A:LYS2–A:SER5": "disulfide"})
    report = analyze_chain(Structure(atoms=moved), "A", decisions=decisions)
    assert any("n'est pas un contact SG–SG" in r for r in report.refusals)


def test_contact_with_a_residue_of_another_chain_is_refused() -> None:
    receptor = [replace(a, z=a.z + 20.0) for a in _chain("GLY", chain="R")]
    receptor = [
        replace(a, x=BELOW_NZ[0], y=BELOW_NZ[1], z=0.0) if a.name == "CA" else a for a in receptor
    ]
    structure = Structure(atoms=(*_chain(side=LYS_NZ), *receptor))
    report = analyze_chain(structure, "A")
    assert any("avec une autre chaîne" in r and "A:LYS2–R:GLY1" in r for r in report.refusals)
    assert not [c for c in report.components if c.residue.startswith("R:")]
    lifted = analyze_chain(
        structure, "A", decisions=ExtractionDecisions(contacts={"A:LYS2–R:GLY1": "not_bonded"})
    )
    assert lifted.accepted


def test_contact_decision_without_object_is_refused() -> None:
    with pytest.raises(ConfigurationError, match="n'est pas un contact"):
        analyze_chain(
            Structure(atoms=tuple(_chain())),
            "A",
            decisions=ExtractionDecisions(contacts={"A:ALA1–A:SER5": "not_bonded"}),
        )


# --- Peptide standard : rien ne change --------------------------------------------------------


def test_standard_peptide_has_no_contact_and_is_extracted_as_before() -> None:
    result = extract_peptides(Structure(atoms=tuple(_chain(side=LYS_NZ))), "A", journal=Journal())
    assert result.report.bonds == () and result.report.accepted
    (peptide,) = result.peptides
    assert STRUCTURE_SET_ASIDE not in peptide.metadata


# --- Décisions : forme de configuration -------------------------------------------------------


def test_component_and_contact_decisions_round_trip_through_a_dictionary() -> None:
    data = {
        "breaks": "reject",
        "unobserved_termini": "reject",
        "components": {"A:DOT101": "exclude_as_unbound"},
        "contacts": {"A:CYS2–A:CYS5": "disulfide"},
    }
    decisions = ExtractionDecisions.from_dict(data)
    assert decisions.to_dict() == data
    assert ExtractionDecisions().to_dict() == {"breaks": "reject", "unobserved_termini": "reject"}


@pytest.mark.parametrize(
    "data",
    [
        {"components": {"A:DOT101": "accept"}},
        {"components": {"A:DOT101": "not_bonded"}},
        {"contacts": {"A:CYS2–A:CYS5": "exclude_as_unbound"}},
        {"components": ["A:DOT101"]},
    ],
)
def test_invalid_component_or_contact_decisions_are_refused(data: dict[str, object]) -> None:
    with pytest.raises(ConfigurationError):
        ExtractionDecisions.from_dict(data)


# --- Rapport : chaque liaison documentée ------------------------------------------------------


def test_each_bond_of_the_report_names_both_ends_distance_source_and_certainty() -> None:
    atoms = _chain(side=LYS_NZ) + _het("DOT", 101, BELOW_NZ, (0, -1, 0))
    report = analyze_chain(Structure(atoms=tuple(atoms)), "A")
    data = report.to_dict()
    yaml.safe_dump(data, allow_unicode=True)  # sérialisable tel quel
    (bond,) = data["bonds"]
    for end in (bond["first"], bond["second"]):
        assert set(end) >= {"component", "chain", "residue", "number", "atom", "element", "alt_loc"}
    assert bond["distance"] == pytest.approx(1.45)
    assert bond["evidence"] == "geometric" and bond["certainty"] == "possible"
    assert bond["status"] == "possible" and bond["reason"]
    assert data["covalent_criteria"]["tolerance"] == 0.45
    (component,) = data["components"]
    assert component["attachment"] == "side_chain"


# --- Prévention des pertes silencieuses -------------------------------------------------------


def test_every_residue_of_the_structure_near_the_chain_is_accounted_for() -> None:
    """Chaque composant au contact de la chaîne est soit rapporté et refusé, soit écarté par
    une décision nommée : aucun ne disparaît entre la structure et le peptide."""
    atoms = (
        _chain(side=LYS_NZ)
        + _het("DOT", 101, BELOW_NZ, (0, -1, 0))
        + _het("DOT", 102, (-1.40, 0.0, 0.0), (-1, 0, 0), chain="Z", serial=930)
        + _het("LIG", 103, (0.0, 30.0, 0.0), (0, 1, 0), serial=960)
    )
    report = analyze_chain(Structure(atoms=tuple(atoms)), "A")
    by_label = {c.residue: c for c in report.components}
    assert set(by_label) == {"A:DOT101", "Z:DOT102", "A:LIG103"}
    for label in ("A:DOT101", "Z:DOT102"):
        assert by_label[label].decision is Decision.REJECT
        assert any(label in reason for reason in report.refusals)
    assert by_label["A:LIG103"].category is ComponentCategory.UNBOUND_LIGAND


DETECTION_MODULES = ("io/connectivity.py", "io/peptide_extraction.py")


@pytest.mark.parametrize("module", DETECTION_MODULES)
def test_no_numeric_threshold_is_written_in_the_detection_code(module: str) -> None:
    """Les seuils viennent des fichiers de données, jamais du code (décision 3)."""
    path = Path(__file__).resolve().parents[1] / "src" / "pepforge" / module
    tree = ast.parse(path.read_text(encoding="utf-8"))
    floats = [
        node.lineno
        for node in ast.walk(tree)
        if isinstance(node, ast.Constant) and isinstance(node.value, float)
    ]
    assert floats == [], f"{module} : constantes décimales lignes {floats}"


def test_a_geometric_contact_is_never_turned_into_a_declared_connection() -> None:
    """La géométrie n'est jamais une preuve : aucun code de détection ne construit une
    `Connection`, qui est réservée à ce que le fichier déclare."""
    root = Path(__file__).resolve().parents[1] / "src" / "pepforge"
    for module in DETECTION_MODULES:
        tree = ast.parse((root / module).read_text(encoding="utf-8"))
        built = [
            node.lineno
            for node in ast.walk(tree)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "Connection"
        ]
        assert built == [], f"{module} construit une Connection lignes {built}"


def test_refused_extraction_names_every_uncertain_bond() -> None:
    atoms = _chain(side=LYS_NZ) + _het("DOT", 101, BELOW_NZ, (0, -1, 0))
    with pytest.raises(ExtractionRefusedError) as info:
        extract_peptides(Structure(atoms=tuple(atoms)), "A", journal=Journal())
    report = info.value.report
    # Les deux extrémités sont rangées par identifiant de résidu (chaîne puis numéro).
    assert [b.label() for b in report.uncertain_bonds()] == ["A:LYS2:NZ–A:DOT101:C1"]  # type: ignore[union-attr]


# --- Tolérance de 0,45 Å : bornes exactes et sensibilité (revue de Firas du 28/11/2026) ---
#
# La tolérance est une convention de détection (Open Babel), validée provisoirement par
# Firas. Ces tests fixent le comportement aux bornes de la fenêtre, montrent comment le
# résultat varie autour de 0,45 Å, et vérifient qu'aucune valeur n'est lue ailleurs que
# dans `covalent_geometry.yaml`.

HEAVY_PAIRS = [("N", "C"), ("C", "C"), ("O", "C"), ("S", "S"), ("S", "C"), ("P", "O"), ("SE", "C")]


@pytest.mark.parametrize(("first", "second"), HEAVY_PAIRS)
def test_window_bounds_are_inclusive_for_every_element_pair(first: str, second: str) -> None:
    """Fenêtre fermée : d = limite est une liaison possible, le flottant suivant ne l'est plus."""
    limit = CRITERIA.max_distance(first, second)
    assert limit == pytest.approx(
        CRITERIA.radii[first] + CRITERIA.radii[second] + CRITERIA.tolerance
    )
    (at_limit,) = _pair(limit, first, second)
    assert at_limit.status is BondStatus.POSSIBLE
    assert at_limit.max_distance == pytest.approx(limit)
    assert _pair(math.nextafter(limit, math.inf), first, second) == []
    (at_minimum,) = _pair(CRITERIA.min_distance, first, second)
    assert at_minimum.status is BondStatus.POSSIBLE
    (below,) = _pair(math.nextafter(CRITERIA.min_distance, 0.0), first, second)
    assert below.certainty is Certainty.OVERLAP


@pytest.mark.parametrize(
    ("offset", "expected"),
    [
        (-0.001, BondStatus.ACCEPTED),
        (0.0, BondStatus.ACCEPTED),
        (0.001, BondStatus.INCOMPATIBLE),
    ],
)
def test_declared_bond_follows_the_same_window(offset: float, expected: BondStatus) -> None:
    """Une liaison déclarée est confirmée dans la fenêtre, contredite juste au-delà."""
    limit = CRITERIA.max_distance("N", "C")
    atoms = _chain(side=LYS_NZ) + _het("DOT", 101, (NZ[0], NZ[1] - limit - offset, 0.0), (0, -1, 0))
    assessment = assess_declared(_nz_link(), _by_residue(atoms), CRITERIA)
    assert assessment.status is expected
    assert assessment.distance == pytest.approx(limit + offset, abs=1e-6)


def test_declared_bond_under_the_minimum_distance_is_contradicted() -> None:
    below = CRITERIA.min_distance - 0.01
    atoms = _chain(side=LYS_NZ) + _het("DOT", 101, (NZ[0], NZ[1] - below, 0.0), (0, -1, 0))
    assessment = assess_declared(_nz_link(), _by_residue(atoms), CRITERIA)
    assert assessment.status is BondStatus.INCOMPATIBLE
    assert assessment.certainty is Certainty.CONTRADICTED


SENSITIVITY = (0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60)


def test_detection_widens_monotonically_with_the_tolerance() -> None:
    """Sensibilité : un contact trouvé à une tolérance l'est encore à toute tolérance plus large."""
    distances = [round(1.40 + 0.02 * k, 2) for k in range(41)]  # 1,40 à 2,20 Å, paire N–C
    found: list[set[float]] = []
    for tolerance in SENSITIVITY:
        criteria = replace(CRITERIA, tolerance=tolerance)
        found.append({d for d in distances if _pair(d, criteria=criteria)})
    for narrower, wider in pairwise(found):
        assert narrower <= wider
    limit_045 = CRITERIA.radii["N"] + CRITERIA.radii["C"] + 0.45
    assert max(found[SENSITIVITY.index(0.45)]) <= limit_045
    assert found[0] != found[-1]  # la tolérance change bien le résultat : elle n'est pas ignorée


@pytest.mark.parametrize("tolerance", SENSITIVITY)
def test_contact_exactly_at_the_reference_limit_depends_only_on_the_tolerance(
    tolerance: float,
) -> None:
    """Un contact à r1 + r2 + 0,45 Å n'est vu que si la tolérance atteint 0,45 Å."""
    distance = CRITERIA.radii["S"] + CRITERIA.radii["S"] + 0.45
    contacts = _pair(distance, "S", "S", criteria=replace(CRITERIA, tolerance=tolerance))
    assert bool(contacts) is (tolerance >= 0.45 - 1e-9)


@pytest.mark.parametrize("tolerance", SENSITIVITY)
def test_a_chelator_at_a_normal_bond_length_is_caught_whatever_the_tolerance(
    tolerance: float,
) -> None:
    """Le chélateur d'essai, à 1,45 Å de NZ, reste possibly_bound sur toute la plage étudiée."""
    atoms = _chain(side=LYS_NZ) + _het("DOT", 101, BELOW_NZ, (0, -1, 0))
    report = analyze_chain(
        Structure(atoms=tuple(atoms)), "A", covalent=replace(CRITERIA, tolerance=tolerance)
    )
    assert _component(report, "A:DOT101").category is ComponentCategory.POSSIBLY_BOUND  # type: ignore[attr-defined]
    assert not report.accepted


@pytest.mark.parametrize(("first", "second"), HEAVY_PAIRS)
def test_a_geometric_contact_is_never_confirmed_at_any_distance(first: str, second: str) -> None:
    """Balayage de 0,05 à 3,00 Å : la géométrie seule ne donne jamais « confirmée »."""
    for step in range(1, 61):
        for contact in _pair(0.05 * step, first, second):
            assert contact.evidence is BondEvidence.GEOMETRIC
            assert contact.certainty in (Certainty.POSSIBLE, Certainty.OVERLAP)
            assert contact.certainty is not Certainty.CONFIRMED


def test_a_disulfide_decision_does_not_make_the_geometric_contact_confirmed() -> None:
    """La décision ajoute un pont au peptide ; le rapport garde la preuve : géométrie seule."""
    structure = _undeclared_disulfide()
    decisions = ExtractionDecisions(contacts={"A:CYS2–A:CYS5": "disulfide"})
    result = extract_peptides(structure, "A", journal=Journal(), decisions=decisions)
    (contact,) = [b for b in result.report.bonds if b.label() == "A:CYS2:SG–A:CYS5:SG"]
    assert contact.evidence is BondEvidence.GEOMETRIC
    assert contact.certainty is Certainty.POSSIBLE
    assert structure.connections == ()  # rien n'est déclaré à la place du fichier


def test_the_tolerance_is_read_from_the_data_file_only() -> None:
    """La valeur vient de `covalent_geometry.yaml` (le code n'a aucune constante décimale)."""
    path = Path(__file__).resolve().parents[1] / "src/pepforge/residues/data/covalent_geometry.yaml"
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert CRITERIA.tolerance == raw["detection"]["tolerance"] == 0.45
    assert CRITERIA.min_distance == raw["detection"]["min_distance"]
    assert "Open Babel" in CRITERIA.tolerance_source
